"""
The canonical data model: one schema for every site.

    sites          id, name, timezone, accreditation scope
    instruments    per site: role, manufacturer, model, serial, firmware
    calibrations   per instrument: record id, type, ports, date, kit serial, kit due, verification history (JSON)
    samples        sample id, kind (production | artefact), part, lot, cable type, nominal length, design
    procedures     id, version, sha256, cable type, text
    runs           one measurement job: site, instrument, calibration, sample, procedure, operator,
                   timestamps, temperature and its method, ambient, trust, verdict, headline margin,
                   attempt count, archive path + manifest hash, software versions (JSON), source
    results        scalar quantities per run: quantity, pair, frequency (NULL for scalars), value, unit,
                   standard uncertainty (NULL until a budget is evaluated), verdict margin if any
    traces         optional full traces per run (quantity, pair, x[], value[] as JSON)
    verifications  calibration verifications: site, instrument, calibration, time, standard, status,
                   max_dev_db, max_dev_deg, min_rl_db
    events         what the platform concluded and when: control-chart alerts, En failures, invalidations

Everything a cross-site comparison, a control chart or an audit needs is a query over these tables.
Rows carry the *source* (archive path and manifest hash) so every number is one hop from its raw file.
"""
from __future__ import annotations

import json
import sqlite3
from pathlib import Path

__all__ = ["PlatformDB"]

_SCHEMA = """
CREATE TABLE IF NOT EXISTS sites (
    id TEXT PRIMARY KEY, name TEXT, timezone TEXT, scope TEXT
);
CREATE TABLE IF NOT EXISTS instruments (
    id INTEGER PRIMARY KEY,
    site TEXT NOT NULL REFERENCES sites(id), role TEXT, manufacturer TEXT, model TEXT, serial TEXT, firmware TEXT,
    UNIQUE(site, role, serial)
);
CREATE TABLE IF NOT EXISTS calibrations (
    id INTEGER PRIMARY KEY,
    site TEXT NOT NULL, instrument_serial TEXT, record_id TEXT, type TEXT, ports TEXT, date TEXT,
    temperature_c REAL, operator TEXT, kit_serial TEXT, kit_due TEXT, verifications_json TEXT,
    UNIQUE(site, instrument_serial, record_id)
);
CREATE TABLE IF NOT EXISTS samples (
    id TEXT PRIMARY KEY, kind TEXT, part_number TEXT, lot TEXT, cable_type TEXT, length_m REAL, design TEXT, notes TEXT
);
CREATE TABLE IF NOT EXISTS procedures (
    sha256 TEXT PRIMARY KEY, id TEXT, version TEXT, kind TEXT, cable_type TEXT, text TEXT
);
CREATE TABLE IF NOT EXISTS runs (
    id TEXT PRIMARY KEY,
    site TEXT NOT NULL REFERENCES sites(id),
    instrument_serial TEXT, calibration_record TEXT, sample_id TEXT REFERENCES samples(id),
    procedure_sha TEXT REFERENCES procedures(sha256), procedure_id TEXT,
    operator TEXT, campaign TEXT, batch TEXT,
    started TEXT, finished TEXT, t_unix REAL,
    temperature_c REAL, temperature_method TEXT, ambient_c REAL, humidity_pct REAL,
    state TEXT, trust TEXT, verdict TEXT, headline TEXT, headline_margin REAL, attempts INTEGER,
    cal_age_days REAL, cal_delta_t_k REAL, verification_status TEXT, verification_dev_db REAL,
    trace_noise_db REAL, ifbw_hz REAL, averages INTEGER, power_dbm REAL, points INTEGER,
    archive_dir TEXT, manifest_sha256 TEXT, software_json TEXT, source TEXT
);
CREATE TABLE IF NOT EXISTS results (
    id INTEGER PRIMARY KEY,
    run_id TEXT NOT NULL REFERENCES runs(id),
    quantity TEXT NOT NULL, pair TEXT, f_hz REAL, value REAL NOT NULL, unit TEXT,
    u_std REAL, u_json TEXT, margin REAL
);
CREATE TABLE IF NOT EXISTS traces (
    id INTEGER PRIMARY KEY,
    run_id TEXT NOT NULL REFERENCES runs(id), quantity TEXT, pair TEXT, unit TEXT, x_json TEXT, value_json TEXT
);
CREATE TABLE IF NOT EXISTS verifications (
    id INTEGER PRIMARY KEY,
    site TEXT NOT NULL, instrument_serial TEXT, calibration_record TEXT, time TEXT, t_unix REAL,
    standard TEXT, status TEXT, max_dev_db REAL, max_dev_deg REAL, min_rl_db REAL,
    UNIQUE(site, instrument_serial, calibration_record, time)
);
CREATE TABLE IF NOT EXISTS events (
    id INTEGER PRIMARY KEY,
    time TEXT, kind TEXT, site TEXT, subject TEXT, severity TEXT, message TEXT, data_json TEXT
);
CREATE INDEX IF NOT EXISTS ix_runs_site ON runs(site);
CREATE INDEX IF NOT EXISTS ix_runs_sample ON runs(sample_id);
CREATE INDEX IF NOT EXISTS ix_results_run ON results(run_id);
CREATE INDEX IF NOT EXISTS ix_results_q ON results(quantity, pair, f_hz);
"""


class PlatformDB:
    def __init__(self, path: str | Path):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.con = sqlite3.connect(self.path)
        self.con.row_factory = sqlite3.Row
        self.con.executescript(_SCHEMA)

    def close(self):
        self.con.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    # ---- writes
    def upsert(self, table: str, row: dict, key: str | list[str] = "id") -> None:
        cols = list(row)
        self.con.execute(f"INSERT OR REPLACE INTO {table} ({','.join(cols)}) VALUES ({','.join('?' * len(cols))})",
                         [row[c] for c in cols])

    def insert_results(self, run_id: str, rows: list[dict]) -> None:
        self.con.execute("DELETE FROM results WHERE run_id = ?", (run_id,))
        self.con.executemany("INSERT INTO results (run_id, quantity, pair, f_hz, value, unit, u_std, u_json, margin) VALUES (?,?,?,?,?,?,?,?,?)",
                             [(run_id, r["quantity"], r.get("pair"), r.get("f_hz"), float(r["value"]), r.get("unit"),
                               r.get("u_std"), json.dumps(r["u_json"]) if r.get("u_json") else None, r.get("margin")) for r in rows])

    def insert_traces(self, run_id: str, rows: list[dict]) -> None:
        self.con.execute("DELETE FROM traces WHERE run_id = ?", (run_id,))
        self.con.executemany("INSERT INTO traces (run_id, quantity, pair, unit, x_json, value_json) VALUES (?,?,?,?,?,?)",
                             [(run_id, r["quantity"], r.get("pair"), r.get("unit"), json.dumps(r["x"]), json.dumps(r["value"])) for r in rows])

    def event(self, time: str, kind: str, site: str | None, subject: str, severity: str, message: str, data: dict | None = None) -> None:
        self.con.execute("INSERT INTO events (time, kind, site, subject, severity, message, data_json) VALUES (?,?,?,?,?,?,?)",
                         (time, kind, site, subject, severity, message, json.dumps(data or {}, default=str)))

    def commit(self):
        self.con.commit()

    # ---- reads
    def q(self, sql: str, args=()) -> list[sqlite3.Row]:
        return list(self.con.execute(sql, args))

    def sites(self) -> list[str]:
        return [r["id"] for r in self.q("SELECT id FROM sites ORDER BY id")]

    def results_matrix(self, quantity: str, pair: str | None, f_hz: float | None, sample_id: str | None = None,
                       trust: tuple[str, ...] = ("trusted", "flagged")) -> list[sqlite3.Row]:
        """Every (site, run, value, u) of one scalar quantity, for comparisons."""
        sql = ("SELECT r.site, r.id AS run_id, r.t_unix, r.started, r.batch, r.sample_id, r.operator, r.instrument_serial, r.calibration_record, "
               "r.temperature_c, r.trust, s.value, s.u_std FROM results s JOIN runs r ON r.id = s.run_id "
               "WHERE s.quantity = ? AND (s.pair IS ? OR s.pair = ?) AND r.trust IN (%s)" % ",".join("?" * len(trust)))
        args: list = [quantity, pair, pair, *trust]
        if f_hz is None:
            sql += " AND s.f_hz IS NULL"
        else:
            sql += " AND ABS(s.f_hz - ?) < 1"
            args.append(f_hz)
        if sample_id:
            sql += " AND r.sample_id = ?"
            args.append(sample_id)
        sql += " ORDER BY r.t_unix"
        return self.q(sql, args)

    def counts(self) -> dict:
        out = {}
        for t in ("sites", "instruments", "calibrations", "samples", "procedures", "runs", "results", "traces", "verifications", "events"):
            out[t] = self.con.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0]
        return out
