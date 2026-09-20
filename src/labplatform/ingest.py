"""
Ingest: labauto archives (the primary source), cablecheck result files, plain CSV.

From a labauto laboratory root the ingester reads ``lab.toml`` (site id and name), every sealed job
under ``archive/`` (manifest, sidecars, cablecheck result), every calibration record, and turns them
into rows of the canonical schema.  Scalar *quantities* are derived from the cablecheck traces:

    insertion_loss[pair] at the platform's comparison frequencies (10, 100, 300, 600 MHz) - dB
    return_loss / lcl / lctl / next : worst value over the limit span (min) and the value at 100 MHz
    impedance_mean, impedance_fitted, impedance_min, impedance_max, nvp, delay_per_metre : scalars
    loop_resistance : from the auxiliary measurement, ohm

Quarantined (rejected) runs are ingested too - trust is a column, and comparisons filter on it - because
a platform that only sees accepted data cannot report rejection rates.
"""
from __future__ import annotations

import json
import tomllib
from datetime import datetime
from pathlib import Path

import numpy as np

from .schema import PlatformDB

__all__ = ["ingest_lab", "ingest_job", "COMPARISON_FREQS_HZ", "derive_quantities"]

COMPARISON_FREQS_HZ = (10e6, 100e6, 300e6, 600e6)
_MIN_QUANTITIES = ("return_loss", "lcl", "lctl", "next", "fext")


def _t_unix(iso: str | None) -> float | None:
    if not iso:
        return None
    try:
        return datetime.fromisoformat(iso.replace("Z", "+00:00")).timestamp()
    except ValueError:
        return None


def derive_quantities(cc_result: dict, aux: dict | None = None) -> tuple[list[dict], list[dict]]:
    """(scalar rows, trace rows) from a cablecheck result JSON (with traces)."""
    rows, traces = [], []
    margins = {(r["quantity"], r["pair"]): r.get("worst_margin") for r in cc_result.get("results", [])}
    for t in cc_result.get("traces", []):
        q, pair, unit = t["quantity"], t["pair"], t.get("unit")
        x, v = np.asarray(t["x"], float), np.asarray(t["value"], float)
        if t.get("xunit", "") == "" or x.size <= 1:
            if v.size:
                rows.append({"quantity": q, "pair": pair, "f_hz": None, "value": float(v[0]), "unit": unit, "margin": margins.get((q, pair))})
            continue
        traces.append({"quantity": q, "pair": pair, "unit": unit, "x": [float(a) for a in x], "value": [float(a) for a in v]})
        if q == "insertion_loss":
            for f in COMPARISON_FREQS_HZ:
                if x[0] <= f <= x[-1]:
                    rows.append({"quantity": q, "pair": pair, "f_hz": f, "value": float(np.interp(f, x, v)), "unit": unit})
            rows.append({"quantity": q, "pair": pair, "f_hz": None, "value": float(v[-1]), "unit": unit, "margin": margins.get((q, pair))})
        elif q in _MIN_QUANTITIES:
            lim = t.get("limit")
            m = np.isfinite(v)
            if lim is not None:
                lm = np.asarray([np.nan if a is None else a for a in lim], float)
                m &= np.isfinite(lm)
            if m.any():
                rows.append({"quantity": q + "_min", "pair": pair, "f_hz": None, "value": float(np.min(v[m])), "unit": unit, "margin": margins.get((q, pair))})
            if x[0] <= 100e6 <= x[-1]:
                rows.append({"quantity": q, "pair": pair, "f_hz": 100e6, "value": float(np.interp(100e6, x, v)), "unit": unit})
        elif q == "impedance_profile":
            pass   # the windowed scalars carry it
        elif q in ("group_delay", "phase_delay"):
            if x[0] <= 100e6 <= x[-1]:
                rows.append({"quantity": q, "pair": pair, "f_hz": 100e6, "value": float(np.interp(100e6, x, v)), "unit": unit})
    if aux:
        for k, v in aux.items():
            if isinstance(v, (int, float)) and np.isfinite(v):
                rows.append({"quantity": k, "pair": None, "f_hz": None, "value": float(v), "unit": "ohm" if "resistance" in k else None})
    return rows, traces


def ingest_job(db: PlatformDB, job_dir: Path, site: str, source: str = "labauto") -> str | None:
    man_p = job_dir / "manifest.json"
    if not man_p.exists():
        return None
    man = json.loads(man_p.read_text(encoding="utf-8"))
    job = man["job"]
    metas = sorted(job_dir.glob("*.meta.json"))
    if not metas:
        return None
    meta = json.loads(metas[0].read_text(encoding="utf-8"))
    sample = job.get("sample") or meta.get("sample") or {}
    vna = meta.get("instrument", {}).get("vna", {})
    ident = vna.get("identity", {})
    cal = meta.get("calibration", {})
    env = meta.get("environment", {})
    sweep = vna.get("sweep", {})
    checks = {c["name"]: c for c in meta.get("validation", {}).get("checks", [])}
    proc = meta.get("procedure", {})
    from hashlib import sha256
    manifest_sha = sha256(man_p.read_bytes()).hexdigest()
    run_id = job["job_id"]
    db.upsert("samples", {"id": sample.get("sample_id"), "kind": "artefact" if str(sample.get("notes", "")).startswith("artefact") or "GOLD" in str(sample.get("sample_id", "")) else "production",
                          "part_number": sample.get("part_number"), "lot": sample.get("lot"), "cable_type": sample.get("cable_type"),
                          "length_m": sample.get("length_m"), "design": sample.get("design"), "notes": sample.get("notes")})
    ptext = (job_dir / "procedure.toml").read_text(encoding="utf-8") if (job_dir / "procedure.toml").exists() else ""
    db.upsert("procedures", {"sha256": proc.get("hash"), "id": proc.get("id"), "version": proc.get("version"), "kind": proc.get("kind"),
                             "cable_type": proc.get("cable_type"), "text": ptext})
    db.upsert("instruments", {"site": site, "role": "vna", "manufacturer": ident.get("manufacturer"), "model": ident.get("model"),
                              "serial": ident.get("serial"), "firmware": ident.get("firmware")})
    if cal.get("id"):
        db.upsert("calibrations", {"site": site, "instrument_serial": ident.get("serial"), "record_id": cal.get("id"), "type": cal.get("type"),
                                   "ports": json.dumps(cal.get("ports")), "date": cal.get("date"), "temperature_c": cal.get("temperature_c"),
                                   "operator": cal.get("operator"), "kit_serial": cal.get("kit_serial"), "kit_due": cal.get("kit_due"),
                                   "verifications_json": json.dumps(cal.get("verification"))})
        v = cal.get("verification") or {}
        if v.get("time"):
            db.con.execute("INSERT OR IGNORE INTO verifications (site, instrument_serial, calibration_record, time, t_unix, standard, status, max_dev_db, max_dev_deg, min_rl_db) "
                           "VALUES (?,?,?,?,?,?,?,?,?,?)", (site, ident.get("serial"), cal.get("id"), v["time"], _t_unix(v["time"]), v.get("standard"),
                                                           v.get("status"), v.get("max_dev_db"), v.get("max_dev_deg"), v.get("min_rl_db")))
    analysis = job.get("analysis", {}).get("cablecheck", {})
    db.upsert("runs", {
        "id": run_id, "site": site, "instrument_serial": ident.get("serial"), "calibration_record": cal.get("id"),
        "sample_id": sample.get("sample_id"), "procedure_sha": proc.get("hash"), "procedure_id": proc.get("id"),
        "operator": job.get("operator"), "campaign": job.get("campaign"), "batch": job.get("batch"),
        "started": job.get("started"), "finished": job.get("finished"), "t_unix": _t_unix(job.get("started")),
        "temperature_c": env.get("sample_temperature_c"), "temperature_method": env.get("sample_temperature_method"),
        "ambient_c": env.get("ambient_c"), "humidity_pct": env.get("humidity_pct"),
        "state": job.get("state"), "trust": job.get("trust"), "verdict": analysis.get("verdict"), "headline": analysis.get("headline"),
        "headline_margin": analysis.get("headline_margin"), "attempts": job.get("attempts"),
        "cal_age_days": cal.get("age_days"), "cal_delta_t_k": cal.get("delta_t_k"),
        "verification_status": (cal.get("verification") or {}).get("status"), "verification_dev_db": (cal.get("verification") or {}).get("max_dev_db"),
        "trace_noise_db": checks.get("trace_noise", {}).get("value"), "ifbw_hz": sweep.get("ifbw_hz"), "averages": sweep.get("averages"),
        "power_dbm": sweep.get("power_dbm"), "points": sweep.get("points"),
        "archive_dir": str(job_dir), "manifest_sha256": manifest_sha, "software_json": json.dumps(meta.get("software", {})), "source": source,
    })
    res_p = job_dir / "analysis" / "cablecheck_result.json"
    aux = json.loads((job_dir / "aux.json").read_text(encoding="utf-8")) if (job_dir / "aux.json").exists() else {}
    aux = {k: v for k, v in aux.items() if not k.endswith("_state")}
    if res_p.exists():
        rows, traces = derive_quantities(json.loads(res_p.read_text(encoding="utf-8")), aux)
        db.insert_results(run_id, rows)
        db.insert_traces(run_id, traces)
    elif aux:
        db.insert_results(run_id, [{"quantity": k, "pair": None, "f_hz": None, "value": float(v), "unit": "ohm"} for k, v in aux.items() if isinstance(v, (int, float))])
    return run_id


def ingest_lab(db: PlatformDB, lab_root: str | Path, site: str | None = None) -> dict:
    """Ingest every sealed job of a labauto laboratory; returns counts."""
    lab_root = Path(lab_root)
    cfg = tomllib.loads((lab_root / "lab.toml").read_text(encoding="utf-8")).get("lab", {})
    site = site or cfg.get("id", lab_root.name)
    db.upsert("sites", {"id": site, "name": cfg.get("name", site), "timezone": cfg.get("timezone", "UTC"), "scope": cfg.get("scope", "")})
    archive = lab_root / cfg.get("archive", "archive")
    n = 0
    for man in sorted(archive.rglob("manifest.json")):
        if ingest_job(db, man.parent, site):
            n += 1
    # jobs that never produced data (refused at a gate, errors): from the laboratory database, so refusal rates are visible
    import sqlite3
    labdb = lab_root / cfg.get("db", "lab.sqlite")
    n_refused = 0
    if labdb.exists():
        con = sqlite3.connect(labdb)
        con.row_factory = sqlite3.Row
        for j in con.execute("SELECT * FROM jobs WHERE state IN ('ABORTED','ERROR')"):
            summ = json.loads(j["summary_json"] or "{}")
            db.upsert("runs", {"id": j["job_id"], "site": site, "instrument_serial": None, "calibration_record": j["calibration_id"], "sample_id": j["sample_id"],
                               "procedure_sha": j["procedure_hash"], "procedure_id": j["procedure_id"], "operator": j["operator"], "campaign": j["campaign"],
                               "batch": j["batch"], "started": j["started"], "finished": j["finished"], "t_unix": _t_unix(j["started"]), "state": j["state"],
                               "trust": None, "verdict": None, "headline": None, "headline_margin": None, "attempts": 0,
                               "archive_dir": None, "manifest_sha256": None, "software_json": json.dumps({"reason": summ.get("reason")}), "source": "labdb"})
            if j["sample_id"]:
                db.con.execute("INSERT OR IGNORE INTO samples (id, kind) VALUES (?, ?)", (j["sample_id"], "artefact" if "GOLD" in j["sample_id"] else "production"))
            n_refused += 1
        con.close()
    # calibration records not (yet) referenced by a run
    cal_dir = lab_root / cfg.get("calibration_dir", "calibration")
    ncal = 0
    for p in sorted(cal_dir.glob("*.json")) if cal_dir.exists() else []:
        rec = json.loads(p.read_text(encoding="utf-8"))
        db.upsert("calibrations", {"site": site, "instrument_serial": rec.get("instrument_serial"), "record_id": rec.get("id"), "type": rec.get("type"),
                                   "ports": json.dumps(rec.get("ports")), "date": rec.get("date"), "temperature_c": rec.get("temperature_c"),
                                   "operator": rec.get("operator"), "kit_serial": rec.get("kit_serial"), "kit_due": rec.get("kit_due"),
                                   "verifications_json": json.dumps(rec.get("verifications", []))})
        for v in rec.get("verifications", []):
            db.con.execute("INSERT OR IGNORE INTO verifications (site, instrument_serial, calibration_record, time, t_unix, standard, status, max_dev_db, max_dev_deg, min_rl_db) "
                           "VALUES (?,?,?,?,?,?,?,?,?,?)", (site, rec.get("instrument_serial"), rec.get("id"), v["time"], _t_unix(v["time"]), v.get("standard"),
                                                           v.get("status"), v.get("max_dev_db"), v.get("max_dev_deg"), v.get("min_rl_db")))
        ncal += 1
    db.commit()
    return {"site": site, "runs": n, "refused": n_refused, "calibrations": ncal}
