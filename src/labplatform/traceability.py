"""
Traceability for accreditation: the chain behind every number, and what a failure invalidates.

    chain(run_id)    result -> run -> raw files and their hashes (from the archive manifest) -> procedure
                     (id, version, sha256) -> instrument (identity) -> calibration record (kit serial and
                     due date, verification against the check standard and its reference file) -> site
                     -> software versions.  Every link is a stored value, not a lookup that could change.

    suspect_runs()   the retrospective question every ISO/IEC 17025 laboratory must be able to answer:
                     when a calibration verification fails, which results since the last *passed*
                     verification on that instrument were produced under a calibration that may already
                     have been wrong?  Those runs are marked as events, not deleted - the decision to
                     re-measure, recall or annotate a certificate is a human one.

    evidence_pack()  a JSON dossier for one run or one sample: the chain, the platform's uncertainty
                     budget, the comparison scores of the site at the time, and the archive integrity
                     state, ready to attach to a test report or hand to an assessor.
"""
from __future__ import annotations

import json
from pathlib import Path

from .schema import PlatformDB

__all__ = ["chain", "suspect_runs", "evidence_pack"]


def chain(db: PlatformDB, run_id: str) -> dict:
    run = db.q("SELECT * FROM runs WHERE id=?", (run_id,))
    if not run:
        raise KeyError(run_id)
    run = dict(run[0])
    site = db.q("SELECT * FROM sites WHERE id=?", (run["site"],))
    inst = db.q("SELECT * FROM instruments WHERE site=? AND serial=?", (run["site"], run["instrument_serial"]))
    cal = db.q("SELECT * FROM calibrations WHERE site=? AND record_id=?", (run["site"], run["calibration_record"]))
    proc = db.q("SELECT id, version, kind, cable_type, sha256 FROM procedures WHERE sha256=?", (run["procedure_sha"],))
    sample = db.q("SELECT * FROM samples WHERE id=?", (run["sample_id"],))
    results = [dict(r) for r in db.q("SELECT quantity, pair, f_hz, value, unit, u_std, u_json, margin FROM results WHERE run_id=? ORDER BY quantity, f_hz", (run_id,))]
    for r in results:
        r["budget"] = json.loads(r.pop("u_json")) if r.get("u_json") else None
    files, integrity = [], None
    man_p = Path(run["archive_dir"]) / "manifest.json" if run.get("archive_dir") else None
    if man_p and man_p.exists():
        man = json.loads(man_p.read_text(encoding="utf-8"))
        files = [{"path": k, **v} for k, v in man.get("files", {}).items()]
        try:
            from labauto.archive import Archive
            integrity = Archive(man_p.parents[2]).verify(man_p.parent)
        except Exception as e:  # pragma: no cover
            integrity = {"ok": None, "error": str(e)}
    verifs = [dict(v) for v in db.q("SELECT time, standard, status, max_dev_db, max_dev_deg, min_rl_db FROM verifications WHERE site=? AND calibration_record=? ORDER BY t_unix",
                                    (run["site"], run["calibration_record"]))]
    return {"run": {k: v for k, v in run.items() if k != "software_json"}, "software": json.loads(run.get("software_json") or "{}"),
            "site": dict(site[0]) if site else None, "instrument": dict(inst[0]) if inst else None,
            "calibration": {k: v for k, v in dict(cal[0]).items() if k != "verifications_json"} | {"verifications": verifs} if cal else None,
            "procedure": dict(proc[0]) if proc else None, "sample": dict(sample[0]) if sample else None,
            "files": files, "archive_integrity": integrity, "results": results}


def suspect_runs(db: PlatformDB, record_events: bool = True) -> list[dict]:
    """Runs whose instrument's next verification (after the run) failed, with no passed verification in between."""
    out = []
    seen = set()
    for inst in db.q("SELECT DISTINCT site, instrument_serial FROM verifications"):
        vs = db.q("SELECT t_unix, time, status, calibration_record, max_dev_db FROM verifications WHERE site=? AND instrument_serial=? ORDER BY t_unix",
                  (inst["site"], inst["instrument_serial"]))
        for k, v in enumerate(vs):
            if v["status"] != "fail":
                continue
            t_prev_pass = max([p["t_unix"] for p in vs[:k] if p["status"] == "pass"], default=0.0)
            # a job starts a few seconds before the verification it triggers, hence the 60 s tolerance
            runs = db.q("SELECT id, sample_id, started, trust, verdict, headline_margin FROM runs WHERE site=? AND instrument_serial=? AND t_unix > ? AND t_unix < ? "
                        "AND trust IN ('trusted','flagged')", (inst["site"], inst["instrument_serial"], t_prev_pass - 60.0, v["t_unix"]))
            for r in runs:
                if r["id"] in seen:
                    continue
                seen.add(r["id"])
                d = {"site": inst["site"], "instrument_serial": inst["instrument_serial"], "run_id": r["id"], "sample_id": r["sample_id"], "started": r["started"],
                     "verdict": r["verdict"], "headline_margin": r["headline_margin"], "failed_verification": v["time"], "failed_calibration": v["calibration_record"],
                     "max_dev_db": v["max_dev_db"]}
                out.append(d)
                if record_events:
                    db.event(v["time"], "suspect-run", inst["site"], r["id"], "action",
                             f"run {r['id']} ({r['sample_id']}) measured between the last passed verification and the failed one at {v['time']} "
                             f"(dA = {v['max_dev_db']:.3f} dB): review", d)
    if record_events:
        db.commit()
    return out


def evidence_pack(db: PlatformDB, run_id: str, out_dir: str | Path, comparisons: list[dict] | None = None) -> Path:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    c = chain(db, run_id)
    site = c["run"]["site"]
    if comparisons:
        c["site_comparison_scores"] = [{"label": cm["label"], **next((s for s in cm["sites"] if s["site"] == site), {})} for cm in comparisons]
    c["events"] = [dict(e) for e in db.q("SELECT time, kind, severity, message FROM events WHERE subject=? OR (site=? AND kind='control-chart') ORDER BY time",
                                          (run_id, site))]
    p = out_dir / f"evidence_{run_id}.json"
    p.write_text(json.dumps(c, indent=1, default=str), encoding="utf-8")
    return p
