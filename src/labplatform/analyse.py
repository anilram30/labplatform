"""One call that runs every analysis of the platform and returns a JSON-serialisable summary."""
from __future__ import annotations

import json
from pathlib import Path

from .control import run_all_charts
from .drift import analyse_drift
from .rnr import repeatability_reproducibility, verdicts_at_risk
from .roundrobin import (
    DEFAULT_QUANTITIES,
    artefact_stability,
    compare_all,
    compare_by_round,
    en_vs_frequency,
)
from .schema import PlatformDB
from .traceability import suspect_runs
from .uncertainty import evaluate_budgets

__all__ = ["analyse_network", "site_overview"]


def site_overview(db: PlatformDB) -> list[dict]:
    out = []
    for s in db.q("SELECT * FROM sites ORDER BY id"):
        site = s["id"]
        n = {r["trust"] or r["state"]: r["c"] for r in db.q("SELECT COALESCE(trust, state) AS trust, COUNT(*) AS c FROM runs WHERE site=? GROUP BY COALESCE(trust, state)", (site,))}
        v = {r["verdict"]: r["c"] for r in db.q("SELECT verdict, COUNT(*) AS c FROM runs WHERE site=? AND trust IN ('trusted','flagged') GROUP BY verdict", (site,))}
        last_ver = db.q("SELECT time, status, max_dev_db FROM verifications WHERE site=? ORDER BY t_unix DESC LIMIT 1", (site,))
        last_cal = db.q("SELECT record_id, date FROM calibrations WHERE site=? ORDER BY date DESC LIMIT 1", (site,))
        inst = db.q("SELECT model, serial FROM instruments WHERE site=? AND role='vna'", (site,))
        ev = db.q("SELECT severity, COUNT(*) AS c FROM events WHERE site=? GROUP BY severity", (site,))
        out.append({"site": site, "name": s["name"], "instrument": f"{inst[0]['model']} {inst[0]['serial']}" if inst else "", "runs": sum(n.values()),
                    "trusted": n.get("trusted", 0), "flagged": n.get("flagged", 0), "rejected": n.get("rejected", 0),
                    "aborted": n.get("ABORTED", 0) + n.get("ERROR", 0), "pass": v.get("PASS", 0), "fail": v.get("FAIL", 0),
                    "last_verification": dict(last_ver[0]) if last_ver else None, "last_calibration": dict(last_cal[0]) if last_cal else None,
                    "events": {e["severity"]: e["c"] for e in ev}})
    return out


def analyse_network(db: PlatformDB, artefact_id: str, reference: str = "consensus", phase1: int = 5, quantities=None) -> dict:
    db.con.execute("DELETE FROM events")
    n_budget = evaluate_budgets(db, artefact_id)
    comparisons = compare_all(db, artefact_id, quantities or DEFAULT_QUANTITIES, reference)
    from datetime import datetime, timezone
    now = datetime.now(timezone.utc).isoformat(timespec="seconds")
    for c in comparisons:
        for s in c.sites:
            if abs(s.en) > 1:
                db.event(now, "En", s.site, c.label, "action" if abs(s.en) > 1.5 else "warning",
                         f"{c.label}: E_n = {s.en:+.2f} (x = {s.x:.4g}, ref = {c.x_ref:.4g} {c.unit or ''})", s.to_dict())
    enf = en_vs_frequency(db, artefact_id, reference=reference)
    by_round = compare_by_round(db, artefact_id, "insertion_loss", "A", 600e6, reference)
    for r in by_round["rounds"]:
        for site, v in r["sites"].items():
            if abs(v["En"]) > 1:
                db.event(datetime.fromtimestamp(r["t_unix"], tz=timezone.utc).isoformat(timespec="seconds"), "En-round", site, f"{r['batch']} IL@600",
                         "action" if abs(v["En"]) > 1.5 else "warning", f"{r['batch']}: E_n(IL @ 600 MHz) = {v['En']:+.2f}, d = {v['d']:+.4f} dB (U = {2 * v['u']:.3f} dB)", v)
    rnr = {}
    for q, pair, f in (quantities or DEFAULT_QUANTITIES):
        r = repeatability_reproducibility(db, artefact_id, q, pair, f)
        if r is not None:
            rnr[r.label] = r.to_dict()
    il600 = rnr.get("insertion_loss[A] @ 600 MHz")
    risk = verdicts_at_risk(db, il600["anova"]["s_R"]) if il600 else None
    charts = run_all_charts(db, artefact_id, phase1)
    stability = {f"{q}@{int(f / 1e6) if f else 'scalar'}": artefact_stability(db, artefact_id, q, pair, f)
                 for q, pair, f in [("insertion_loss", "A", 600e6), ("insertion_loss", "A", 100e6), ("impedance_fitted", "A", None), ("loop_resistance", None, None)]}
    suspects = suspect_runs(db)
    drift = analyse_drift(db, artefact_id, by_round, phase1)
    for f in drift:
        if f.drifting:
            k = f.change.get("index")
            when = f.labels[k] if (k is not None and f.change.get("significant")) else (f.labels[f.ewma_first] if f.ewma_first is not None else f.labels[-1])
            db.event(now, "drift", f.site, f.series, "action" if f.change.get("significant") else "warning", f.conclusion, {"round": when})
    # trust card of one artefact result per site (the latest trusted round) for the dashboard
    from .trust import trust_card
    cards = {}
    for site in db.sites():
        r = db.q("SELECT id FROM runs WHERE site=? AND sample_id=? AND trust='trusted' ORDER BY t_unix DESC LIMIT 1", (site, artefact_id))
        if r:
            try:
                cards[site] = trust_card(db, r[0]["id"], "insertion_loss", "A", 600e6, artefact_id, reference, suspects)
            except KeyError:
                pass
    db.commit()
    return {"artefact": artefact_id, "reference": reference, "budgets_evaluated": n_budget, "overview": site_overview(db),
            "comparisons": [c.to_dict() for c in comparisons], "en_vs_frequency": enf, "by_round": by_round, "rnr": rnr, "verdicts_at_risk": risk,
            "charts": [c.to_dict() for c in charts], "stability": stability, "suspect_runs": suspects,
            "drift": [f.to_dict() for f in drift], "trust_cards": cards,
            "events": [dict(e) for e in db.q("SELECT time, kind, site, subject, severity, message FROM events ORDER BY time")],
            "counts": db.counts()}


def write_summary(summary: dict, path: str | Path) -> Path:
    path = Path(path)
    slim = dict(summary)
    slim["charts"] = [{k: v for k, v in c.items() if k not in ("t", "x", "labels", "z")} for c in summary["charts"]]
    path.write_text(json.dumps(slim, indent=1, default=float), encoding="utf-8")
    return path
