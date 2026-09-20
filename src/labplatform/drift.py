"""
Longitudinal analysis: is the measurement system itself drifting, and why?

For each site and each artefact series (round means of IL at the comparison frequencies, worst-case
return loss, fitted impedance; the verification deviation dA; the degree of equivalence d_i per round)
three detectors run:

    trend         ordinary least squares slope with its p-value (a slow drift)
    change point  single change point by maximum two-sample t statistic over all splits (binary
                  segmentation, first level) with a Bonferroni-corrected threshold: the round at which
                  the mean shifted, with the shift and its significance
    EWMA          the first round at which the EWMA of :mod:`control` signals

A series is *drifting* when any detector fires.  Then comes the question a laboratory manager actually
asks - what changed? - and the answer comes from project B's metadata, which the platform holds per
round for the same site:

    calibration   the recalibration date (every round in the demonstration), the verification dA and
                  return loss, the kit serial
    environment   ambient temperature and humidity at the sample, temperature label method
    people/proc   operator, procedure hash, instrument firmware
    product       the site's *production* results over the same rounds

Root-cause attribution is evidence, not proof: for each candidate the platform reports (i) whether it
changed at the same round as the detected change point, (ii) the Pearson correlation between the series
and the candidate's numeric series over the rounds, and (iii) a plain-language hypothesis.  The
artefact/production discrimination is the important one: if the artefact shifted *and* the production
samples shifted together, the measurement system did; if only production shifted, the product did -
which is where project E's manufacturing correlation takes over.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field

import numpy as np
from scipy import stats as sps

from .control import ewma
from .schema import PlatformDB
from .stats import linear_trend

__all__ = ["change_point", "DriftFinding", "site_series", "analyse_drift"]


def change_point(x, min_seg: int = 2, alpha: float = 0.05) -> dict:
    """Single change point by the maximum |t| over all splits with at least ``min_seg`` points per side."""
    x = np.asarray(x, float)
    n = x.size
    if n < 2 * min_seg + 1:
        return {"index": None, "shift": 0.0, "t": 0.0, "p": 1.0, "significant": False}
    best = (None, 0.0, 0.0)
    for k in range(min_seg, n - min_seg + 1):
        a, b = x[:k], x[k:]
        sp = np.sqrt(((a.size - 1) * a.var(ddof=1) + (b.size - 1) * b.var(ddof=1)) / (n - 2)) if n > 2 else 0.0
        if sp <= 0:
            continue
        t = (b.mean() - a.mean()) / (sp * np.sqrt(1 / a.size + 1 / b.size))
        if abs(t) > abs(best[1]):
            best = (k, float(t), float(b.mean() - a.mean()))
    if best[0] is None:
        return {"index": None, "shift": 0.0, "t": 0.0, "p": 1.0, "significant": False}
    p = float(2 * sps.t.sf(abs(best[1]), n - 2))
    n_tests = n - 2 * min_seg + 1
    return {"index": best[0], "shift": best[2], "t": best[1], "p": p, "p_bonferroni": min(1.0, p * n_tests),
            "significant": bool(p * n_tests < alpha)}


@dataclass
class DriftFinding:
    site: str
    series: str
    unit: str | None
    n: int
    labels: list[str]
    x: list[float]
    trend: dict
    change: dict
    ewma_first: int | None
    drifting: bool
    causes: list[dict] = field(default_factory=list)
    conclusion: str = ""

    def to_dict(self) -> dict:
        return {"site": self.site, "series": self.series, "unit": self.unit, "n": self.n, "labels": self.labels, "x": self.x,
                "trend": self.trend, "change_point": self.change, "ewma_first_signal": self.ewma_first, "drifting": self.drifting,
                "causes": self.causes, "conclusion": self.conclusion}


def site_series(db: PlatformDB, site: str, artefact_id: str, by_round: dict | None = None) -> tuple[list[str], dict, dict]:
    """(round labels, {series name: (unit, values by round)}, {metadata name: values by round}) for one site."""
    rows = db.q("SELECT r.batch, r.t_unix, r.id, r.operator, r.calibration_record, r.verification_dev_db, r.ambient_c, r.humidity_pct, r.temperature_c, "
                "r.procedure_sha, r.software_json FROM runs r WHERE r.site=? AND r.sample_id=? AND r.trust IN ('trusted','flagged') AND r.batch LIKE 'RR-%' ORDER BY r.t_unix",
                (site, artefact_id))
    if not rows:
        return [], {}, {}
    batches = sorted({r["batch"] for r in rows})
    series: dict[str, tuple] = {}
    for q, pair, f, name in [("insertion_loss", "A", 100e6, "IL @ 100 MHz"), ("insertion_loss", "A", 600e6, "IL @ 600 MHz"),
                             ("return_loss_min", "A", None, "worst-case RL"), ("impedance_fitted", "A", None, "fitted impedance")]:
        vals = []
        for b in batches:
            v = db.q("SELECT s.value FROM results s JOIN runs r ON r.id=s.run_id WHERE r.site=? AND r.sample_id=? AND r.batch=? AND s.quantity=? AND (s.pair IS ? OR s.pair=?) "
                     "AND (s.f_hz IS ? OR ABS(s.f_hz-?)<1) AND r.trust IN ('trusted','flagged')", (site, artefact_id, b, q, pair, pair, f, f if f is not None else 0))
            vals.append(float(np.mean([x["value"] for x in v])) if v else np.nan)
        unit = db.q("SELECT unit FROM results WHERE quantity=? LIMIT 1", (q,))
        series[name] = (unit[0]["unit"] if unit else None, vals)
    ver = []
    for b in batches:
        r = [x for x in rows if x["batch"] == b]
        ver.append(float(np.nanmean([x["verification_dev_db"] for x in r if x["verification_dev_db"] is not None])) if r else np.nan)
    series["verification dA"] = ("dB", ver)
    if by_round:
        d = [next((rr["sites"][site]["d"] for rr in by_round["rounds"] if rr["batch"] == b and site in rr["sites"]), np.nan) for b in batches]
        series["degree of equivalence d_i (IL @ 600 MHz)"] = ("dB", d)
    # metadata per round
    meta: dict[str, list] = {"ambient_c": [], "humidity_pct": [], "sample_temperature_c": [], "operator": [], "calibration": [], "procedure": [],
                             "firmware": [], "kit_serial": [], "verification_rl_db": []}
    for b in batches:
        r = [x for x in rows if x["batch"] == b]
        meta["ambient_c"].append(float(np.nanmean([x["ambient_c"] for x in r if x["ambient_c"] is not None])) if r else np.nan)
        meta["humidity_pct"].append(float(np.nanmean([x["humidity_pct"] for x in r if x["humidity_pct"] is not None])) if r else np.nan)
        meta["sample_temperature_c"].append(float(np.nanmean([x["temperature_c"] for x in r if x["temperature_c"] is not None])) if r else np.nan)
        meta["operator"].append(r[0]["operator"] if r else None)
        meta["calibration"].append(r[0]["calibration_record"] if r else None)
        meta["procedure"].append((r[0]["procedure_sha"] or "")[:12] if r else None)
        sw = json.loads(r[0]["software_json"] or "{}") if r else {}
        meta["firmware"].append(sw.get("labauto"))
        cal = db.q("SELECT kit_serial FROM calibrations WHERE site=? AND record_id=?", (site, r[0]["calibration_record"])) if r else []
        meta["kit_serial"].append(cal[0]["kit_serial"] if cal else None)
        rl = db.q("SELECT min_rl_db FROM verifications WHERE site=? AND calibration_record=? ORDER BY t_unix LIMIT 1", (site, r[0]["calibration_record"])) if r else []
        meta["verification_rl_db"].append(rl[0]["min_rl_db"] if rl else np.nan)
    # production samples of the site per round: mean IL @ 600 MHz of its production runs in the same fortnight
    prod = []
    for b in batches:
        rr = [x for x in rows if x["batch"] == b]
        t0 = min(x["t_unix"] for x in rr)
        v = db.q("SELECT s.value FROM results s JOIN runs r ON r.id=s.run_id JOIN samples sm ON sm.id=r.sample_id WHERE r.site=? AND sm.kind='production' "
                 "AND r.t_unix BETWEEN ? AND ? AND s.quantity='insertion_loss' AND ABS(s.f_hz-600e6)<1 AND r.trust IN ('trusted','flagged')", (site, t0 - 3600, t0 + 7 * 86400))
        prod.append(float(np.mean([x["value"] for x in v])) if v else np.nan)
    meta["production IL @ 600 MHz"] = prod
    return batches, series, meta


def _corr(a, b) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 4 or np.std(a[m]) == 0 or np.std(b[m]) == 0:
        return float("nan")
    return float(np.corrcoef(a[m], b[m])[0, 1])


def attribute(finding: DriftFinding, meta: dict, series: dict) -> None:
    k = finding.change.get("index")
    causes = []
    # categorical changes at the change point
    for name in ("operator", "procedure", "firmware", "kit_serial"):
        vals = meta.get(name, [])
        if k is not None and 0 < k < len(vals) and vals[k] != vals[k - 1]:
            causes.append({"candidate": name, "kind": "categorical change at the change point", "from": vals[k - 1], "to": vals[k], "score": 1.0})
    cal = meta.get("calibration", [])
    n_recal = sum(1 for a, b in zip(cal[:-1], cal[1:]) if a != b)
    routine = len(cal) > 2 and n_recal >= 0.5 * (len(cal) - 1)
    if k is not None and 0 < k < len(cal) and cal[k] != cal[k - 1]:
        causes.append({"candidate": "recalibration", "kind": "routine recalibration every round: uninformative" if routine else "the change point coincides with a recalibration",
                       "from": cal[k - 1], "to": cal[k], "score": 0.0 if routine else 0.7})
    # numeric correlations over rounds
    for name in ("verification dA",):
        if name in series and finding.series != name:
            r = _corr(finding.x, series[name][1])
            if np.isfinite(r):
                causes.append({"candidate": "calibration / test cable (verification dA)", "kind": "correlation over rounds", "r": r, "score": abs(r)})
    for name, label in (("verification_rl_db", "calibration (verification return loss)"), ("ambient_c", "ambient temperature"),
                        ("humidity_pct", "humidity"), ("sample_temperature_c", "sample temperature")):
        r = _corr(finding.x, meta.get(name, []))
        if np.isfinite(r):
            causes.append({"candidate": label, "kind": "correlation over rounds", "r": r, "score": abs(r)})
    # measurement system vs product (only meaningful for an artefact quantity, not for the check-standard series)
    artefact_series = finding.series != "verification dA"
    prod = meta.get("production IL @ 600 MHz", [])
    r = _corr(finding.x, prod) if artefact_series else float("nan")
    cp_prod = change_point([v for v in prod if np.isfinite(v)]) if sum(np.isfinite(prod)) >= 5 else {"significant": False}
    if artefact_series:
        causes.append({"candidate": "production samples move with the artefact", "kind": "artefact vs production", "r": r if np.isfinite(r) else None,
                       "production_change_point": cp_prod.get("significant", False), "score": abs(r) if np.isfinite(r) else 0.0})
    causes.sort(key=lambda c: -c["score"])
    finding.causes = causes
    # conclusion
    if not finding.drifting:
        finding.conclusion = "no drift detected"
        return
    top = [c for c in causes if c["score"] >= 0.6 and c["candidate"] != "production samples move with the artefact"]
    lead = top[0]["candidate"] if top else ("no metadata candidate (|r| < 0.6, no categorical change" + ("; recalibration is routine" if routine else "") + ")")
    when = f"at round {finding.labels[k]} (shift {finding.change['shift']:+.4g} {finding.unit or ''})" if k is not None and finding.change.get("significant") else \
        (f"as a trend of {finding.trend.get('slope_per_30d', float('nan')):+.4g} {finding.unit or ''}/30 d" if finding.trend.get("significant") else "flagged by EWMA only")
    if not artefact_series:
        system = "instrument-side evidence from the check standard, independent of the artefact: test cable / connectors / calibration kit"
    else:
        prod_c = next((c for c in causes if c["kind"] == "artefact vs production"), {})
        if prod_c.get("r") is not None and prod_c["r"] > 0.6:
            system = "the production results moved with the artefact: the measurement system changed, not the product"
        elif prod_c.get("r") is not None and prod_c["r"] < 0.2:
            system = "the production results did not follow the artefact: a measurement-system change confined to the artefact hook-up, or the artefact itself"
        else:
            system = "artefact/production relation inconclusive"
    finding.conclusion = f"{finding.series} at {finding.site} drifts {when}; leading candidate: {lead}; {system}."


def analyse_drift(db: PlatformDB, artefact_id: str, by_round: dict | None = None, phase1: int = 5) -> list[DriftFinding]:
    out = []
    for site in db.sites():
        batches, series, meta = site_series(db, site, artefact_id, by_round)
        if len(batches) < 5:
            continue
        for name, (unit, vals) in series.items():
            x = np.asarray(vals, float)
            m = np.isfinite(x)
            if m.sum() < 5:
                continue
            xs, labels = x[m], [b for b, ok in zip(batches, m) if ok]
            t = np.arange(xs.size) * 14.0
            tr = linear_trend(t, xs)
            tr["slope_per_30d"] = tr.pop("slope") * 30 if np.isfinite(tr.get("slope", np.nan)) else float("nan")
            tr["significant"] = bool(np.isfinite(tr["p"]) and tr["p"] < 0.05)
            cp = change_point(xs)
            base = xs[:min(phase1, xs.size)]
            sigma = float(np.mean(np.abs(np.diff(base))) / 1.128) if base.size > 1 else 1e-9
            sigma = sigma if sigma > 0 else (float(np.std(base, ddof=1)) if base.size > 1 else 1e-9)
            e = ewma(xs, float(base.mean()), sigma, start=min(phase1, xs.size))
            drifting = bool(tr["significant"] or cp["significant"] or e["first_signal"] is not None)
            f = DriftFinding(site, name, unit, int(xs.size), labels, xs.tolist(), tr, cp, e["first_signal"], drifting)
            attribute(f, {k: [v for v, ok in zip(vs, m) if ok] for k, vs in meta.items()},
                      {k: (u, [v for v, ok in zip(vs, m) if ok]) for k, (u, vs) in series.items()})
            out.append(f)
    return out
