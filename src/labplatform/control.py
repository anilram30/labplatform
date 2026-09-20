"""
Control charts on what each site measures repeatedly: the artefact, the check standard, the noise.

Individuals chart (Shewhart I-MR) on a series x_1..x_N with a phase-I window of the first m points:

    centre  CL = mean(x_1..x_m)
    sigma   sigma_hat = mean(|x_k - x_{k-1}|) / d2,   d2 = 1.128 for a moving range of 2
    limits  CL +- 3 sigma_hat (action), CL +- 2 sigma_hat (warning)

with the Western Electric rules evaluated over the whole series: (1) one point beyond 3 sigma,
(2) two of three beyond 2 sigma on the same side, (3) four of five beyond 1 sigma on the same side,
(4) eight in a row on one side.  An EWMA chart on the same series,

    z_k = lambda x_k + (1 - lambda) z_{k-1},  z_0 = CL,  lambda = 0.2
    limits  CL +- L sigma_hat sqrt( lambda / (2 - lambda) (1 - (1 - lambda)^{2k}) ),  L = 2.7

detects small sustained shifts that Shewhart misses.  Phase-I points set the limits and are not
judged; rule violations and EWMA signals are evaluated from the first phase-II point on, and each
becomes an event with the site, the series and the point.

Rational subgroups matter: repeats within a round share a calibration, so charting single runs
would set limits from the re-connection scatter alone (a few mdB) and flag every recalibration.
Series charted per site are therefore: the artefact's insertion loss at the comparison frequencies
and its fitted impedance as *round means* (one point per round, temperature-corrected; the within-round
scatter is the repeatability of the R&R analysis); the calibration verification deviation dA and return
loss (one point per verification); and the trace-noise estimate of every artefact run.  Rules 1 and 2
and the EWMA signal become events; rules 3 and 4 are annotated on the chart only.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .schema import PlatformDB

__all__ = ["Chart", "individuals_chart", "ewma", "chart_series", "run_all_charts"]

D2 = 1.128


@dataclass
class Chart:
    site: str
    series: str
    unit: str | None
    t: list[float]
    x: list[float]
    labels: list[str]
    cl: float
    sigma: float
    phase1: int
    ucl: float
    lcl: float
    violations: list[dict] = field(default_factory=list)
    ewma: dict = field(default_factory=dict)

    @property
    def n(self) -> int:
        return len(self.x)

    def to_dict(self) -> dict:
        return {"site": self.site, "series": self.series, "unit": self.unit, "n": self.n, "cl": self.cl, "sigma": self.sigma, "phase1": self.phase1,
                "ucl": self.ucl, "lcl": self.lcl, "violations": self.violations, "ewma": {k: v for k, v in self.ewma.items() if k != "z"},
                "t": self.t, "x": self.x, "labels": self.labels, "z": self.ewma.get("z")}


def western_electric(x: np.ndarray, cl: float, sigma: float, start: int = 0) -> list[dict]:
    """Rule violations at indices >= start (phase II); phase-I points set the limits and are not judged."""
    out = []
    if sigma <= 0 or x.size == 0:
        return out
    z = (x - cl) / sigma
    n = x.size
    for k in range(start, n):
        if abs(z[k]) > 3:
            out.append({"rule": 1, "index": k, "message": f"point beyond 3 sigma (z = {z[k]:+.2f})"})
    for k in range(max(2, start), n):
        w = z[k - 2:k + 1]
        for sgn in (1, -1):
            if np.sum(sgn * w > 2) >= 2:
                out.append({"rule": 2, "index": k, "message": "two of three beyond 2 sigma" + (" (high)" if sgn > 0 else " (low)")})
                break
    for k in range(max(4, start), n):
        w = z[k - 4:k + 1]
        for sgn in (1, -1):
            if np.sum(sgn * w > 1) >= 4:
                out.append({"rule": 3, "index": k, "message": "four of five beyond 1 sigma" + (" (high)" if sgn > 0 else " (low)")})
                break
    for k in range(max(7, start), n):
        w = z[k - 7:k + 1]
        if np.all(w > 0) or np.all(w < 0):
            out.append({"rule": 4, "index": k, "message": "eight in a row on one side"})
    # keep the first index at which each rule fires plus later distinct events spaced by >= 3 points
    dedup, last = [], {}
    for v in sorted(out, key=lambda d: (d["index"], d["rule"])):
        if v["rule"] not in last or v["index"] - last[v["rule"]] >= 3:
            dedup.append(v)
            last[v["rule"]] = v["index"]
    return dedup


def ewma(x: np.ndarray, cl: float, sigma: float, lam: float = 0.2, L: float = 2.7, start: int = 0) -> dict:
    z = np.empty(x.size)
    prev = cl
    for k in range(x.size):
        prev = lam * x[k] + (1 - lam) * prev
        z[k] = prev
    k = np.arange(1, x.size + 1)
    half = L * sigma * np.sqrt(lam / (2 - lam) * (1 - (1 - lam) ** (2 * k)))
    out = [int(i) for i in np.where(np.abs(z - cl) > half)[0] if i >= start]
    return {"lambda": lam, "L": L, "z": z.tolist(), "ucl": (cl + half).tolist(), "lcl": (cl - half).tolist(),
            "first_signal": out[0] if out else None, "signals": out}


def individuals_chart(site: str, series: str, unit: str | None, t, x, labels, phase1: int = 6) -> Chart:
    x = np.asarray(x, float)
    m = min(max(phase1, 3), x.size)
    base = x[:m]
    cl = float(base.mean())
    mr = np.abs(np.diff(base))
    sigma = float(mr.mean() / D2) if mr.size else float("nan")
    if not np.isfinite(sigma) or sigma <= 0:
        sigma = float(np.std(base, ddof=1)) if base.size > 1 else 1e-9
    ch = Chart(site, series, unit, list(map(float, t)), x.tolist(), list(labels), cl, sigma, m, cl + 3 * sigma, cl - 3 * sigma)
    ch.violations = western_electric(x, cl, sigma, start=m)
    ch.ewma = ewma(x, cl, sigma, start=m)
    return ch


def chart_series(db: PlatformDB, site: str, artefact_id: str) -> list[tuple[str, str | None, list, list, list]]:
    """(series name, unit, t, x, labels) for everything worth charting at a site."""
    out = []
    for q, pair, f, name in [("insertion_loss", "A", 100e6, "artefact IL @ 100 MHz"), ("insertion_loss", "A", 600e6, "artefact IL @ 600 MHz"),
                             ("impedance_fitted", "A", None, "artefact fitted impedance"), ("loop_resistance", None, None, "artefact loop resistance")]:
        rows = [r for r in db.results_matrix(q, pair, f, artefact_id) if r["site"] == site]
        cells: dict[str, list] = {}
        for r in rows:
            cells.setdefault(r["batch"], []).append(r)
        if len(cells) >= 4:
            unit = db.q("SELECT unit FROM results WHERE quantity=? LIMIT 1", (q,))[0]["unit"]
            batches = sorted(cells, key=lambda b: min(r["t_unix"] for r in cells[b]))
            out.append((name, unit, [float(np.mean([r["t_unix"] for r in cells[b]])) for b in batches],
                        [float(np.mean([r["value"] for r in cells[b]])) for b in batches], batches))

    ver = db.q("SELECT t_unix, max_dev_db, min_rl_db, calibration_record FROM verifications WHERE site=? ORDER BY t_unix", (site,))
    if len(ver) >= 4:
        out.append(("verification |S21| deviation", "dB", [v["t_unix"] for v in ver], [v["max_dev_db"] for v in ver], [v["calibration_record"] for v in ver]))
        out.append(("verification return loss", "dB", [v["t_unix"] for v in ver], [v["min_rl_db"] for v in ver], [v["calibration_record"] for v in ver]))
    noise = db.q("SELECT t_unix, trace_noise_db, id FROM runs WHERE site=? AND sample_id=? AND trust IN ('trusted','flagged') AND trace_noise_db IS NOT NULL ORDER BY t_unix",
                 (site, artefact_id))
    if len(noise) >= 4:
        out.append(("trace noise (artefact runs)", "dB", [r["t_unix"] for r in noise], [r["trace_noise_db"] for r in noise], [r["id"] for r in noise]))
    return out


def run_all_charts(db: PlatformDB, artefact_id: str, phase1: int = 5, record_events: bool = True) -> list[Chart]:
    from datetime import datetime, timezone
    charts = []
    for site in db.sites():
        for name, unit, t, x, labels in chart_series(db, site, artefact_id):
            ch = individuals_chart(site, name, unit, t, x, labels, phase1 if "trace noise" not in name else max(phase1, 8))
            charts.append(ch)
            if record_events:
                for v in ch.violations:
                    if v["rule"] > 2:
                        continue
                    k = v["index"]
                    db.event(datetime.fromtimestamp(ch.t[k], tz=timezone.utc).isoformat(timespec="seconds"), "control-chart", site, name,
                             "action" if v["rule"] == 1 else "warning", f"{name}: {v['message']} at {labels[k]} (x = {ch.x[k]:.4g} {unit or ''})",
                             {"rule": v["rule"], "index": k, "value": ch.x[k], "cl": ch.cl, "sigma": ch.sigma})
                if ch.ewma.get("first_signal") is not None:
                    k = ch.ewma["first_signal"]
                    db.event(datetime.fromtimestamp(ch.t[k], tz=timezone.utc).isoformat(timespec="seconds"), "ewma", site, name, "warning",
                             f"{name}: EWMA signal at {labels[k]} (z = {ch.ewma['z'][k]:.4g}, CL {ch.cl:.4g})", {"index": k})
    db.commit()
    return charts
