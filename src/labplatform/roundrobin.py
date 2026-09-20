"""
Interlaboratory comparison of the round-robin artefact.

For a quantity measured by every site on the same artefact:

    site value       x_i = mean of the site's temperature-corrected results
    site uncertainty u_i^2 = u_sys,i^2 + s_i^2 / n_i
                     u_sys,i : the systematic part of the budget (calibration + temperature), averaged
                     s_i     : observed standard deviation of the site's values across repeats and
                               rounds (covers connector repeatability and calibration-to-calibration
                               variation), or the pooled repeatability when n_i < 3
    reference value  either a designated reference site (x_ref, u_ref) or the consensus value of the
                     largest consistent subset (Cox 2002): the inverse-variance weighted mean of the
                     sites that pass the chi-squared consistency test
    degree of equivalence  d_i = x_i - x_ref,   u(d_i)^2 = u_i^2 - u_ref^2  (member of the consensus)
                                                u(d_i)^2 = u_i^2 + u_ref^2  (otherwise / designated ref)
    scores           E_n = d_i / (2 u(d_i))              ISO 13528 / ISO 17043: |E_n| <= 1 satisfactory
                     zeta = d_i / sqrt(u_i^2 + u_ref^2)  (same, without the k = 2)
                     z    = (x_i - x*) / s*               robust mean and SD of the site values (Algorithm A)

The same machinery run at 30 frequencies over the insertion-loss traces gives E_n(f) per site.
Artefact stability is the linear trend of the reference (or consensus) value over the rounds.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field

import numpy as np

from .schema import PlatformDB
from .stats import algorithm_a, largest_consistent_subset, linear_trend, weighted_mean
from .uncertainty import T_REF_C, TEMPERATURE_COEFFICIENTS

__all__ = ["SiteValue", "Comparison", "compare", "compare_all", "compare_by_round", "en_vs_frequency", "artefact_stability", "DEFAULT_QUANTITIES"]

DEFAULT_QUANTITIES = [("insertion_loss", "A", 10e6), ("insertion_loss", "A", 100e6), ("insertion_loss", "A", 300e6), ("insertion_loss", "A", 600e6),
                      ("impedance_fitted", "A", None), ("impedance_mean", "A", None), ("return_loss_min", "A", None),
                      ("nvp", "A", None), ("loop_resistance", None, None)]


@dataclass
class SiteValue:
    site: str
    n: int
    x: float
    s: float
    u_sys: float
    u: float
    d: float = float("nan")
    u_d: float = float("nan")
    en: float = float("nan")
    zeta: float = float("nan")
    z: float = float("nan")
    member: bool = False

    @property
    def status(self) -> str:
        a = abs(self.en)
        return "compatible" if a <= 1 else ("questionable" if a <= 1.5 else "incompatible")

    def statement(self, unit: str | None, x_ref: float) -> str:
        """The sentence a metrologist would write: the difference, the uncertainty of the difference, the verdict."""
        u = unit or ""
        U_d = 2 * self.u_d if np.isfinite(self.u_d) else float("nan")
        if not np.isfinite(U_d):
            return f"{self.site} is the reference ({self.x:.4g} {u}, U = {2 * self.u:.2g} {u})."
        if abs(self.en) <= 1:
            return (f"{self.site} differs from the reference by {self.d:+.4g} {u}; the expanded uncertainty of that difference is "
                    f"{U_d:.2g} {u} (E_n = {self.en:+.2f}): the difference is not statistically distinguishable from zero - compatible.")
        if abs(self.en) <= 1.5:
            return (f"{self.site} differs from the reference by {self.d:+.4g} {u} against U(d) = {U_d:.2g} {u} (E_n = {self.en:+.2f}): "
                    f"questionable - the difference exceeds its uncertainty but not by a factor 1.5.")
        return (f"{self.site} differs from the reference by {self.d:+.4g} {u} against U(d) = {U_d:.2g} {u} (E_n = {self.en:+.2f}): "
                f"the difference is statistically significant - incompatible; the site's uncertainty budget or its result is wrong.")

    def to_dict(self) -> dict:
        return {"site": self.site, "n": self.n, "x": self.x, "s": self.s, "u_sys": self.u_sys, "u": self.u, "U_k2": 2 * self.u,
                "d": self.d, "u_d": self.u_d, "U_d_k2": 2 * self.u_d if np.isfinite(self.u_d) else None, "En": self.en, "zeta": self.zeta, "z": self.z,
                "member": self.member, "status": self.status,
                "verdict": "satisfactory" if abs(self.en) <= 1 else ("questionable" if abs(self.en) <= 1.5 else "unsatisfactory")}


@dataclass
class Comparison:
    quantity: str
    pair: str | None
    f_hz: float | None
    unit: str | None
    reference: str
    x_ref: float
    u_ref: float
    sites: list[SiteValue]
    consensus: dict = field(default_factory=dict)
    robust: dict = field(default_factory=dict)

    @property
    def label(self) -> str:
        return f"{self.quantity}[{self.pair}]" + (f" @ {self.f_hz / 1e6:.0f} MHz" if self.f_hz else "") if self.pair else self.quantity

    def to_dict(self) -> dict:
        return {"quantity": self.quantity, "pair": self.pair, "f_hz": self.f_hz, "unit": self.unit, "label": self.label, "reference": self.reference,
                "x_ref": self.x_ref, "u_ref": self.u_ref, "U_ref_k2": 2 * self.u_ref,
                "sites": [s.to_dict() | {"statement": s.statement(self.unit, self.x_ref)} for s in self.sites],
                "consensus": {k: v for k, v in self.consensus.items() if k != "chi2"} | ({"chi2": self.consensus["chi2"]} if self.consensus.get("chi2") else {}),
                "robust": self.robust}


def _site_values(db: PlatformDB, artefact_id: str, quantity: str, pair: str | None, f_hz: float | None,
                 t_min: float | None = None, t_max: float | None = None) -> tuple[list[SiteValue], str | None]:
    rows = db.results_matrix(quantity, pair, f_hz, artefact_id)
    if t_min is not None:
        rows = [r for r in rows if r["t_unix"] >= t_min]
    if t_max is not None:
        rows = [r for r in rows if r["t_unix"] <= t_max]
    by: dict[str, list] = {}
    for r in rows:
        by.setdefault(r["site"], []).append(r)
    unit = None
    u_rows = db.q("SELECT s.unit FROM results s WHERE s.quantity=? LIMIT 1", (quantity,))
    if u_rows:
        unit = u_rows[0]["unit"]
    out = []
    for site, rs in sorted(by.items()):
        x = np.array([r["value"] for r in rs], float)
        budgets = []
        for r in rs:
            b = db.q("SELECT u_json FROM results WHERE run_id=? AND quantity=? AND (pair IS ? OR pair=?) AND (f_hz IS ? OR ABS(f_hz-?)<1)",
                     (r["run_id"], quantity, pair, pair, f_hz, f_hz if f_hz is not None else 0))
            if b and b[0]["u_json"]:
                budgets.append(json.loads(b[0]["u_json"]))
        if budgets:
            u_sys = float(np.sqrt(np.mean([b.get("u_sys", np.sqrt(b["u_cal"] ** 2 + b["u_T"] ** 2 + b["u_noise"] ** 2)) ** 2 for b in budgets])))
            u_rep_pooled = float(np.mean([b["u_rep"] for b in budgets]))
        else:
            u_sys, u_rep_pooled = 0.0, 0.01
        n = x.size
        s = float(np.std(x, ddof=1)) if n >= 3 else u_rep_pooled
        u = float(np.sqrt(u_sys ** 2 + s ** 2 / n))
        out.append(SiteValue(site, n, float(x.mean()), s, u_sys, u))
    return out, unit


def compare(db: PlatformDB, artefact_id: str, quantity: str, pair: str | None = "A", f_hz: float | None = None,
            reference: str = "consensus", t_min: float | None = None, t_max: float | None = None) -> Comparison | None:
    sites, unit = _site_values(db, artefact_id, quantity, pair, f_hz, t_min, t_max)
    if len(sites) < 2:
        return None
    x = np.array([s.x for s in sites])
    u = np.array([s.u for s in sites])
    labels = [s.site for s in sites]
    consensus = largest_consistent_subset(x, u, labels)
    robust = algorithm_a(x)
    if reference == "consensus":
        x_ref, u_ref = consensus["value"], consensus["u"]
        members = set(consensus["members"])
    else:
        ref = next((s for s in sites if s.site == reference), None)
        if ref is None:
            raise KeyError(f"reference site {reference!r} has no data for {quantity}")
        x_ref, u_ref, members = ref.x, ref.u, set()
    for s in sites:
        s.d = s.x - x_ref
        if reference == "consensus" and s.site in members:
            s.member = True
            s.u_d = float(np.sqrt(max(s.u ** 2 - u_ref ** 2, 1e-30)))
        elif reference != "consensus" and s.site == reference:
            s.member, s.u_d = True, float("nan")
        else:
            s.u_d = float(np.sqrt(s.u ** 2 + u_ref ** 2))
        s.en = s.d / (2 * s.u_d) if np.isfinite(s.u_d) and s.u_d > 0 else 0.0
        s.zeta = s.d / np.sqrt(s.u ** 2 + u_ref ** 2) if reference != "consensus" or not s.member else s.d / s.u_d
        s.z = (s.x - robust["mean"]) / robust["std"] if robust["std"] > 0 else 0.0
    return Comparison(quantity, pair, f_hz, unit, reference, float(x_ref), float(u_ref), sites, consensus, robust)


def compare_by_round(db: PlatformDB, artefact_id: str, quantity: str = "insertion_loss", pair: str | None = "A", f_hz: float | None = 600e6,
                     reference: str = "consensus") -> dict:
    """The comparison repeated per round (batch RR-xx): a time series of d_i, U_i and E_n per site.
    Within a round a site's n is the number of re-connections, so s_i/sqrt(n) covers only repeatability and
    the calibration of that round enters through u_cal - which is exactly what makes a degrading site visible
    round by round while the all-rounds comparison averages it away."""
    batches = [b["batch"] for b in db.q("SELECT DISTINCT batch FROM runs WHERE sample_id=? AND batch LIKE 'RR-%' ORDER BY batch", (artefact_id,))]
    rounds = []
    for b in batches:
        rows = db.q("SELECT MIN(t_unix) AS t0, MAX(t_unix) AS t1 FROM runs WHERE sample_id=? AND batch=?", (artefact_id, b))
        c = compare(db, artefact_id, quantity, pair, f_hz, reference, t_min=rows[0]["t0"] - 1, t_max=rows[0]["t1"] + 1)
        if c is None:
            continue
        rounds.append({"batch": b, "t_unix": float(rows[0]["t0"]), "x_ref": c.x_ref, "u_ref": c.u_ref, "members": c.consensus.get("members", []),
                       "sites": {s.site: {"x": s.x, "u": s.u, "d": s.d, "u_d": s.u_d, "En": s.en, "z": s.z, "n": s.n} for s in c.sites}})
    sites = sorted({s for r in rounds for s in r["sites"]})
    return {"quantity": quantity, "pair": pair, "f_hz": f_hz, "reference": reference, "sites": sites, "rounds": rounds}


def compare_all(db: PlatformDB, artefact_id: str, quantities=None, reference: str = "consensus", **kw) -> list[Comparison]:
    out = []
    for q, pair, f in (quantities or DEFAULT_QUANTITIES):
        c = compare(db, artefact_id, q, pair, f, reference, **kw)
        if c is not None:
            out.append(c)
    return out


def en_vs_frequency(db: PlatformDB, artefact_id: str, n_freq: int = 30, reference: str = "consensus", pair: str = "A") -> dict:
    """E_n(f) per site from the insertion-loss traces (temperature-corrected, budget as for the scalars)."""
    rows = db.q("SELECT r.site, r.id AS run_id, r.temperature_c, r.verification_dev_db, r.trace_noise_db, t.x_json, t.value_json "
                "FROM traces t JOIN runs r ON r.id = t.run_id WHERE r.sample_id=? AND t.quantity='insertion_loss' AND t.pair=? AND r.trust IN ('trusted','flagged')",
                (artefact_id, pair))
    if not rows:
        return {}
    x0 = np.asarray(json.loads(rows[0]["x_json"]), float)
    freqs = np.linspace(max(x0[0], 10e6), x0[-1], n_freq)
    c = TEMPERATURE_COEFFICIENTS["insertion_loss"]
    by: dict[str, dict] = {}
    for r in rows:
        x = np.asarray(json.loads(r["x_json"]), float)
        v = np.asarray(json.loads(r["value_json"]), float)
        T = r["temperature_c"] if r["temperature_c"] is not None else T_REF_C
        il = np.interp(freqs, x, v) / (1 + c * (T - T_REF_C))
        d = by.setdefault(r["site"], {"il": [], "u_cal": [], "u_T": []})
        d["il"].append(il)
        d["u_cal"].append((r["verification_dev_db"] or 0.1) / np.sqrt(3))
        d["u_T"].append(np.abs(il) * np.sqrt((c * 0.5) ** 2 + (0.3 * c * (T - T_REF_C)) ** 2))
    sites = sorted(by)
    X = np.array([np.mean(by[s]["il"], axis=0) for s in sites])                       # (sites, freqs)
    S = np.array([np.std(by[s]["il"], axis=0, ddof=1) if len(by[s]["il"]) >= 3 else np.full(freqs.size, 0.01) for s in sites])
    N = np.array([len(by[s]["il"]) for s in sites])[:, None]
    Usys = np.array([np.sqrt(np.mean(by[s]["u_cal"]) ** 2 + np.mean(by[s]["u_T"], axis=0) ** 2) for s in sites])
    U = np.sqrt(Usys ** 2 + S ** 2 / N)
    en = np.zeros_like(X)
    xref = np.zeros(freqs.size)
    uref = np.zeros(freqs.size)
    members = []
    for k in range(freqs.size):
        if reference == "consensus":
            lcs = largest_consistent_subset(X[:, k], U[:, k], sites)
            xref[k], uref[k] = lcs["value"], lcs["u"]
            mem = set(lcs["members"])
        else:
            i = sites.index(reference)
            xref[k], uref[k], mem = X[i, k], U[i, k], set()
        members.append(sorted(mem))
        for i, s in enumerate(sites):
            d = X[i, k] - xref[k]
            ud = np.sqrt(max(U[i, k] ** 2 - uref[k] ** 2, 1e-30)) if s in mem else np.sqrt(U[i, k] ** 2 + uref[k] ** 2)
            en[i, k] = d / (2 * ud) if not (reference != "consensus" and s == reference) else 0.0
    return {"f_hz": freqs.tolist(), "sites": sites, "il": X.tolist(), "u": U.tolist(), "x_ref": xref.tolist(), "u_ref": uref.tolist(),
            "En": en.tolist(), "members": members, "reference": reference}


def artefact_stability(db: PlatformDB, artefact_id: str, quantity: str = "insertion_loss", pair: str | None = "A", f_hz: float | None = 600e6,
                       site: str | None = None) -> dict:
    """Linear trend of the artefact's value over time: per site, and of the per-batch consensus."""
    rows = db.results_matrix(quantity, pair, f_hz, artefact_id)
    out = {"quantity": quantity, "pair": pair, "f_hz": f_hz, "per_site": {}, "consensus": None}
    by: dict[str, dict] = {}
    for r in rows:
        by.setdefault(r["site"], {}).setdefault(r["batch"], []).append((r["t_unix"], r["value"]))
    for s, cells in sorted(by.items()):
        # one point per round (repeats within a round share a calibration and are not independent)
        pts = [(float(np.mean([p[0] for p in v])), float(np.mean([p[1] for p in v]))) for v in cells.values()]
        t = np.array([p[0] for p in pts]) / 86400.0
        y = np.array([p[1] for p in pts])
        tr = linear_trend(t, y)
        tr["slope_per_day"] = tr.pop("slope")
        tr["slope_per_30d"] = tr["slope_per_day"] * 30 if np.isfinite(tr["slope_per_day"]) else float("nan")
        tr["significant"] = bool(np.isfinite(tr["p"]) and tr["p"] < 0.05)
        out["per_site"][s] = tr
    # consensus per batch (round)
    batches = db.q("SELECT DISTINCT batch FROM runs WHERE sample_id=? AND batch LIKE 'RR-%' ORDER BY batch", (artefact_id,))
    pts = []
    for b in batches:
        rs = db.q("SELECT r.t_unix, s.value, s.u_std FROM results s JOIN runs r ON r.id=s.run_id WHERE r.sample_id=? AND r.batch=? AND s.quantity=? "
                  "AND (s.pair IS ? OR s.pair=?) AND (s.f_hz IS ? OR ABS(s.f_hz-?)<1) AND r.trust IN ('trusted','flagged')",
                  (artefact_id, b["batch"], quantity, pair, pair, f_hz, f_hz if f_hz is not None else 0))
        if len(rs) >= 2:
            y, uy = weighted_mean([r["value"] for r in rs], [max(r["u_std"] or 0.01, 1e-6) for r in rs])
            pts.append((float(np.mean([r["t_unix"] for r in rs])) / 86400.0, y, uy, b["batch"]))
    if len(pts) >= 3:
        tr = linear_trend([p[0] for p in pts], [p[1] for p in pts])
        tr["slope_per_day"] = tr.pop("slope")
        tr["slope_per_30d"] = tr["slope_per_day"] * 30
        tr["significant"] = bool(tr["p"] < 0.05)
        tr["points"] = [{"t_days": p[0], "value": p[1], "u": p[2], "batch": p[3]} for p in pts]
        out["consensus"] = tr
    return out
