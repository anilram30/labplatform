"""
Statistical building blocks used by the comparison, precision and control modules.

    weighted_mean        inverse-variance weighted mean and its standard uncertainty
    chi2_consistency     Birge / chi-squared consistency test of a set of (x_i, u_i)
    largest_consistent_subset   Cox's procedure for a consensus value (Metrologia 39, 2002)
    algorithm_a          ISO 13528 Annex C robust mean and standard deviation
    cochran              Cochran's test for a within-cell variance outlier (ISO 5725-2 §7.3.3)
    grubbs               Grubbs' test for one outlying cell mean (ISO 5725-2 §7.3.4)
    nested_anova         balanced-approximation nested ANOVA: sites / conditions within sites / repeats
    linear_trend         slope, its standard error and p-value (artefact stability)
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import stats as sps

__all__ = ["weighted_mean", "chi2_consistency", "largest_consistent_subset", "algorithm_a", "cochran", "grubbs",
           "nested_anova", "NestedAnova", "linear_trend"]


def weighted_mean(x, u) -> tuple[float, float]:
    x, u = np.asarray(x, float), np.asarray(u, float)
    w = 1.0 / u ** 2
    y = float(np.sum(w * x) / np.sum(w))
    return y, float(1.0 / np.sqrt(np.sum(w)))


def chi2_consistency(x, u, y: float | None = None, alpha: float = 0.05) -> dict:
    """chi2_obs = sum((x_i - y)^2 / u_i^2) against chi2_{1-alpha}(N-1)."""
    x, u = np.asarray(x, float), np.asarray(u, float)
    if y is None:
        y, _ = weighted_mean(x, u)
    n = x.size
    chi2 = float(np.sum((x - y) ** 2 / u ** 2))
    crit = float(sps.chi2.ppf(1 - alpha, max(n - 1, 1)))
    return {"chi2": chi2, "dof": n - 1, "critical": crit, "consistent": bool(chi2 <= crit), "birge": float(np.sqrt(chi2 / max(n - 1, 1)))}


def largest_consistent_subset(x, u, labels=None, alpha: float = 0.05) -> dict:
    """Cox (2002) procedure: weighted mean of all; if inconsistent, drop the member with the largest |E_n|
    and repeat until the remaining set is consistent.  Returns the reference value, its uncertainty,
    the members and the excluded labels."""
    x, u = np.asarray(x, float), np.asarray(u, float)
    labels = list(labels) if labels is not None else list(range(x.size))
    idx = list(range(x.size))
    excluded = []
    while len(idx) >= 2:
        xs, us = x[idx], u[idx]
        y, uy = weighted_mean(xs, us)
        c = chi2_consistency(xs, us, y, alpha)
        if c["consistent"]:
            return {"value": y, "u": uy, "members": [labels[i] for i in idx], "excluded": excluded, "chi2": c}
        # degrees of equivalence for members: u(d_i)^2 = u_i^2 - u_y^2
        en = np.abs(xs - y) / (2 * np.sqrt(np.maximum(us ** 2 - uy ** 2, 1e-30)))
        k = int(np.argmax(en))
        excluded.append(labels[idx[k]])
        idx.pop(k)
    y, uy = weighted_mean(x[idx], u[idx]) if idx else (float("nan"), float("nan"))
    return {"value": y, "u": uy, "members": [labels[i] for i in idx], "excluded": excluded, "chi2": None}


def algorithm_a(x, tol: float = 1e-6, max_iter: int = 100) -> dict:
    """ISO 13528:2015 Annex C.3.1 robust mean x* and standard deviation s*."""
    x = np.asarray(x, float)
    n = x.size
    if n < 2:
        return {"mean": float(x.mean()) if n else float("nan"), "std": float("nan"), "iterations": 0}
    xs = float(np.median(x))
    ss = float(1.483 * np.median(np.abs(x - xs)))
    if ss == 0:
        ss = float(np.std(x, ddof=1)) or 1e-12
    for it in range(max_iter):
        delta = 1.5 * ss
        xi = np.clip(x, xs - delta, xs + delta)
        xs_new = float(np.mean(xi))
        ss_new = float(1.134 * np.sqrt(np.sum((xi - xs_new) ** 2) / (n - 1)))
        if abs(xs_new - xs) < tol * max(abs(xs), 1e-12) and abs(ss_new - ss) < tol * max(ss, 1e-12):
            xs, ss = xs_new, ss_new
            break
        xs, ss = xs_new, ss_new
    return {"mean": xs, "std": ss, "iterations": it + 1}


def cochran(variances, n: int, alpha: float = 0.05) -> dict:
    """C = s_max^2 / sum s_i^2; critical value from the F distribution:
    C_crit = 1 / (1 + (p - 1) / F_{alpha/p}(nu, (p - 1) nu)),  nu = n - 1."""
    s2 = np.asarray(variances, float)
    p = s2.size
    if p < 2 or n < 2 or s2.sum() <= 0:
        return {"C": float("nan"), "critical": float("nan"), "outlier": None, "index": None}
    nu = n - 1
    F = sps.f.ppf(1 - alpha / p, nu, (p - 1) * nu)
    crit = 1.0 / (1.0 + (p - 1) / F)
    k = int(np.argmax(s2))
    C = float(s2[k] / s2.sum())
    return {"C": C, "critical": float(crit), "outlier": bool(C > crit), "index": k}


def grubbs(values, alpha: float = 0.05) -> dict:
    """Two-sided Grubbs test for a single outlier among cell means."""
    x = np.asarray(values, float)
    n = x.size
    if n < 3:
        return {"G": float("nan"), "critical": float("nan"), "outlier": None, "index": None}
    m, s = x.mean(), x.std(ddof=1)
    if s == 0:
        return {"G": 0.0, "critical": float("nan"), "outlier": False, "index": None}
    k = int(np.argmax(np.abs(x - m)))
    G = float(abs(x[k] - m) / s)
    t = sps.t.ppf(1 - alpha / (2 * n), n - 2)
    crit = float((n - 1) / np.sqrt(n) * np.sqrt(t ** 2 / (n - 2 + t ** 2)))
    return {"G": G, "critical": crit, "outlier": bool(G > crit), "index": k}


@dataclass
class NestedAnova:
    a: int                          # sites
    b_mean: float                   # conditions (rounds) per site, mean
    n_mean: float                   # repeats per cell, mean
    grand_mean: float
    ms_site: float
    ms_cond: float
    ms_err: float
    var_site: float
    var_cond: float
    var_err: float
    site_means: dict = field(default_factory=dict)

    @property
    def s_r(self) -> float:         # repeatability
        return float(np.sqrt(self.var_err))

    @property
    def s_i(self) -> float:         # intermediate precision (condition-to-condition within a site)
        return float(np.sqrt(self.var_err + self.var_cond))

    @property
    def s_l(self) -> float:         # between-site
        return float(np.sqrt(self.var_site))

    @property
    def s_R(self) -> float:         # reproducibility
        return float(np.sqrt(self.var_err + self.var_cond + self.var_site))

    def to_dict(self) -> dict:
        return {"sites": self.a, "conditions_per_site": self.b_mean, "repeats_per_cell": self.n_mean, "grand_mean": self.grand_mean,
                "ms_site": self.ms_site, "ms_condition": self.ms_cond, "ms_error": self.ms_err,
                "var_site": self.var_site, "var_condition": self.var_cond, "var_error": self.var_err,
                "s_r": self.s_r, "s_I": self.s_i, "s_L": self.s_l, "s_R": self.s_R,
                "r_limit": 2.8 * self.s_r, "R_limit": 2.8 * self.s_R, "site_means": self.site_means}


def nested_anova(cells: dict[tuple[str, str], list[float]]) -> NestedAnova:
    """cells[(site, condition)] = repeat values.  Two-stage nested ANOVA with expected mean squares
    E[MS_err] = s_e^2,  E[MS_cond] = s_e^2 + n s_c^2,  E[MS_site] = s_e^2 + n s_c^2 + n b s_a^2
    (balanced formulae with mean n and b; negative components are set to zero)."""
    sites = sorted({s for s, _ in cells})
    cell_means, cell_n, ss_err, n_cells = {}, {}, 0.0, 0
    for k, v in cells.items():
        v = np.asarray(v, float)
        if v.size == 0:
            continue
        cell_means[k] = float(v.mean())
        cell_n[k] = v.size
        ss_err += float(np.sum((v - v.mean()) ** 2))
        n_cells += 1
    n_mean = float(np.mean(list(cell_n.values())))
    b_mean = float(n_cells / len(sites))
    site_means = {s: float(np.mean([m for (ss, _), m in cell_means.items() if ss == s])) for s in sites}
    grand = float(np.mean(list(cell_means.values())))
    ss_cond = sum(cell_n[k] * (m - site_means[k[0]]) ** 2 for k, m in cell_means.items())
    ss_site = sum(sum(cell_n[k] for k in cell_means if k[0] == s) * (site_means[s] - grand) ** 2 for s in sites)
    df_err = max(sum(nk - 1 for nk in cell_n.values()), 1)
    df_cond = max(n_cells - len(sites), 1)
    df_site = max(len(sites) - 1, 1)
    ms_err, ms_cond, ms_site = ss_err / df_err, ss_cond / df_cond, ss_site / df_site
    var_err = ms_err
    var_cond = max((ms_cond - ms_err) / n_mean, 0.0)
    var_site = max((ms_site - ms_cond) / (n_mean * b_mean), 0.0)
    return NestedAnova(len(sites), b_mean, n_mean, grand, ms_site, ms_cond, ms_err, var_site, var_cond, var_err, site_means)


def linear_trend(t, y) -> dict:
    """Ordinary least squares y = a + b t with the standard error of b and the two-sided p-value of b = 0."""
    t, y = np.asarray(t, float), np.asarray(y, float)
    n = t.size
    if n < 3:
        return {"slope": float("nan"), "se": float("nan"), "p": float("nan"), "intercept": float("nan"), "n": n}
    res = sps.linregress(t, y)
    return {"slope": float(res.slope), "se": float(res.stderr), "p": float(res.pvalue), "intercept": float(res.intercept), "n": n,
            "r2": float(res.rvalue ** 2)}
