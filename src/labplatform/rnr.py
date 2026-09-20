"""
Repeatability and reproducibility of the network (ISO 5725-2 / -6 style).

The round-robin artefact gives a nested design: sites (p) x rounds within sites (each round on a fresh
calibration, a different day) x repeats within a round (re-connections).  The nested ANOVA of
:mod:`stats` yields the variance components

    s_r   repeatability            (re-connection, trace noise)          within a round
    s_I   intermediate precision   (+ calibration-to-calibration, day)   within a site
    s_R   reproducibility          (+ site-to-site)                      across the network

and the limits  r = 2.8 s_r,  R = 2.8 s_R  (the difference between two single results that is
exceeded with 5 % probability under repeatability / reproducibility conditions).  Before the ANOVA,
Cochran's test looks for a site whose within-round variance is out of line, and Grubbs' test for a
site whose mean is; both are reported, neither removes data automatically.

For production samples (each measured once per site, at most) only s_r from the artefact applies; the
platform reports how many production margins are within 2 s_R of their limit, i.e. verdicts that a
different site might have reversed.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .schema import PlatformDB
from .stats import NestedAnova, cochran, grubbs, nested_anova

__all__ = ["RnR", "repeatability_reproducibility", "verdicts_at_risk"]


@dataclass
class RnR:
    quantity: str
    pair: str | None
    f_hz: float | None
    unit: str | None
    anova: NestedAnova
    cochran: dict
    grubbs: dict
    cells: dict = field(default_factory=dict)
    per_site: dict = field(default_factory=dict)

    @property
    def label(self) -> str:
        return f"{self.quantity}[{self.pair}]" + (f" @ {self.f_hz / 1e6:.0f} MHz" if self.f_hz else "") if self.pair else self.quantity

    def to_dict(self) -> dict:
        return {"quantity": self.quantity, "pair": self.pair, "f_hz": self.f_hz, "unit": self.unit, "label": self.label,
                "anova": self.anova.to_dict(), "cochran": self.cochran, "grubbs": self.grubbs, "per_site": self.per_site}


def repeatability_reproducibility(db: PlatformDB, artefact_id: str, quantity: str, pair: str | None = "A", f_hz: float | None = None) -> RnR | None:
    rows = db.results_matrix(quantity, pair, f_hz, artefact_id)
    if not rows:
        return None
    cells: dict[tuple[str, str], list[float]] = {}
    for r in rows:
        cells.setdefault((r["site"], r["batch"]), []).append(float(r["value"]))
    if len({s for s, _ in cells}) < 2:
        return None
    an = nested_anova(cells)
    # Cochran on per-site pooled within-round variances (n = mean repeats), Grubbs on site means
    sites = sorted({s for s, _ in cells})
    per_site = {}
    var_site, n_rep = [], []
    for s in sites:
        vs = [np.asarray(v, float) for (ss, _), v in cells.items() if ss == s and len(v) >= 2]
        ss_ = sum(float(np.sum((v - v.mean()) ** 2)) for v in vs)
        dof = sum(v.size - 1 for v in vs)
        s_r = float(np.sqrt(ss_ / dof)) if dof > 0 else float("nan")
        means = [float(np.mean(v)) for (ss, _), v in cells.items() if ss == s]
        per_site[s] = {"s_r": s_r, "rounds": len(means), "s_between_rounds": float(np.std(means, ddof=1)) if len(means) >= 2 else float("nan"),
                       "mean": float(np.mean(means))}
        var_site.append(s_r ** 2 if np.isfinite(s_r) else 0.0)
        n_rep.append(int(round(np.mean([v.size for v in vs]))) if vs else 1)
    coch = cochran(var_site, max(min(n_rep), 2)) if len(sites) >= 2 else {}
    if coch.get("index") is not None:
        coch["site"] = sites[coch["index"]]
    gr = grubbs([per_site[s]["mean"] for s in sites])
    if gr.get("index") is not None:
        gr["site"] = sites[gr["index"]]
    unit_row = db.q("SELECT unit FROM results WHERE quantity=? LIMIT 1", (quantity,))
    return RnR(quantity, pair, f_hz, unit_row[0]["unit"] if unit_row else None, an, coch, gr,
               {f"{s}|{b}": v for (s, b), v in cells.items()}, per_site)


def verdicts_at_risk(db: PlatformDB, s_R_db: float, quantity: str = "insertion_loss") -> dict:
    """Production runs whose headline insertion-loss margin lies within 2 s_R of the limit."""
    rows = db.q("SELECT r.id, r.site, r.sample_id, r.verdict, r.headline, r.headline_margin FROM runs r JOIN samples s ON s.id = r.sample_id "
                "WHERE s.kind = 'production' AND r.trust IN ('trusted','flagged') AND r.headline LIKE ?", (quantity + "%",))
    at_risk = [dict(r) for r in rows if r["headline_margin"] is not None and abs(r["headline_margin"]) < 2 * s_R_db]
    return {"n_production": len(rows), "n_at_risk": len(at_risk), "band_db": 2 * s_R_db, "runs": at_risk}
