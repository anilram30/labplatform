"""
Uncertainty budgets per run and quantity, from the ingredients the sidecars carry.

Every scalar result carries a GUM-style budget [JCGM 100]

    u_c^2 = u_repeat^2 + u_reprod^2 + u_cal^2 + u_fixture^2 + u_instrument^2 + u_environment^2 + u_noise^2 + u_res^2

with the components (all in the quantity's unit):

    u_repeat      Type A repeatability of the *site* for this quantity: the pooled within-round standard
                  deviation of the round-robin artefact re-connections; a single-run value carries s, not
                  s/sqrt(n).  Before any repeats exist a declared default is used and flagged.
    u_reprod      the site's calibration-to-calibration (round-to-round) variation that the calibration
                  bound does not already cover: with sigma_c the site's between-round standard deviation of
                  the artefact round means (repeatability part removed),  u_reprod = sqrt(max(sigma_c^2 - u_cal^2, 0)).
                  Counting the full sigma_c as well as u_cal would count the same effect twice.
    u_cal         systematic instrument error bounded by the last calibration verification: the maximum
                  |S21| deviation dA from the check standard's certificate is the half-width of a rectangular
                  distribution, u_cal = dA / sqrt(3) (transmission quantities in dB); for reflection- and
                  leakage-type levels L (return loss, LCL, LCTL, NEXT) the residual floor is bounded by the
                  return loss RL seen on the check standard, u_Gamma = 10^(-RL/20)/sqrt(3), and
                  u_L = 8.686 u_Gamma / 10^(-L/20)  - a -27 dB return loss against a -42 dB floor is +-1 dB,
                  which is the honest state of such a comparison; for impedance quantities u_Z = 2 Z_ref u_Gamma;
                  for NVP the artefact's physical length enters at a declared +-0.3 %.
    u_fixture     declared per fixture method of the procedure: none 0; port extension 0.01 dB; 2x-thru
                  de-embedding 0.02 dB (residual of the bisection); impedance 0.2 ohm when a fixture is removed.
    u_instrument  declared receiver linearity / dynamic accuracy: 0.01 dB, 0.05 ohm, 0.02 % for NVP.
    u_environment temperature: the value is corrected to T_ref = 23 C with a stated relative coefficient c
                  per quantity kind, x_23 = x / (1 + c (T - 23)); the correction's uncertainty combines
                  u(T) = 0.5 K (label uncertainty) and a 30 % relative uncertainty of c,
                  u_T = |x| sqrt( (c u(T))^2 + (0.3 c (T - 23))^2 ); plus the ambient window of the
                  procedure as a rectangular term on quantities with a coefficient.
    u_noise       the trace-noise estimate of the validation check (dB); negligible for a healthy sweep
    u_res         resolution: half the last stored digit (negligible)

Not evaluated (declared): mismatch of the test-port match with the cable's impedance, cable movement
between calibration and measurement beyond what the verification saw.  Expanded uncertainty
U = k u_c with k = 2 (~95 %).  The budget is stored with every result, so a comparison can say not only
"site A differs from the reference by 0.8 dB" but "by 0.8 dB against a combined U(d) of 1.2 dB: not
distinguishable" - the compatibility statement of :mod:`roundrobin`.
"""
from __future__ import annotations

import json
from dataclasses import dataclass

import numpy as np

from .schema import PlatformDB

__all__ = ["Budget", "evaluate_budgets", "TEMPERATURE_COEFFICIENTS", "T_REF_C", "pooled_repeatability", "site_reproducibility"]

T_REF_C = 23.0
TEMPERATURE_COEFFICIENTS = {
    "insertion_loss": 0.0025,      # ~0.5 alpha_rho on the conductor term + beta_d on the dielectric term, blended
    "loop_resistance": 0.00393,    # copper
    "impedance_mean": -6e-5,       # C'(T) grows ~1.2e-4/K -> Z falls ~0.6e-4/K
    "impedance_fitted": -6e-5,
    "impedance_min": -6e-5,
    "impedance_max": -6e-5,
    "nvp": -6e-5,
}
DEFAULT_REPEATABILITY = {"dB": 0.01, "ohm": 0.3, "": 0.001, None: 0.01}
U_TEMPERATURE_K = 0.5
REL_U_COEFF = 0.3
FIXTURE_U_DB = {"none": 0.0, "port-extension": 0.01, "2xthru": 0.02, "files": 0.02}
FIXTURE_U_OHM = {"none": 0.0, "port-extension": 0.1, "2xthru": 0.2, "files": 0.2}
INSTRUMENT_U = {"dB": 0.01, "ohm": 0.05, "": 0.0002, None: 0.0}


@dataclass
class Budget:
    value: float                # temperature-corrected value
    value_raw: float
    u_repeat: float
    u_reprod: float
    u_cal: float
    u_fixture: float
    u_instrument: float
    u_environment: float
    u_noise: float
    u_res: float
    rep_source: str
    reprod_source: str

    @property
    def components(self) -> dict:
        return {"u_repeat": self.u_repeat, "u_reprod": self.u_reprod, "u_cal": self.u_cal, "u_fixture": self.u_fixture,
                "u_instrument": self.u_instrument, "u_environment": self.u_environment, "u_noise": self.u_noise, "u_res": self.u_res}

    @property
    def u_c(self) -> float:
        return float(np.sqrt(sum(v ** 2 for v in self.components.values())))

    @property
    def u_systematic(self) -> float:
        """The part that does not average down with repeats at one site."""
        return float(np.sqrt(self.u_cal ** 2 + self.u_fixture ** 2 + self.u_instrument ** 2 + self.u_environment ** 2 + self.u_noise ** 2))

    @property
    def U(self) -> float:
        return 2.0 * self.u_c

    def to_dict(self) -> dict:
        d = {"value": self.value, "value_raw": self.value_raw, **self.components, "u_T": self.u_environment,
             "u_rep": self.u_repeat, "u_sys": self.u_systematic, "u_c": self.u_c, "U_k2": self.U,
             "rep_source": self.rep_source, "reprod_source": self.reprod_source}
        return d


def _kind(quantity: str) -> str:
    return quantity.replace("_min", "") if quantity.endswith("_min") else quantity


def pooled_repeatability(db: PlatformDB, artefact_id: str) -> dict[tuple, float]:
    """Pooled within-(site, batch) standard deviation of the artefact repeats, per site and quantity."""
    rows = db.q("SELECT r.site, r.batch, s.quantity, s.pair, s.f_hz, s.value FROM results s JOIN runs r ON r.id = s.run_id "
                "WHERE r.sample_id = ? AND r.trust IN ('trusted','flagged')", (artefact_id,))
    groups: dict = {}
    for r in rows:
        groups.setdefault((r["site"], r["quantity"], r["pair"], r["f_hz"]), {}).setdefault(r["batch"], []).append(r["value"])
    out = {}
    for key, cells in groups.items():
        ss, dof = 0.0, 0
        for v in cells.values():
            v = np.asarray(v, float)
            if v.size >= 2:
                ss += float(np.sum((v - v.mean()) ** 2))
                dof += v.size - 1
        if dof >= 2:
            out[key] = float(np.sqrt(ss / dof))
    return out


def site_reproducibility(db: PlatformDB, artefact_id: str, rep: dict[tuple, float]) -> dict[tuple, float]:
    """Between-round standard deviation sigma_c of the artefact round means per site and quantity, with the
    repeatability contribution s_r^2 / n removed (negative -> 0).  Needs >= 3 rounds."""
    rows = db.q("SELECT r.site, r.batch, s.quantity, s.pair, s.f_hz, s.value FROM results s JOIN runs r ON r.id = s.run_id "
                "WHERE r.sample_id = ? AND r.trust IN ('trusted','flagged')", (artefact_id,))
    groups: dict = {}
    for r in rows:
        groups.setdefault((r["site"], r["quantity"], r["pair"], r["f_hz"]), {}).setdefault(r["batch"], []).append(r["value"])
    out = {}
    for key, cells in groups.items():
        means = [float(np.mean(v)) for v in cells.values()]
        n_mean = float(np.mean([len(v) for v in cells.values()]))
        if len(means) >= 3:
            s_between = float(np.var(means, ddof=1))
            s_r = rep.get(key, 0.0)
            out[key] = float(np.sqrt(max(s_between - s_r ** 2 / n_mean, 0.0)))
    return out


def evaluate_budgets(db: PlatformDB, artefact_id: str | None = None) -> int:
    """Fill results.u_std / u_json (and replace value by the temperature-corrected value, keeping the raw one) for every result."""
    rep = pooled_repeatability(db, artefact_id) if artefact_id else {}
    reprod = site_reproducibility(db, artefact_id, rep) if artefact_id else {}
    runs = {r["id"]: r for r in db.q("SELECT * FROM runs")}
    fixtures = {}
    for p in db.q("SELECT sha256, text FROM procedures"):
        m = "none"
        for line in (p["text"] or "").splitlines():
            if line.strip().startswith("method"):
                m = line.split("=")[1].strip().strip('"')
                break
        fixtures[p["sha256"]] = m
    n = 0
    for s in db.q("SELECT * FROM results"):
        r = runs.get(s["run_id"])
        if r is None:
            continue
        q, unit = s["quantity"], s["unit"]
        kind = _kind(q)
        raw = float(s["value"]) if s["u_json"] is None else json.loads(s["u_json"]).get("value_raw", float(s["value"]))
        # temperature correction and environment
        T = r["temperature_c"] if r["temperature_c"] is not None else T_REF_C
        c = TEMPERATURE_COEFFICIENTS.get(kind, 0.0)
        corr = 1.0 + c * (T - T_REF_C)
        value = raw / corr if corr > 0 else raw
        u_T = abs(value) * float(np.sqrt((c * U_TEMPERATURE_K) ** 2 + (REL_U_COEFF * c * (T - T_REF_C)) ** 2))
        u_env = u_T
        # calibration
        dA = r["verification_dev_db"]
        rl = None
        v = db.q("SELECT min_rl_db FROM verifications WHERE site=? AND calibration_record=? ORDER BY t_unix DESC LIMIT 1", (r["site"], r["calibration_record"]))
        if v:
            rl = v[0]["min_rl_db"]
        fx = fixtures.get(r["procedure_sha"], "none")
        if unit == "dB" and kind in ("return_loss", "lcl", "lctl", "next", "fext"):
            # a reflection- or leakage-type level L (dB) measured against a residual floor bounded by the check
            # standard's return loss: u_Gamma = 10^(-RL/20)/sqrt(3) adds to |Gamma| = 10^(-L/20), so
            # u_L = 8.686 u_Gamma / |Gamma|  (dominant for weak reflections: a -27 dB return loss against a -42 dB floor is +-1 dB)
            u_gamma = 10 ** (-(rl if rl is not None else 30.0) / 20) / np.sqrt(3)
            gamma = 10 ** (-abs(value) / 20)
            u_cal = float(8.686 * u_gamma / max(gamma, 1e-6))
            u_fix, u_inst = FIXTURE_U_DB.get(fx, 0.02), INSTRUMENT_U["dB"]
        elif unit == "dB":
            u_cal = (dA if dA is not None else 0.1) / np.sqrt(3)
            u_fix, u_inst = FIXTURE_U_DB.get(fx, 0.02), INSTRUMENT_U["dB"]
        elif unit == "ohm" and kind.startswith("impedance"):
            u_gamma = 10 ** (-(rl if rl is not None else 30.0) / 20) / np.sqrt(3)
            u_cal = 2 * 100.0 * u_gamma
            u_fix, u_inst = FIXTURE_U_OHM.get(fx, 0.2), INSTRUMENT_U["ohm"]
        elif unit == "ohm":
            u_cal, u_fix, u_inst = 0.002, 0.0, 0.001          # 4-wire ohmmeter: contact, no fixture, linearity
        elif kind in ("nvp", "delay_per_metre"):
            u_cal, u_fix, u_inst = 0.003 * abs(value), 0.0, INSTRUMENT_U[""] * abs(value)
        else:
            u_cal, u_fix, u_inst = 0.0, 0.0, 0.0
        u_noise = float(r["trace_noise_db"] or 0.0) if unit == "dB" else 0.0
        key = (r["site"], q, s["pair"], s["f_hz"])
        if key in rep:
            u_rep, src = rep[key], "pooled artefact repeats"
        else:
            u_rep, src = DEFAULT_REPEATABILITY.get(unit, 0.01), "declared default (no repeats yet)"
        if key in reprod:
            u_rp, src2 = float(np.sqrt(max(reprod[key] ** 2 - u_cal ** 2, 0.0))), "site round-to-round sigma_c beyond the calibration bound"
        else:
            u_rp, src2 = 0.0, "not available (< 3 rounds)"
        u_res = 0.5e-4 if unit == "dB" else 0.5e-3
        b = Budget(value, raw, float(u_rep), u_rp, float(u_cal), float(u_fix), float(u_inst), float(u_env), u_noise, u_res, src, src2)
        db.con.execute("UPDATE results SET value = ?, u_std = ?, u_json = ? WHERE id = ?", (value, b.u_c, json.dumps(b.to_dict()), s["id"]))
        n += 1
    db.commit()
    return n
