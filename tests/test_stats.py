"""The statistical building blocks on synthetic data with known answers."""
import numpy as np
import pytest

from labplatform.control import ewma, individuals_chart
from labplatform.stats import (
    algorithm_a,
    chi2_consistency,
    cochran,
    grubbs,
    largest_consistent_subset,
    linear_trend,
    nested_anova,
    weighted_mean,
)


def test_weighted_mean_and_consistency():
    x, u = [1.0, 1.1, 0.9], [0.1, 0.1, 0.1]
    y, uy = weighted_mean(x, u)
    assert y == pytest.approx(1.0) and uy == pytest.approx(0.1 / np.sqrt(3))
    c = chi2_consistency(x, u)
    assert c["consistent"] and c["chi2"] == pytest.approx(2.0)
    c2 = chi2_consistency([1.0, 1.1, 2.0], [0.1, 0.1, 0.1])
    assert not c2["consistent"]


def test_largest_consistent_subset_drops_the_outlier():
    x = [10.00, 10.02, 9.99, 10.30, 10.01]
    u = [0.02, 0.02, 0.02, 0.02, 0.02]
    r = largest_consistent_subset(x, u, ["A", "B", "C", "D", "E"])
    assert r["excluded"] == ["D"] and set(r["members"]) == {"A", "B", "C", "E"}
    assert r["value"] == pytest.approx(10.005, abs=1e-3) and r["chi2"]["consistent"]
    r2 = largest_consistent_subset([1.0, 1.0, 1.0], [0.1, 0.1, 0.1])
    assert r2["excluded"] == [] and r2["u"] == pytest.approx(0.1 / np.sqrt(3))


def test_algorithm_a_is_robust():
    rng = np.random.default_rng(0)
    x = rng.normal(5.0, 0.1, 40)
    x[3] = 9.0                                  # one gross outlier
    r = algorithm_a(x)
    assert abs(r["mean"] - 5.0) < 0.05 and 0.06 < r["std"] < 0.16
    assert abs(np.mean(x) - 5.0) > 0.08          # the plain mean is pulled


def test_cochran_and_grubbs():
    c = cochran([0.01, 0.012, 0.011, 0.2, 0.01], n=3)
    assert c["outlier"] and c["index"] == 3
    # critical values reproduce ISO 5725-2 Table 4 (p = 5): 0.684 for n = 3, 0.544 for n = 5
    assert c["critical"] == pytest.approx(0.684, abs=0.002)
    assert cochran([1, 1, 1, 1, 1], n=5)["critical"] == pytest.approx(0.544, abs=0.002)
    c2 = cochran([0.01, 0.012, 0.011, 0.013, 0.01], n=3)
    assert not c2["outlier"]
    g = grubbs([10.0, 10.1, 9.9, 10.05, 12.0])
    assert g["outlier"] and g["index"] == 4
    assert not grubbs([10.0, 10.1, 9.9, 10.05, 10.02])["outlier"]


def test_nested_anova_recovers_variance_components():
    rng = np.random.default_rng(1)
    s_site, s_cond, s_err = 0.05, 0.02, 0.005
    cells = {}
    for i in range(6):
        b_site = rng.normal(0, s_site)
        for j in range(8):
            b_cond = rng.normal(0, s_cond)
            cells[(f"S{i}", f"R{j}")] = list(10.0 + b_site + b_cond + rng.normal(0, s_err, 4))
    a = nested_anova(cells)
    assert a.s_r == pytest.approx(s_err, rel=0.25)
    assert np.sqrt(a.var_cond) == pytest.approx(s_cond, rel=0.4)
    assert a.s_l == pytest.approx(s_site, rel=0.6)
    assert a.s_R > a.s_i > a.s_r
    d = a.to_dict()
    assert d["r_limit"] == pytest.approx(2.8 * a.s_r) and d["sites"] == 6


def test_linear_trend():
    t = np.arange(10.0)
    y = 1.0 + 0.01 * t + np.array([0.001, -0.001] * 5)
    tr = linear_trend(t, y)
    assert tr["slope"] == pytest.approx(0.01, abs=5e-4) and tr["p"] < 1e-6
    flat = linear_trend(t, 1.0 + np.array([0.001, -0.001] * 5))
    assert flat["p"] > 0.05


def test_control_chart_rules_and_ewma():
    rng = np.random.default_rng(2)
    x = 1.0 + rng.normal(0, 0.01, 30)
    x[20:] += 0.05                                       # a shift of 5 sigma after point 20
    ch = individuals_chart("S", "series", "dB", np.arange(30) * 86400.0, x, [str(i) for i in range(30)], phase1=10)
    assert abs(ch.cl - 1.0) < 0.01 and 0.005 < ch.sigma < 0.02
    rules = {v["rule"] for v in ch.violations}
    assert 1 in rules and min(v["index"] for v in ch.violations) >= 20
    assert ch.ewma["first_signal"] is not None and 20 <= ch.ewma["first_signal"] <= 22
    # a small sustained shift (1.5 sigma) is caught by EWMA and rule 4 but usually not rule 1
    y = 1.0 + rng.normal(0, 0.01, 40)
    y[20:] += 0.015
    ch2 = individuals_chart("S", "s", "dB", np.arange(40) * 86400.0, y, [str(i) for i in range(40)], phase1=10)
    assert ch2.ewma["first_signal"] is not None
    assert any(v["rule"] in (3, 4) for v in ch2.violations)
    # nothing in phase I is ever judged
    assert all(v["index"] >= 10 for v in ch2.violations)
    e = ewma(np.full(10, 1.0), 1.0, 0.01, start=0)
    assert e["first_signal"] is None


def test_change_point_detection():
    from labplatform.drift import change_point
    rng = np.random.default_rng(3)
    x = 1.0 + rng.normal(0, 0.005, 12)
    x[7:] += 0.04
    cp = change_point(x)
    assert cp["significant"] and cp["index"] == 7 and cp["shift"] == pytest.approx(0.04, abs=0.01)
    flat = change_point(1.0 + rng.normal(0, 0.005, 12))
    assert not flat["significant"]
    assert change_point([1.0, 1.1, 1.2])["index"] is None
