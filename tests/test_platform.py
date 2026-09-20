"""The platform end to end on a small simulated network: ingest, budgets, comparison, R&R, charts, traceability, dashboard, CLI."""
import json
import logging
from pathlib import Path

import numpy as np
import pytest

from labplatform.analyse import analyse_network
from labplatform.cli import main
from labplatform.dashboard import build_dashboard
from labplatform.ingest import ingest_lab
from labplatform.rnr import repeatability_reproducibility
from labplatform.roundrobin import compare, compare_by_round, en_vs_frequency
from labplatform.schema import PlatformDB
from labplatform.simulate import ARTEFACT_ID, SITES, build_network
from labplatform.traceability import chain, evidence_pack, suspect_runs
from labplatform.uncertainty import evaluate_budgets

logging.getLogger("labauto").setLevel(logging.WARNING)

# two sites, five rounds, two re-connections: LABA is the reference-quality site, LABB degrades and finally fails verification
TEST_SITES = {
    "LABA": dict(SITES["LAB1"], name="A"),
    "LABB": dict(SITES["LAB4"], name="B", degrade_round=4, degrade_quality=1.8, fail_round=5, fail_quality=4.5),
}


@pytest.fixture(scope="module")
def network(tmp_path_factory):
    root = tmp_path_factory.mktemp("net")
    out = build_network(root, rounds=5, repeats=2, production_per_round=1, sites=TEST_SITES)
    db = PlatformDB(root / "platform.sqlite")
    for site, info in out["sites"].items():
        ingest_lab(db, info["root"])
    return root, out, db


def test_ingest_is_complete(network):
    root, out, db = network
    c = db.counts()
    assert c["sites"] == 2 and c["instruments"] == 2
    # LABA: 5 rounds x (2 artefact + 1 production) = 15 runs; LABB: round 5 refused (3 jobs) but recorded from the lab DB
    runs = {r["site"]: r["n"] for r in db.q("SELECT site, COUNT(*) AS n FROM runs GROUP BY site")}
    assert runs["LABA"] == 15 and runs["LABB"] == 15
    refused = db.q("SELECT COUNT(*) AS n FROM runs WHERE state='ABORTED'")[0]["n"]
    assert refused == 3
    assert db.q("SELECT COUNT(*) AS n FROM runs WHERE trust='trusted'")[0]["n"] == 27
    # quantities derived from the traces, with the raw archive linked by manifest hash
    q = {(r["quantity"], r["f_hz"]) for r in db.q("SELECT DISTINCT quantity, f_hz FROM results")}
    assert ("insertion_loss", 600e6) in q and ("impedance_fitted", None) in q and ("loop_resistance", None) in q
    r = db.q("SELECT archive_dir, manifest_sha256, verification_dev_db, cal_age_days FROM runs WHERE trust='trusted' LIMIT 1")[0]
    assert Path(r["archive_dir"]).joinpath("manifest.json").exists() and len(r["manifest_sha256"]) == 64
    assert 0 < r["verification_dev_db"] < 0.1 and 0 < r["cal_age_days"] < 3
    assert db.q("SELECT COUNT(*) AS n FROM verifications")[0]["n"] >= 10
    assert db.q("SELECT kind FROM samples WHERE id=?", (ARTEFACT_ID,))[0]["kind"] == "artefact"
    # idempotent
    ingest_lab(db, out["sites"]["LABA"]["root"])
    assert db.counts()["runs"] == c["runs"] and db.counts()["results"] == c["results"]


def test_uncertainty_budgets(network):
    root, out, db = network
    n = evaluate_budgets(db, ARTEFACT_ID)
    assert n == db.counts()["results"]
    rows = db.q("SELECT s.value, s.u_std, s.u_json, r.temperature_c FROM results s JOIN runs r ON r.id=s.run_id "
                "WHERE s.quantity='insertion_loss' AND ABS(s.f_hz-600e6)<1 AND r.site='LABA' AND r.sample_id=?", (ARTEFACT_ID,))
    assert rows
    b = json.loads(rows[0]["u_json"])
    comps = ("u_repeat", "u_reprod", "u_cal", "u_fixture", "u_instrument", "u_environment", "u_noise", "u_res")
    assert b["u_c"] == pytest.approx(np.sqrt(sum(b[k] ** 2 for k in comps)))
    assert b["u_fixture"] == 0.0 and b["u_instrument"] == 0.01 and b["reprod_source"]
    assert b["u_cal"] == pytest.approx(0.1 / np.sqrt(3), abs=0.06) and b["u_cal"] > 0.005
    assert b["rep_source"] == "pooled artefact repeats" and 0 < b["u_rep"] < 0.03
    assert rows[0]["temperature_c"] == 23.0 and b["value"] == pytest.approx(b["value_raw"])   # no correction at 23 C
    # a warm site is corrected downwards and carries a temperature term
    ohm = json.loads(db.q("SELECT u_json FROM results WHERE quantity='impedance_fitted' LIMIT 1")[0]["u_json"])
    assert ohm["u_cal"] > 0.3   # bounded by the residual directivity seen on the check standard


def test_comparison_scores_and_reference(network):
    root, out, db = network
    evaluate_budgets(db, ARTEFACT_ID)
    c = compare(db, ARTEFACT_ID, "insertion_loss", "A", 600e6, "consensus")
    assert {s.site for s in c.sites} == {"LABA", "LABB"}
    for s in c.sites:
        assert s.n >= 8 and s.u > 0 and np.isfinite(s.en)
        if s.member:
            assert s.u_d == pytest.approx(np.sqrt(s.u ** 2 - c.u_ref ** 2))
    d = c.to_dict()
    assert d["sites"][0]["verdict"] in ("satisfactory", "questionable", "unsatisfactory")
    # designated reference: the reference site itself has En = 0 and is a member
    c2 = compare(db, ARTEFACT_ID, "insertion_loss", "A", 600e6, "LABA")
    ref = next(s for s in c2.sites if s.site == "LABA")
    assert ref.member and ref.en == 0.0
    other = next(s for s in c2.sites if s.site == "LABB")
    assert other.u_d == pytest.approx(np.sqrt(other.u ** 2 + ref.u ** 2))
    with pytest.raises(KeyError):
        compare(db, ARTEFACT_ID, "insertion_loss", "A", 600e6, "NOPE")
    enf = en_vs_frequency(db, ARTEFACT_ID, n_freq=8)
    assert len(enf["f_hz"]) == 8 and np.array(enf["En"]).shape == (2, 8)
    br = compare_by_round(db, ARTEFACT_ID)
    # round 5 has one site only (LABB was refused) -> no comparison possible for that round
    assert [r["batch"] for r in br["rounds"]] == ["RR-01", "RR-02", "RR-03", "RR-04"]
    assert set(br["rounds"][0]["sites"]) == {"LABA", "LABB"}


def test_rnr(network):
    root, out, db = network
    evaluate_budgets(db, ARTEFACT_ID)
    r = repeatability_reproducibility(db, ARTEFACT_ID, "insertion_loss", "A", 600e6)
    a = r.anova
    assert 0.001 < a.s_r < 0.02 and a.s_R >= a.s_i >= a.s_r
    assert a.a == 2 and a.n_mean == pytest.approx(2.0)
    assert set(r.per_site) == {"LABA", "LABB"} and r.per_site["LABA"]["rounds"] == 5
    assert r.cochran["C"] == r.cochran["C"]            # computed (not NaN); Grubbs needs >= 3 sites and is NaN here
    assert r.grubbs["outlier"] is None


def test_analysis_charts_suspects_and_dashboard(network):
    root, out, db = network
    s = analyse_network(db, ARTEFACT_ID, phase1=3)
    assert s["counts"]["runs"] == 30
    ov = {o["site"]: o for o in s["overview"]}
    assert ov["LABB"]["aborted"] == 3 and ov["LABB"]["last_verification"]["status"] == "fail"
    assert ov["LABA"]["last_verification"]["status"] == "pass"
    # the degradation of LABB shows on its verification chart
    chart = next(c for c in s["charts"] if c["site"] == "LABB" and c["series"] == "verification |S21| deviation")
    assert chart["x"][-1] > 0.1 and any(v["rule"] == 1 for v in chart["violations"])
    # retrospective invalidation: LABB's round-4 runs sit between its last passed and its failed verification
    sus = s["suspect_runs"]
    assert sus and all(x["site"] == "LABB" for x in sus)
    assert {x["sample_id"] for x in sus} >= {ARTEFACT_ID}
    assert all("RR-04" in x["run_id"] or "P-04" in x["run_id"] or x["started"] >= "2026-04-13" for x in sus)
    assert len({x["run_id"] for x in sus}) == len(sus)      # no duplicates across repeated failed verifications
    assert any(e["kind"] == "suspect-run" for e in s["events"])
    # traceability chain and evidence pack
    run_id = sus[0]["run_id"]
    ch = chain(db, run_id)
    assert ch["archive_integrity"]["ok"] and ch["calibration"]["kit_serial"].startswith("85052D")
    assert any(f["path"].endswith(".s4p") for f in ch["files"]) and ch["procedure"]["sha256"]
    assert ch["results"][0]["budget"]["u_c"] > 0
    p = evidence_pack(db, run_id, root / "evidence", s["comparisons"])
    assert p.exists() and json.loads(p.read_text())["run"]["id"] == run_id
    # dashboard is self-contained
    html = build_dashboard(db, s, root / "dashboard.html")
    text = html.read_text()
    assert "data:image/png;base64" in text and "http" not in text.split("<body>")[1].replace("http-equiv", "")
    assert "LABB" in text and "Suspect runs" in text


def test_cli(network, capsys):
    root, out, db = network
    dbp = str(root / "platform.sqlite")
    assert main(["analyse", "--db", dbp, "--artefact", ARTEFACT_ID, "--out", str(root / "s.json")]) == 0
    text = capsys.readouterr().out
    assert "E_n scores" in text and "repeatability" in text and (root / "s.json").exists()
    assert main(["compare", "--db", dbp, "--artefact", ARTEFACT_ID, "--quantity", "insertion_loss", "--pair", "A", "--f", "100e6"]) == 0
    assert '"En"' in capsys.readouterr().out
    assert main(["suspects", "--db", dbp]) == 1
    run_id = json.loads(capsys.readouterr().out)[0]["run_id"]
    assert main(["chain", "--db", dbp, run_id]) == 0
    assert '"archive_integrity"' in capsys.readouterr().out
    assert main(["dashboard", "--db", dbp, "--artefact", ARTEFACT_ID, "--out", str(root / "d2.html")]) == 0
    assert (root / "d2.html").exists()


def test_drift_and_root_cause(network):
    root, out, db = network
    from labplatform.drift import analyse_drift, site_series
    from labplatform.roundrobin import compare_by_round
    evaluate_budgets(db, ARTEFACT_ID)
    br = compare_by_round(db, ARTEFACT_ID)
    batches, series, meta = site_series(db, "LABB", ARTEFACT_ID, br)
    assert batches == ["RR-01", "RR-02", "RR-03", "RR-04"] and "verification dA" in series and "production IL @ 600 MHz" in meta
    assert len(meta["calibration"]) == 4 and len(set(meta["calibration"])) == 4        # recalibrated every round
    findings = analyse_drift(db, ARTEFACT_ID, br, phase1=3)
    # with four rounds only the verification series of LABB can show its step; the finding carries the machinery
    by = {(f.site, f.series): f for f in findings}
    assert by == {} or all(f.n >= 5 for f in findings)
    # on a longer synthetic series the attribution names the check-standard evidence and dismisses routine recalibration
    from labplatform.drift import DriftFinding, attribute, change_point
    x = [0.03, 0.031, 0.029, 0.03, 0.032, 0.09, 0.095, 0.1]
    f = DriftFinding("S", "verification dA", "dB", 8, [f"RR-{i:02d}" for i in range(1, 9)], x, {"significant": False, "p": 0.5},
                     change_point(x), None, True)
    meta2 = {"operator": ["a"] * 8, "procedure": ["p"] * 8, "firmware": ["1"] * 8, "kit_serial": ["k"] * 8,
             "calibration": [f"C{i}" for i in range(8)], "verification_rl_db": [41, 42, 41, 42, 41, 42, 41, 42],
             "ambient_c": [23] * 8, "humidity_pct": [45] * 8, "sample_temperature_c": [23] * 8, "production IL @ 600 MHz": [8.0] * 8}
    attribute(f, meta2, {"verification dA": ("dB", x)})
    assert f.change["significant"] and f.change["index"] == 5
    recal = next(c for c in f.causes if c["candidate"] == "recalibration")
    assert recal["score"] == 0.0 and "routine" in recal["kind"]
    assert "instrument-side evidence" in f.conclusion and "RR-06" in f.conclusion
    # an operator change at the change point is named
    meta3 = dict(meta2, operator=["a"] * 5 + ["b"] * 3)
    g = DriftFinding("S", "IL @ 600 MHz", "dB", 8, f.labels, [7.9] * 5 + [7.95] * 3, {"significant": False, "p": 0.5}, change_point([7.9] * 5 + [7.95, 7.951, 7.949]), None, True)
    attribute(g, meta3, {"verification dA": ("dB", x)})
    assert g.causes[0]["candidate"] == "operator" and "operator" in g.conclusion


def test_trust_card(network):
    root, out, db = network
    from labplatform.trust import format_card, trust_card
    evaluate_budgets(db, ARTEFACT_ID)
    sus = suspect_runs(db, record_events=False)
    run = db.q("SELECT id FROM runs WHERE site='LABB' AND sample_id=? AND trust='trusted' ORDER BY t_unix DESC LIMIT 1", (ARTEFACT_ID,))[0]["id"]
    card = trust_card(db, run, "insertion_loss", "A", 600e6, ARTEFACT_ID, "consensus", sus)
    assert card["value"] > 7 and card["U_k2"] > 0 and card["raw_data"]["integrity_ok"]
    assert card["calibration"]["verification"]["status"] == "pass" and card["instrument"]["serial"] == "VNA-4004"
    assert set(card["uncertainty"]) >= {"u_repeat", "u_reprod", "u_cal", "u_fixture", "u_instrument", "u_environment", "u_noise"}
    assert card["doe"]["round"] == "RR-04" and card["doe"]["status"] in ("compatible", "questionable", "incompatible")
    assert "SUSPECT" in card["status"]
    text = format_card(card)
    for key in ("Sample", "Measurement", "Site", "Instrument", "Calibration", "Fixture", "Temperature", "Procedure", "Raw data", "Software", "Uncertainty", "Round", "Status"):
        assert key in text
    with pytest.raises(KeyError):
        trust_card(db, run, "no_such_quantity")
    # compatibility statements
    c = compare(db, ARTEFACT_ID, "insertion_loss", "A", 600e6)
    d = c.to_dict()
    assert all("differs from the reference by" in s["statement"] or "is the reference" in s["statement"] for s in d["sites"])
    assert all(s["status"] in ("compatible", "questionable", "incompatible") for s in d["sites"])
