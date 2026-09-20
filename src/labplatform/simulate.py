"""
A simulated multi-site laboratory network, built from labauto's simulators.

Four sites, each a complete labauto laboratory with its own analyser (serial, residual-error
signature, calibration quality), ambient temperature and operator, measure the same round-robin
artefact - a 15 m 1000BASE-T1 pair, identical physics at every site because the model is seeded
from the sample id - in periodic rounds, three re-connections per round, recalibrating before each
round.  Site 4's test cable degrades: from round 7 its calibration 'quality' factor is 2.2 (a
transmission bias of ~0.04 dB and a larger residual ripple; the verification still passes the
0.1 dB tolerance), and in the last round 4.5 (the verification fails, the engine refuses to
measure, and the platform must flag the previous round's results as suspect).  Each site also
measures a few of its own production samples per round.

Everything the platform then does - ingest, uncertainty, En scores, R&R, control charts,
traceability - reads the sites' archives exactly as it would read real ones.
"""
from __future__ import annotations

import json
import logging
import shutil
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
from cablecheck.io.touchstone import write_touchstone
from labauto.barcode import make_internal_barcode
from labauto.calibration import CalRecord, CalStore
from labauto.clock import SimClock
from labauto.engine import Engine
from labauto.lab import Lab, LabConfig
from labauto.sim.bench import check_standard_network

__all__ = ["SITES", "ARTEFACT_ID", "build_network", "PROCEDURE_TEXT"]

log = logging.getLogger("labplatform.simulate")

ARTEFACT_ID = "GOLD-RR26-0001"
ARTEFACT_BARCODE = make_internal_barcode("GOLD", "RR26", "0001")

SITES = {
    "LAB1": dict(name="Site 1 - Roth (reference laboratory)", serial="VNA-1001", seed=101, ambient=23.0, quality=1.0, operator="a.schmidt",
                 dmm_serial="DMM-1001", degrade_round=None, degrade_quality=1.0, sample_profiles=["good", "good", "marginal", "good"]),
    "LAB2": dict(name="Site 2 - Kitzingen", serial="VNA-2002", seed=202, ambient=24.5, quality=1.6, operator="m.rossi",
                 dmm_serial="DMM-2002", degrade_round=None, degrade_quality=1.0, sample_profiles=["good", "lossy", "good", "good"]),
    "LAB3": dict(name="Site 3 - Jelenia Gora", serial="VNA-3003", seed=303, ambient=22.0, quality=1.3, operator="j.kowalski",
                 dmm_serial="DMM-3003", degrade_round=None, degrade_quality=1.0, sample_profiles=["good", "defect", "good", "marginal"]),
    "LAB4": dict(name="Site 4 - Changzhou", serial="VNA-4004", seed=404, ambient=23.5, quality=1.2, operator="l.wang",
                 dmm_serial="DMM-4004", degrade_round=7, degrade_quality=2.2, fail_round=10, fail_quality=4.5,
                 sample_profiles=["good", "good", "good", "lossy"]),
}

PROCEDURE_TEXT = '''# Round-robin measurement of the reference artefact (identical at every site).
[procedure]
id = "roundrobin-1000base-t1-artefact"
version = "1.0"
title = "Round-robin artefact, 1000BASE-T1 single pair, 1-600 MHz"
kind = "s-parameters"
cable_type = "1000base-t1-link-segment"
description = "Same sweep, same calibration policy and same checks at every site; three re-connections per round."

[requires.vna]
ports = 4
fmin_hz = 1e6
fmax_hz = 600e6
features = ["s-parameters", "calsets"]

[requires.dmm]
optional = true
features = ["4-wire-resistance"]

[sweep]
start_hz = 1e6
stop_hz = 600e6
points = 1200
ifbw_hz = 1000
power_dbm = 0
averages = 1

[calibration]
type = "SOLT"
ports = [1, 2, 3, 4]
max_age_days = 7
max_delta_t_k = 3.0
verify = true
verify_standard = "CHECK-ATT20"
verify_reference = "check_att20.s2p"
verify_port_pairs = [[1, 2], [3, 4]]
verify_tol_db = 0.10
verify_tol_deg = 3.0
verify_rl_min_db = 30
verify_every_hours = 8

[fixture]
method = "none"

[environment]
ambient_min_c = 18
ambient_max_c = 28

[[files]]
name = "pairA"
ports = [1, 2, 3, 4]
port_map = "A+near,A-near,A+far,A-far"
prompt = "Artefact pair A: near +/- to ports 1/2, far +/- to ports 3/4"

[[aux]]
name = "loop_resistance"
instrument = "dmm"
method = "resistance_4w"

[checks]
max_repeats = 2

[metadata]
mandatory = ["operator", "sample.sample_id", "instrument.vna.identity.serial", "calibration.id",
             "calibration.verification.status", "environment.ambient_c", "procedure.hash"]

[downstream]
analysers = ["cablecheck"]
report = false
'''


def _write_site(root: Path, site: str, cfg: dict, round_no: int, t_round: float) -> Path:
    """(Re)write a site's lab.toml, registry and the calibration record for this round."""
    root.mkdir(parents=True, exist_ok=True)
    (root / "procedures").mkdir(exist_ok=True)
    (root / "procedures" / "roundrobin-1000base-t1-artefact.toml").write_text(PROCEDURE_TEXT, encoding="utf-8")
    # site-local copy of the production procedure with only the mandatory analyser (keeps the simulation fast)
    from labauto.procedure import load_procedure
    ptxt = load_procedure("sparam-1000base-t1-pair").text.replace('analysers = ["cablecheck", "zprofile", "cableanalytics"]', 'analysers = ["cablecheck"]').replace("report = true", "report = false")
    (root / "procedures" / "sparam-1000base-t1-pair.toml").write_text(ptxt, encoding="utf-8")
    quality = cfg["quality"]
    if cfg.get("degrade_round") is not None and round_no >= cfg["degrade_round"]:
        quality = cfg["degrade_quality"]
    if cfg.get("fail_round") is not None and round_no >= cfg["fail_round"]:
        quality = cfg["fail_quality"]
    (root / "lab.toml").write_text(f'''[lab]
id = "{site}"
name = "{cfg['name']}"
archive = "archive"
db = "lab.sqlite"
results_db = "results.sqlite"
registry = "registry.csv"
calibration_dir = "calibration"
procedures_dir = "procedures"

[instruments.vna]
driver = "vna.sim"
address = "sim"
serial = "{cfg['serial']}"
seed = {cfg['seed']}

[instruments.dmm]
driver = "dmm.sim"
address = "sim"

[simulation]
ambient_c = {cfg['ambient']}
fixture = false
seed = {cfg['seed']}
cal_quality = {quality}
''', encoding="utf-8")
    rows = ["barcode,sample_id,part_number,lot,cable_type,length_m,design,notes,sim_profile"]
    rows.append(f"{ARTEFACT_BARCODE},{ARTEFACT_ID},GOLD,RR26,1000base-t1-link-segment,15,D100-PE-035,artefact: round-robin reference cable,good")
    n = site[-1]
    for k, prof in enumerate(cfg["sample_profiles"]):
        for r in range(1, 9):
            bc = make_internal_barcode("C1000T1A", f"L{n}{r:02d}", f"{k + 1:04d}")
            rows.append(f"{bc},C1000T1A-L{n}{r:02d}-{k + 1:04d},C1000T1A,L{n}{r:02d},1000base-t1-link-segment,15,D100-PE-035,,{prof}")
    (root / "registry.csv").write_text("\n".join(rows) + "\n", encoding="utf-8")
    store = CalStore(root / "calibration")
    if not store.reference_path("check_att20.s2p").exists():
        write_touchstone(check_standard_network(np.linspace(1e6, 1e9, 1000)), store.reference_path("check_att20.s2p"))
    t_cal = datetime.fromtimestamp(t_round, tz=timezone.utc) - timedelta(days=1, hours=2)
    rec = CalRecord(id=f"CAL-{site}-{t_cal.strftime('%Y%m%d')}", instrument_serial=cfg["serial"], instrument_calset=f"SOLT4_{t_cal.strftime('%Y%m%d')}",
                    type="SOLT", ports=[1, 2, 3, 4], date=t_cal.isoformat(timespec="seconds"), temperature_c=cfg["ambient"] - 0.3,
                    operator=cfg["operator"], kit_serial=f"85052D-{site}", kit_due=(t_cal + timedelta(days=300)).isoformat(timespec="seconds"),
                    notes=f"weekly recalibration before round {round_no}")
    store.save(rec)
    return root


def build_network(root: str | Path, rounds: int = 10, days_between: float = 14.0, repeats: int = 3, production_per_round: int = 2,
                  start: str = "2026-03-02T08:00:00+00:00", sites: dict | None = None, quiet: bool = True) -> dict:
    """Run every round at every site; returns per-site job counts and the site roots."""
    root = Path(root)
    sites = sites or SITES
    if quiet:
        logging.getLogger("labauto").setLevel(logging.WARNING)
    t0 = datetime.fromisoformat(start).timestamp()
    out = {"sites": {}, "rounds": rounds, "artefact": ARTEFACT_ID}
    for site, cfg in sites.items():
        sroot = root / site
        shutil.rmtree(sroot, ignore_errors=True)
        counts = {"artefact": 0, "production": 0, "aborted": 0, "rejected": 0}
        for rnd in range(1, rounds + 1):
            t_round = t0 + (rnd - 1) * days_between * 86400 + list(sites).index(site) * 3600.0
            _write_site(sroot, site, cfg, rnd, t_round)
            with Lab(LabConfig.load(sroot / "lab.toml"), clock=SimClock(t_round)) as lab:
                eng = Engine(lab)
                for k in range(repeats):
                    job = eng.new_job(ARTEFACT_BARCODE, "roundrobin-1000base-t1-artefact", cfg["operator"], site, f"RR-{rnd:02d}", "roundrobin-2026")
                    eng.run_job(job)
                    counts["artefact"] += job.state == "DONE"
                    counts["aborted"] += job.state == "ABORTED"
                    counts["rejected"] += job.state == "REJECTED"
                    lab.clock.sleep(900)
                reg = list(lab.registry.entries.values())
                prod = [e for e in reg if e.sample_id != ARTEFACT_ID]
                for k in range(production_per_round):
                    e = prod[((rnd - 1) * production_per_round + k) % len(prod)]
                    job = eng.new_job(e.barcode, "sparam-1000base-t1-pair", cfg["operator"], site, f"P-{rnd:02d}", "production-2026")
                    eng.run_job(job)
                    counts["production"] += job.state == "DONE"
                    counts["aborted"] += job.state == "ABORTED"
                    counts["rejected"] += job.state == "REJECTED"
                    lab.clock.sleep(600)
        out["sites"][site] = {"root": str(sroot), **counts, **{k: v for k, v in cfg.items() if k != "sample_profiles"}}
        log.info("%s: %s", site, counts)
    (root / "network.json").write_text(json.dumps(out, indent=1))
    return out
