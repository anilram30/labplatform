"""
The trust layer: a result is not a number, it is a number with its provenance, its uncertainty and its
standing in the network.

    RL_min = 26.85 dB   U = 0.18 dB (k = 2)
    Sample ........ GOLD-RR26-0001 (artefact, 15 m, 1000base-t1-link-segment)
    Measurement ... 2026-06-08T11:15 UTC, procedure roundrobin-1000base-t1-artefact v1.0, attempt 1
    Site .......... LAB4 - Site 4 - Changzhou, operator l.wang
    Instrument .... labauto-sim VNA4-8500 VNA-4004, firmware A.01.20
    Calibration ... CAL-LAB4-20260607 (SOLT, 2026-06-07, kit 85052D-LAB4 due 2027-04-03), verified 2026-06-08 pass dA = 0.056 dB
    Fixture ....... none
    Temperature ... 23.5 C (ambient), corrected to 23 C
    Procedure ..... sha256 3c1f...
    Raw data ...... pairA.s4p sha256 9a2b...   manifest 41c0...   archive integrity OK
    Software ...... labauto 0.1.1, cablecheck 0.1.1
    Uncertainty ... u_repeat 4.9 mdB, u_reprod 0, u_cal 32 mdB, u_fixture 0, u_instrument 10 mdB, u_environment 9.9 mdB, u_noise 0.2 mdB
    Round ......... RR-08; DoE d = -0.041 dB, U(d) = 0.17 dB, E_n = -0.48
    Status ........ compatible; SUSPECT - a later verification on this instrument failed (2026-06-22)

``trust_card`` builds this record for any (run, quantity) as a dict and as text; the dashboard shows one
and the evidence pack embeds them all.
"""
from __future__ import annotations

import numpy as np

from .roundrobin import compare
from .schema import PlatformDB
from .traceability import chain

__all__ = ["trust_card", "format_card"]


def trust_card(db: PlatformDB, run_id: str, quantity: str = "insertion_loss", pair: str | None = "A", f_hz: float | None = 600e6,
               artefact_id: str | None = None, reference: str = "consensus", suspects: list[dict] | None = None) -> dict:
    c = chain(db, run_id)
    run = c["run"]
    res = next((r for r in c["results"] if r["quantity"] == quantity and (r["pair"] == pair or pair is None)
                and ((f_hz is None and r["f_hz"] is None) or (f_hz is not None and r["f_hz"] is not None and abs(r["f_hz"] - f_hz) < 1))), None)
    if res is None:
        raise KeyError(f"{run_id} has no result {quantity}[{pair}] @ {f_hz}")
    b = res.get("budget") or {}
    cal = c.get("calibration") or {}
    ver = (cal.get("verifications") or [{}])[-1] if cal else {}
    inst = c.get("instrument") or {}
    site = c.get("site") or {}
    sample = c.get("sample") or {}
    files = [f for f in c["files"] if f["path"].endswith((".s4p", ".s2p", ".s8p"))]
    doe = None
    art = artefact_id or (run["sample_id"] if sample.get("kind") == "artefact" else None)
    if art and run["sample_id"] == art and run.get("batch"):
        rows = db.q("SELECT MIN(t_unix) AS t0, MAX(t_unix) AS t1 FROM runs WHERE sample_id=? AND batch=?", (art, run["batch"]))
        try:
            cm = compare(db, art, quantity, pair, f_hz, reference, t_min=rows[0]["t0"] - 1, t_max=rows[0]["t1"] + 1)
        except KeyError:
            cm = None
        if cm is not None:
            sv = next((s for s in cm.sites if s.site == run["site"]), None)
            if sv is not None:
                doe = {"round": run["batch"], "d": sv.d, "U_d": 2 * sv.u_d if np.isfinite(sv.u_d) else None, "En": sv.en, "status": sv.status,
                       "x_ref": cm.x_ref, "U_ref": 2 * cm.u_ref, "reference": reference}
    suspect = None
    if suspects is not None:
        suspect = next((s for s in suspects if s["run_id"] == run_id), None)
    status = (doe["status"] if doe else ("trusted" if run.get("trust") == "trusted" else run.get("trust")))
    if suspect:
        status += f"; SUSPECT - a later verification on this instrument failed ({suspect['failed_verification'][:10]})"
    integrity = c.get("archive_integrity") or {}
    fixture = "none"
    proc_text = db.q("SELECT text FROM procedures WHERE sha256=?", (run.get("procedure_sha"),))
    if proc_text and proc_text[0]["text"]:
        for line in proc_text[0]["text"].splitlines():
            if line.strip().startswith("method"):
                fixture = line.split("=")[1].strip().strip('"')
                break
    return {
        "quantity": quantity, "pair": pair, "f_hz": f_hz, "value": res["value"], "value_raw": b.get("value_raw"), "unit": res["unit"],
        "U_k2": b.get("U_k2"), "u_c": b.get("u_c"),
        "sample": {"id": run["sample_id"], "kind": sample.get("kind"), "length_m": sample.get("length_m"), "cable_type": sample.get("cable_type"),
                   "part_number": sample.get("part_number"), "lot": sample.get("lot")},
        "measurement": {"started": run.get("started"), "procedure_id": run.get("procedure_id"), "procedure_sha256": run.get("procedure_sha"),
                        "attempts": run.get("attempts"), "trust": run.get("trust"), "verdict": run.get("verdict"), "sweep": {"points": run.get("points"), "ifbw_hz": run.get("ifbw_hz"), "power_dbm": run.get("power_dbm"), "averages": run.get("averages")}},
        "site": {"id": run["site"], "name": site.get("name"), "operator": run.get("operator")},
        "instrument": {"manufacturer": inst.get("manufacturer"), "model": inst.get("model"), "serial": inst.get("serial"), "firmware": inst.get("firmware")},
        "calibration": {"record": cal.get("record_id"), "type": cal.get("type"), "date": cal.get("date"), "kit_serial": cal.get("kit_serial"), "kit_due": cal.get("kit_due"),
                        "age_days": run.get("cal_age_days"), "delta_t_k": run.get("cal_delta_t_k"),
                        "verification": {"time": ver.get("time"), "status": ver.get("status"), "max_dev_db": ver.get("max_dev_db"), "min_rl_db": ver.get("min_rl_db"), "standard": ver.get("standard")}},
        "fixture": fixture,
        "temperature": {"sample_c": run.get("temperature_c"), "method": run.get("temperature_method"), "ambient_c": run.get("ambient_c"), "humidity_pct": run.get("humidity_pct"), "corrected_to_c": 23.0},
        "raw_data": {"files": files, "manifest_sha256": run.get("manifest_sha256"), "archive_dir": run.get("archive_dir"), "integrity_ok": integrity.get("ok")},
        "software": c.get("software"),
        "uncertainty": {k: b.get(k) for k in ("u_repeat", "u_reprod", "u_cal", "u_fixture", "u_instrument", "u_environment", "u_noise", "u_res")} | {"rep_source": b.get("rep_source"), "reprod_source": b.get("reprod_source")},
        "doe": doe, "status": status,
    }


def _fmt(v, nd=4):
    return "—" if v is None else (f"{v:.{nd}g}" if isinstance(v, (int, float)) else str(v))


def format_card(card: dict) -> str:
    u = card["unit"] or ""
    q = f"{card['quantity']}[{card['pair']}]" + (f" @ {card['f_hz'] / 1e6:.0f} MHz" if card["f_hz"] else "") if card["pair"] else card["quantity"]
    lines = [f"{q} = {card['value']:.4g} {u}   U = {_fmt(card['U_k2'], 2)} {u} (k = 2)"]
    s, m, si, i, c, t, r, sw, ub, d = (card["sample"], card["measurement"], card["site"], card["instrument"], card["calibration"], card["temperature"],
                                       card["raw_data"], card["software"] or {}, card["uncertainty"], card["doe"])
    lines.append(f"Sample ........ {s['id']} ({s.get('kind')}, {_fmt(s.get('length_m'))} m, {s.get('cable_type')})")
    lines.append(f"Measurement ... {(m.get('started') or '')[:16]} UTC, procedure {m.get('procedure_id')}, attempt {m.get('attempts')}, trust {m.get('trust')}, verdict {m.get('verdict')}")
    lines.append(f"Site .......... {si['id']} - {si.get('name')}, operator {si.get('operator')}")
    lines.append(f"Instrument .... {i.get('manufacturer')} {i.get('model')} {i.get('serial')}, firmware {i.get('firmware')}")
    v = c.get("verification") or {}
    lines.append(f"Calibration ... {c.get('record')} ({c.get('type')}, {(c.get('date') or '')[:10]}, kit {c.get('kit_serial')} due {(c.get('kit_due') or '')[:10]}), "
                 f"verified {(v.get('time') or '')[:10]} {v.get('status')} dA = {_fmt(v.get('max_dev_db'), 3)} dB, RL = {_fmt(v.get('min_rl_db'), 3)} dB")
    lines.append(f"Fixture ....... {card['fixture']}")
    lines.append(f"Temperature ... {_fmt(t.get('sample_c'), 3)} C ({t.get('method')}), corrected to {t.get('corrected_to_c'):.0f} C; ambient {_fmt(t.get('ambient_c'), 3)} C")
    lines.append(f"Procedure ..... sha256 {(m.get('procedure_sha256') or '')[:16]}…")
    files = ", ".join(f"{f['path']} sha256 {f['sha256'][:12]}…" for f in r["files"][:2])
    lines.append(f"Raw data ...... {files}; manifest {(r.get('manifest_sha256') or '')[:12]}…; archive integrity {'OK' if r.get('integrity_ok') else ('FAILED' if r.get('integrity_ok') is False else 'not checked')}")
    lines.append(f"Software ...... labauto {sw.get('labauto')}, cablecheck {sw.get('cablecheck')}, numpy {sw.get('numpy')}")
    scale = 1000 if u == "dB" else 1
    unit_u = "mdB" if u == "dB" else u
    lines.append("Uncertainty ... " + ", ".join(f"{k[2:]} {_fmt((ub.get(k) or 0) * scale, 2)} {unit_u}" for k in ("u_repeat", "u_reprod", "u_cal", "u_fixture", "u_instrument", "u_environment", "u_noise")))
    if d:
        lines.append(f"Round ......... {d['round']}; DoE d = {d['d']:+.4g} {u}, U(d) = {_fmt(d['U_d'], 2)} {u}, E_n = {d['En']:+.2f} (reference {d['reference']}: {d['x_ref']:.4g} ± {d['U_ref']:.2g} {u})")
    lines.append(f"Status ........ {card['status']}")
    return "\n".join(lines)
