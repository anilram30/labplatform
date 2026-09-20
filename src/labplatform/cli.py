"""
Command line::

    labplatform simulate DIR [--rounds 10 --repeats 3]        build the simulated four-site network (labauto)
    labplatform ingest --db platform.sqlite LAB_ROOT...      ingest labauto laboratories (site id from lab.toml)
    labplatform analyse --db platform.sqlite --artefact ID [--reference consensus|LAB1] [--out summary.json]
    labplatform dashboard --db platform.sqlite --artefact ID --out dashboard.html
    labplatform compare --db platform.sqlite --artefact ID [--quantity insertion_loss --pair A --f 600e6]
    labplatform chain --db platform.sqlite RUN_ID             the traceability chain of one run (JSON)
    labplatform evidence --db platform.sqlite RUN_ID --out DIR
    labplatform suspects --db platform.sqlite                 runs under a calibration whose next verification failed
    labplatform card --db platform.sqlite RUN_ID [--quantity insertion_loss --pair A --f 600e6]   the trust card of one result
    labplatform drift --db platform.sqlite --artefact ID      drift detection + root-cause attribution per site
    labplatform demo DIR                                      simulate + ingest + analyse + dashboard
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

from . import __version__


def _db(args):
    from .schema import PlatformDB
    return PlatformDB(args.db)


def cmd_simulate(args):
    from .simulate import build_network
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    out = build_network(args.dir, rounds=args.rounds, repeats=args.repeats, production_per_round=args.production)
    print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk in ("artefact", "production", "aborted", "rejected")} for k, v in out["sites"].items()}, indent=1))
    return 0


def cmd_ingest(args):
    from .ingest import ingest_lab
    with _db(args) as db:
        for root in args.labs:
            print(json.dumps(ingest_lab(db, root)))
        print(json.dumps(db.counts()))
    return 0


def cmd_analyse(args):
    from .analyse import analyse_network, write_summary
    with _db(args) as db:
        s = analyse_network(db, args.artefact, args.reference, args.phase1)
    _print_summary(s)
    if args.out:
        write_summary(s, args.out)
        print("summary:", args.out)
    return 0


def _print_summary(s):
    print(f"{'site':<6} {'runs':>5} {'trusted':>8} {'rejected':>9} {'PASS':>5} {'FAIL':>5}  last verification")
    for o in s["overview"]:
        lv = o["last_verification"] or {}
        print(f"{o['site']:<6} {o['runs']:>5} {o['trusted']:>8} {o['rejected']:>9} {o['pass']:>5} {o['fail']:>5}  {lv.get('time', '')[:10]} {lv.get('status', '')} "
              f"{('dA=%.3f dB' % lv['max_dev_db']) if lv.get('max_dev_db') is not None else ''}")
    sites = sorted({x["site"] for c in s["comparisons"] for x in c["sites"]})
    print("\nE_n scores (reference: %s)" % s["reference"])
    print(f"{'quantity':<34} {'reference':>22} " + " ".join(f"{x:>8}" for x in sites))
    for c in s["comparisons"]:
        by = {x["site"]: x for x in c["sites"]}
        print(f"{c['label']:<34} {c['x_ref']:>12.5g} ± {c['U_ref_k2']:<7.2g} " + " ".join(f"{by[x]['En']:>+7.2f}{'' if by[x]['member'] else '×'}" if x in by else f"{'—':>8}" for x in sites))
    print("\nrepeatability / reproducibility")
    for label, r in s["rnr"].items():
        a = r["anova"]
        print(f"  {label:<34} s_r={a['s_r']:.4g}  s_I={a['s_I']:.4g}  s_R={a['s_R']:.4g}  R={a['R_limit']:.4g}"
              + (f"  Cochran outlier: {r['cochran'].get('site')}" if r["cochran"].get("outlier") else "") + (f"  Grubbs outlier: {r['grubbs'].get('site')}" if r["grubbs"].get("outlier") else ""))
    if s.get("verdicts_at_risk"):
        v = s["verdicts_at_risk"]
        print(f"  production verdicts within 2 s_R ({v['band_db']:.3f} dB) of the limit: {v['n_at_risk']} of {v['n_production']}")
    print(f"\nevents: {len(s['events'])}  (action: {sum(1 for e in s['events'] if e['severity'] == 'action')})")
    for e in s["events"][:40]:
        print(f"  {e['time'][:16]} {e['severity']:<7} {e['site'] or '':<5} {e['message'][:110]}")
    if s["suspect_runs"]:
        print(f"\nsuspect runs: {len(s['suspect_runs'])}")


def cmd_dashboard(args):
    from .analyse import analyse_network, write_summary
    from .dashboard import build_dashboard
    with _db(args) as db:
        s = analyse_network(db, args.artefact, args.reference, args.phase1)
        p = build_dashboard(db, s, args.out, args.title)
        if args.summary:
            write_summary(s, args.summary)
    print("dashboard:", p)
    return 0


def cmd_compare(args):
    from .roundrobin import compare
    from .uncertainty import evaluate_budgets
    with _db(args) as db:
        evaluate_budgets(db, args.artefact)
        c = compare(db, args.artefact, args.quantity, args.pair, args.f, args.reference)
    print(json.dumps(c.to_dict() if c else None, indent=1, default=float))
    return 0


def cmd_chain(args):
    from .traceability import chain
    with _db(args) as db:
        print(json.dumps(chain(db, args.run), indent=1, default=str))
    return 0


def cmd_evidence(args):
    from .traceability import evidence_pack
    with _db(args) as db:
        p = evidence_pack(db, args.run, args.out)
    print("evidence pack:", p)
    return 0


def cmd_suspects(args):
    from .traceability import suspect_runs
    with _db(args) as db:
        s = suspect_runs(db, record_events=False)
    print(json.dumps(s, indent=1, default=str))
    return 0 if not s else 1


def cmd_card(args):
    from .traceability import suspect_runs
    from .trust import format_card, trust_card
    with _db(args) as db:
        sus = suspect_runs(db, record_events=False)
        card = trust_card(db, args.run, args.quantity, args.pair or None, args.f if args.f else None, args.artefact, args.reference, sus)
    if args.json:
        print(json.dumps(card, indent=1, default=str))
    else:
        print(format_card(card))
    return 0


def cmd_drift(args):
    from .drift import analyse_drift
    from .roundrobin import compare_by_round
    from .uncertainty import evaluate_budgets
    with _db(args) as db:
        evaluate_budgets(db, args.artefact)
        br = compare_by_round(db, args.artefact, "insertion_loss", "A", 600e6, args.reference)
        fs = analyse_drift(db, args.artefact, br)
    for f in fs:
        if f.drifting or args.all:
            print(f"{f.site:<6} {f.series:<42} {'DRIFT' if f.drifting else 'ok   '}  {f.conclusion}")
            for c in f.causes[:3]:
                print(f"       - {c['candidate']}: {c['kind']}" + (f", r = {c['r']:+.2f}" if c.get('r') is not None else "") + (f" ({c.get('from')} -> {c.get('to')})" if c.get('from') else ""))
    return 0


def cmd_demo(args):
    from .analyse import analyse_network, write_summary
    from .dashboard import build_dashboard
    from .ingest import ingest_lab
    from .schema import PlatformDB
    from .simulate import ARTEFACT_ID, build_network
    root = Path(args.dir)
    out = build_network(root, rounds=args.rounds, repeats=args.repeats, production_per_round=args.production)
    with PlatformDB(root / "platform.sqlite") as db:
        for site, info in out["sites"].items():
            print(json.dumps(ingest_lab(db, info["root"])))
        s = analyse_network(db, ARTEFACT_ID, args.reference)
        _print_summary(s)
        write_summary(s, root / "summary.json")
        build_dashboard(db, s, root / "dashboard.html")
    print("dashboard:", root / "dashboard.html")
    return 0


def build_parser():
    p = argparse.ArgumentParser(prog="labplatform", description="Multi-site metrology and data platform")
    p.add_argument("--version", action="version", version=f"labplatform {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("simulate"); s.add_argument("dir"); s.add_argument("--rounds", type=int, default=10); s.add_argument("--repeats", type=int, default=3)
    s.add_argument("--production", type=int, default=2); s.set_defaults(func=cmd_simulate)
    i = sub.add_parser("ingest"); i.add_argument("--db", required=True); i.add_argument("labs", nargs="+"); i.set_defaults(func=cmd_ingest)
    a = sub.add_parser("analyse"); a.add_argument("--db", required=True); a.add_argument("--artefact", required=True); a.add_argument("--reference", default="consensus")
    a.add_argument("--phase1", type=int, default=5); a.add_argument("--out"); a.set_defaults(func=cmd_analyse)
    d = sub.add_parser("dashboard"); d.add_argument("--db", required=True); d.add_argument("--artefact", required=True); d.add_argument("--reference", default="consensus")
    d.add_argument("--phase1", type=int, default=5); d.add_argument("--out", default="dashboard.html"); d.add_argument("--summary"); d.add_argument("--title", default="Cable laboratory network")
    d.set_defaults(func=cmd_dashboard)
    c = sub.add_parser("compare"); c.add_argument("--db", required=True); c.add_argument("--artefact", required=True); c.add_argument("--quantity", default="insertion_loss")
    c.add_argument("--pair", default="A"); c.add_argument("--f", type=float, default=600e6); c.add_argument("--reference", default="consensus"); c.set_defaults(func=cmd_compare)
    ch = sub.add_parser("chain"); ch.add_argument("--db", required=True); ch.add_argument("run"); ch.set_defaults(func=cmd_chain)
    e = sub.add_parser("evidence"); e.add_argument("--db", required=True); e.add_argument("run"); e.add_argument("--out", default="evidence"); e.set_defaults(func=cmd_evidence)
    su = sub.add_parser("suspects"); su.add_argument("--db", required=True); su.set_defaults(func=cmd_suspects)
    cd = sub.add_parser("card"); cd.add_argument("--db", required=True); cd.add_argument("run"); cd.add_argument("--quantity", default="insertion_loss")
    cd.add_argument("--pair", default="A"); cd.add_argument("--f", type=float, default=600e6); cd.add_argument("--artefact"); cd.add_argument("--reference", default="consensus")
    cd.add_argument("--json", action="store_true"); cd.set_defaults(func=cmd_card)
    dr = sub.add_parser("drift"); dr.add_argument("--db", required=True); dr.add_argument("--artefact", required=True); dr.add_argument("--reference", default="consensus")
    dr.add_argument("--all", action="store_true"); dr.set_defaults(func=cmd_drift)
    m = sub.add_parser("demo"); m.add_argument("dir"); m.add_argument("--rounds", type=int, default=10); m.add_argument("--repeats", type=int, default=3)
    m.add_argument("--production", type=int, default=2); m.add_argument("--reference", default="consensus"); m.set_defaults(func=cmd_demo)
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
