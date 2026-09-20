"""Figures for the report from the simulated four-site network.  ``python docs/make_figures.py [--reuse DIR]``"""
from __future__ import annotations

import base64
import json
import shutil
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from labplatform import dashboard as dash
from labplatform.analyse import analyse_network, write_summary
from labplatform.ingest import ingest_lab
from labplatform.schema import PlatformDB
from labplatform.simulate import ARTEFACT_ID, build_network

OUT = Path(__file__).parent / "figures"
OUT.mkdir(exist_ok=True)
WORK = Path(__file__).parent / "_figwork"
C1, C2, C3, C4, CK, CG = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#0b0b0b", "#9a9a96"


def save_data_uri(uri: str, path: Path):
    path.write_bytes(base64.b64decode(uri.split(",", 1)[1]))


def style(ax):
    ax.grid(True, color="#e6e6e3", lw=0.7)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def fig_budget(db, summary):
    """Stacked uncertainty components per site for IL @ 600 MHz (mean over the site's artefact runs)."""
    sites = db.sites()
    comps = ["u_cal", "u_repeat", "u_reprod", "u_environment", "u_instrument", "u_noise"]
    vals = {c: [] for c in comps}
    for s in sites:
        rows = db.q("SELECT s.u_json FROM results s JOIN runs r ON r.id=s.run_id WHERE r.site=? AND r.sample_id=? AND s.quantity='insertion_loss' AND ABS(s.f_hz-600e6)<1 AND r.trust='trusted'",
                    (s, ARTEFACT_ID))
        bs = [json.loads(r["u_json"]) for r in rows if r["u_json"]]
        for c in comps:
            vals[c].append(float(np.sqrt(np.mean([b[c] ** 2 for b in bs]))) * 1000 if bs else 0.0)
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.6))
    ax = axes[0]
    x = np.arange(len(sites))
    bottom = np.zeros(len(sites))
    for c, col, lab in zip(comps, [C1, C3, "#7b5ea7", C4, CG, C2], ["calibration (verification ΔA/√3)", "repeatability (pooled repeats)", "reproducibility beyond the cal bound",
                                                                   "environment (temperature)", "instrument (declared)", "trace noise"]):
        v = np.array(vals[c]) ** 2
        ax.bar(x, v, bottom=bottom, color=col, width=0.6, label=lab)
        bottom += v
    ax.set_xticks(x), ax.set_xticklabels(sites)
    ax.set_ylabel("variance contribution / mdB²")
    ax.set_title("uncertainty budget of a single IL @ 600 MHz result (u_c² split)", loc="left", fontsize=10)
    ax.legend(fontsize=7.5, frameon=False)
    style(ax)
    ax = axes[1]
    c = next(c for c in summary["comparisons"] if c["label"].endswith("600 MHz"))
    for i, s in enumerate(c["sites"]):
        ax.errorbar(i, (s["x"] - c["x_ref"]) * 1000, yerr=2 * s["u"] * 1000, fmt="o", color=[C1, C2, C3, C4][i % 4], capsize=4, ms=6)
        ax.text(i + 0.12, (s["x"] - c["x_ref"]) * 1000, f"E_n {s['En']:+.2f}\nn = {s['n']}", fontsize=8, va="center")
    ax.axhspan(-2 * c["u_ref"] * 1000, 2 * c["u_ref"] * 1000, color=CG, alpha=0.2, label="U of the consensus")
    ax.axhline(0, color=CK, lw=0.8)
    ax.set_xticks(range(len(c["sites"]))), ax.set_xticklabels([s["site"] for s in c["sites"]])
    ax.set_ylabel("x_i − x_ref / mdB")
    ax.set_title("all-rounds comparison, IL @ 600 MHz (bars: U, k = 2)", loc="left", fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    style(ax)
    fig.tight_layout()
    fig.savefig(OUT / "budget_and_comparison.png", dpi=140)
    plt.close(fig)


def fig_rnr(summary):
    labels = [k for k in summary["rnr"] if k.startswith("insertion_loss")]
    fig, ax = plt.subplots(figsize=(8, 3.4))
    x = np.arange(len(labels))
    for k, (comp, col, lab) in enumerate([("var_error", C1, "repeatability s_r²"), ("var_condition", C3, "round-to-round (calibration) s_c²"), ("var_site", C2, "between-site s_L²")]):
        v = np.array([summary["rnr"][l]["anova"][comp] for l in labels]) * 1e6
        ax.bar(x + (k - 1) * 0.25, v, width=0.24, color=col, label=lab)
    ax.set_xticks(x), ax.set_xticklabels([l.replace("insertion_loss[A] @ ", "IL @ ") for l in labels], fontsize=9)
    ax.set_ylabel("variance component / mdB²")
    ax.set_title("nested ANOVA of the artefact: where the network's spread comes from", loc="left", fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    style(ax)
    fig.tight_layout()
    fig.savefig(OUT / "rnr.png", dpi=140)
    plt.close(fig)


def main():
    reuse = None
    if "--reuse" in sys.argv:
        reuse = Path(sys.argv[sys.argv.index("--reuse") + 1])
    root = reuse or WORK
    if reuse is None:
        shutil.rmtree(WORK, ignore_errors=True)
        out = build_network(WORK)
        sites = {k: v["root"] for k, v in out["sites"].items()}
    else:
        sites = {p.name: str(p) for p in sorted(root.iterdir()) if (p / "lab.toml").exists()}
    dbp = root / "platform_fig.sqlite"
    if dbp.exists():
        dbp.unlink()
    with PlatformDB(dbp) as db:
        for site, r in sites.items():
            ingest_lab(db, r)
        summary = analyse_network(db, ARTEFACT_ID)
        write_summary(summary, OUT / "summary.json")
        dash.build_dashboard(db, summary, OUT / "dashboard.html")
        save_data_uri(dash.fig_en_freq(summary["en_vs_frequency"]), OUT / "en_vs_frequency.png")
        save_data_uri(dash.fig_by_round(summary["by_round"]), OUT / "by_round.png")
        save_data_uri(dash.fig_charts(summary["charts"]), OUT / "control_charts.png")
        save_data_uri(dash.fig_stability(summary["stability"]), OUT / "stability.png")
        save_data_uri(dash.fig_drift(summary["drift"]), OUT / "drift.png")
        from labplatform.trust import format_card
        (OUT / "trust_card_example.txt").write_text("\n\n".join(format_card(c) for c in summary["trust_cards"].values()))
        (OUT / "drift_findings.json").write_text(json.dumps([f for f in summary["drift"] if f["drifting"]], indent=1, default=float))
        fig_budget(db, summary)
        fig_rnr(summary)
        # tables for the report
        (OUT / "en_table.json").write_text(json.dumps(summary["comparisons"], indent=1, default=float))
        (OUT / "rnr_table.json").write_text(json.dumps({k: v["anova"] for k, v in summary["rnr"].items()}, indent=1, default=float))
        (OUT / "events.json").write_text(json.dumps(summary["events"], indent=1, default=str))
        (OUT / "suspects.json").write_text(json.dumps(summary["suspect_runs"], indent=1, default=str))
        (OUT / "overview.json").write_text(json.dumps(summary["overview"], indent=1, default=str))
        from labplatform.traceability import chain
        run = summary["suspect_runs"][0]["run_id"] if summary["suspect_runs"] else db.q("SELECT id FROM runs WHERE trust='trusted' LIMIT 1")[0]["id"]
        c = chain(db, run)
        c["results"] = c["results"][:4]
        (OUT / "chain_example.json").write_text(json.dumps(c, indent=1, default=str)[:12000])
    if reuse is None:
        shutil.rmtree(WORK, ignore_errors=True)
    print("done:", sorted(p.name for p in OUT.iterdir()))


if __name__ == "__main__":
    main()
