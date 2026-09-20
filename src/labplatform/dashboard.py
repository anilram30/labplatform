"""
The dashboard: one self-contained HTML file built from the platform database.

Sections: network overview, round-robin scores (E_n table and E_n(f) curves), repeatability and
reproducibility, control charts per site (individuals + EWMA), artefact stability, traceability
(suspect runs, run lookup) and the event log.  Figures are matplotlib PNGs embedded as data URIs;
no external resources, so the file can be mailed, archived or opened from a USB stick at an audit.
"""
from __future__ import annotations

import base64
import html
import io
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .schema import PlatformDB

__all__ = ["build_dashboard"]

C1, C2, C3, C4, CK, CG, CPASS, CFAIL = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#0b0b0b", "#9a9a96", "#008300", "#e34948"
SITE_COLOURS = [C1, C2, C3, C4, "#7b5ea7", "#d63384"]

CSS = """
:root{--ink:#0b0b0b;--muted:#5f5f5b;--line:#e6e6e3;--blue:#2a78d6;--orange:#eb6834;--green:#008300;--red:#e34948;--yellow:#eda100;--bg:#ffffff;--card:#fafaf8}
@media (prefers-color-scheme: dark){:root:not([data-theme=light]){--ink:#ececea;--muted:#a8a8a4;--line:#33332f;--bg:#141413;--card:#1d1d1b}}
:root[data-theme=dark]{--ink:#ececea;--muted:#a8a8a4;--line:#33332f;--bg:#141413;--card:#1d1d1b}
body{font-family:system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;color:var(--ink);background:var(--bg);margin:0;padding:0 16px 48px;line-height:1.45}
main{max-width:1180px;margin:0 auto}
h1{font-size:1.5rem;margin:24px 0 4px}h2{font-size:1.15rem;margin:36px 0 8px;border-bottom:1px solid var(--line);padding-bottom:4px}
h3{font-size:1rem;margin:20px 0 6px}
p.sub{color:var(--muted);margin:0 0 12px}
table{border-collapse:collapse;width:100%;font-size:.86rem;margin:8px 0 12px}th,td{padding:5px 8px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}
th{color:var(--muted);font-weight:600}td.num{text-align:right;font-variant-numeric:tabular-nums}
.ok{background:rgba(0,131,0,.12)}.warn{background:rgba(237,161,0,.18)}.bad{background:rgba(227,73,72,.16)}
.badge{display:inline-block;padding:1px 7px;border-radius:9px;font-size:.78rem;border:1px solid var(--line)}
.card{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:12px 14px;margin:10px 0}
img{max-width:100%;height:auto;display:block;margin:6px 0}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(280px,1fr));gap:12px}
.kpi{font-size:1.4rem;font-weight:600}.kpi small{font-size:.8rem;color:var(--muted);font-weight:400;display:block}
details summary{cursor:pointer;color:var(--muted)}
code{font-size:.82rem}
.tblwrap{overflow-x:auto}
"""


def _png(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=120, bbox_inches="tight")
    plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def _style(ax):
    ax.grid(True, color="#e6e6e3", lw=0.7)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def _esc(x) -> str:
    return html.escape("" if x is None else str(x))


def _en_class(en: float) -> str:
    a = abs(en)
    return "ok" if a <= 1 else ("warn" if a <= 1.5 else "bad")


def fig_en_table(comparisons: list[dict]) -> str:
    sites = sorted({s["site"] for c in comparisons for s in c["sites"]})
    rows = ["<table><tr><th>quantity</th><th>reference (U, k=2)</th>" + "".join(f"<th>{_esc(s)}</th>" for s in sites) + "</tr>"]
    for c in comparisons:
        by = {s["site"]: s for s in c["sites"]}
        cells = ""
        for s in sites:
            v = by.get(s)
            if v is None:
                cells += "<td>—</td>"
            else:
                cells += (f"<td class='num {_en_class(v['En'])}' title='x = {v['x']:.5g}, U = {v['U_k2']:.3g}, n = {v['n']}, z = {v['z']:+.2f}'>"
                          f"{v['En']:+.2f}{'' if v['member'] else ' ×'}</td>")
        unit = c.get("unit") or ""
        rows.append(f"<tr><td>{_esc(c['label'])}</td><td class='num'>{c['x_ref']:.5g} ± {c['U_ref_k2']:.2g} {_esc(unit)}</td>{cells}</tr>")
    rows.append("</table><p class='sub'>E<sub>n</sub> = (x<sub>i</sub> − x<sub>ref</sub>) / (2 u(d<sub>i</sub>)); |E<sub>n</sub>| ≤ 1 satisfactory, ≤ 1.5 questionable (amber), > 1.5 unsatisfactory (red). "
                "× = site excluded from the consensus (largest consistent subset). Hover a cell for x, U, n and the robust z-score.</p>")
    return "".join(rows)


def fig_en_freq(enf: dict) -> str:
    if not enf:
        return ""
    f = np.array(enf["f_hz"]) / 1e6
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.6))
    ax = axes[0]
    for i, s in enumerate(enf["sites"]):
        ax.plot(f, enf["En"][i], color=SITE_COLOURS[i % len(SITE_COLOURS)], lw=1.5, label=s)
    ax.axhspan(-1, 1, color="#008300", alpha=0.07)
    ax.axhline(1, color=CK, lw=0.8, ls="--"), ax.axhline(-1, color=CK, lw=0.8, ls="--")
    ax.set_xlabel("frequency / MHz"), ax.set_ylabel("E_n")
    ax.set_title("E_n of the artefact's insertion loss vs frequency", loc="left", fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    _style(ax)
    ax = axes[1]
    xref = np.array(enf["x_ref"])
    for i, s in enumerate(enf["sites"]):
        ax.plot(f, np.array(enf["il"][i]) - xref, color=SITE_COLOURS[i % len(SITE_COLOURS)], lw=1.5, label=s)
        ax.fill_between(f, np.array(enf["il"][i]) - xref - 2 * np.array(enf["u"][i]), np.array(enf["il"][i]) - xref + 2 * np.array(enf["u"][i]),
                        color=SITE_COLOURS[i % len(SITE_COLOURS)], alpha=0.08)
    ax.plot(f, 2 * np.array(enf["u_ref"]), color=CK, lw=0.8, ls=":"), ax.plot(f, -2 * np.array(enf["u_ref"]), color=CK, lw=0.8, ls=":")
    ax.axhline(0, color=CK, lw=0.8)
    ax.set_xlabel("frequency / MHz"), ax.set_ylabel("IL − consensus / dB")
    ax.set_title("degrees of equivalence with U (k = 2) bands; dotted: U of the reference", loc="left", fontsize=10)
    _style(ax)
    fig.tight_layout()
    return _png(fig)


def fig_by_round(br: dict) -> str:
    if not br or not br.get("rounds"):
        return ""
    sites = br["sites"]
    rounds = br["rounds"]
    t = np.array([r["t_unix"] for r in rounds]); t = (t - t[0]) / 86400.0
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.6))
    ax = axes[0]
    for i, s in enumerate(sites):
        d = np.array([r["sites"].get(s, {}).get("d", np.nan) for r in rounds])
        u = np.array([r["sites"].get(s, {}).get("u", np.nan) for r in rounds])
        col = SITE_COLOURS[i % len(SITE_COLOURS)]
        ax.errorbar(t + 0.6 * i, d * 1000, yerr=2 * u * 1000, fmt="o-", color=col, ms=4, lw=1, capsize=2, label=s)
    ax.axhline(0, color=CK, lw=0.8)
    ax.set_xlabel("days since the first round"), ax.set_ylabel("d_i = x_i − x_ref / mdB")
    ax.set_title("degree of equivalence per round, IL @ 600 MHz (bars: U, k = 2)", loc="left", fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    _style(ax)
    ax = axes[1]
    for i, s in enumerate(sites):
        en = np.array([r["sites"].get(s, {}).get("En", np.nan) for r in rounds])
        ax.plot(t, en, "o-", color=SITE_COLOURS[i % len(SITE_COLOURS)], ms=4, lw=1, label=s)
    ax.axhspan(-1, 1, color="#008300", alpha=0.07)
    ax.axhline(1, color=CK, lw=0.8, ls="--"), ax.axhline(-1, color=CK, lw=0.8, ls="--")
    ax.set_xlabel("days since the first round"), ax.set_ylabel("E_n")
    ax.set_title("E_n per round", loc="left", fontsize=10)
    _style(ax)
    fig.tight_layout()
    return _png(fig)


def fig_drift(drift: list[dict]) -> str:
    sel = [f for f in drift if f["series"] in ("IL @ 600 MHz", "verification dA", "degree of equivalence d_i (IL @ 600 MHz)")]
    if not sel:
        return ""
    sites = sorted({f["site"] for f in sel})
    names = ["IL @ 600 MHz", "verification dA", "degree of equivalence d_i (IL @ 600 MHz)"]
    fig, axes = plt.subplots(len(sites), 3, figsize=(15, 2.3 * len(sites)), squeeze=False)
    for i, site in enumerate(sites):
        for j, name in enumerate(names):
            ax = axes[i, j]
            f = next((f for f in sel if f["site"] == site and f["series"] == name), None)
            if f is None:
                ax.axis("off")
                continue
            x = np.array(f["x"]); k = np.arange(x.size)
            ax.plot(k, x, "o-", color=C1, lw=1, ms=3.5)
            cp = f["change_point"]
            if cp.get("index") is not None:
                m = cp["index"]
                ax.plot([0, m - 1], [x[:m].mean()] * 2, color=C2 if cp["significant"] else CG, lw=1.5)
                ax.plot([m, x.size - 1], [x[m:].mean()] * 2, color=C2 if cp["significant"] else CG, lw=1.5)
                if cp["significant"]:
                    ax.axvline(m - 0.5, color=C2, lw=0.8, ls="--")
            if f["trend"].get("significant"):
                ax.plot(k, f["trend"]["intercept"] + f["trend"]["slope_per_30d"] / 30 * k * 14, color=C4, lw=1, ls=":")
            ax.set_title(f"{site}: {name}" + ("  — DRIFT" if f["drifting"] else ""), loc="left", fontsize=8.5, color=CFAIL if f["drifting"] else CK)
            ax.set_xticks(k), ax.set_xticklabels([l.replace("RR-", "") for l in f["labels"]], fontsize=7)
            ax.set_ylabel(f["unit"] or "")
            _style(ax)
    fig.tight_layout()
    return _png(fig)


def fig_charts(charts: list[dict], series_filter: tuple[str, ...] = ("artefact IL @ 600 MHz", "verification |S21| deviation")) -> str:
    sel = [c for c in charts if c["series"] in series_filter]
    if not sel:
        return ""
    sites = sorted({c["site"] for c in sel})
    ncol = len(series_filter)
    fig, axes = plt.subplots(len(sites), ncol, figsize=(5.4 * ncol, 2.5 * len(sites)), squeeze=False)
    for i, site in enumerate(sites):
        for j, series in enumerate(series_filter):
            ax = axes[i, j]
            c = next((c for c in sel if c["site"] == site and c["series"] == series), None)
            if c is None:
                ax.axis("off")
                continue
            t = (np.array(c["t"]) - c["t"][0]) / 86400.0
            x = np.array(c["x"])
            ax.plot(t, x, "o-", color=C1, lw=1, ms=3.5, label="x")
            ax.plot(t, c["z"], color=C2, lw=1.2, label="EWMA")
            ax.axhline(c["cl"], color=CK, lw=0.8)
            ax.axhline(c["ucl"], color=CFAIL, lw=0.8, ls="--"), ax.axhline(c["lcl"], color=CFAIL, lw=0.8, ls="--")
            ax.plot(t, c["ewma"]["ucl"] if "ucl" in c["ewma"] else [np.nan] * len(t), color=C2, lw=0.6, ls=":")
            ax.plot(t, c["ewma"]["lcl"] if "lcl" in c["ewma"] else [np.nan] * len(t), color=C2, lw=0.6, ls=":")
            ax.axvspan(t[0], t[min(c["phase1"], len(t)) - 1], color=CG, alpha=0.1)
            for v in c["violations"]:
                k = v["index"]
                ax.plot(t[k], x[k], "o", color=CFAIL, ms=8, mfc="none", mew=1.5)
            if c["ewma"].get("first_signal") is not None:
                k = c["ewma"]["first_signal"]
                ax.plot(t[k], c["z"][k], "s", color=C2, ms=8, mfc="none", mew=1.5)
            ax.set_title(f"{site}: {series}", loc="left", fontsize=9)
            ax.set_xlabel("days" if i == len(sites) - 1 else ""), ax.set_ylabel(c["unit"] or "")
            _style(ax)
            if i == 0 and j == 0:
                ax.legend(fontsize=7, frameon=False, loc="upper left")
    fig.tight_layout()
    return _png(fig)


def fig_stability(stab: dict) -> str:
    s = stab.get("insertion_loss@600") or {}
    cons = s.get("consensus")
    if not cons:
        return ""
    fig, ax = plt.subplots(figsize=(7, 3))
    pts = cons["points"]
    t = np.array([p["t_days"] for p in pts]); t = t - t[0]
    y = np.array([p["value"] for p in pts]); u = np.array([p["u"] for p in pts])
    ax.errorbar(t, y, yerr=2 * u, fmt="o", color=C1, ms=4, capsize=2, label="per-round consensus ± U")
    ax.plot(t, cons["intercept"] + cons["slope_per_day"] * (t + pts[0]["t_days"]), color=C2, lw=1.2,
            label=f"trend {cons['slope_per_30d'] * 1000:+.1f} mdB / 30 d (p = {cons['p']:.2f}){'  — significant' if cons['significant'] else ''}")
    ax.set_xlabel("days since the first round"), ax.set_ylabel("IL @ 600 MHz / dB")
    ax.set_title("artefact stability over the campaign", loc="left", fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    _style(ax)
    fig.tight_layout()
    return _png(fig)


def build_dashboard(db: PlatformDB, summary: dict, out_path: str | Path, title: str = "Cable laboratory network") -> Path:
    out_path = Path(out_path)
    ov = summary["overview"]
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    parts = [f"<!DOCTYPE html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width, initial-scale=1'>"
             f"<title>{_esc(title)}</title><style>{CSS}</style></head><body><main>",
             f"<h1>{_esc(title)}</h1><p class='sub'>built {now} from {summary['counts']['runs']} runs at {len(ov)} sites · artefact {_esc(summary['artefact'])} · "
             f"reference: {_esc(summary['reference'])}</p>"]
    # KPIs
    n_trusted = sum(o["trusted"] for o in ov); n_rej = sum(o["rejected"] for o in ov); n_ab = sum(o["aborted"] for o in ov)
    n_en_bad = sum(1 for c in summary["comparisons"] for s in c["sites"] if abs(s["En"]) > 1)
    n_alerts = sum(1 for e in summary["events"] if e["severity"] == "action")
    parts.append("<div class='grid'>" + "".join(
        f"<div class='card'><div class='kpi'>{v}<small>{k}</small></div></div>" for k, v in [
            ("trusted runs", n_trusted), ("quarantined", n_rej), ("refused at a gate", n_ab), ("E_n scores outside ±1", n_en_bad),
            ("action-level events", n_alerts), ("suspect runs (failed verification)", len(summary["suspect_runs"]))]) + "</div>")
    # overview table
    parts.append("<h2>Sites</h2><div class='tblwrap'><table><tr><th>site</th><th>analyser</th><th>runs</th><th>trusted</th><th>flagged</th><th>rejected</th>"
                 "<th>PASS</th><th>FAIL</th><th>last verification</th><th>last calibration</th><th>events</th></tr>")
    for o in ov:
        lv = o["last_verification"] or {}
        lc = o["last_calibration"] or {}
        parts.append(f"<tr><td><b>{_esc(o['site'])}</b><br><span class='sub'>{_esc(o['name'])}</span></td><td>{_esc(o['instrument'])}</td>"
                     f"<td class='num'>{o['runs']}</td><td class='num'>{o['trusted']}</td><td class='num'>{o['flagged']}</td><td class='num'>{o['rejected']}</td>"
                     f"<td class='num'>{o['pass']}</td><td class='num'>{o['fail']}</td>"
                     f"<td class='{'ok' if lv.get('status') == 'pass' else 'bad'}'>{_esc(lv.get('time', '')[:10])} {_esc(lv.get('status', ''))} "
                     f"{('ΔA %.3f dB' % lv['max_dev_db']) if lv.get('max_dev_db') is not None else ''}</td>"
                     f"<td>{_esc(lc.get('record_id', ''))}<br><span class='sub'>{_esc(lc.get('date', '')[:10])}</span></td>"
                     f"<td>{' '.join(f'<span class=badge>{_esc(k)}: {v}</span>' for k, v in o['events'].items()) or '—'}</td></tr>")
    parts.append("</table></div>")
    # round robin
    parts.append("<h2>Round-robin comparison</h2><p class='sub'>Each site's value of the circulating artefact against the reference value, with its uncertainty budget "
                 "(calibration verification, trace noise, repeatability, temperature correction).</p>")
    parts.append("<div class='tblwrap'>" + fig_en_table(summary["comparisons"]) + "</div>")
    c600 = next((c for c in summary["comparisons"] if c["label"].endswith("600 MHz")), None)
    if c600:
        parts.append("<div class='card'><b>Compatibility statements, IL @ 600 MHz</b><ul>" + "".join(f"<li>{_esc(s['statement'])}</li>" for s in c600["sites"]) + "</ul></div>")
    img = fig_en_freq(summary["en_vs_frequency"])
    if img:
        parts.append(f"<img src='{img}' alt='E_n versus frequency'>")
    img = fig_by_round(summary.get("by_round", {}))
    if img:
        parts.append("<h3>Round by round</h3><p class='sub'>The all-rounds comparison averages a site's history; per round, a site's uncertainty is that round's "
                     "calibration and re-connection scatter only, so a site that starts drifting shows as a growing degree of equivalence.</p>")
        parts.append(f"<img src='{img}' alt='per-round comparison'>")
    # R&R
    parts.append("<h2>Repeatability and reproducibility</h2><div class='tblwrap'><table><tr><th>quantity</th><th>s<sub>r</sub> (repeat)</th><th>s<sub>I</sub> (site, day-to-day)</th>"
                 "<th>s<sub>R</sub> (network)</th><th>r = 2.8 s<sub>r</sub></th><th>R = 2.8 s<sub>R</sub></th><th>Cochran</th><th>Grubbs</th></tr>")
    for label, r in summary["rnr"].items():
        a = r["anova"]
        co, gr = r["cochran"], r["grubbs"]
        parts.append(f"<tr><td>{_esc(label)}</td><td class='num'>{a['s_r']:.4g}</td><td class='num'>{a['s_I']:.4g}</td><td class='num'>{a['s_R']:.4g}</td>"
                     f"<td class='num'>{a['r_limit']:.4g}</td><td class='num'>{a['R_limit']:.4g}</td>"
                     f"<td class='{'bad' if co.get('outlier') else 'ok'}'>{('C = %.2f (crit %.2f) %s' % (co['C'], co['critical'], co.get('site', ''))) if co.get('C') == co.get('C') else '—'}</td>"
                     f"<td class='{'bad' if gr.get('outlier') else 'ok'}'>{('G = %.2f (crit %.2f) %s' % (gr['G'], gr['critical'], gr.get('site', ''))) if gr.get('G') == gr.get('G') else '—'}</td></tr>")
    parts.append("</table></div>")
    risk = summary.get("verdicts_at_risk")
    if risk:
        parts.append(f"<p class='sub'>Production verdicts within 2 s<sub>R</sub> = {risk['band_db']:.3f} dB of the insertion-loss limit (a different site might reverse them): "
                     f"<b>{risk['n_at_risk']}</b> of {risk['n_production']}.</p>")
    # control charts
    parts.append("<h2>Control charts</h2><p class='sub'>Individuals chart (blue, 3σ limits dashed, phase-I window shaded) and EWMA (orange, λ = 0.2, L = 2.7) per site; "
                 "red rings mark Western-Electric rule violations, orange squares the first EWMA signal.</p>")
    img = fig_charts(summary["charts"])
    if img:
        parts.append(f"<img src='{img}' alt='control charts'>")
    img = fig_charts(summary["charts"], ("artefact fitted impedance", "trace noise (artefact runs)"))
    if img:
        parts.append(f"<details><summary>impedance and trace-noise charts</summary><img src='{img}' alt='more charts'></details>")
    # drift and root cause
    parts.append("<h2>Drift and root cause</h2><p class='sub'>Per site and series: trend, single change point (binary segmentation with a Bonferroni-corrected t-test), "
                 "EWMA. A drifting series is attributed against project B's metadata for the same rounds — recalibration, verification deviation, ambient, "
                 "operator, procedure, firmware — and against the site's own production results (did the product move with the artefact?).</p>")
    img = fig_drift(summary.get("drift", []))
    if img:
        parts.append(f"<img src='{img}' alt='drift'>")
    drifting = [f for f in summary.get("drift", []) if f["drifting"]]
    if drifting:
        parts.append("<div class='tblwrap'><table><tr><th>site</th><th>series</th><th>detectors</th><th>conclusion</th><th>top candidates</th></tr>")
        for f in drifting:
            det = []
            if f["change_point"].get("significant"):
                det.append(f"change point at {f['labels'][f['change_point']['index']]} (shift {f['change_point']['shift']:+.4g} {f['unit'] or ''}, p = {f['change_point']['p_bonferroni']:.3f})")
            if f["trend"].get("significant"):
                det.append(f"trend {f['trend']['slope_per_30d']:+.4g} {f['unit'] or ''}/30 d (p = {f['trend']['p']:.3f})")
            if f["ewma_first_signal"] is not None:
                det.append(f"EWMA at {f['labels'][f['ewma_first_signal']]}")
            top = "; ".join(f"{c['candidate']}" + (f" (r = {c['r']:+.2f})" if c.get("r") is not None else "") for c in f["causes"][:3])
            parts.append(f"<tr><td>{_esc(f['site'])}</td><td>{_esc(f['series'])}</td><td>{_esc('; '.join(det))}</td><td>{_esc(f['conclusion'])}</td><td>{_esc(top)}</td></tr>")
        parts.append("</table></div>")
    else:
        parts.append("<p class='sub'>No series drifts.</p>")
    # trust cards
    parts.append("<h2>Trust layer</h2><p class='sub'>A result as the platform hands it on: value, uncertainty budget, provenance, standing in the round, status. "
                 "One card per site for the artefact's latest trusted IL @ 600 MHz.</p>")
    from .trust import format_card
    for site, card in (summary.get("trust_cards") or {}).items():
        parts.append(f"<details><summary>{_esc(site)}: {_esc(card['quantity'])} = {card['value']:.4g} {_esc(card['unit'] or '')} — {_esc(card['status'])}</summary>"
                     f"<pre class='card' style='font-size:.8rem;overflow-x:auto'>{_esc(format_card(card))}</pre></details>")
    # stability
    img = fig_stability(summary["stability"])
    parts.append("<h2>Artefact stability</h2>")
    if img:
        parts.append(f"<img src='{img}' alt='artefact stability'>")
    parts.append("<div class='tblwrap'><table><tr><th>quantity</th><th>site</th><th>rounds</th><th>slope / 30 d</th><th>p</th><th>significant</th></tr>")
    for key, st in summary["stability"].items():
        for site, tr in st["per_site"].items():
            parts.append(f"<tr><td>{_esc(key)}</td><td>{_esc(site)}</td><td class='num'>{tr['n']}</td><td class='num'>{tr['slope_per_30d']:+.4g}</td>"
                         f"<td class='num'>{tr['p']:.3f}</td><td class='{'warn' if tr['significant'] else 'ok'}'>{'yes' if tr['significant'] else 'no'}</td></tr>")
    parts.append("</table></div>")
    # traceability
    parts.append("<h2>Traceability</h2>")
    if summary["suspect_runs"]:
        parts.append("<div class='card bad'><b>Suspect runs</b> — measured between the last passed and a failed calibration verification:<ul>" +
                     "".join(f"<li>{_esc(s['site'])} {_esc(s['run_id'])} ({_esc(s['sample_id'])}, {_esc(s['verdict'])}) — verification failed {_esc(s['failed_verification'][:16])}, ΔA = {s['max_dev_db']:.3f} dB</li>"
                             for s in summary["suspect_runs"]) + "</ul></div>")
    else:
        parts.append("<p class='sub'>No run is under a calibration whose next verification failed.</p>")
    runs = db.q("SELECT id, site, sample_id, started, trust, verdict, headline, headline_margin, calibration_record, verification_dev_db, archive_dir, manifest_sha256 "
                "FROM runs ORDER BY t_unix DESC LIMIT 400")
    parts.append("<details><summary>run register (latest 400)</summary><div class='tblwrap'><table><tr><th>started</th><th>site</th><th>sample</th><th>trust</th><th>verdict</th>"
                 "<th>headline</th><th>margin</th><th>calibration</th><th>ΔA</th><th>manifest sha256</th></tr>")
    for r in runs:
        margin = "" if r["headline_margin"] is None else "%+.2f" % r["headline_margin"]
        dev = "" if r["verification_dev_db"] is None else "%.3f" % r["verification_dev_db"]
        tcls = "ok" if r["trust"] == "trusted" else ("warn" if r["trust"] == "flagged" else "bad")
        vcls = "ok" if r["verdict"] == "PASS" else ("bad" if r["verdict"] == "FAIL" else "")
        parts.append(f"<tr><td>{_esc((r['started'] or '')[:16])}</td><td>{_esc(r['site'])}</td><td>{_esc(r['sample_id'])}</td>"
                     f"<td class='{tcls}'>{_esc(r['trust'])}</td><td class='{vcls}'>{_esc(r['verdict'] or '')}</td>"
                     f"<td>{_esc(r['headline'] or '')}</td><td class='num'>{margin}</td>"
                     f"<td>{_esc(r['calibration_record'])}</td><td class='num'>{dev}</td>"
                     f"<td><code title='{_esc(r['archive_dir'])}'>{_esc((r['manifest_sha256'] or '')[:12])}</code></td></tr>")
    parts.append("</table></div></details>")
    # events
    parts.append("<h2>Events</h2><div class='tblwrap'><table><tr><th>time</th><th>kind</th><th>site</th><th>severity</th><th>message</th></tr>")
    for e in summary["events"]:
        parts.append(f"<tr><td>{_esc(e['time'][:16])}</td><td>{_esc(e['kind'])}</td><td>{_esc(e['site'])}</td>"
                     f"<td class='{'bad' if e['severity'] == 'action' else 'warn'}'>{_esc(e['severity'])}</td><td>{_esc(e['message'])}</td></tr>")
    if not summary["events"]:
        parts.append("<tr><td colspan=5>none</td></tr>")
    parts.append("</table></div>")
    parts.append("<p class='sub'>labplatform · every number links to a sealed archive directory by its manifest hash; uncertainties are GUM budgets from the sidecar ingredients; "
                 "E<sub>n</sub> per ISO 13528 with a Cox largest-consistent-subset reference; precision per ISO 5725-2 nested ANOVA.</p></main></body></html>")
    out_path.write_text("".join(parts), encoding="utf-8")
    return out_path
