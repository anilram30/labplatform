# labplatform — multi-site metrology platform

**Compares laboratories against each other, quantifies the uncertainty of every result, and detects measurement-system drift.**

[![CI](https://github.com/anilram30/labplatform/actions/workflows/ci.yml/badge.svg)](https://github.com/anilram30/labplatform/actions/workflows/ci.yml)
[![Python 3.11+](https://img.shields.io/badge/python-3.11%20%7C%203.12%20%7C%203.13-blue)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](LICENSE)
[![Tests](https://img.shields.io/badge/tests-16%20passing-brightgreen)](tests/)
[![Report](https://img.shields.io/badge/report-11%20pages-informational)](docs/report.pdf)

> **Part of the [HF cable toolchain](https://github.com/anilram30/hf-cable-toolchain)** — seven packages that take a high-frequency cable from a raw measurement to a predicted Ethernet link.
> 
> [A · cablecheck](https://github.com/anilram30/cablecheck)  ·  [B · labauto](https://github.com/anilram30/labauto)  ·  [C · shieldeval](https://github.com/anilram30/shieldeval)  ·  [D · zprofile](https://github.com/anilram30/zprofile)  ·  [E · cableanalytics](https://github.com/anilram30/cableanalytics)  ·  **F · labplatform**  ·  [G · linktwin](https://github.com/anilram30/linktwin)

---

## The problem it solves

When several laboratories measure the same product, the hard question is not what each one measured but whether they agree — and how you would know. `labplatform` gives every single result a full uncertainty budget, so a disagreement of 0.8 dB against a combined uncertainty of 1.2 dB is correctly reported as *not distinguishable* rather than as a problem. On top of that it runs the interlaboratory comparison, estimates repeatability and reproducibility, keeps control charts, detects drift over months and attributes it to the calibration, fixture or environment change that caused it.

## At a glance

|  |  |
|---|---|
| **Takes** | Sealed archives and laboratory databases from any number of `labauto` sites |
| **Produces** | An uncertainty budget per result, E_n comparison scores, repeatability and reproducibility, drift alerts, trust cards and a dashboard |
| **Checked against** | Worked values from ISO 13528 and ISO 5725-2, and a simulated four-site network with a planted fault |
| **Technical report** | [`docs/report.pdf`](docs/report.pdf) — 11 pages, 42 references, every method stated with its mathematics and its limitations |
| **Tests** | 16, run against Python 3.11, 3.12 and 3.13 on every push |
| **Data** | Entirely synthetic. No proprietary or customer measurements are used anywhere in this toolchain. |

## Install

Python 3.11 or newer.

```sh
pip install "git+https://github.com/anilram30/cablecheck.git" \
            "git+https://github.com/anilram30/labauto.git" \
            "git+https://github.com/anilram30/labplatform.git"
```

---

## What it does

The layer above `labauto` (project B) — the trust layer of the network: one canonical data model for
every site, GUM uncertainty budgets (repeat, reprod, cal, fixture, instrument, environment, noise) stored
with every result, the round-robin comparison with E_n / ζ / z scores and compatibility statements
against a designated reference or a Cox largest-consistent-subset consensus, repeatability–intermediate
precision–reproducibility by nested ANOVA (ISO 5725-2), Shewhart + EWMA control charts, artefact
stability, longitudinal drift detection (trend, change point, EWMA) with root-cause attribution against
B's calibration/environment/operator/procedure metadata and the site's production results, traceability
chains, evidence packs, retrospective invalidation after a failed verification, trust cards, and a
self-contained HTML dashboard.

It adds no measurements of its own: it reads the sealed archives and laboratory databases that
`labauto` laboratories produce.

## Quick start

```bash
pip install cablecheck_projectA.zip labauto_projectB.zip labplatform_projectF.zip

labplatform demo network            # four simulated sites, 10 fortnightly rounds (~2.5 min): simulate → ingest → analyse → dashboard
open network/dashboard.html         # or any browser; no external resources

# on real laboratories
labplatform ingest --db platform.sqlite /path/to/LAB1 /path/to/LAB2 ...
labplatform analyse --db platform.sqlite --artefact GOLD-RR26-0001 [--reference consensus|LAB1] --out summary.json
labplatform dashboard --db platform.sqlite --artefact GOLD-RR26-0001 --out dashboard.html
labplatform compare  --db platform.sqlite --artefact GOLD-RR26-0001 --quantity insertion_loss --pair A --f 600e6
labplatform chain    --db platform.sqlite <RUN_ID>          # the traceability chain of one result, JSON
labplatform evidence --db platform.sqlite <RUN_ID> --out evidence/
labplatform suspects --db platform.sqlite                   # runs under a calibration whose next verification failed
labplatform drift    --db platform.sqlite --artefact GOLD-RR26-0001   # drift + root cause per site
labplatform card     --db platform.sqlite <RUN_ID> --artefact GOLD-RR26-0001   # the trust card of one result
```

A trust card:

```
insertion_loss[A] @ 600 MHz = 7.918 dB   U = 0.071 dB (k = 2)
Sample ........ GOLD-RR26-0001 (artefact, 15 m, 1000base-t1-link-segment)
Measurement ... 2026-06-08T11:00 UTC, procedure roundrobin-1000base-t1-artefact, attempt 1, trust trusted, verdict PASS
Site .......... LAB4 - Site 4 - Changzhou, operator l.wang
Instrument .... labauto-sim VNA4-8500 VNA-4004, firmware A.01.20
Calibration ... CAL-LAB4-20260607 (SOLT, 2026-06-07, kit 85052D-LAB4 due 2027-04-03), verified 2026-06-08 pass dA = 0.0556 dB, RL = 41.7 dB
Fixture ....... none
Temperature ... 23.5 C (ambient (no chamber)), corrected to 23 C; ambient 23.5 C
Procedure ..... sha256 71412b7884673a1a…
Raw data ...... pairA.s4p sha256 6795bccbe89d…; manifest 06c5806d4290…; archive integrity OK
Software ...... labauto 0.1.1, cablecheck 0.1.1, numpy 2.4.4
Uncertainty ... repeat 4.7 mdB, reprod 0 mdB, cal 32 mdB, fixture 0 mdB, instrument 10 mdB, environment 10 mdB, noise 0.18 mdB
Round ......... RR-08; DoE d = +0.002164 dB, U(d) = 0.063 dB, E_n = +0.03 (reference consensus: 7.919 ± 0.031 dB)
Status ........ compatible; SUSPECT - a later verification on this instrument failed (2026-06-22)
```

## What the demonstration shows

Four sites (own analyser serial and residual-error signature, calibration quality, ambient, operator)
measure the same 15 m 1000BASE-T1 artefact three times per round with re-connection, recalibrating before
each round, plus two production samples per round. Site 4's test cable degrades from round 7 and fails
verification in the last rounds. The platform:

- puts every site's insertion loss within |E_n| ≤ 0.4 of the consensus (7.9149 ± 0.030 dB at 600 MHz) and states, per site, whether a difference is statistically distinguishable;
- compares the worst-case return loss with the ±0.9 dB it deserves against a −42 dB residual floor (an |S21|-based bound would have manufactured an E_n of −2);
- finds s_r = 4.6 mdB, s_I = 23.6 mdB (calibration to calibration), s_R = 23.6 mdB, R = 0.066 dB at 600 MHz, and counts production verdicts within 2 s_R of a limit (none);
- shows the artefact stable (−0.6 mdB / 30 d, p = 0.82);
- sees site 4's verification deviation leave its control limits at round 7 (and finds the change point there, attributed to instrument-side evidence from the check standard), three rounds before labauto's gate refuses to measure, and marks the five results of round 8 as suspect after the failed verification;
- renders it all in `dashboard.html` and writes `summary.json`.

## Layout

```
src/labplatform/
  schema.py        the canonical SQLite schema and PlatformDB
  ingest.py        labauto archives + laboratory DB → runs, results (scalar quantities from the traces), traces, verifications
  uncertainty.py   GUM budgets per result: calibration verification, trace noise, pooled repeatability, temperature correction
  stats.py         weighted mean, chi², Cox LCS, ISO 13528 Algorithm A, Cochran, Grubbs, nested ANOVA, trend
  roundrobin.py    site values, reference value, E_n / ζ / z, E_n(f), per-round comparison, artefact stability
  rnr.py           repeatability / intermediate precision / reproducibility, verdicts at risk
  control.py       individuals + EWMA charts with Western Electric rules, rational subgroups, events
  drift.py         trend / change point / EWMA per site and series, root-cause attribution against B's metadata and production results
  trust.py         trust_card(run, quantity): value, U, provenance, budget, DoE, status
  traceability.py  chain(run), suspect_runs(), evidence_pack()
  analyse.py       one call that runs everything → summary dict
  dashboard.py     self-contained HTML (matplotlib PNGs as data URIs)
  simulate.py      the four-site network built from labauto's simulators
  cli.py
tests/             16 tests: statistics against known answers; the platform end to end on a two-site network with a degrading site
docs/report.md     the report (docs/build.sh → report.pdf); docs/make_figures.py regenerates the figures (add --reuse DIR to skip the simulation)
```

## Honesty notes

All sites are simulated (project B's simulators); site differences are the simulator's. Temperature
coefficients, the 0.5 K label uncertainty and the ±0.3 % artefact length are declared constants. The
budget omits mismatch and cable-movement terms. The ANOVA uses balanced formulae on a nearly balanced
design. Statistics are the standard ones (ISO 13528, ISO 5725-2/6, Cox 2002, GUM), implemented with every
formula in the open so a laboratory can run its own data through them.

---

## Contributing

Bug reports, questions about the methods, and pull requests are all welcome — see
[CONTRIBUTING.md](CONTRIBUTING.md). Numerical changes need a numerical test, and a change to a method
is also a change to `docs/report.md`.

## Licence and attribution

MIT — see [LICENSE](LICENSE). Author: Sreeram Anil.

Part of the **[HF cable toolchain](https://github.com/anilram30/hf-cable-toolchain)** · [Report an issue](https://github.com/anilram30/labplatform/issues) ·
[Changelog](CHANGELOG.md)
