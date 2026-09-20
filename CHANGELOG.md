# Changelog

All notable changes to this project are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] — 2026-09-18

### Added
- First release: the multi-site trust layer.
- Canonical schema and ingest from `labauto` archives and laboratory databases, including jobs that
  were refused at a gate — a network that cannot see its refusal rate cannot manage it.
- GUM uncertainty budgets stored with every result: repeatability, reproducibility, calibration,
  fixture, instrument, environment, noise and resolution, with the temperature correction to 23 °C.
- Round-robin comparison with E_n, ζ and z scores against a designated reference or a Cox
  largest-consistent-subset consensus, per frequency and per round, with compatibility statements.
- Repeatability, intermediate precision and reproducibility by nested ANOVA (ISO 5725-2), with
  Cochran and Grubbs screening and verdicts-at-risk.
- Shewhart and EWMA control charts with Western Electric rules on rational subgroups; artefact
  stability testing.
- Longitudinal drift detection (trend, change point, EWMA) with root-cause attribution against
  calibration, environment, operator and procedure metadata.
- Traceability chain, evidence packs, suspect-run identification, trust cards and a self-contained
  dashboard.
- Four-site network simulation, CLI, 16 tests and an 11-page technical report.
