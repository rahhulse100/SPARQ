# SPARQ news

## Development

- Added `sparq_assess_procurement()` and `sparq_report_procurement()` for one
  shared-reference assessment across explicit uniform, local-hole, edge, and
  quality-weighted support-loss scenarios. The workflow records scenario-level
  support, loss audits, quality summaries, and per-scenario report bundles;
  it explicitly does not perform batch correction or imputation.
- Reports now write a reproducibility manifest, loaded-package versions,
  session information, settings, input data/row-ID fingerprints, and direct
  input row identifiers when an output directory is supplied.

## 0.2.0 — 2026-10-04

- Added workflow, script, spatial-object, and adaptive assessment interfaces.
- Added `quick`, `standard`, and `thorough` presets.
- Added conventional assessment `summary()`, `plot()`, and `as.data.frame()` methods.
- Added a cohort workflow runner with aggregate progress, failure handling, and audit manifest export.
- Added R help pages, a citation record, continuous package checks, and an updated project README.
- Added direct tolerance-based reliability certificates and an anytime-valid
  beta-binomial mixture confidence sequence for adaptive stopping.
- Added explicit quality-weighted, edge-weighted, contiguous, uniform, and
  custom spatial loss mechanisms with retained-spot audit tables.
- Added discrete reliability-frontier estimation and empirical-Bayes
  beta-binomial cohort reliability inference.

## 0.0.3.9000

- Development version used for the initial support-quality, known-truth benchmark,
  reliability, and ablation implementations.
