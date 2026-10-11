<img src="docs/assets/sparq-hero.png" align="right" width="150" alt="SPARQ spatial-transcriptomics perturbation illustration">

# SPARQ

**Spatial Perturbation Analysis for Reliability Quantification**

[![R-CMD-check](https://github.com/rahhulse100/SPARQ1/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/rahhulse100/SPARQ1/actions/workflows/R-CMD-check.yaml)
[![R version](https://img.shields.io/badge/R-%3E%3D%204.1-276DC3)](https://www.r-project.org/)

<br clear="right">

SPARQ asks whether an analysis result remains consistent when plausible spot loss is introduced into a spatial-transcriptomic section. It does not denoise the data, replace an analysis method, or establish biological truth. It quantifies how sensitive a specified result is to a defined perturbation model.

SPARQ supports four built-in output types:

- **scalar**: one numeric estimate, score, or statistic;
- **ranked**: an ordered feature, pathway, or interaction list;
- **partition**: spot-level domain or cluster labels;
- **graph**: canonical edge IDs describing a network or neighborhood graph.

## Preflight an analysis before the full run

Before launching an expensive assessment, use one reference run plus one small,
deterministic spot-removal run to verify that the analysis survives missing
spots. The preflight checks input size, coordinate columns, output type, and,
for partitions, exact alignment between returned and retained spot IDs.

```r
sparq_validate_workflow(
  data = spot_data,
  analysis_function = function(data) mean(data$signature_score),
  output_type = "scalar"
)
```

## Installation

Install the development package from GitHub:

```r
install.packages("remotes")
remotes::install_github("rahhulse100/SPARQ1")
library(SPARQ)
```

## Standard workflow

Give SPARQ a spot-level data frame and a function that reruns the analysis after rows have been removed. For a scalar result:

```r
set.seed(1)
spot_data <- data.frame(
  signature_score = rnorm(100),
  x = runif(100),
  y = runif(100)
)

assessment <- sparq_assess_workflow(
  data = spot_data,
  analysis_function = function(data) mean(data$signature_score),
  output_type = "scalar",
  stress_model = "uniform_random",
  retention = 0.75,
  n_iterations = 100,
  result_id = "section_001",
  reference_scale = 1,
  seed = 1
)

assessment$summary
assessment$comparisons
```

## Presets and familiar result methods

Start with a named protocol instead of choosing every setting. SPARQ prints the
resolved settings at the start of an assessment and records them in the output.
Explicit arguments always override one preset value.

```r
assessment <- sparq_assess_workflow(
  data = spot_data,
  analysis_function = function(data) mean(data$signature_score),
  output_type = "scalar",
  preset = "standard"
)

summary(assessment)
as.data.frame(assessment)
plot(assessment)
```

| Preset | Perturbations | Minimum completed | Instability bootstrap draws |
| --- | ---: | ---: | ---: |
| `quick` | 25 | 15 | 200 |
| `standard` | 100 | 50 | 2,000 |
| `thorough` | 500 | 100 | 10,000 |

All presets retain 75% of observations by default. Use `sparq_preset()` to
inspect the complete setting list, or override any setting directly.

## Tolerance-based reliability decisions

Support quality remains a continuous score that can be calibrated on known
truth. Separately, SPARQ can make a direct three-way decision from the observed
perturbations when the analyst prespecifies an acceptable normalized change.
The decision is `supported`, `not_supported`, or `inconclusive`; failed reruns
count conservatively as not preserving the result.

```r
assessment <- sparq_assess_workflow(
  data = spot_data,
  analysis_function = function(data) mean(data$signature_score),
  output_type = "scalar",
  reliability_tolerance = 0.10,
  required_reliability = 0.90
)

assessment$reliability_decision
```

SPARQ prints a compact execution log by default: the reference analysis,
perturbation progress, and the final support result. Printing the assessment
itself gives the same decision-relevant summary without requiring users to
search through a wide table. Set `verbose = FALSE` for notebooks, batch jobs,
or programmatic runs where console output is not wanted.

```r
assessment

# <SPARQ workflow assessment: section_001>
#   Output type: scalar
#   Support quality: 0.91 (adequate support)
#   Completed perturbations: 100/100; failed: 0
```

For domain analysis, the supplied analysis function must return a named vector of spot IDs and domain labels. For ranked and graph analyses, it must return an ordered unique feature vector or canonical edge-ID vector, respectively.

## Reusing an existing R script

`sparq_assess_script()` allows a user to reuse a script rather than writing an R function. SPARQ injects the current perturbed data under `input_object`; the script must create `result_object`.

```r
assessment <- sparq_assess_script(
  data = spot_data,
  script = "banksy_analysis.R",
  input_object = "seurat_obj",
  result_object = "banksy_domains",
  output_type = "partition",
  retention = 0.75,
  n_iterations = 100
)
```

The script must use the injected object rather than repeatedly reading an unchanged input file.

## Spatial-object adapters

`sparq_assess_spatial()` accepts a Seurat object, SpatialExperiment object,
Giotto object, or a coordinate-bearing data frame. SPARQ retains the original
object and subsets it for every perturbation, so the supplied analysis function
receives the object type it already expects.

```r
assessment <- sparq_assess_spatial(
  object = seurat_spatial_object,
  analysis_function = function(x) {
    # Rerun BANKSY or another domain method on x.
    # Return named spot labels, e.g. setNames(new_domains, colnames(x)).
  },
  output_type = "partition",
  retention = 0.75,
  n_iterations = 100
)
```

For custom spatial classes, provide an adapter with an `extract()` function
that returns spot IDs and named coordinates, plus a `subset()` function that
returns the object restricted to retained spot IDs. See
`?sparq_prepare_spatial` for the exact contract.

## Adaptive stopping

For expensive analyses, `sparq_assess_adaptive()` evaluates the support score in batches and stops when a bootstrap interval lies entirely above or below a user-specified decision threshold. It does not infer reliability from a single baseline analysis.

```r
adaptive_assessment <- sparq_assess_adaptive(
  data = spot_data,
  analysis_function = function(data) mean(data$signature_score),
  comparator = sparq_compare_scalar,
  retention = 0.75,
  support_threshold = 0.80,
  batch_size = 10,
  min_iterations = 50,
  max_iterations = 100
)

adaptive_assessment$adaptive_history
adaptive_assessment$adaptive_decision
```

For tolerance-based stopping, use `reliability_tolerance` instead of
`support_threshold`. This avoids the support-score bootstrap at each adaptive
check and stops as soon as the sequential reliability interval supports a final
decision.

### Sequentially valid stopping

Adaptive reliability stopping uses an anytime-valid beta-binomial mixture
confidence sequence, rather than repeatedly checking an ordinary fixed-time
interval. Its bounds remain valid across the accumulating perturbation stream.
This means the stopping decision remains tied to the stated confidence level
when SPARQ stops early.

## Spatial loss mechanisms

Every SPARQ claim is conditional on a stated loss mechanism. In addition to
uniform loss and a contiguous spatial hole, define quality-weighted or
edge-weighted loss explicitly. The retained spot set and its per-spot weights
are audited for every application.

```r
quality_loss <- sparq_define_loss_model(
  type = "quality_weighted",
  quality_col = "library_size",
  quality_direction = "higher_quality_retained",
  quality_strength = 1
)

assessment <- sparq_assess_workflow(
  data = spot_data,
  analysis_function = function(data) mean(data$signature_score),
  output_type = "scalar",
  loss_model = quality_loss,
  reliability_tolerance = 0.10,
  required_reliability = 0.90
)
```

`edge_weighted` loss requires coordinates and preferentially removes spots
farther from the coordinate-wise tissue center. These are specified stress
distributions, not assertions that one mechanism is the true missingness model.

## Procurement-aware sensitivity analysis

When a specimen cannot be recollected, the relevant question is not whether
SPARQ can reconstruct tissue that was never captured. It cannot. Instead,
`sparq_assess_procurement()` runs one shared reference analysis and evaluates
the same conclusion under every declared, data-supported loss scenario:
uniform loss, a contiguous missing region, edge-biased loss, and loss weighted
by a user-supplied spot-level quality metric. The output names the scenario
under which a conclusion becomes fragile and preserves every loss audit.

```r
spot_data$library_size <- rpois(nrow(spot_data), lambda = 1500) + 1L

procurement <- sparq_assess_procurement(
  data = spot_data,
  analysis_function = function(data) mean(data$signature_score),
  output_type = "scalar",
  x_col = "x",
  y_col = "y",
  quality_col = "library_size",
  reliability_tolerance = 0.10,
  required_reliability = 0.90
)

summary(procurement)
plot(procurement)

sparq_report_procurement(
  procurement,
  output_dir = "sparq_procurement_report"
)
```

This is a sensitivity analysis, not batch correction, denoising, imputation,
or a claim that any single loss model is the true procurement process. Batch
correction, when appropriate, is an upstream analysis choice; SPARQ can assess
the stability of its downstream result.

## Reliability frontiers

One retention value answers whether a conclusion survives one sampling regime.
A frontier answers how much sampling can be lost before the conclusion is no
longer certified on a prespecified grid. SPARQ reports the lowest *evaluated*
retention that remains supported at every higher evaluated retention; it never
interpolates a critical value between grid points.

```r
frontier <- sparq_assess_frontier(
  data = spot_data,
  analysis_function = function(data) mean(data$signature_score),
  output_type = "scalar",
  retentions = c(0.50, 0.60, 0.75, 0.90, 1.00),
  reliability_tolerance = 0.10,
  required_reliability = 0.90
)

frontier$frontier
frontier$critical_retention
plot(frontier)
```

## Cohort analyses

For many independent sections, use `sparq_assess_cohort_workflow()` rather than
printing every perturbation from every section. It reports aggregate completion
progress and an estimated remaining time. With `output_dir`, it exports one
manifest containing section IDs, status, failures, elapsed time, seeds, and
support summaries, plus the cohort settings.

```r
cohort <- sparq_assess_cohort_workflow(
  sample_data = sections,
  analysis_function = run_analysis,
  output_type = "partition",
  preset = "standard",
  output_dir = "sparq_cohort_output"
)

cohort$manifest
cohort$failures
```

When every section uses the same tolerance and reliability requirement, cohort
assessment also returns `cohort$reliability`. It uses an empirical-Bayes
beta-binomial model to distinguish finite perturbation uncertainty within each
section from reliability variation between sections. It does not treat spots as
biological replicates.

```r
cohort$reliability$section_reliability
cohort$reliability$cohort_summary
```

## Reports and audit bundles

`sparq_report()` returns an inspectable report object. With an `output_dir`, it
writes a self-contained audit bundle: the support summary, perturbation table,
failure ledger, settings, reference result, any retained perturbation results,
fragility tables, and a PDF overview. SPARQ never creates an output directory
unless one is explicitly supplied.

```r
report <- sparq_report(
  assessment,
  output_dir = "sparq_section_001_report"
)

report
report$output_files
```

## Precomputed outputs

If reference and perturbed outputs have already been generated elsewhere, use `sparq_assess_precomputed()` or `sparq_assess_precomputed_files()`. A single baseline output alone cannot be assessed for perturbation robustness.

## Interpretation

SPARQ returns a continuous support-quality score, instability summaries, failed-iteration records, and when perturbed outputs are retained—fragility localization. A high score supports consistency under the specified perturbation model; it does not prove a result is biologically correct or unchanged under every possible source of variation.

## Development status

SPARQ is an installable R source package. The current development focus is broader workflow adapters, adaptive stopping validation for every output type, documentation, and external biological validation.
