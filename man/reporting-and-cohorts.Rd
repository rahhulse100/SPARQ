\name{sparq_report}
\alias{sparq_report}
\alias{sparq_localize_fragility}
\alias{sparq_assess_precomputed}
\alias{sparq_assess_precomputed_files}
\alias{sparq_read_result_file}
\alias{sparq_assess_cohort}
\alias{sparq_assess_cohort_parallel}
\alias{sparq_assess_cohort_workflow}
\alias{summary.sparq_cohort_assessment}
\alias{as.data.frame.sparq_cohort_assessment}
\title{Report, localize, and aggregate SPARQ assessments}
\description{
Create audit-ready reports, identify unstable features, spots, or graph edges,
or assess results that were produced outside SPARQ.
}
\section{Calling conventions}{\preformatted{
sparq_report(assessment, output_dir = NULL, include_plots = TRUE, verbose = TRUE)

sparq_localize_fragility(full_result, perturbed_results,
  output_type = c("scalar", "ranked", "partition", "graph"), ...)

sparq_assess_precomputed(full_results, perturbed_results, comparator, ...)
sparq_assess_precomputed_files(full_file, perturbed_file, comparator = sparq_compare_scalar, ...)
sparq_assess_cohort(sample_data, analysis_function, comparator, ...)
sparq_assess_cohort_parallel(sample_data, analysis_function, comparator, workers = 1, ...)
sparq_assess_cohort_workflow(sample_data, analysis_function, output_type,
  preset = "standard", output_dir = NULL, ...)
}}
\arguments{
  \item{assessment}{A completed SPARQ assessment.}
  \item{output_dir}{An optional new directory for TSV, RDS, and PDF report files.}
  \item{include_plots}{Whether to create the base-R summary plot.}
  \item{verbose}{Whether report progress is printed.}
  \item{full_result, perturbed_results}{Reference output and perturbation outputs.}
  \item{output_type}{Output representation used for localization.}
  \item{full_results}{Named reference outputs for multiple samples.}
  \item{comparator}{Comparator appropriate to the supplied result type.}
  \item{full_file, perturbed_file}{TSV, CSV, TXT, or RDS inputs with results.}
  \item{sample_data}{Named list of sample-level data frames. Each sample is
  assessed independently.}
  \item{analysis_function}{Function rerun on each sample or perturbation.}
  \item{workers}{Number of local workers for cohort assessment.}
  \item{preset}{Assessment preset passed to every section in a workflow cohort.}
  \item{...}{Assessment, localization, and file-reader controls.}
}
\value{
\code{sparq_report} returns a \code{sparq_report} object and writes nothing
unless \code{output_dir} is supplied. When present in the assessment, reports
retain the per-iteration spot-loss audit, the fixed-time reliability decision,
and the anytime-valid reliability confidence sequence. These are written as
separate TSV files when \code{output_dir} is supplied. Every written report
also includes a reproducibility manifest, package-version table, session
information, and, for direct spot-level assessments, the input row identifiers
plus data and identifier fingerprints. Localization returns a
\code{sparq_fragility} object. Cohort functions return one summary row per
sample plus a separate failure ledger. \code{sparq_assess_cohort_workflow}
returns a \code{sparq_cohort_assessment} with summaries, failures, a per-section
manifest, settings, and optional output paths.
}
\examples{
comparisons <- data.frame(delta = c(0.1, -0.1, 0.05),
  instability = c(0.1, 0.1, 0.05), iteration = 1:3,
  status = "ok", error_message = NA_character_)
assessment <- list(result_id = "example", full_result = 1,
  comparisons = comparisons,
  summary = sparq_support_from_comparisons(comparisons, 3,
    min_iterations = 3, instability_bootstrap_B = 20),
  failures = comparisons[0, ], settings = list(seed = 1))
sparq_report(assessment, include_plots = FALSE, verbose = FALSE)
}
