\name{sparq_assess_procurement}
\alias{sparq_assess_procurement}
\alias{sparq_report_procurement}
\alias{print.sparq_procurement_assessment}
\alias{summary.sparq_procurement_assessment}
\alias{as.data.frame.sparq_procurement_assessment}
\alias{plot.sparq_procurement_assessment}
\title{Assess sensitivity to declared procurement-related spot-loss scenarios}
\description{
Runs one reference analysis and evaluates the same result under declared,
plausible support-loss scenarios. This interface is intended to make sensitivity
to spatial coverage, local tissue loss, edge loss, and spot-level quality
variation transparent. It does not correct batch effects, impute unobserved
tissue, or recover biological information absent from a specimen.
}
\usage{
sparq_assess_procurement(
  data, analysis_function,
  output_type = c("scalar", "ranked", "partition", "graph"),
  x_col = NULL, y_col = NULL, quality_col = NULL, scenarios = NULL,
  retention = NULL, n_iterations = NULL, preset = "standard",
  quality_direction = c("higher_quality_retained", "lower_quality_retained"),
  quality_strength = 1, edge_strength = 1, top_k = 100,
  result_id = "result", reference_scale = 1, seed = 1,
  preflight = TRUE, verbose = NULL, ...
)

sparq_report_procurement(assessment, output_dir = NULL,
  include_plots = TRUE, verbose = TRUE)
}
\arguments{
  \item{data}{A nonempty spot-level data frame.}
  \item{analysis_function}{A function that reruns the selected analysis on a
  complete or perturbed spot-level data frame.}
  \item{output_type}{One of \code{"scalar"}, \code{"ranked"},
  \code{"partition"}, or \code{"graph"}.}
  \item{x_col, y_col}{Optional finite numeric coordinate columns. They are
  required for \code{"contiguous_hole"} and \code{"edge_weighted"} scenarios.}
  \item{quality_col}{Optional finite numeric spot-level quality metric, such as
  library size or detected-gene count. It is required for
  \code{"quality_weighted"}.}
  \item{scenarios}{Unique values among \code{"uniform_random"},
  \code{"contiguous_hole"}, \code{"edge_weighted"}, and
  \code{"quality_weighted"}. When \code{NULL}, SPARQ evaluates every scenario
  supported by the supplied columns.}
  \item{retention, n_iterations, preset}{Perturbation settings passed to
  \code{sparq_assess_workflow}.}
  \item{quality_direction}{Whether higher or lower quality values receive
  greater retention probability in quality-weighted loss.}
  \item{quality_strength, edge_strength}{Nonnegative loss-model strengths.}
  \item{top_k}{Ranking depth for a ranked result.}
  \item{result_id}{Identifier for the result or section.}
  \item{reference_scale}{Comparator reference scale.}
  \item{seed}{Master random seed. Scenario-specific seeds are deterministically
  derived from this seed.}
  \item{preflight}{Whether to run one preflight reference analysis and small
  spot-removal rerun before the scenario assessments. The preflight reference
  is reused by all scenarios.}
  \item{verbose}{Whether to print progress.}
  \item{...}{Additional arguments to \code{sparq_assess_workflow}, including
  \code{reliability_tolerance}, \code{required_reliability}, or
  \code{cache_dir}.}
  \item{assessment}{A \code{sparq_procurement_assessment}.}
  \item{output_dir}{Optional new directory in which to write an audit bundle.}
  \item{include_plots}{Whether per-scenario reports include PDF plots.}
}
\value{
\code{sparq_assess_procurement} returns a
\code{sparq_procurement_assessment}: one shared reference result,
scenario-specific SPARQ assessments, a scenario summary, retained-spot audits,
a quality profile when available, settings, and a non-corrective interpretation.
\code{summary()} and \code{as.data.frame()} return the scenario summary.
\code{sparq_report_procurement} returns in-memory scenario reports and, when
\code{output_dir} is supplied, writes the shared reference, scenario summary,
quality profile, settings, full assessment, and one complete report bundle per
scenario. Every written bundle includes a reproducibility manifest,
package-version table, session information, and input row-ID/data fingerprints.
}
\details{
Every conclusion is conditional on the declared scenarios. A quality metric is
used only to define a sensitivity analysis; it is not treated as evidence of a
true missingness mechanism. The function never adjusts the original result.
If a batch-correction workflow is appropriate, it is an upstream analysis step;
this function can then assess the stability of that downstream result.
}
\examples{
spots <- data.frame(
  value = rnorm(60), x = runif(60), y = runif(60),
  library_size = rexp(60, rate = 1 / 1000)
)

assessment <- sparq_assess_procurement(
  data = spots,
  analysis_function = function(x) mean(x$value),
  output_type = "scalar",
  x_col = "x", y_col = "y", quality_col = "library_size",
  preset = "quick", n_iterations = 10, min_iterations = 5,
  instability_bootstrap_B = 20, verbose = FALSE
)

summary(assessment)
}
