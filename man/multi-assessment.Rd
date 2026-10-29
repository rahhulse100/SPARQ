\name{sparq_assess_many}
\alias{sparq_assess_many}
\alias{sparq_report_many}
\alias{summary.sparq_multi_assessment}
\alias{print.sparq_multi_assessment}
\title{Assess multiple analyses using shared spatial perturbations}
\usage{
sparq_assess_many(data, analyses, stress_model = "uniform_random",
  retention = NULL, n_iterations = NULL, x_col = NULL, y_col = NULL,
  loss_model = NULL, seed = 1, failure_action = c("record", "stop"),
  keep_perturbed_results = FALSE, preset = "standard", verbose = NULL, ...)
}
\section{Audit bundles}{
\code{sparq_report_many()} writes a top-level summary, failure ledger, saved
assessment, report manifest, and one complete \code{sparq_report()} bundle for
each analysis. The output directory must not already exist.
}
\arguments{
  \item{data}{A nonempty spot-level data frame.}
  \item{analyses}{A named list whose elements provide \code{analysis_function}
  and \code{output_type}; optional \code{top_k} and \code{reference_scale}
  values are supported.}
  \item{stress_model, retention, n_iterations, x_col, y_col, loss_model, seed}{SPARQ perturbation controls.}
  \item{failure_action}{Whether failed method iterations are recorded or stop the assessment.}
  \item{keep_perturbed_results}{Whether to retain every method output.}
  \item{preset, verbose}{Standard SPARQ settings.}
  \item{...}{Arguments passed to \code{sparq_support_from_comparisons}.}
}
\value{A \code{sparq_multi_assessment} containing one assessment per named analysis,
shared loss audits, and an error ledger.}
\description{Runs multiple user-defined analyses on one shared reference dataset and
the same retained spot set at every perturbation iteration.}
\examples{
spots <- data.frame(value = rnorm(30), x = runif(30), y = runif(30))
fit <- sparq_assess_many(
  spots,
  analyses = list(mean_value = list(
    analysis_function = function(x) mean(x$value), output_type = "scalar"
  )),
  preset = "quick", n_iterations = 5, min_iterations = 3,
  instability_bootstrap_B = 20, verbose = FALSE
)
summary(fit)
}
