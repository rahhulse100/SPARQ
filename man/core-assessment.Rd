\name{sparq_run}
\alias{sparq_run}
\alias{sparq_assess}
\alias{sparq_run_with_stress_model}
\alias{sparq_stress_model}
\alias{sparq_compare_scalar}
\alias{sparq_compare_ranked}
\alias{sparq_compare_graph}
\alias{sparq_support_from_comparisons}
\alias{sparq_bound01}
\alias{sparq_bootstrap_ci}
\alias{sparq_effect_scale}
\alias{sparq_support_reason}
\title{Core perturbation assessment tools}
\description{
Low-level SPARQ functions generate perturbations, compare a reference output
with perturbed outputs, and convert the resulting comparison table into a
support-quality summary.
}
\section{Calling conventions}{\preformatted{
sparq_run(data, analysis_function, comparator, stress_model = "uniform_random",
  retention = 0.75, n_iterations = 100, ...)

sparq_assess(full_result, perturbed_results, comparator, ...)

sparq_compare_scalar(full_result, perturbed_result)
sparq_compare_ranked(full_rank, perturbed_rank, k = 100)
sparq_compare_graph(full_graph, perturbed_graph)

sparq_support_from_comparisons(comparison_table, n_iterations = nrow(comparison_table), ...)
}}
\arguments{
  \item{data}{A nonempty data frame with one row per spatial observation.}
  \item{analysis_function}{A function run on the full and perturbed data.}
  \item{comparator}{A function returning a nonnegative instability value and,
  when meaningful, a signed displacement.}
  \item{stress_model}{\code{"uniform_random"}, \code{"contiguous_hole"},
  \code{"none"}, or \code{"custom"}.}
  \item{loss_model}{Optional \code{sparq_loss_model}. When supplied, this
  explicit loss distribution determines retained spots in place of the built-in
  stress-model selection.}
  \item{retention}{Fraction of observations retained in each perturbation.}
  \item{n_iterations}{Number of requested perturbation runs.}
  \item{full_result, perturbed_results}{Previously computed reference and
  perturbed outputs for \code{sparq_assess}.}
  \item{full_rank, perturbed_rank, perturbed_result}{Ordered unique feature IDs
  or a perturbed output passed to a comparator.}
  \item{k}{Size of the ranked top-feature set.}
  \item{full_partition, perturbed_partition}{Named vectors of observation IDs
  and labels.}
  \item{full_graph, perturbed_graph}{Canonical graph edge-ID vectors.}
  \item{comparison_table}{Rows of comparator outputs, including
  \code{instability} and optionally \code{delta}.}
  \item{...}{Support parameters, coordinate columns, seed, cache controls,
  error policy, and progress controls.}
}
\value{
\code{sparq_run} and \code{sparq_assess} return an assessment list with
\code{summary}, \code{comparisons}, \code{failures}, and the reference result.
Comparators return one comparison row. \code{sparq_support_from_comparisons}
returns a one-row support-quality table.
}
\examples{
sparq_compare_scalar(1, 0.8)
sparq_compare_ranked(c("A", "B", "C"), c("B", "A", "D"), k = 2)
sparq_compare_graph(c("a--b", "b--c"), c("a--b"))
}
