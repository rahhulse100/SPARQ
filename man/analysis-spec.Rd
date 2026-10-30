\name{sparq_analysis_spec}
\alias{sparq_analysis_spec}
\alias{sparq_partition_analysis}
\alias{sparq_ranked_analysis}
\alias{sparq_graph_analysis}
\alias{sparq_scalar_analysis}
\alias{sparq_custom_analysis}
\title{Declare an output-specific SPARQ analysis contract}
\usage{
sparq_analysis_spec(analysis_function, output_type, comparator,
  reliability_tolerance, required_reliability = 0.90,
  reliability_confidence = 0.95, reference_scale = 1,
  top_k = NULL, rule_name = "custom")

sparq_partition_analysis(analysis_function, min_ari,
  required_reliability = 0.90, reliability_confidence = 0.95)

sparq_ranked_analysis(analysis_function, top_k, min_jaccard,
  required_reliability = 0.90, reliability_confidence = 0.95)

sparq_graph_analysis(analysis_function, min_edge_jaccard,
  required_reliability = 0.90, reliability_confidence = 0.95)

sparq_scalar_analysis(analysis_function, max_normalized_change,
  reference_scale = 1, required_reliability = 0.90,
  reliability_confidence = 0.95)

sparq_custom_analysis(analysis_function, comparator, reliability_tolerance,
  reference_scale = 1, required_reliability = 0.90,
  reliability_confidence = 0.95,
  rule_name = "custom instability rule")
}
\arguments{
  \item{analysis_function}{Function accepting retained input data and returning one result.}
  \item{output_type}{One of \code{scalar}, \code{ranked}, \code{partition}, \code{graph}, or \code{custom}.}
  \item{comparator}{Function comparing the baseline and perturbed results.}
  \item{reliability_tolerance}{Largest normalized instability that counts as preservation.}
  \item{required_reliability}{Required probability of preservation.}
  \item{reliability_confidence}{Confidence level for the reliability decision.}
  \item{reference_scale}{Positive scale used for scalar or custom instability.}
  \item{top_k}{Number of leading ranked features to compare.}
  \item{rule_name}{Human-readable description of a custom preservation rule.}
  \item{min_ari}{Minimum adjusted Rand index that counts as preservation.}
  \item{min_jaccard}{Minimum top-k Jaccard overlap that counts as preservation.}
  \item{min_edge_jaccard}{Minimum edge Jaccard overlap that counts as preservation.}
  \item{max_normalized_change}{Largest normalized scalar change that counts as preservation.}
}
\value{An object of class \code{sparq_analysis_spec}, accepted by \code{sparq_assess_many}.}
\description{
Constructors that pair a native rerun function with an output-specific
comparator and an explicit preservation rule. They do not select a loss model
or run an assessment.
}
\details{
The built-in result classes are partitions, ranked results, graphs, and scalar
estimates. Use \code{sparq_custom_analysis} for a different result type; its
comparator must return a non-negative \code{instability} column.
}
\examples{
domain_spec <- sparq_partition_analysis(
  function(x) setNames(ifelse(x$value > 0, "high", "low"), x$id),
  min_ari = 0.60
)
}
