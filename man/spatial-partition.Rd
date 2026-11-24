\name{sparq_compare_spatial_partition}
\alias{sparq_compare_spatial_partition}
\title{Compare spatial partitions using boundaries and domain adjacency}
\usage{
sparq_compare_spatial_partition(
  full_labels,
  perturbed_labels,
  edges,
  weights = c(ari = 1, boundary_f1 = 1, adjacency_jaccard = 1)
)
}
\arguments{
  \item{full_labels}{Named reference partition labels.}
  \item{perturbed_labels}{Named native-rerun partition labels.}
  \item{edges}{A two-column spot-neighborhood edge table or matrix.}
  \item{weights}{Named nonnegative weights for ARI, tolerant boundary F1, and domain-adjacency Jaccard. They are normalized internally.}
}
\value{A one-row data frame with ARI, tolerant boundary precision, recall and F1, domain-adjacency Jaccard, the normalized component weights, and combined instability.}
\description{
Provides an optional output-specific comparator for spatial domain analyses. Boundary agreement permits one neighborhood hop of displacement. Domain adjacency is calculated after a one-to-one maximum-overlap matching of native-rerun labels to reference labels. Unmatched native-rerun domains remain distinct so that domain splitting creates a detectable new adjacency. No pathology labels or clean-reference results are used.
}
\details{
The supplied graph is induced on spot IDs shared by both partitions. The default combined instability gives ARI, boundary F1, and adjacency Jaccard equal weight. A study should declare any alternative weights before held-out evaluation.
}
\examples{
edges <- data.frame(source = c("a", "b", "c"), target = c("b", "c", "d"))
full <- c(a = "A", b = "A", c = "B", d = "B")
rerun <- c(a = "x", b = "x", c = "y", d = "y")
sparq_compare_spatial_partition(full, rerun, edges)
}
