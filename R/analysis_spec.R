#' Declare an output-specific SPARQ analysis contract
#'
#' These constructors make the output contract explicit before an assessment is
#' launched.  They do not run an analysis or choose a loss model.  Instead,
#' they pair a user-supplied rerun function with the appropriate comparator and
#' a documented preservation rule.  The resulting object can be supplied to
#' [sparq_assess_many()].
#'
#' SPARQ provides constructors for four common result classes: partitions,
#' ranked results, graphs, and scalar estimates.  Other result classes can use
#' [sparq_custom_analysis()] with a comparator that returns a non-negative
#' `instability` column.
#'
#' @param analysis_function Function accepting the retained input data and
#'   returning one analysis result.
#' @param output_type Output type: `"partition"`, `"ranked"`, `"graph"`, or
#'   `"scalar"`.
#' @param comparator Function comparing the baseline and perturbed results.
#' @param reliability_tolerance Largest normalized instability that counts as
#'   preservation.
#' @param required_reliability Required probability of preservation.
#' @param reliability_confidence Confidence level for a reliability decision.
#' @param reference_scale Positive scale for scalar or custom instability.
#' @param top_k Number of leading ranked features to compare.
#' @param min_ari Minimum adjusted Rand index that counts as preservation.
#' @param min_jaccard Minimum top-k Jaccard overlap that counts as preservation.
#' @param min_edge_jaccard Minimum edge Jaccard overlap that counts as
#'   preservation.
#' @param max_normalized_change Largest normalized scalar change that counts as
#'   preservation.
#'
#' @return An object of class `sparq_analysis_spec`.
#' @export
sparq_analysis_spec <- function(
    analysis_function, output_type,
    comparator, reliability_tolerance,
    required_reliability = 0.90,
    reliability_confidence = 0.95,
    reference_scale = 1,
    top_k = NULL,
    rule_name = "custom") {

  if (!is.function(analysis_function)) {
    stop("analysis_function must be a function.", call. = FALSE)
  }
  if (!is.character(output_type) || length(output_type) != 1L ||
      is.na(output_type) || !output_type %in% c("scalar", "ranked", "partition", "graph", "custom")) {
    stop("output_type must be scalar, ranked, partition, graph, or custom.", call. = FALSE)
  }
  if (!is.function(comparator)) {
    stop("comparator must be a function.", call. = FALSE)
  }
  sparq_check_number(reliability_tolerance, "reliability_tolerance", 0)
  sparq_check_number(required_reliability, "required_reliability",
    .Machine$double.eps, 1 - .Machine$double.eps)
  sparq_check_number(reliability_confidence, "reliability_confidence",
    .Machine$double.eps, 1 - .Machine$double.eps)
  sparq_check_number(reference_scale, "reference_scale")
  if (!is.null(top_k)) {
    sparq_check_number(top_k, "top_k", 1, integer = TRUE)
  }
  if (!is.character(rule_name) || length(rule_name) != 1L ||
      is.na(rule_name) || !nzchar(rule_name)) {
    stop("rule_name must be one nonempty character value.", call. = FALSE)
  }

  out <- list(
    analysis_function = analysis_function,
    output_type = output_type,
    comparator = comparator,
    top_k = top_k,
    reference_scale = reference_scale,
    decision_rule = list(
      name = rule_name,
      reliability_tolerance = reliability_tolerance,
      required_reliability = required_reliability,
      reliability_confidence = reliability_confidence,
      reference_scale = reference_scale
    )
  )
  class(out) <- "sparq_analysis_spec"
  out
}

#' @rdname sparq_analysis_spec
#' @export
sparq_partition_analysis <- function(
    analysis_function, min_ari,
    required_reliability = 0.90,
    reliability_confidence = 0.95) {
  sparq_check_number(min_ari, "min_ari", 0, 1)
  sparq_analysis_spec(
    analysis_function = analysis_function,
    output_type = "partition",
    comparator = sparq_compare_partition,
    reliability_tolerance = 1 - min_ari,
    required_reliability = required_reliability,
    reliability_confidence = reliability_confidence,
    rule_name = paste0("adjusted_Rand_index >= ", format(min_ari, trim = TRUE))
  )
}

#' @rdname sparq_analysis_spec
#' @export
sparq_ranked_analysis <- function(
    analysis_function, top_k, min_jaccard,
    required_reliability = 0.90,
    reliability_confidence = 0.95) {
  sparq_check_number(top_k, "top_k", 1, integer = TRUE)
  sparq_check_number(min_jaccard, "min_jaccard", 0, 1)
  comparator <- function(full_result, perturbed_result) {
    sparq_compare_ranked(full_result, perturbed_result, k = top_k)
  }
  sparq_analysis_spec(
    analysis_function = analysis_function,
    output_type = "ranked",
    comparator = comparator,
    reliability_tolerance = 1 - min_jaccard,
    required_reliability = required_reliability,
    reliability_confidence = reliability_confidence,
    top_k = top_k,
    rule_name = paste0("top_", top_k, "_Jaccard >= ", format(min_jaccard, trim = TRUE))
  )
}

#' @rdname sparq_analysis_spec
#' @export
sparq_graph_analysis <- function(
    analysis_function, min_edge_jaccard,
    required_reliability = 0.90,
    reliability_confidence = 0.95) {
  sparq_check_number(min_edge_jaccard, "min_edge_jaccard", 0, 1)
  sparq_analysis_spec(
    analysis_function = analysis_function,
    output_type = "graph",
    comparator = sparq_compare_graph,
    reliability_tolerance = 1 - min_edge_jaccard,
    required_reliability = required_reliability,
    reliability_confidence = reliability_confidence,
    rule_name = paste0("edge_Jaccard >= ", format(min_edge_jaccard, trim = TRUE))
  )
}

#' @rdname sparq_analysis_spec
#' @export
sparq_scalar_analysis <- function(
    analysis_function, max_normalized_change,
    reference_scale = 1,
    required_reliability = 0.90,
    reliability_confidence = 0.95) {
  sparq_check_number(max_normalized_change, "max_normalized_change", 0)
  sparq_check_number(reference_scale, "reference_scale")
  sparq_analysis_spec(
    analysis_function = analysis_function,
    output_type = "scalar",
    comparator = sparq_compare_scalar,
    reliability_tolerance = max_normalized_change,
    required_reliability = required_reliability,
    reliability_confidence = reliability_confidence,
    reference_scale = reference_scale,
    rule_name = paste0("normalized_absolute_change <= ",
      format(max_normalized_change, trim = TRUE))
  )
}

#' @rdname sparq_analysis_spec
#' @export
sparq_custom_analysis <- function(
    analysis_function, comparator, reliability_tolerance,
    reference_scale = 1,
    required_reliability = 0.90,
    reliability_confidence = 0.95,
    rule_name = "custom instability rule") {
  sparq_analysis_spec(
    analysis_function = analysis_function,
    output_type = "custom",
    comparator = comparator,
    reliability_tolerance = reliability_tolerance,
    required_reliability = required_reliability,
    reliability_confidence = reliability_confidence,
    reference_scale = reference_scale,
    rule_name = rule_name
  )
}

#' @export
print.sparq_analysis_spec <- function(x, ...) {
  cat("<SPARQ analysis specification>\n")
  cat("  Output type: ", x$output_type, "\n", sep = "")
  cat("  Preservation rule: ", x$decision_rule$name, "\n", sep = "")
  cat("  Required reliability: ", x$decision_rule$required_reliability,
    " at confidence ", x$decision_rule$reliability_confidence, "\n", sep = "")
  invisible(x)
}
