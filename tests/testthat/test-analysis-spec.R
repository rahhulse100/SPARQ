test_that("output-specific analysis specifications record explicit rules", {
  partition <- sparq_partition_analysis(
    analysis_function = function(x) setNames(ifelse(x$value > 0, "high", "low"), x$id),
    min_ari = 0.60
  )
  ranked <- sparq_ranked_analysis(
    analysis_function = function(x) c("g1", "g2", "g3"),
    top_k = 2,
    min_jaccard = 0.70
  )
  graph <- sparq_graph_analysis(
    analysis_function = function(x) c("a--b", "b--c"),
    min_edge_jaccard = 0.80
  )
  scalar <- sparq_scalar_analysis(
    analysis_function = function(x) mean(x$value),
    max_normalized_change = 0.15,
    reference_scale = 2
  )

  expect_s3_class(partition, "sparq_analysis_spec")
  expect_equal(partition$output_type, "partition")
  expect_equal(partition$decision_rule$reliability_tolerance, 0.40)
  expect_equal(ranked$top_k, 2L)
  expect_equal(ranked$decision_rule$reliability_tolerance, 0.30)
  expect_equal(graph$decision_rule$reliability_tolerance, 0.20)
  expect_equal(scalar$decision_rule$reliability_tolerance, 0.15)
  expect_equal(scalar$reference_scale, 2)
  expect_output(print(partition), "adjusted_Rand_index")
})

test_that("analysis specifications require explicit valid preservation rules", {
  expect_error(
    sparq_partition_analysis(function(x) x, min_ari = 1.1),
    "min_ari"
  )
  expect_error(
    sparq_ranked_analysis(function(x) x, top_k = 0, min_jaccard = 0.6),
    "top_k"
  )
  expect_error(
    sparq_scalar_analysis(function(x) 1, max_normalized_change = -0.1),
    "max_normalized_change"
  )
  expect_error(
    sparq_custom_analysis(function(x) x, comparator = 1, reliability_tolerance = 0.1),
    "comparator"
  )
})

test_that("sparq_assess_many applies each analysis specification separately", {
  set.seed(19)
  spots <- data.frame(
    id = paste0("spot_", seq_len(30)),
    value = rnorm(30)
  )

  analyses <- list(
    domains = sparq_partition_analysis(
      function(x) setNames(ifelse(x$value > 0, "high", "low"), x$id),
      min_ari = 0.60
    ),
    genes = sparq_ranked_analysis(
      function(x) c("g1", "g2", "g3"),
      top_k = 2,
      min_jaccard = 0.60
    ),
    graph = sparq_graph_analysis(
      function(x) c("a--b", "b--c"),
      min_edge_jaccard = 0.60
    ),
    score = sparq_scalar_analysis(
      function(x) mean(x$value),
      max_normalized_change = 1,
      reference_scale = 1
    )
  )

  fit <- sparq_assess_many(
    data = spots,
    analyses = analyses,
    retention = 0.75,
    n_iterations = 4,
    min_iterations = 2,
    instability_bootstrap_B = 20,
    verbose = FALSE,
    seed = 42
  )

  expect_s3_class(fit, "sparq_multi_assessment")
  expect_equal(names(fit$assessments), names(analyses))
  expect_true(all(vapply(fit$assessments, function(x) {
    is.data.frame(x$reliability_decision) && nrow(x$reliability_decision) == 1L
  }, logical(1))))
  expect_true(all(vapply(fit$assessments, function(x) {
    !is.null(x$settings$decision_rule)
  }, logical(1))))
  expect_true(all(c("decision_rule", "reliability_decision") %in%
    names(summary(fit))))
})
