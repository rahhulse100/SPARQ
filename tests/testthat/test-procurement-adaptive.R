test_that("adaptive procurement keeps scenario decisions independent", {
  set.seed(8)
  data <- data.frame(
    value = rnorm(40),
    x = runif(40),
    y = runif(40),
    quality = rexp(40)
  )

  fit <- sparq_assess_procurement_adaptive(
    data = data,
    analysis_function = function(x) mean(x$value),
    output_type = "scalar",
    x_col = "x", y_col = "y", quality_col = "quality",
    retention = 0.75,
    reliability_tolerance = 1,
    required_reliability = 0.5,
    reliability_confidence = 0.8,
    batch_size = 2,
    min_iterations = 2,
    max_iterations = 4,
    instability_bootstrap_B = 20,
    support_bootstrap_B = 20,
    verbose = FALSE,
    seed = 8
  )

  expect_s3_class(fit, "sparq_procurement_adaptive_assessment")
  expect_equal(names(fit$assessments), c(
    "uniform_random", "contiguous_hole", "edge_weighted", "quality_weighted"
  ))
  expect_equal(nrow(fit$scenario_summary), 4)
  expect_true(all(vapply(fit$assessments, function(x) {
    is.data.frame(x$adaptive_decision) && nrow(x$adaptive_decision) == 1L
  }, logical(1))))
})
