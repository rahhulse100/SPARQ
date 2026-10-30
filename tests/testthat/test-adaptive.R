test_that("adaptive assessment stops after a conclusive support interval", {
  set.seed(2)
  data <- data.frame(value = rnorm(60), x = runif(60), y = runif(60))

  fit <- sparq_assess_adaptive(
    data = data,
    analysis_function = function(x) mean(x$value),
    comparator = sparq_compare_scalar,
    retention = 0.75,
    support_threshold = 0.80,
    batch_size = 5,
    min_iterations = 10,
    max_iterations = 20,
    support_bootstrap_B = 20,
    instability_bootstrap_B = 20,
    seed = 3
  )

  expect_true(nrow(fit$adaptive_history) >= 1)
  expect_true(fit$adaptive_decision$n_attempted <= 20)
  expect_equal(nrow(fit$comparisons), fit$adaptive_decision$n_attempted)
  expect_true(fit$adaptive_decision$decision %in%
    c("above_threshold", "below_threshold", "inconclusive"))
})

test_that("adaptive assessment requires a user-specified threshold", {
  data <- data.frame(value = 1:10)
  expect_error(
    sparq_assess_adaptive(
      data = data,
      analysis_function = function(x) mean(x$value),
      comparator = sparq_compare_scalar
    ),
    "support_threshold"
  )
})

test_that("adaptive assessments provide a compact S3 summary", {
  data <- data.frame(value = 1:30)
  fit <- sparq_assess_adaptive(
    data = data,
    analysis_function = function(x) mean(x$value),
    comparator = sparq_compare_scalar,
    retention = 0.75,
    support_threshold = 0.8,
    batch_size = 2,
    min_iterations = 2,
    max_iterations = 4,
    support_bootstrap_B = 10,
    instability_bootstrap_B = 10,
    verbose = FALSE
  )

  expect_output(print(fit), "Adaptive decision")
})

test_that("adaptive and fixed retention-one assessments agree on support", {
  data <- data.frame(value = seq_len(30))
  fixed <- sparq_assess_workflow(
    data = data, analysis_function = function(x) mean(x$value),
    output_type = "scalar", retention = 1, n_iterations = 40,
    min_iterations = 20, instability_bootstrap_B = 20,
    reliability_tolerance = 0.1, required_reliability = 0.9,
    seed = 22, verbose = FALSE
  )
  adaptive <- sparq_assess_adaptive(
    data = data, analysis_function = function(x) mean(x$value),
    comparator = sparq_compare_scalar, retention = 1,
    reliability_tolerance = 0.1, required_reliability = 0.9,
    batch_size = 10, min_iterations = 20, max_iterations = 100,
    support_bootstrap_B = 20, instability_bootstrap_B = 20,
    seed = 22, verbose = FALSE
  )
  expect_equal(fixed$reliability_decision$decision, "supported")
  expect_equal(adaptive$adaptive_decision$decision, "supported")
  expect_lte(adaptive$adaptive_decision$n_attempted, 100)
})
