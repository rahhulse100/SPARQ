test_that("multi-analysis assessment shares each retained spot set", {
  set.seed(41)
  data <- data.frame(value = rnorm(40), x = runif(40), y = runif(40))
  seen <- list(left = list(), right = list())
  analyses <- list(
    left = list(
      output_type = "scalar",
      analysis_function = function(x) {
        seen$left[[length(seen$left) + 1L]] <<- sort(x$value)
        mean(x$value)
      }
    ),
    right = list(
      output_type = "scalar",
      analysis_function = function(x) {
        seen$right[[length(seen$right) + 1L]] <<- sort(x$value)
        stats::median(x$value)
      }
    )
  )
  fit <- sparq_assess_many(data, analyses, retention = 0.75,
    n_iterations = 5, min_iterations = 3, instability_bootstrap_B = 20,
    seed = 42, verbose = FALSE)
  expect_s3_class(fit, "sparq_multi_assessment")
  expect_equal(names(fit$assessments), c("left", "right"))
  expect_equal(nrow(summary(fit)), 2)
  expect_equal(length(fit$loss_audits), 5)
  expect_equal(length(seen$left), 6)
  expect_equal(length(seen$right), 6)
  expect_equal(seen$left[-1], seen$right[-1])
})

test_that("multi-analysis report creates top-level and per-analysis audits", {
  data <- data.frame(value = seq_len(20), x = seq_len(20), y = seq_len(20))
  fit <- sparq_assess_many(data, list(mean_value = list(
    output_type = "scalar", analysis_function = function(x) mean(x$value)
  )), retention = 1, n_iterations = 3, min_iterations = 2,
  instability_bootstrap_B = 20, verbose = FALSE)
  output <- tempfile()
  report <- sparq_report_many(fit, output, include_plots = FALSE, verbose = FALSE)
  expect_true(all(file.exists(unlist(report[c("summary", "errors", "manifest")] ))))
  expect_true(file.exists(file.path(output, "mean_value", "reproducibility_manifest.tsv")))
})
