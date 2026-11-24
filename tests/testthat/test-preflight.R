test_that("workflow preflight validates scalar reruns without launching a full assessment", {
  set.seed(61)
  spots <- data.frame(value = rnorm(40), x = runif(40), y = runif(40))
  fit <- sparq_validate_workflow(
    data = spots,
    analysis_function = function(x) mean(x$value),
    output_type = "scalar",
    verbose = FALSE,
    seed = 5
  )
  expect_s3_class(fit, "sparq_workflow_preflight")
  expect_true(fit$passed)
  expect_lt(fit$retained_size, fit$input_size)
  expect_true(all(fit$checks$status %in% c("passed", "not required")))
})

test_that("workflow preflight rejects partition labels that do not align to retained spots", {
  spots <- data.frame(value = 1:30, x = 1:30, y = 30:1)
  edges <- data.frame(source = as.character(1:29), target = as.character(2:30))
  expect_error(
    sparq_validate_workflow(
      data = spots,
      analysis_function = function(x) setNames(rep("domain", nrow(x)), paste0("wrong_", seq_len(nrow(x)))),
      output_type = "partition",
      partition_graph = edges,
      verbose = FALSE
    ),
    "spot IDs do not match retained data"
  )
})
