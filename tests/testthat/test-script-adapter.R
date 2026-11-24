test_that("script adapter injects input and retrieves the declared result", {
  script <- tempfile(fileext = ".R")
  writeLines(
    "banksy_domains <- setNames(ifelse(seurat_obj$value > 0, 'tumor', 'stroma'), rownames(seurat_obj))",
    script
  )
  on.exit(unlink(script), add = TRUE)

  set.seed(3)
  data <- data.frame(value = rnorm(40), x = runif(40), y = runif(40))
  rownames(data) <- paste0("spot_", seq_len(nrow(data)))
  edges <- data.frame(source = rownames(data)[-nrow(data)], target = rownames(data)[-1L])

  fit <- sparq_assess_script(
    data = data,
    script = script,
    input_object = "seurat_obj",
    result_object = "banksy_domains",
    output_type = "partition",
    partition_graph = edges,
    retention = 0.75,
    n_iterations = 5,
    min_iterations = 3,
    x_col = "x",
    y_col = "y",
    seed = 4,
    instability_bootstrap_B = 20
  )

  expect_s3_class(fit, "sparq_script_assessment")
  expect_equal(nrow(fit$comparisons), 5)
  expect_equal(fit$script$result_object, "banksy_domains")
})

test_that("script adapter rejects a missing declared reference result", {
  script <- tempfile(fileext = ".R")
  writeLines("other_result <- 1", script)
  on.exit(unlink(script), add = TRUE)

  expect_error(
    sparq_assess_script(
      data = data.frame(value = 1:10),
      script = script,
      input_object = "input_data",
      result_object = "expected_result",
      output_type = "scalar",
      n_iterations = 2,
      min_iterations = 2,
      instability_bootstrap_B = 2
    ),
    "did not create"
  )
})
