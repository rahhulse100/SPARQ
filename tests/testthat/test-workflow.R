test_that("workflow assessment selects the scalar comparator and records outputs", {
  set.seed(1)
  data <- data.frame(value = rnorm(40), x = runif(40), y = runif(40))

  fit <- sparq_assess_workflow(
    data = data,
    analysis_function = function(x) mean(x$value),
    output_type = "scalar",
    retention = 0.75,
    n_iterations = 6,
    seed = 12,
    min_iterations = 3,
    instability_bootstrap_B = 20,
    keep_perturbed_results = TRUE
  )

  expect_s3_class(fit, "sparq_workflow_assessment")
  expect_equal(nrow(fit$comparisons), 6)
  expect_equal(length(fit$perturbed_results), 6)
  expect_true(all(c("support_quality", "instability_median") %in% names(fit$summary)))
  expect_true(is.list(fit$input_provenance))
  expect_equal(fit$input_provenance$n_observations, nrow(data))
  expect_equal(length(fit$input_provenance$row_ids), nrow(data))
})

test_that("workflow assessments provide a compact S3 summary", {
  data <- data.frame(value = 1:20)
  fit <- sparq_assess_workflow(
    data = data,
    analysis_function = function(x) mean(x$value),
    output_type = "scalar",
    n_iterations = 3,
    min_iterations = 2,
    instability_bootstrap_B = 10,
    verbose = FALSE
  )

  expect_output(print(fit), "SPARQ workflow assessment")
  expect_output(print(fit), "Support quality")
})

test_that("fragility localization summarizes all built-in output types", {
  scalar <- sparq_localize_fragility(1, c(0.8, 1.1, 0.9), "scalar")
  expect_null(scalar$localization)
  expect_equal(scalar$global$n_perturbations, 3)

  ranked <- sparq_localize_fragility(
    c("a", "b", "c"),
    list(c("a", "b", "c"), c("b", "a", "d")),
    "ranked",
    top_k = 2
  )
  expect_true("top_k_selection_frequency" %in% names(ranked$features))

  partition <- sparq_localize_fragility(
    c(a = "x", b = "x", c = "y"),
    list(c(a = "x", b = "x"), c(a = "x", c = "y")),
    "partition"
  )
  expect_equal(nrow(partition$spots), 3)

  graph <- sparq_localize_fragility(
    c("a--b", "b--c"),
    list(c("a--b"), c("a--b", "c--d")),
    "graph"
  )
  expect_true("selection_frequency" %in% names(graph$edges))
})

test_that("reports can be returned without writing files", {
  comparison <- data.frame(
    delta = c(0.1, -0.1, 0.05),
    instability = c(0.1, 0.1, 0.05),
    iteration = 1:3,
    status = "ok",
    error_message = NA_character_,
    n_retained = 10L
  )
  assessment <- list(
    result_id = "example",
    full_result = 1,
    comparisons = comparison,
    summary = sparq_support_from_comparisons(
      comparison, n_iterations = 3, min_iterations = 3,
      instability_bootstrap_B = 20
    ),
    failures = comparison[0, , drop = FALSE],
    settings = list(seed = 1)
  )

  report <- sparq_report(assessment, include_plots = FALSE)
  expect_s3_class(report, "sparq_report")
  expect_equal(report$result_id, "example")
})

test_that("reports write a self-contained audit bundle when requested", {
  comparison <- data.frame(
    delta = c(0.1, -0.1, 0.05),
    instability = c(0.1, 0.1, 0.05),
    iteration = 1:3,
    status = "ok",
    error_message = NA_character_,
    n_retained = 10L
  )
  assessment <- list(
    result_id = "example",
    full_result = 1,
    comparisons = comparison,
    summary = sparq_support_from_comparisons(
      comparison, n_iterations = 3, min_iterations = 3,
      instability_bootstrap_B = 20
    ),
    failures = comparison[0, , drop = FALSE],
    settings = list(seed = 1)
  )
  output_dir <- tempfile("sparq-report-")
  on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)

  report <- sparq_report(assessment, output_dir = output_dir,
    include_plots = FALSE)

  expect_true(file.exists(report$output_files$summary))
  expect_true(file.exists(report$output_files$comparisons))
  expect_true(file.exists(report$output_files$report))
  expect_true(file.exists(report$output_files$reproducibility_manifest))
  expect_true(file.exists(report$output_files$package_versions))
  expect_true(file.exists(report$output_files$session_info))
})

test_that("reports export direct-input fingerprints and row identifiers", {
  input <- data.frame(value = rnorm(25))
  assessment <- sparq_assess_workflow(
    data = input,
    analysis_function = function(x) mean(x$value),
    output_type = "scalar",
    n_iterations = 4,
    min_iterations = 2,
    instability_bootstrap_B = 10,
    verbose = FALSE
  )
  output_dir <- tempfile("sparq-provenance-")
  on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)
  report <- sparq_report(assessment, output_dir = output_dir,
    include_plots = FALSE, verbose = FALSE)

  expect_true(file.exists(report$output_files$input_row_ids))
  manifest <- utils::read.delim(report$output_files$reproducibility_manifest,
    stringsAsFactors = FALSE)
  expect_true(all(c("package.git_commit", "input.data_fingerprint",
    "input.row_id_fingerprint", "setting.seed") %in% manifest$key))
})
