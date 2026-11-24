test_that("spatial partition comparator is label invariant", {
  edges <- data.frame(source = c("a", "b", "c"), target = c("b", "c", "d"))
  full <- c(a = "tumor", b = "tumor", c = "stroma", d = "stroma")
  rerun <- c(a = "two", b = "two", c = "one", d = "one")
  x <- sparq_compare_spatial_partition(full, rerun, edges)
  y <- sparq_compare_partition(full, rerun, edges)

  expect_equal(x$ari, 1)
  expect_equal(x$boundary_f1, 1)
  expect_equal(x$adjacency_jaccard, 1)
  expect_equal(x$instability, 0)
  expect_equal(y, x)
})

test_that("boundary agreement permits one-hop displacement", {
  edges <- data.frame(
    source = c("a", "b", "c", "d"),
    target = c("b", "c", "d", "e")
  )
  full <- c(a = "A", b = "A", c = "B", d = "B", e = "B")
  rerun <- c(a = "x", b = "x", c = "x", d = "y", e = "y")
  x <- sparq_compare_spatial_partition(full, rerun, edges)

  # The reference interface is b/c and the rerun interface is c/d. These
  # interface sets are one spatial hop apart, so tolerant agreement is exact.
  expect_equal(x$boundary_f1, 1)
  expect_equal(x$adjacency_jaccard, 1)
  expect_lt(x$ari, 1)
})

test_that("spatial partition comparator aligns a retained rerun by spot ID", {
  edges <- data.frame(source = c("a", "b", "c"), target = c("b", "c", "d"))
  full <- c(a = "A", b = "A", c = "B", d = "B")
  rerun <- c(b = "first", c = "second", d = "second")
  x <- sparq_compare_spatial_partition(full, rerun, edges)

  expect_equal(x$n_aligned, 3L)
  expect_equal(x$ari, 1)
  expect_equal(x$boundary_f1, 1)
  expect_equal(x$adjacency_jaccard, 1)
})

test_that("domain adjacency detects an unmatched split domain", {
  edges <- data.frame(source = c("a", "b", "c", "d"), target = c("b", "c", "d", "e"))
  full <- c(a = "A", b = "A", c = "A", d = "B", e = "B")
  rerun <- c(a = "x", b = "x", c = "y", d = "z", e = "z")
  x <- sparq_compare_spatial_partition(full, rerun, edges)

  expect_lt(x$adjacency_jaccard, 1)
})

test_that("spatial partition comparator validates weights and graph alignment", {
  edges <- data.frame(source = c("a", "b"), target = c("b", "c"))
  labels <- c(a = "A", b = "A", c = "B")
  expect_error(
    sparq_compare_spatial_partition(labels, labels, edges, weights = c(ari = 1)),
    "weights must be named"
  )
  expect_error(
    sparq_compare_spatial_partition(labels[c("a", "c")], labels[c("a", "c")], edges),
    "No supplied neighborhood edges"
  )
})

test_that("spatial partition comparator can drive an adaptive assessment", {
  spots <- data.frame(
    spot_id = letters[1:8],
    x = rep(1:4, each = 2),
    y = rep(1:2, 4),
    value = c(0, 0, 0, 0, 1, 1, 1, 1)
  )
  edges <- data.frame(source = letters[1:7], target = letters[2:8])
  analysis <- function(x) {
    stats::setNames(ifelse(x$value < 0.5, "low", "high"), x$spot_id)
  }
  comparator <- function(full, perturbed) {
    sparq_compare_spatial_partition(full, perturbed, edges)
  }
  x <- sparq_assess_procurement_adaptive(
    data = spots, analysis_function = analysis, output_type = "partition",
    comparator = comparator, partition_graph = edges,
    x_col = "x", y_col = "y", retention = 0.75,
    scenarios = "uniform_random", reliability_tolerance = 0.5,
    min_iterations = 2, max_iterations = 2, batch_size = 1, verbose = FALSE
  )
  expect_true("boundary_f1" %in% names(x$assessments$uniform_random$comparisons))
  expect_true("adjacency_jaccard" %in% names(x$assessments$uniform_random$comparisons))
})

test_that("spatial partition comparator can drive a fixed procurement assessment", {
  spots <- data.frame(
    spot_id = letters[1:8],
    x = rep(1:4, each = 2),
    y = rep(1:2, 4),
    value = c(0, 0, 0, 0, 1, 1, 1, 1)
  )
  edges <- data.frame(source = letters[1:7], target = letters[2:8])
  analysis <- function(x) {
    stats::setNames(ifelse(x$value < 0.5, "low", "high"), x$spot_id)
  }
  comparator <- function(full, perturbed) {
    sparq_compare_spatial_partition(full, perturbed, edges)
  }
  x <- sparq_assess_procurement(
    data = spots, analysis_function = analysis, output_type = "partition",
    comparator = comparator, partition_graph = edges,
    x_col = "x", y_col = "y", retention = 0.75,
    n_iterations = 2, scenarios = "uniform_random", preflight = FALSE,
    min_iterations = 2, instability_bootstrap_B = 20, verbose = FALSE
  )
  expect_true("boundary_f1" %in% names(x$assessments$uniform_random$comparisons))
  expect_true("adjacency_jaccard" %in% names(x$assessments$uniform_random$comparisons))
})
