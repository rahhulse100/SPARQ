library(SPARQ)
sparq_test_methods <-
function (full_size = FALSE) 
{
    passed <- 0L
    check <- function(ok, label) {
        if (!isTRUE(ok)) 
            stop("FAILED: ", label, call. = FALSE)
        passed <<- passed + 1L
        cat("PASS:", label, "\n")
    }
    errors <- function(code) inherits(tryCatch({
        force(code)
        NULL
    }, error = function(e) e), "error")
    check(all(sparq_bound01(c(-1, 0.40000000000000002, 2)) == 
        c(0, 0.40000000000000002, 1)), "bounding")
    set.seed(97)
    before <- .Random.seed
    ci <- sparq_bootstrap_ci(1:9, B = 300)
    check(identical(.Random.seed, before), "bootstrap preserves caller RNG")
    check(sparq_seed(sparq_seed(1, 1, 12), 2) != sparq_seed(sparq_seed(1, 
        2, 12), 1), "different sample/iteration combinations do not reuse adjacent seeds")
    q1 <- sparq_bootstrap_ci(c(0, 2), B = 20000, conf = 0.5, 
        seed = 51)
    check(all(q1 == c(estimate = 1, lower = 0, upper = 2)) || 
        all(q1 == c(estimate = 1, lower = 1, upper = 2)) || all(q1 == 
        c(estimate = 1, lower = 0, upper = 1)), "two-point median bootstrap support")
    for (x in list(c(0, 0, 1, 3), c(0, 1, 3, 8, 12))) {
        fast <- sparq_bootstrap_ci(x, B = 20000, conf = 0.80000000000000004, 
            seed = 52)
        brute <- sparq_bootstrap_ci(x, statistic = function(z) stats::median(z), 
            B = 20000, conf = 0.80000000000000004, seed = 53)
        check(max(abs(fast - brute)) < 1e-10, "optimized median bootstrap agrees with ordinary resampling quantiles")
    }
    dat <- data.frame(value = seq_len(50), x = seq_len(50), y = 0)
    pt <- sparq_stress_model(dat, 0.75, seed = 1)
    check(nrow(pt) == 37L && !anyDuplicated(pt$value), "75% of 50 retains 37 unique observations")
    check(identical(pt, sparq_stress_model(dat, 0.75, seed = 1)), 
        "retention reproducibility")
    bad <- dat
    bad$x[1] <- NA_real_
    check(errors(sparq_stress_model(bad, 0.75, "contiguous_hole", 
        "x", "y")), "invalid spatial coordinates rejected")
    center <- sparq_with_seed(sparq_seed(17, 1), sample.int(nrow(dat), 
        1))
    expected_hole <- dat[order(-(dat$x - dat$x[center])^2, seq_len(nrow(dat)))[1:25], 
        , drop = FALSE]
    check(identical(sparq_stress_model(dat, 0.5, "contiguous_hole", 
        "x", "y", seed = 17), expected_hole), "spatial hole retains the farthest observations")
    check(identical(sparq_stress_model(dat, 1, "none"), dat), 
        "none preserves input")
    check(sparq_compare_ranked(c("a", "b", "c"), c("c", "b", 
        "a"), 2)$spearman == -1, "rank correlation aligns feature identity")
    check(abs(sparq_compare_ranked(c("a", "b", "c"), c("b", "c", 
        "d"), 2)$replacement - 0.5) < 9.9999999999999998e-13, 
        "top-feature replacement fraction")
    partition_edges <- data.frame(source = c("a", "b", "c"), target = c("b", "c", "d"))
    check(sparq_compare_partition(c(a = 1, b = 1, c = 2, d = 2), 
        c(d = 9, c = 9, a = 8), partition_edges)$ari == 1, "partition IDs align after removal/reordering")
    check(errors(sparq_compare_partition(1:4, 1:3, partition_edges)), "unnamed unequal partitions rejected")
    check(abs(sparq_compare_graph(c("a-b", "b-c"), c("b-c", "c-d"))$instability - 
        2/3) < 9.9999999999999998e-13, "graph Jaccard")
    cmp <- data.frame(delta = c(-0.10000000000000001, 0.10000000000000001), 
        instability = c(0.10000000000000001, 0.10000000000000001))
    s <- sparq_support_from_comparisons(cmp, 2, min_iterations = 2, 
        instability_bootstrap_B = 50)
    check(s$bias == 0 && s$uncertainty == 0 && s$support_quality == 
        1, "stated score formula retained, including its plateau")
    cmp2 <- rbind(cmp, data.frame(delta = NA_real_, instability = NA_real_))
    s2 <- sparq_support_from_comparisons(cmp2, 3, min_iterations = 3, 
        instability_bootstrap_B = 50)
    check(s2$n_iterations == 2 && s2$n_failed == 1 && abs(s2$iteration_support - 
        2/3) < 9.9999999999999998e-13, "only successful iterations count")
    check(errors(sparq_support_from_comparisons(cmp, 2, reference_scale = NA_real_)), 
        "invalid scale rejected")
    zero <- sparq_support_from_comparisons(cmp, 2, reference_scale = 0, 
        instability_bootstrap_B = 50)
    check(zero$reference_scale == 9.9999999999999995e-07, "explicit scale floor")
    base <- sparq_run(dat, sparq_benchmark_mean, sparq_compare_scalar, 
        retention = 1, n_iterations = 5, min_iterations = 2, 
        instability_bootstrap_B = 20)
    check(all(base$comparisons$instability == 0) && base$summary$support_quality == 
        1, "full retention deterministic baseline")
    failure <- function(data, retention, iteration) {
        if (iteration == 2) 
            stop("intentional test failure")
        data[1:10, , drop = FALSE]
    }
    f <- sparq_run(dat, sparq_benchmark_mean, sparq_compare_scalar, 
        "custom", 0.5, 3, custom_function = failure, min_iterations = 3, 
        instability_bootstrap_B = 20)
    check(nrow(f$failures) == 1 && f$summary$n_iterations == 
        2 && f$summary$n_failed == 1, "failure ledger")
    cache <- tempfile()
    dir.create(cache)
    on.exit(unlink(cache, recursive = TRUE), add = TRUE)
    c1 <- sparq_run(dat, sparq_benchmark_mean, sparq_compare_scalar, 
        n_iterations = 3, cache_dir = cache, cache_key = "test-mean-v1", 
        instability_bootstrap_B = 20)
    c2 <- sparq_run(dat, sparq_benchmark_mean, sparq_compare_scalar, 
        n_iterations = 3, cache_dir = cache, cache_key = "test-mean-v1", 
        instability_bootstrap_B = 20)
    check(identical(c1, c2), "cache resume")
    c3 <- sparq_run(dat, sparq_benchmark_mean, sparq_compare_scalar, 
        retention = 0.25, n_iterations = 3, cache_dir = cache, 
        cache_key = "test-mean-v1", instability_bootstrap_B = 20)
    check(!identical(c1$comparisons, c3$comparisons), "cache cannot reuse another retention")
    perfect <- sparq_binary_metrics(c(1, 1, 0, 0), c(0.90000000000000002, 
        0.80000000000000004, 0.20000000000000001, 0.10000000000000001), 
        0.5)
    check(all(unlist(perfect[c("sensitivity", "specificity", 
        "balanced_accuracy", "auroc", "auprc")]) == 1), "perfect prediction metric reference")
    tied <- sparq_binary_metrics(c(1, 0, 0, 0), rep(0.5, 4), 
        0.5)
    check(tied$auroc == 0.5 && tied$auprc == 0.75, "tied-score AUROC/AP equals chance/prevalence")
    check(identical(tied, sparq_binary_metrics(c(0, 0, 1, 0), 
        rep(0.5, 4), 0.5)), "AP is invariant to tied-row ordering")
    scores <- c(0.10000000000000001, 0.20000000000000001, 0.20000000000000001, 
        0.59999999999999998, 0.80000000000000004, 0.90000000000000002)
    truth <- c(0, 1, 0, 1, 1, 0)
    tab <- sparq_threshold_candidates(scores, truth)
    brute <- vapply(tab$threshold, function(t) sparq_binary_metrics(truth, 
        scores, t)$balanced_accuracy, numeric(1))
    check(max(abs(tab$balanced_accuracy - brute)) < 9.9999999999999998e-13, 
        "threshold objective independently recomputed")
    a <- sparq_calibrate_threshold(scores, truth, c(0.10000000000000001, 
        0.90000000000000002), c(0, 1), B = 20)
    b <- sparq_calibrate_threshold(scores, truth, c(0.10000000000000001, 
        0.90000000000000002), c(1, 0), B = 20)
    check(identical(a$threshold, b$threshold) && identical(a$bootstrap_thresholds, 
        b$bootstrap_thresholds), "held-out labels cannot affect calibration")
    classes <- c("Stable", "Mild instability", "Moderate instability", 
        "Severe instability")
    nclass <- if (full_size) 
        300L
    else 10L
    design <- sparq_methods_design()
    design$n_per_class <- nclass
    design$n_samples <- 4L * nclass
    design$n_replicates <- 2L
    design$n_iterations <- 5L
    design$min_iterations <- 3L
    for (k in c("threshold_bootstrap_B", "instability_bootstrap_B", 
        "reliability_bootstrap_B", "reliability_permutation_B", 
        "stress_bootstrap_B")) design[[k]] <- 20L
    design$min_observations_per_bin <- 2L
    design$reliability_bins <- 4L
    test_data <- sparq_with_seed(4321, data.frame(benchmark_id = paste0("id", 
        seq_len(4 * nclass)), difficulty = rep(classes, each = nclass), 
        true_instability = pmax(0, stats::rnorm(4 * nclass, rep(c(0.02, 
            0.070000000000000007, 0.13, 0.23999999999999999), 
            each = nclass), rep(c(0.029999999999999999, 0.059999999999999998, 
            0.080000000000000002, 0.12), each = nclass)))))
    d1 <- sparq_known_truth_benchmark(test_data, design, n_workers = 1)
    d2 <- sparq_known_truth_benchmark(test_data, design, n_workers = 2)
    check(identical(d1$replicate_results, d2$replicate_results) && 
        identical(d1$comparisons, d2$comparisons), "serial and parallel benchmark equality")
    check(nrow(d1$replicate_results) == 2 * (4 * nclass - (floor(nclass * 
        0.69999999999999996) + floor(3 * nclass * 0.69999999999999996))), 
        "held-out row audit")
    check(all(d1$comparisons$instability == abs(d1$comparisons$delta)), 
        "instability calculated from actual scalar comparison")
    check(all(d1$comparisons$n_retained == 37L), "benchmark actually applies 75% retention")
    check(!anyDuplicated(paste(d1$all_results$replicate, d1$all_results$benchmark_id)), 
        "unique sample-replicate records")
    check(!any(paste(d1$training_results$replicate, d1$training_results$benchmark_id) %in% 
        paste(d1$replicate_results$replicate, d1$replicate_results$benchmark_id)), 
        "train/held-out sample separation")
    latent <- setNames(test_data$true_instability, test_data$benchmark_id)
    contaminated <- test_data
    contaminated$support_quality <- as.numeric(contaminated$difficulty == 
        "Stable")
    check(errors(sparq_validate_benchmark(contaminated, design)), 
        "precomputed truth-derived scores cannot enter benchmark inputs")
    check(all(d1$replicate_results$true_instability == latent[d1$replicate_results$benchmark_id]), 
        "exact input latent values preserved")
    obs <- sparq_make_observations(test_data, design, 1)
    ob2 <- sparq_make_observations(test_data, design, 2)
    check(!identical(obs, ob2), "replicate observation generation varies")
    shuffled <- test_data
    shuffled$difficulty <- rev(shuffled$difficulty)
    check(identical(obs, sparq_make_observations(shuffled, design, 
        1)), "class labels do not enter observation generator")
    one <- sparq_run(obs[[1]], sparq_benchmark_mean, sparq_compare_scalar, 
        n_iterations = 5, min_iterations = 3, instability_bootstrap_B = 20)
    obs[[1]]$truth <- 1
    obs[[1]]$true_instability <- 100
    two <- sparq_run(obs[[1]], sparq_benchmark_mean, sparq_compare_scalar, 
        n_iterations = 5, min_iterations = 3, instability_bootstrap_B = 20)
    check(identical(one$comparisons, two$comparisons), "unused truth columns cannot change analysis/comparison")
    r <- d1$reliability
    check(all(r$curve$n >= 2) && sum(r$curve$n) == nrow(d1$replicate_results), 
        "adaptive bin membership")
    check(r$n_samples == length(unique(d1$replicate_results$benchmark_id)), 
        "reliability resampling unit is benchmark ID")
    calc <- abs(stats::median(d1$replicate_results$true_instability) - 
        (1 - stats::median(d1$replicate_results$support_quality)))
    check(identical(r$median_scale_discrepancy, calc), "requested median discrepancy formula")
    check(is.na(r$monotonicity_p_value) || r$monotonicity_p_value >= 
        1/21, "Monte Carlo P-values cannot be zero")
    for (x in list(c(1, 2), c(1, 5, 7, 9))) {
        w <- seq_along(x)
        check(sparq_weighted_median(x, w) == stats::median(rep(x, 
            w)), "cluster weighted median")
    }
    sp <- sparq_known_truth_stress(test_data, design, calibration_results = d1, 
        n_workers = 2)
    wrong_input <- test_data
    wrong_input$true_instability[1] <- wrong_input$true_instability[1] + 
        0.10000000000000001
    check(errors(sparq_known_truth_stress(wrong_input, design, 
        calibration_results = d1, n_workers = 1)), "stress refuses a different calibration input table")
    check(nrow(sp$sample_results) == 4 * nclass * 2 * 5, "stress result count")
    for (rid in 1:2) {
        z <- sp$sample_results[sp$sample_results$replicate == 
            rid, , drop = FALSE]
        check(length(unique(z$threshold)) == 1 && unique(z$threshold) == 
            d1$replicate_metrics$threshold[rid], "stress threshold transferred unchanged")
        z75 <- z[z$retention == 0.75, , drop = FALSE]
        orig <- d1$all_results[d1$all_results$replicate == rid, 
            , drop = FALSE]
        check(identical(z75$support_quality, orig$support_quality), 
            "stress 75% reuses identical dataset/seeds and score")
    }
    check(all(sp$sample_results$support_quality[sp$sample_results$retention == 
        1] == 1), "stress 100% baseline")
    check(all(sp$comparisons$n_retained[sp$comparisons$retention == 
        0.10000000000000001] == 5L), "stress 10% actually retains five observations")
    check(nrow(sp$failures) == 0, "stress failure ledger empty only when all iterations succeed")
    spatial_design <- design
    spatial_design$n_replicates <- 1L
    spatial_design$stress_retentions <- c(1, 0.5)
    spatial <- sparq_stress_response(list(a = dat, b = dat), 
        sparq_benchmark_mean, sparq_compare_scalar, design = spatial_design, 
        stress_models = "contiguous_hole", x_col = "x", y_col = "y", 
        n_workers = 1)
    check(nrow(spatial$sample_results) == 4L && all(is.na(spatial$sample_results$threshold)), 
        "unlabeled spatial stress reports no invented classifier")
    check(all(spatial$comparisons$n_retained[spatial$comparisons$retention == 
        0.5] == 25), "spatial stress actually removes half the observations")
    dfile <- tempfile()
    on.exit(unlink(dfile), add = TRUE)
    sparq_write_design(design, dfile)
    check(isTRUE(all.equal(design, sparq_read_design(dfile))), 
        "versioned YAML round trip")
    output <- tempfile()
    on.exit(unlink(output, recursive = TRUE), add = TRUE)
    sparq_export_benchmark(d1, output)
    check(all(file.exists(file.path(output, c("replicate_results.tsv", 
        "replicate_metrics.tsv", "threshold_summary.tsv", "reliability_bins.tsv", 
        "reliability_summary.tsv", "figure3_benchmark.rds", "design.yml")))), 
        "C/D/E exports contain all requested outputs")
    check(nrow(utils::read.delim(file.path(output, "replicate_results.tsv"))) == 
        nrow(d1$replicate_results), "TSV row counts survive export")
    check(errors(sparq_export_benchmark(d1, output)), "existing benchmark outputs cannot be overwritten")
    stress_output <- tempfile()
    on.exit(unlink(stress_output, recursive = TRUE), add = TRUE)
    sparq_export_stress(sp, stress_output)
    check(all(file.exists(file.path(stress_output, c("sample_results.tsv", 
        "replicate_metrics.tsv", "summary.tsv", "paired_differences.tsv", 
        "failures.tsv", "stress_response.rds")))), "stress exports contain summaries and failure ledger")
    cat("All", passed, "checks passed. These are software tests, not publication results.\n")
    invisible(list(n_checks = passed, benchmark = d1, stress = sp))
}
environment(sparq_test_methods) <- asNamespace("SPARQ")
sparq_test_methods()
