sparq_assess <-
function (full_result, perturbed_results, comparator, result_id = "result",
    reference_scale = 1, seed = 1, ...)
{
    if (!length(perturbed_results))
        stop("perturbed_results is empty.")
    rows <- lapply(seq_along(perturbed_results), function(i) sparq_comparison_row(comparator,
        full_result, perturbed_results[[i]], i))
    comparisons <- sparq_bind_rows(rows)
    summary <- sparq_support_from_comparisons(comparisons, length(perturbed_results),
        reference_scale, seed = seed, ...)
    summary$result_id <- result_id
    out <- list(result_id = result_id, full_result = full_result,
        perturbed_results = perturbed_results, comparisons = comparisons,
        summary = summary,
        failures = comparisons[comparisons$status != "ok", , drop = FALSE],
        settings = list(reference_scale = reference_scale, seed = seed,
            n_iterations = length(perturbed_results)))
    class(out) <- "sparq_assessment"
    out
}
sparq_assess_cohort <-
function (sample_data, analysis_function, comparator, ...)
{
    out <- sparq_assess_cohort_parallel(sample_data, analysis_function,
        comparator, workers = 1, ...)
    if (nrow(out$errors))
        warning("Some cohort assessments failed; inspect attr(result, 'errors').")
    ans <- out$results
    attr(ans, "errors") <- out$errors
    ans
}
sparq_assess_cohort_parallel <-
function (sample_data, analysis_function, comparator, stress_model = "uniform_random",
    retention = 0.75, n_iterations = 100, workers = 1, checkpoint_dir = NULL,
    resume = TRUE, error_log = NULL, seed = 1, cache_key = NULL,
    ...)
{
    if (!length(sample_data))
        stop("sample_data is empty.")
    if (is.null(names(sample_data)))
        names(sample_data) <- paste0("sample_", seq_along(sample_data))
    if (anyDuplicated(names(sample_data)) || anyNA(names(sample_data)))
        stop("Sample IDs must be unique.")
    out <- sparq_parallel_map(seq_along(sample_data), function(i) {
        tryCatch(list(ok = TRUE, value = sparq_run(sample_data[[i]],
            analysis_function, comparator, stress_model, retention,
            n_iterations, names(sample_data)[i], seed = sparq_seed(seed,
                i), cache_dir = checkpoint_dir, resume = resume,
            cache_key = cache_key, ...)), error = function(e) list(ok = FALSE,
            error = conditionMessage(e)))
    }, workers)
    rows <- lapply(seq_along(out), function(i) if (out[[i]]$ok)
        out[[i]]$value$summary)
    errors <- sparq_bind_rows(lapply(seq_along(out), function(i) {
        z <- out[[i]]
        if (!z$ok)
            return(data.frame(result_id = names(sample_data)[i],
                error_message = z$error))
        f <- z$value$failures
        if (nrow(f)) {
            f$result_id <- names(sample_data)[i]
            f
        }
        else NULL
    }))
    if (nrow(errors) && !is.null(error_log))
        utils::write.table(errors, error_log, sep = "\t", row.names = FALSE,
            quote = FALSE)
    list(results = sparq_bind_rows(rows), errors = errors, expected_samples = names(sample_data),
        completed_samples = names(sample_data)[vapply(out, function(z) z$ok &&
            z$value$summary$n_failed == 0L, logical(1))], assessments = out)
}
sparq_atomic_rds <-
function (object, path)
{
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    f <- tempfile(tmpdir = dirname(path))
    on.exit(unlink(f))
    saveRDS(object, f)
    if (!file.rename(f, path))
        stop("Could not write checkpoint: ", path)
    invisible(path)
}
sparq_benchmark_mean <-
function (data)
mean(data$value)
sparq_benchmark_replicate <-
function (b, d, r, keep_comparisons = TRUE)
{
    split <- sparq_split_benchmark(b, d$train_fraction, sparq_seed(d$seed +
        r, 0, 11))
    datasets <- sparq_make_observations(b, d, r)
    sample_ids <- sort(b$benchmark_id)
    fits <- lapply(seq_len(nrow(b)), function(i) {
        fit <- sparq_run(datasets[[i]], sparq_benchmark_mean,
            sparq_compare_scalar, stress_model = "uniform_random",
            retention = d$retention, n_iterations = d$n_iterations,
            result_id = b$benchmark_id[i], reference_scale = d$reference_scale,
            seed = sparq_seed(d$seed + r, match(b$benchmark_id[i],
                sample_ids), 12), min_iterations = d$min_iterations,
            min_reproducibility = d$min_reproducibility, max_relative_bias = d$max_relative_bias,
            max_relative_uncertainty = d$max_relative_uncertainty,
            instability_bootstrap_B = d$instability_bootstrap_B,
            confidence = d$confidence)
        row <- cbind(b[i, , drop = FALSE], fit$summary[, setdiff(names(fit$summary),
            "result_id"), drop = FALSE])
        row$reference_result <- fit$full_result
        row$replicate <- r
        row$seed <- d$seed + r
        row$observation_fingerprint <- sparq_fingerprint(datasets[[i]])
        row$split <- if (i %in% split$training)
            "training"
        else "held_out"
        comparison <- fit$comparisons
        comparison$benchmark_id <- b$benchmark_id[i]
        comparison$replicate <- r
        list(row = row, comparisons = if (keep_comparisons) comparison else NULL)
    })
    rows <- sparq_bind_rows(lapply(fits, `[[`, "row"))
    comparisons <- sparq_bind_rows(lapply(fits, `[[`, "comparisons"))
    if (any(rows$n_failed > 0L) || any(!is.finite(rows$support_quality)))
        stop("Replicate ", r, " contains failed/invalid assessments; it is not a complete benchmark.")
    training <- rows[split$training, , drop = FALSE]
    held <- rows[split$held_out, , drop = FALSE]
    calibration <- sparq_calibrate_threshold(training$support_quality,
        training$truth, held$support_quality, held$truth, B = d$threshold_bootstrap_B,
        confidence = d$confidence, seed = sparq_seed(d$seed +
            r, 0, 13))
    rows$threshold <- calibration$threshold
    rows$threshold_ci_low <- calibration$threshold_interval$lower
    rows$threshold_ci_high <- calibration$threshold_interval$upper
    rows$predicted_stable <- as.integer(rows$support_quality >=
        calibration$threshold)
    metrics <- calibration$held_out_metrics
    metrics$replicate <- r
    metrics$seed <- d$seed + r
    metrics$threshold <- calibration$threshold
    metrics$threshold_ci_low <- calibration$threshold_interval$lower
    metrics$threshold_ci_high <- calibration$threshold_interval$upper
    metrics$n_training <- nrow(training)
    metrics$n_held_out <- nrow(held)
    metrics$informative_training_scores <- calibration$informative
    list(rows = rows, metrics = metrics, comparisons = comparisons,
        threshold_bootstrap = data.frame(replicate = r, draw = seq_len(d$threshold_bootstrap_B),
            threshold = calibration$bootstrap_thresholds))
}
sparq_binary_metrics <-
function (truth, score, threshold)
{
    sparq_validate_scores(score, truth)
    sparq_check_number(threshold, "threshold")
    if (length(unique(truth)) != 2L)
        stop("Both truth classes are required for performance metrics.")
    predicted <- score >= threshold
    positive <- truth == 1
    tp <- sum(predicted & positive)
    fn <- sum(!predicted & positive)
    tn <- sum(!predicted & !positive)
    fp <- sum(predicted & !positive)
    n1 <- sum(positive)
    n0 <- sum(!positive)
    auroc <- (sum(rank(score, ties.method = "average")[positive]) -
        n1 * (n1 + 1)/2)/(n1 * n0)
    o <- order(1 - score, decreasing = TRUE)
    s <- (1 - score)[o]
    y <- (1 - truth)[o]
    ends <- c(which(diff(s) != 0), length(s))
    hits <- cumsum(y)[ends]
    recall <- hits/sum(y)
    precision <- hits/ends
    ap <- sum(diff(c(0, recall)) * precision)
    data.frame(sensitivity = tp/n1, specificity = tn/n0, balanced_accuracy = (tp/n1 +
        tn/n0)/2, auroc = auroc, auprc = ap, tp = tp, fn = fn,
        tn = tn, fp = fp, unstable_prevalence = n0/(n1 + n0))
}
sparq_bind_rows <-
function (rows)
{
    rows <- Filter(function(x) is.data.frame(x) && nrow(x) >
        0L, rows)
    if (!length(rows))
        return(data.frame())
    cols <- unique(unlist(lapply(rows, names)))
    out <- do.call(rbind, lapply(rows, function(x) {
        for (nm in setdiff(cols, names(x))) x[[nm]] <- NA
        x[, cols, drop = FALSE]
    }))
    rownames(out) <- NULL
    out
}
sparq_bootstrap_ci <-
function (x, statistic = stats::median, B = 2000, conf = 0.94999999999999996,
    seed = 1)
{
    sparq_check_number(B, "B", 1, integer = TRUE)
    sparq_check_number(conf, "conf", .Machine$double.eps, 1 -
        .Machine$double.eps)
    x <- x[is.finite(x)]
    if (length(x) < 2L)
        return(c(estimate = if (length(x)) statistic(x) else NA_real_,
            lower = NA_real_, upper = NA_real_))
    values <- sparq_with_seed(seed, {
        if (identical(statistic, stats::median)) {
            ordered <- sort(x)
            n <- length(x)
            if (n%%2L) {
                k <- (n + 1)/2
                u <- stats::rbeta(B, k, n + 1 - k)
                ordered[pmin(n, pmax(1, ceiling(n * u)))]
            }
            else {
                k <- n/2
                u <- stats::rbeta(B, k, n + 1 - k)
                v <- u + (1 - u) * stats::rbeta(B, 1, n - k)
                (ordered[pmin(n, pmax(1, ceiling(n * u)))] +
                  ordered[pmin(n, pmax(1, ceiling(n * v)))])/2
            }
        }
        else replicate(B, statistic(x[sample.int(length(x), length(x),
            TRUE)]))
    })
    c(estimate = statistic(x), lower = unname(stats::quantile(values,
        (1 - conf)/2)), upper = unname(stats::quantile(values,
        1 - (1 - conf)/2)))
}
sparq_bound01 <-
function (x)
pmax(0, pmin(1, x))
sparq_calibrate_threshold <-
function (support_quality, truth, held_out_support_quality = NULL,
    held_out_truth = NULL, folds = 0, B = 2000, seed = 1, min_class_n = 2,
    confidence = 0.94999999999999996)
{
    sparq_validate_scores(support_quality, truth)
    sparq_check_number(B, "B", 1, integer = TRUE)
    sparq_check_number(confidence, "confidence", .Machine$double.eps,
        1 - .Machine$double.eps)
    sparq_check_number(folds, "folds", 0, integer = TRUE)
    if (sum(truth == 1) < min_class_n || sum(truth == 0) < min_class_n)
        stop("Too few training cases in one or both classes.")
    best <- sparq_choose_threshold(support_quality, truth)
    draws <- sparq_with_seed(seed, {
        groups <- split(seq_along(truth), truth)
        replicate(B, {
            ix <- unlist(lapply(groups, function(g) g[sample.int(length(g),
                length(g), TRUE)]), use.names = FALSE)
            sparq_choose_threshold(support_quality[ix], truth[ix])$threshold
        })
    })
    ci <- unname(stats::quantile(draws, c((1 - confidence)/2,
        1 - (1 - confidence)/2)))
    held <- NULL
    if (xor(is.null(held_out_support_quality), is.null(held_out_truth)))
        stop("Supply both held-out vectors.")
    if (!is.null(held_out_truth))
        held <- sparq_binary_metrics(held_out_truth, held_out_support_quality,
            best$threshold)
    cv <- NULL
    if (folds > 0L) {
        if (folds < 2L || folds > min(table(truth)))
            stop("folds must be 2 through the smallest class size.")
        fid <- sparq_with_seed(sparq_seed(seed, 0, 6), {
            z <- integer(length(truth))
            for (g in split(seq_along(truth), truth)) z[g] <- sample(rep(seq_len(folds),
                length.out = length(g)))
            z
        })
        cv <- sparq_bind_rows(lapply(seq_len(folds), function(f) {
            th <- sparq_choose_threshold(support_quality[fid !=
                f], truth[fid != f])$threshold
            d <- sparq_binary_metrics(truth[fid == f], support_quality[fid ==
                f], th)
            d$fold <- f
            d$threshold <- th
            d
        }))
    }
    list(threshold = best$threshold, training_metrics = best,
        held_out_metrics = held, threshold_interval = data.frame(threshold = best$threshold,
            lower = ci[1], upper = ci[2]), bootstrap_thresholds = draws,
        cross_validation = cv, informative = length(unique(support_quality)) >
            1L, bootstrap_unit = "training sample; class-stratified",
        confidence = confidence, tie_rule = "central maximizing candidate (lower central candidate for an even tie)")
}
sparq_check_number <-
function (x, name, lower = -Inf, upper = Inf, integer = FALSE)
{
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x < lower || x > upper || (integer && x != floor(x)))
        stop(name, " is invalid.", call. = FALSE)
    invisible(x)
}
sparq_choose_threshold <-
function (score, truth)
{
    tab <- sparq_threshold_candidates(score, truth)
    best <- which(abs(tab$balanced_accuracy - max(tab$balanced_accuracy)) <=
        9.9999999999999998e-13)
    tab[best[ceiling(length(best)/2)], , drop = FALSE]
}
sparq_compare_graph <-
function (full_edges, perturbed_edges)
{
    if (is.data.frame(full_edges) || is.data.frame(perturbed_edges))
        stop("Supply canonical edge-ID vectors, not edge data frames.")
    a <- unique(as.character(full_edges))
    b <- unique(as.character(perturbed_edges))
    if (anyNA(c(a, b)))
        stop("Edge IDs cannot contain missing values.")
    u <- union(a, b)
    if (!length(u))
        stop("Two empty graphs do not provide assessable edge agreement.")
    j <- length(intersect(a, b))/length(u)
    data.frame(edge_jaccard = j, instability = 1 - j)
}
sparq_compare_partition <-
function (full_labels, perturbed_labels)
{
    a_names <- names(full_labels)
    b_names <- names(perturbed_labels)
    if (xor(is.null(a_names), is.null(b_names)))
        stop("Name both partitions by observation ID, or neither.")
    if (!is.null(a_names)) {
        if (anyNA(c(a_names, b_names)) || anyDuplicated(a_names) ||
            anyDuplicated(b_names))
            stop("Partition IDs must be unique and nonmissing.")
        common <- intersect(a_names, b_names)
        full_labels <- full_labels[match(common, a_names)]
        perturbed_labels <- perturbed_labels[match(common, b_names)]
    }
    else if (length(full_labels) != length(perturbed_labels))
        stop("Unequal partitions need observation-ID names for alignment.")
    if (length(full_labels) < 2L || anyNA(full_labels) || anyNA(perturbed_labels))
        stop("At least two nonmissing aligned partition labels are required.")
    tab <- table(as.character(full_labels), as.character(perturbed_labels))
    ch2 <- function(x) x * (x - 1)/2
    nij <- sum(ch2(tab))
    ai <- sum(ch2(rowSums(tab)))
    bj <- sum(ch2(colSums(tab)))
    expected <- ai * bj/ch2(sum(tab))
    denom <- (ai + bj)/2 - expected
    ari <- if (abs(denom) < .Machine$double.eps)
        1
    else (nij - expected)/denom
    data.frame(ari = ari, instability = 1 - ari, n_aligned = sum(tab))
}
sparq_compare_ranked <-
function (full_rank, perturbed_rank, k = 100)
{
    sparq_check_number(k, "k", 1, integer = TRUE)
    full_rank <- as.character(full_rank)
    perturbed_rank <- as.character(perturbed_rank)
    if (anyNA(c(full_rank, perturbed_rank)) || anyDuplicated(full_rank) ||
        anyDuplicated(perturbed_rank))
        stop("Ranked feature identifiers must be unique and nonmissing.")
    a <- head(full_rank, k)
    b <- head(perturbed_rank, k)
    if (!length(a) || !length(b))
        stop("Ranked feature sets must not be empty.")
    common <- intersect(full_rank, perturbed_rank)
    rho <- if (length(common) >= 2L)
        stats::cor(match(common, full_rank), match(common, perturbed_rank),
            method = "spearman")
    else NA_real_
    j <- length(intersect(a, b))/length(union(a, b))
    data.frame(spearman = rho, jaccard = j, replacement = 1 -
        length(intersect(a, b))/length(a), instability = 1 -
        j, n_rank_common = length(common))
}
sparq_compare_scalar <-
function (full_result, perturbed_result)
{
    sparq_check_number(full_result, "full_result")
    sparq_check_number(perturbed_result, "perturbed_result")
    delta <- perturbed_result - full_result
    data.frame(delta = delta, instability = abs(delta))
}
sparq_comparison_row <-
function (comparator, full_result, perturbed_result, iteration)
{
    d <- as.data.frame(comparator(full_result, perturbed_result))
    if (nrow(d) != 1L || !"instability" %in% names(d) || !is.numeric(d$instability) ||
        !is.finite(d$instability) || d$instability < 0 || ("delta" %in%
        names(d) && (!is.numeric(d$delta) || !is.finite(d$delta))))
        stop("Comparator must return one row with finite nonnegative instability and finite delta if provided.")
    d$iteration <- iteration
    d$status <- "ok"
    d$error_message <- ""
    d
}
sparq_curves <-
function (truth, score)
{
    sparq_validate_scores(score, truth)
    if (length(unique(truth)) != 2L)
        stop("Both classes required.")
    one <- function(y, s, kind) {
        o <- order(s, decreasing = TRUE)
        y <- y[o]
        s <- s[o]
        ends <- c(which(diff(s) != 0), length(s))
        tp <- cumsum(y)[ends]
        if (kind == "roc")
            data.frame(threshold = c(Inf, s[ends]), false_positive_rate = c(0,
                (ends - tp)/sum(y == 0)), sensitivity = c(0,
                tp/sum(y)))
        else data.frame(threshold = c(Inf, s[ends]), recall = c(0,
            tp/sum(y)), precision = c(1, tp/ends))
    }
    list(roc = one(truth, score, "roc"), precision_recall = one(1 -
        truth, 1 - score, "pr"))
}
sparq_engine_signature <-
function ()
{
    e <- environment(sparq_engine_signature)
    nms <- ls(e, all.names = TRUE)
    nms <- sort(nms[startsWith(nms, "sparq_") & !nms %in% c("sparq_test_methods",
        "sparq_install_methods")])
    defs <- lapply(nms, function(n) {
        f <- get(n, e, inherits = FALSE)
        if (is.function(f))
            list(name = n, arguments = formals(f), body = body(f))
        else NULL
    })
    sparq_fingerprint(defs)
}
sparq_export_benchmark <-
function (result, output_dir)
{
    if (dir.exists(output_dir) && length(list.files(output_dir,
        all.files = TRUE, no.. = TRUE)))
        stop("Use a new, empty output directory; existing results will not be overwritten.")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    keys <- c("replicate_results", "replicate_metrics", "threshold_summary",
        "metric_summary", "training_results", "all_results",
        "comparisons", "threshold_bootstrap", "roc", "precision_recall",
        "benchmark_data", "difficulty_metrics")
    for (k in keys) if (is.data.frame(result[[k]]))
        utils::write.table(result[[k]], file.path(output_dir,
            paste0(k, ".tsv")), sep = "\t", row.names = FALSE,
            quote = FALSE)
    if (!is.null(result$reliability)) {
        r <- result$reliability
        utils::write.table(r$curve, file.path(output_dir, "reliability_bins.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
        summary <- data.frame(spearman = r$monotonicity_spearman,
            p_value = r$monotonicity_p_value, median_scale_discrepancy = r$median_scale_discrepancy,
            n_samples = r$n_samples, n_observations = r$n_observations,
            bootstrap_B = r$bootstrap_B, permutation_B = r$permutation_B)
        utils::write.table(summary, file.path(output_dir, "reliability_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
    }
    sparq_write_design(result$design, file.path(output_dir, "design.yml"))
    saveRDS(result, file.path(output_dir, "figure3_benchmark.rds"))
    writeLines(capture.output(result$session_info), file.path(output_dir,
        "sessionInfo.txt"))
    writeLines(result$input_fingerprint, file.path(output_dir,
        "input_fingerprint.txt"))
    writeLines(result$implementation_fingerprint, file.path(output_dir,
        "implementation_fingerprint.txt"))
    invisible(normalizePath(output_dir))
}
sparq_export_stress <-
function (result, output_dir)
{
    if (dir.exists(output_dir) && length(list.files(output_dir,
        all.files = TRUE, no.. = TRUE)))
        stop("Use a new, empty output directory.")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    for (k in c("sample_results", "replicate_metrics", "summary",
        "paired_differences", "thresholds", "comparisons", "failures")) utils::write.table(result[[k]],
        file.path(output_dir, paste0(k, ".tsv")), sep = "\t",
        row.names = FALSE, quote = FALSE)
    sparq_write_design(result$design, file.path(output_dir, "design.yml"))
    saveRDS(result, file.path(output_dir, "stress_response.rds"))
    writeLines(capture.output(result$session_info), file.path(output_dir,
        "sessionInfo.txt"))
    writeLines(c(result$input_fingerprint, result$implementation_fingerprint),
        file.path(output_dir, "fingerprints.txt"))
    invisible(normalizePath(output_dir))
}
sparq_finalize_support <-
function (support_table, support_threshold)
{
    sparq_check_number(support_threshold, "support_threshold")
    if (!is.data.frame(support_table) || !"support_quality" %in%
        names(support_table))
        stop("support_table must contain support_quality.")
    out <- support_table
    out$threshold <- support_threshold
    out$predicted_stable <- ifelse(is.finite(out$support_quality),
        as.integer(out$support_quality >= support_threshold),
        NA_integer_)
    out$decision_rule <- "support_quality >= training threshold; instability intervals are not support-score intervals"
    out
}
sparq_fingerprint <-
function (object)
{
    f <- tempfile()
    on.exit(unlink(f))
    saveRDS(object, f, version = 2)
    unname(tools::md5sum(f))
}
sparq_known_truth_benchmark <-
function (benchmark_data, design = sparq_methods_design(), n_workers = 1,
    checkpoint_dir = NULL, resume = TRUE, keep_comparisons = TRUE,
    compute_reliability = TRUE)
{
    b <- sparq_validate_benchmark(benchmark_data, design)
    d <- design
    expected_held <- sum(vapply(split(b$truth, b$truth), function(x) length(x) -
        floor(length(x) * d$train_fraction), numeric(1)))
    cat("Expected:", nrow(b), "samples; 1 scenario;", d$n_iterations,
        "iterations/sample;", d$n_replicates, "replicates;",
        expected_held * d$n_replicates, "held-out rows;", nrow(b) *
            d$n_replicates * d$n_iterations, "comparison rows.\n")
    key <- sparq_fingerprint(list(b, d, keep_comparisons, implementation = sparq_engine_signature()))
    outputs <- sparq_parallel_map(seq_len(d$n_replicates), function(r) {
        path <- if (!is.null(checkpoint_dir))
            file.path(checkpoint_dir, key, paste0("replicate_",
                r, ".rds"))
        else NULL
        if (resume && !is.null(path) && file.exists(path))
            return(readRDS(path))
        z <- sparq_benchmark_replicate(b, d, r, keep_comparisons)
        if (!is.null(path))
            sparq_atomic_rds(z, path)
        z
    }, n_workers)
    all <- sparq_bind_rows(lapply(outputs, `[[`, "rows"))
    met <- sparq_bind_rows(lapply(outputs, `[[`, "metrics"))
    held <- all[all$split == "held_out", , drop = FALSE]
    if (nrow(held) != expected_held * d$n_replicates || nrow(met) !=
        d$n_replicates || anyDuplicated(paste(all$replicate,
        all$benchmark_id)))
        stop("Benchmark output count audit failed.")
    curves <- lapply(split(held, held$replicate), function(z) {
        a <- sparq_curves(z$truth, z$support_quality)
        a$roc$replicate <- z$replicate[1]
        a$precision_recall$replicate <- z$replicate[1]
        a
    })
    difficulty_metrics <- sparq_bind_rows(lapply(split(held,
        interaction(held$replicate, held$difficulty, drop = TRUE)),
        function(z) data.frame(replicate = z$replicate[1], difficulty = z$difficulty[1],
            n = nrow(z), correctly_classified = sum(z$predicted_stable ==
                z$truth), class_accuracy = mean(z$predicted_stable ==
                z$truth), predicted_unstable_fraction = mean(z$predicted_stable ==
                0))))
    reliability <- if (compute_reliability)
        sparq_reliability_curve(held$support_quality, held$true_instability,
            held$benchmark_id, bins = d$reliability_bins, min_observations_per_bin = d$min_observations_per_bin,
            B = d$reliability_bootstrap_B, permutation_B = d$reliability_permutation_B,
            confidence = d$confidence, seed = sparq_seed(d$seed,
                0, 14))
    else NULL
    cat("Actual:", length(unique(all$benchmark_id)), "samples; 1 scenario;",
        paste(range(all$n_iterations), collapse = "-"), "completed iterations/sample;",
        nrow(held), "held-out rows; missing samples:", length(setdiff(b$benchmark_id,
            all$benchmark_id)), "\n")
    list(replicate_results = held, replicate_metrics = met, threshold_summary = sparq_threshold_summary(met$threshold,
        d$threshold_bootstrap_B, d$confidence, d$seed), metric_summary = sparq_metric_summary(met,
        d$confidence), difficulty_metrics = difficulty_metrics,
        training_results = all[all$split == "training", , drop = FALSE],
        all_results = all, comparisons = sparq_bind_rows(lapply(outputs,
            `[[`, "comparisons")), threshold_bootstrap = sparq_bind_rows(lapply(outputs,
            `[[`, "threshold_bootstrap")), held_out_results = held,
        held_out_metrics = met, benchmark_data = b, reliability = reliability,
        roc = sparq_bind_rows(lapply(curves, `[[`, "roc")), precision_recall = sparq_bind_rows(lapply(curves,
            `[[`, "precision_recall")), design = d, input_fingerprint = key,
        implementation_fingerprint = sparq_engine_signature(),
        session_info = utils::sessionInfo())
}
sparq_known_truth_benchmark_parallel <-
function (benchmark_data, design = sparq_methods_design(), n_workers = 1,
    ...)
{
    sparq_known_truth_benchmark(benchmark_data, design, n_workers = n_workers,
        ...)
}
sparq_known_truth_stress <-
function (benchmark_data, design = sparq_methods_design(), calibration_results = NULL,
    n_workers = 1, checkpoint_dir = NULL)
{
    b <- sparq_validate_benchmark(benchmark_data, design)
    if (!is.null(calibration_results)) {
        original <- calibration_results$benchmark_data
        fields <- c("benchmark_id", "difficulty", "true_instability",
            "truth")
        if (!is.data.frame(original) || !all(fields %in% names(original)) ||
            !isTRUE(all.equal(b[, fields], original[match(b$benchmark_id,
                original$benchmark_id), fields], check.attributes = FALSE)))
            stop("The calibration benchmark table differs from the stress benchmark table.")
    }
    factory <- function(r) setNames(sparq_make_observations(b,
        design, r), b$benchmark_id)
    sparq_stress_response(factory(1), sparq_benchmark_mean, sparq_compare_scalar,
        truth = setNames(b$truth, b$benchmark_id), design = design,
        n_workers = n_workers, checkpoint_dir = checkpoint_dir,
        cache_key = sparq_fingerprint(list(b, design)), calibration_results = calibration_results,
        replicate_factory = factory)
}
sparq_make_observations <-
function (benchmark_data, design, replicate_id)
{
    b <- sparq_validate_benchmark(benchmark_data, design)
    ids <- sort(b$benchmark_id)
    lapply(seq_len(nrow(b)), function(i) {
        idx <- match(b$benchmark_id[i], ids)
        sample_seed <- sparq_seed(design$seed + replicate_id,
            idx, 10)
        sd <- max(design$observation_sd_floor, design$observation_sd_multiplier *
            b$true_instability[i])
        data.frame(observation_id = seq_len(design$n_observations),
            value = sparq_with_seed(sample_seed, stats::rnorm(design$n_observations,
                design$reference_mean, sd)))
    })
}
sparq_methods_design <-
function ()
{
    list(protocol_version = "methods-v1", seed = 20260912L, n_samples = 1200L,
        n_per_class = 300L, n_observations = 50L, reference_mean = 1,
        observation_sd_multiplier = 1, observation_sd_floor = 0,
        n_iterations = 100L, retention = 0.75, n_replicates = 20L,
        train_fraction = 0.69999999999999996, reference_scale = 1,
        min_iterations = 50L, min_reproducibility = 0.5, max_relative_bias = 0.5,
        max_relative_uncertainty = 0.25, threshold_bootstrap_B = 2000L,
        instability_bootstrap_B = 2000L, reliability_bootstrap_B = 2000L,
        reliability_permutation_B = 2000L, confidence = 0.94999999999999996,
        reliability_bins = 10L, min_observations_per_bin = 10L,
        stress_retentions = c(1, 0.75, 0.5, 0.25, 0.10000000000000001),
        stress_calibration_retention = 0.75, stress_bootstrap_B = 2000L)
}
sparq_metric_summary <-
function (metrics, confidence = 0.94999999999999996)
{
    keys <- intersect(c("sensitivity", "specificity", "balanced_accuracy",
        "auroc", "auprc"), names(metrics))
    sparq_bind_rows(lapply(keys, function(k) {
        x <- metrics[[k]]
        if (any(!is.finite(x)))
            stop("Nonfinite metric: ", k)
        q <- stats::quantile(x, c((1 - confidence)/2, 1 - (1 -
            confidence)/2))
        data.frame(metric = k, n_replicates = length(x), mean = mean(x),
            sd = stats::sd(x), median = stats::median(x), replicate_interval_low = unname(q[1]),
            replicate_interval_high = unname(q[2]), interval_type = "central empirical replicate interval; not a confidence interval for the mean")
    }))
}
sparq_parallel_map <-
function (x, fun, n_workers = 1)
{
    sparq_check_number(n_workers, "n_workers", 1, integer = TRUE)
    safe <- function(i) tryCatch(list(ok = TRUE, value = fun(i)),
        error = function(e) list(ok = FALSE, error = conditionMessage(e)))
    if (n_workers > 1L && .Platform$OS.type != "unix")
        stop("Multiple workers require Unix/Linux/macOS; use n_workers=1 on Windows.")
    out <- if (n_workers == 1L)
        lapply(x, safe)
    else parallel::mclapply(x, safe, mc.cores = min(n_workers,
        length(x)), mc.preschedule = FALSE, mc.set.seed = FALSE)
    bad <- which(!vapply(out, function(z) is.list(z) && isTRUE(z$ok),
        logical(1)))
    if (length(bad)) {
        messages <- vapply(bad, function(i) {
            z <- out[[i]]
            detail <- if (is.list(z) && !is.null(z$error))
                z$error
            else paste(z, collapse = " ")
            if (!nzchar(detail))
                detail <- "worker terminated without returning a result"
            paste0("Task ", x[[i]], ": ", detail)
        }, character(1))
        stop(paste(messages, collapse = "\n"), call. = FALSE)
    }
    lapply(out, `[[`, "value")
}
sparq_read_design <-
function (path)
{
    lines <- readLines(path, warn = FALSE)
    lines <- trimws(lines)
    lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
    decode <- function(s) {
        if (grepl("^\"[^\"\\\\]*\"$", s))
            return(sub("\"$", "", sub("^\"", "", s)))
        if (grepl("^\\[.*\\]$", s)) {
            vals <- strsplit(sub("\\]$", "", sub("^\\[", "",
                s)), ",", fixed = TRUE)[[1]]
            ans <- suppressWarnings(as.numeric(trimws(vals)))
            if (!length(ans) || anyNA(ans))
                stop("Invalid numeric list in design.")
            return(ans)
        }
        v <- suppressWarnings(as.numeric(s))
        if (length(v) != 1L || !is.finite(v))
            stop("Unsupported design value: ", s)
        v
    }
    keys <- sub(":.*$", "", lines)
    if (anyDuplicated(keys) || any(!grepl("^[a-zA-Z_][a-zA-Z_0-9]*:",
        lines)))
        stop("Invalid flat YAML design.")
    d <- setNames(lapply(sub("^[^:]+:[[:space:]]*", "", lines),
        decode), keys)
    sparq_validate_design(d)
    d
}
sparq_reliability_curve <-
function (support_quality, true_error, benchmark_id = NULL, bins = 10,
    min_observations_per_bin = 10, B = 2000, seed = 1, confidence = 0.94999999999999996,
    permutation_B = 2000)
{
    if (length(support_quality) != length(true_error) || any(!is.finite(support_quality)) ||
        any(support_quality < 0 | support_quality > 1) || any(!is.finite(true_error)) ||
        any(true_error < 0))
        stop("Reliability requires paired finite support scores and nonnegative latent instability.")
    if (is.null(benchmark_id))
        stop("benchmark_id is required for sample-level resampling.")
    if (length(benchmark_id) != length(true_error) || anyNA(benchmark_id))
        stop("Invalid benchmark_id.")
    sparq_check_number(bins, "bins", 1, integer = TRUE)
    sparq_check_number(min_observations_per_bin, "min_observations_per_bin",
        2, integer = TRUE)
    sparq_check_number(B, "B", 1, integer = TRUE)
    sparq_check_number(permutation_B, "permutation_B", 1, integer = TRUE)
    sparq_check_number(confidence, "confidence", .Machine$double.eps,
        1 - .Machine$double.eps)
    ids <- sort(unique(as.character(benchmark_id)))
    cluster <- match(benchmark_id, ids)
    truth_by_id <- vapply(seq_along(ids), function(i) {
        z <- unique(true_error[cluster == i])
        if (length(z) != 1L)
            stop("Latent instability must be constant within benchmark_id.")
        z
    }, numeric(1))
    nb <- min(bins, floor(length(support_quality)/min_observations_per_bin))
    if (nb < 1L || length(ids) < 3L)
        stop("Too few samples for reliability analysis.")
    o <- order(support_quality, as.character(benchmark_id), seq_along(support_quality))
    bin_id <- integer(length(o))
    bin_id[o] <- pmin(nb, ceiling(seq_along(o) * nb/length(o)))
    curve <- sparq_bind_rows(lapply(seq_len(nb), function(b) {
        ix <- which(bin_id == b)
        e <- true_error[ix]
        s <- support_quality[ix]
        data.frame(bin = b, n = length(ix), n_samples = length(unique(cluster[ix])),
            support_min = min(s), support_max = max(s), support_median = stats::median(s),
            true_error_median = stats::median(e), q25_true_error = unname(stats::quantile(e,
                0.25)), q75_true_error = unname(stats::quantile(e,
                0.75)))
    }))
    draws <- sparq_with_seed(seed, replicate(B, {
        w <- tabulate(sample.int(length(ids), length(ids), TRUE),
            nbins = length(ids))[cluster]
        vapply(seq_len(nb), function(b) sparq_weighted_median(true_error[bin_id ==
            b], w[bin_id == b]), numeric(1))
    }))
    draws <- matrix(draws, nrow = nb)
    qs <- t(vapply(seq_len(nb), function(b) {
        z <- draws[b, ]
        z <- z[is.finite(z)]
        if (!length(z))
            c(NA_real_, NA_real_)
        else unname(stats::quantile(z, c((1 - confidence)/2,
            1 - (1 - confidence)/2)))
    }, numeric(2)))
    curve$true_error_ci_low <- qs[, 1]
    curve$true_error_ci_high <- qs[, 2]
    curve$n_valid_bootstrap <- rowSums(is.finite(draws))
    variable <- length(unique(support_quality)) > 1L && length(unique(true_error)) >
        1L
    rho <- if (variable)
        stats::cor(support_quality, true_error, method = "spearman")
    else NA_real_
    p <- NA_real_
    if (variable) {
        perm <- sparq_with_seed(sparq_seed(seed, 0, 7), replicate(permutation_B,
            stats::cor(support_quality, truth_by_id[sample.int(length(ids))][cluster],
                method = "spearman")))
        p <- (1 + sum(abs(perm) >= abs(rho) - 9.9999999999999998e-13))/(permutation_B +
            1)
    }
    discrepancy <- abs(stats::median(true_error) - (1 - stats::median(support_quality)))
    list(curve = curve, assignments = data.frame(row = seq_along(bin_id),
        benchmark_id = benchmark_id, bin = bin_id), monotonicity_spearman = rho,
        monotonicity_p_value = p, overall_calibration_error = discrepancy,
        median_scale_discrepancy = discrepancy, n_observations = length(support_quality),
        n_samples = length(ids), n_bins = nb, bootstrap_B = B,
        confidence = confidence, permutation_B = permutation_B,
        bootstrap_unit = "benchmark_id; all held-out rows of a resampled ID travel together; fixed original bins",
        p_value_method = "two-sided Monte Carlo permutation of latent instability across benchmark IDs",
        interpretation = "Median scale discrepancy is descriptive, not probability calibration error.")
}
sparq_reliability_decision <-
function (comparison_table, tolerance, required_reliability = 0.9,
    confidence = 0.95, reference_scale = 1)
{
    if (!is.data.frame(comparison_table) || !nrow(comparison_table) ||
        !all(c("status", "instability") %in% names(comparison_table)))
        stop("comparison_table must contain attempted perturbations with status and instability.",
            call. = FALSE)
    sparq_check_number(tolerance, "tolerance", 0)
    sparq_check_number(required_reliability, "required_reliability",
        .Machine$double.eps, 1 - .Machine$double.eps)
    sparq_check_number(confidence, "confidence", .Machine$double.eps,
        1 - .Machine$double.eps)
    sparq_check_number(reference_scale, "reference_scale")
    scale <- max(abs(reference_scale), 9.9999999999999995e-07)
    n_attempted <- nrow(comparison_table)
    distance <- suppressWarnings(as.numeric(comparison_table$instability))
    completed <- comparison_table$status == "ok" & is.finite(distance)
    preserved <- completed & distance/scale <= tolerance
    n_preserved <- sum(preserved)
    alpha <- 1 - confidence
    lower <- if (n_preserved == 0L)
        0
    else stats::qbeta(alpha/2, n_preserved, n_attempted - n_preserved +
        1L)
    upper <- if (n_preserved == n_attempted)
        1
    else stats::qbeta(1 - alpha/2, n_preserved + 1L, n_attempted -
        n_preserved)
    decision <- if (lower >= required_reliability)
        "supported"
    else if (upper < required_reliability)
        "not_supported"
    else "inconclusive"
    data.frame(n_attempted = n_attempted, n_completed = sum(completed),
        n_failed = n_attempted - sum(completed), n_preserved = n_preserved,
        preservation_rate = n_preserved/n_attempted,
        lower_confidence_bound = lower, upper_confidence_bound = upper,
        tolerance = tolerance, required_reliability = required_reliability,
        confidence = confidence, reference_scale = scale,
        decision = decision, stringsAsFactors = FALSE)
}
sparq_run <-
function (data, analysis_function, comparator, stress_model = "uniform_random",
    retention = NULL, n_iterations = NULL, result_id = "result",
    x_col = NULL, y_col = NULL, custom_function = NULL, reference_scale = 1,
    seed = 1, cache_dir = NULL, resume = TRUE, cache_key = NULL,
    failure_action = c("record", "stop"), full_result = NULL,
    keep_perturbed_results = FALSE, verbose = NULL, progress_every = NULL,
    preset = "standard", min_iterations = NULL,
    instability_bootstrap_B = NULL, reliability_tolerance = NULL,
    required_reliability = 0.9, reliability_confidence = 0.95,
    loss_model = NULL, ...)
{
    failure_action <- match.arg(failure_action)
    preset_values <- sparq_resolve_preset(
        preset = preset, retention = retention, n_iterations = n_iterations,
        min_iterations = min_iterations,
        instability_bootstrap_B = instability_bootstrap_B,
        progress_every = progress_every, verbose = verbose
    )
    retention <- preset_values$retention
    n_iterations <- preset_values$n_iterations
    min_iterations <- preset_values$min_iterations
    instability_bootstrap_B <- preset_values$instability_bootstrap_B
    progress_every <- preset_values$progress_every
    verbose <- preset_values$verbose
    sparq_inform_preset(preset_values)
    if (!is.function(analysis_function) || !is.function(comparator)) {
        stop("analysis_function and comparator must be functions.",
            call. = FALSE)
    }
    if (!is.data.frame(data) || !nrow(data)) {
        stop("data must be a nonempty data.frame.", call. = FALSE)
    }
    input_provenance <- sparq_input_provenance(data)
    if (!is.null(loss_model) && !inherits(loss_model, "sparq_loss_model")) {
        stop("loss_model must be created by sparq_define_loss_model().",
            call. = FALSE)
    }
    sparq_inform(
        paste0(
            "Starting ", result_id, ": ", nrow(data), " observations; ",
            n_iterations, " perturbations; ", round(retention * 100),
            "% retention; ", stress_model, " stress."
        ),
        verbose
    )
    if (is.null(full_result)) {
        sparq_inform("Running reference analysis.", verbose)
        full_result <- sparq_with_seed(sparq_seed(seed, 0, 1),
            analysis_function(data))
    }
    cache_prefix <- NULL
    if (!is.null(cache_dir)) {
        if (length(cache_key) != 1L || is.na(cache_key) || !nzchar(cache_key)) {
            stop(paste("Caching requires cache_key identifying analysis code,",
                "parameters, and external dependencies."), call. = FALSE)
        }
        cache_identity <- list(data = data, full_result = full_result,
            analysis = body(analysis_function), comparator = body(comparator),
            custom = if (is.function(custom_function)) {
                body(custom_function)
            } else {
                NULL
            }, stress_model = stress_model, retention = retention,
            seed = seed, x_col = x_col, y_col = y_col, cache_key = cache_key,
            loss_model = loss_model,
            keep_perturbed_results = keep_perturbed_results,
            implementation = sparq_engine_signature())
        cache_prefix <- file.path(cache_dir, sparq_fingerprint(cache_identity))
        dir.create(cache_prefix, recursive = TRUE, showWarnings = FALSE)
    }
    iteration_output <- lapply(seq_len(n_iterations), function(iteration) {
        comparison_path <- if (!is.null(cache_prefix)) {
            file.path(cache_prefix, paste0(iteration, ".rds"))
        }
        else {
            NULL
        }
        result_path <- if (!is.null(cache_prefix)) {
            file.path(cache_prefix, paste0(iteration, "_perturbed_result.rds"))
        }
        else {
            NULL
        }
        loss_audit_path <- if (!is.null(cache_prefix) && !is.null(loss_model)) {
            file.path(cache_prefix, paste0(iteration, "_loss_audit.rds"))
        }
        else {
            NULL
        }
        if (resume && !is.null(comparison_path) && file.exists(comparison_path)) {
            cached_comparison <- readRDS(comparison_path)
            cached_ok <- (is.data.frame(cached_comparison) &&
                nrow(cached_comparison) == 1L && identical(cached_comparison$status,
                "ok"))
            cached_result_available <- (!keep_perturbed_results ||
                file.exists(result_path))
            cached_audit_available <- (is.null(loss_audit_path) ||
                file.exists(loss_audit_path))
            if (cached_ok && cached_result_available && cached_audit_available) {
                if (iteration %% progress_every == 0L || iteration == n_iterations) {
                    sparq_progress(iteration, n_iterations, verbose)
                }
                return(list(comparison = cached_comparison, perturbed_result = if (keep_perturbed_results) {
                  readRDS(result_path)
                } else {
                  NULL
                }, loss_audit = if (!is.null(loss_audit_path)) {
                    readRDS(loss_audit_path)
                } else {
                    NULL
                }))
            }
        }
        output <- tryCatch(sparq_with_seed(sparq_seed(seed, iteration,
            2), {
            perturbed_data <- sparq_stress_model(data = data,
                retention = retention, model = stress_model,
                x_col = x_col, y_col = y_col, custom_function = custom_function,
                iteration = iteration, seed = seed, loss_model = loss_model)
            perturbed_result <- analysis_function(perturbed_data)
            comparison <- sparq_comparison_row(comparator = comparator,
                full_result = full_result, perturbed_result = perturbed_result,
                iteration = iteration)
                comparison$n_retained <- nrow(perturbed_data)
            comparison$loss_model <- if (is.null(loss_model)) stress_model else loss_model$type
            list(comparison = comparison, perturbed_result = if (keep_perturbed_results) {
                perturbed_result
            } else {
                NULL
            }, loss_audit = attr(perturbed_data, "sparq_loss_audit"))
        }), error = function(error) {
            if (failure_action == "stop") {
                stop(paste0("Iteration ", iteration, ": ", conditionMessage(error)),
                  call. = FALSE)
            }
            list(comparison = data.frame(instability = NA_real_,
                iteration = iteration, status = "failed", error_message = conditionMessage(error),
                n_retained = NA_integer_, stringsAsFactors = FALSE),
                perturbed_result = NULL)
        })
        if (!is.null(comparison_path)) {
            sparq_atomic_rds(output$comparison, comparison_path)
            if (keep_perturbed_results && identical(output$comparison$status,
                "ok")) {
                sparq_atomic_rds(output$perturbed_result, result_path)
            }
            if (!is.null(loss_audit_path) && !is.null(output$loss_audit)) {
                sparq_atomic_rds(output$loss_audit, loss_audit_path)
            }
        }
        if (iteration %% progress_every == 0L || iteration == n_iterations) {
            sparq_progress(iteration, n_iterations, verbose)
        }
        output
    })
    comparisons <- sparq_bind_rows(lapply(iteration_output, function(x) {
        x$comparison
    }))
    loss_audits <- lapply(iteration_output, function(x) x$loss_audit)
    if (all(vapply(loss_audits, is.null, logical(1)))) {
        loss_audits <- NULL
    } else {
        names(loss_audits) <- paste0("iteration_", seq_len(n_iterations))
    }
    summary <- sparq_support_from_comparisons(comparison_table = comparisons,
        n_iterations = n_iterations, reference_scale = reference_scale,
        min_iterations = min_iterations,
        instability_bootstrap_B = instability_bootstrap_B, seed = seed, ...)
    summary$result_id <- result_id
    summary$stress_model <- stress_model
    summary$retention <- retention
    summary$assessment_mode <- if (stress_model == "none") {
        "no_perturbation_check"
    }
    else {
        "automatic_rerun"
    }
    perturbed_results <- NULL
    if (keep_perturbed_results) {
        perturbed_results <- lapply(iteration_output, function(x) {
            x$perturbed_result
        })
        names(perturbed_results) <- paste0("iteration_", seq_len(n_iterations))
    }
    out <- list(result_id = result_id, full_result = full_result, perturbed_results = perturbed_results,
        loss_audits = loss_audits,
        comparisons = comparisons, summary = summary, settings = list(stress_model = stress_model,
            retention = retention, n_iterations = n_iterations,
            reference_scale = reference_scale, seed = seed, keep_perturbed_results = keep_perturbed_results,
            progress_every = progress_every, preset = preset_values$preset,
            min_iterations = min_iterations,
            instability_bootstrap_B = instability_bootstrap_B,
            loss_model = loss_model),
        input_provenance = input_provenance,
        failures = comparisons[comparisons$status != "ok", ,
            drop = FALSE])
    if (!is.null(reliability_tolerance)) {
        reliability <- sparq_reliability_decision(comparisons,
            tolerance = reliability_tolerance,
            required_reliability = required_reliability,
            confidence = reliability_confidence,
            reference_scale = reference_scale)
        out$reliability_decision <- reliability
        for (name in names(reliability)) {
            out$summary[[paste0("reliability_", name)]] <- reliability[[name]][1]
        }
        out$settings$reliability_tolerance <- reliability_tolerance
        out$settings$required_reliability <- required_reliability
        out$settings$reliability_confidence <- reliability_confidence
    }
    class(out) <- unique(c("sparq_assessment", class(out)))
    sparq_inform(
        paste0(
            "Completed ", result_id, ": support quality ",
            format(summary$support_quality[1], digits = 4), "; ",
            sparq_assessment_status(summary), "."
        ),
        verbose
    )
    out
}
sparq_run_with_stress_model <-
function (data, analysis_function, comparator, stress_model = "uniform_random",
    retention = 0.75, n_iterations = 100, x_col = NULL, y_col = NULL,
    custom_function = NULL, result_id = "result", reference_scale = 1,
    seed = 1, ...)
{
    sparq_run(data, analysis_function, comparator, stress_model,
        retention, n_iterations, result_id, x_col, y_col, custom_function,
        reference_scale, seed, ...)$summary
}
sparq_seed <-
function (seed, index, stream = 0)
{
    as.integer((48271 * as.double(seed) + 104729 * as.double(index) +
        1000003 * stream)%%2147483647)
}
sparq_split_benchmark <-
function (b, train_fraction, seed)
{
    sparq_with_seed(seed, {
        train <- unlist(lapply(split(seq_len(nrow(b)), b$truth),
            function(ix) {
                ix <- ix[order(b$benchmark_id[ix])]
                n <- floor(length(ix) * train_fraction)
                if (n < 2L || n >= length(ix))
                  stop("Split must leave at least two training cases and one held-out case per class.")
                ix[sample.int(length(ix), n)]
            }), use.names = FALSE)
        list(training = sort(train), held_out = setdiff(seq_len(nrow(b)),
            train))
    })
}
sparq_stress_model <-
function (data, retention = 0.75, model = c("uniform_random",
    "contiguous_hole", "none", "custom"), x_col = NULL, y_col = NULL,
    custom_function = NULL, iteration = 1, seed = NULL, loss_model = NULL)
{
    model <- match.arg(model)
    if (!is.data.frame(data) || !nrow(data))
        stop("data must be a nonempty data.frame.")
    sparq_check_number(retention, "retention", .Machine$double.eps,
        1)
    sparq_check_number(iteration, "iteration", 1, integer = TRUE)
    if (!is.null(loss_model)) {
        if (!inherits(loss_model, "sparq_loss_model"))
            stop("loss_model must be created by sparq_define_loss_model().",
                call. = FALSE)
        applied <- sparq_apply_loss_model(data = data,
            loss_model = loss_model, retention = retention, x_col = x_col,
            y_col = y_col, iteration = iteration, seed = seed)
        perturbed <- applied$data
        attr(perturbed, "sparq_loss_audit") <- applied$audit
        attr(perturbed, "sparq_loss_model") <- loss_model
        return(perturbed)
    }
    generate <- function() {
        if (model == "none")
            return(data)
        if (model == "custom") {
            if (!is.function(custom_function))
                stop("custom_function must be a function.")
            ans <- custom_function(data = data, retention = retention,
                iteration = iteration)
            if (!is.data.frame(ans) || !nrow(ans))
                stop("custom_function returned no data.")
            return(ans)
        }
        n <- nrow(data)
        n_keep <- min(n, max(2L, floor(n * retention)))
        if (model == "uniform_random") {
            if (n_keep == n)
                return(data)
            return(data[sample.int(n, n_keep), , drop = FALSE])
        }
        if (length(x_col) != 1L || length(y_col) != 1L || !all(c(x_col,
            y_col) %in% names(data)))
            stop("Valid x_col and y_col are required.")
        xy <- data[, c(x_col, y_col), drop = FALSE]
        if (!all(vapply(xy, is.numeric, logical(1))) || any(!is.finite(as.matrix(xy))))
            stop("All spatial coordinates must be finite numeric values; none are silently removed.")
        if (n_keep == n)
            return(data)
        center <- sample.int(n, 1L)
        d2 <- (xy[[1]] - xy[[1]][center])^2 + (xy[[2]] - xy[[2]][center])^2
        data[order(-d2, seq_len(n))[seq_len(n_keep)], , drop = FALSE]
    }
    if (is.null(seed))
        generate()
    else sparq_with_seed(sparq_seed(seed, iteration), generate())
}
sparq_stress_response <-
function (sample_data, analysis_function, comparator, truth = NULL,
    design = sparq_methods_design(), stress_models = "uniform_random",
    x_col = NULL, y_col = NULL, n_workers = 1, checkpoint_dir = NULL,
    resume = TRUE, cache_key = NULL, calibration_results = NULL,
    replicate_factory = NULL)
{
    sparq_validate_design(design)
    d <- design
    if (!is.list(sample_data) || !length(sample_data) || is.null(names(sample_data)) ||
        anyNA(names(sample_data)) || anyDuplicated(names(sample_data)))
        stop("sample_data must be a named list with unique sample IDs.")
    ids <- names(sample_data)
    if (!all(stress_models %in% c("uniform_random", "contiguous_hole")) ||
        anyDuplicated(stress_models))
        stop("Supported stress models are uniform_random and contiguous_hole.")
    if ("contiguous_hole" %in% stress_models)
        for (dat in sample_data) sparq_stress_model(dat, retention = d$retention,
            model = "contiguous_hole", x_col = x_col, y_col = y_col,
            seed = d$seed)
    labels <- !is.null(truth)
    if (labels) {
        if (is.null(names(truth)) || anyDuplicated(names(truth)) ||
            !setequal(names(truth), ids))
            stop("truth must be named by sample ID.")
        truth <- truth[ids]
        if (anyNA(truth) || !all(truth %in% c(0, 1)) || length(unique(truth)) !=
            2L)
            stop("Invalid truth labels.")
    }
    if (!is.null(checkpoint_dir) && (length(cache_key) != 1L ||
        !nzchar(cache_key)))
        stop("Checkpointing stress runs requires cache_key identifying the complete analysis and inputs.")
    key <- sparq_fingerprint(list(sample_data, truth, d, stress_models,
        x_col, y_col, cache_key, calibration_results, body(analysis_function),
        body(comparator), if (is.function(replicate_factory)) body(replicate_factory),
        sparq_engine_signature()))
    nr <- length(d$stress_retentions)
    nm <- length(stress_models)
    cat("Expected:", length(ids), "samples;", nr * nm, "stress conditions;",
        d$n_iterations, "iterations/sample;", d$n_replicates,
        "paired replicates;", length(ids) * nr * nm * d$n_replicates,
        "sample-condition summaries.\n")
    outputs <- sparq_parallel_map(seq_len(d$n_replicates), function(r) {
        path <- if (!is.null(checkpoint_dir))
            file.path(checkpoint_dir, key, paste0("stress_",
                r, ".rds"))
        else NULL
        if (resume && !is.null(path) && file.exists(path))
            return(readRDS(path))
        dats <- if (is.null(replicate_factory))
            sample_data
        else replicate_factory(r)
        if (!identical(names(dats), ids))
            stop("replicate_factory changed sample IDs or their order.")
        split <- if (labels)
            sparq_split_benchmark(data.frame(benchmark_id = ids,
                truth = truth), d$train_fraction, sparq_seed(d$seed +
                r, 0, 11))
        else list(training = integer(), held_out = seq_along(ids))
        summaries <- list()
        comparisons <- list()
        thresholds <- list()
        references <- lapply(seq_along(ids), function(i) sparq_with_seed(sparq_seed(sparq_seed(d$seed +
            r, match(ids[i], sort(ids)), 12), 0, 1), analysis_function(dats[[i]])))
        for (model in stress_models) {
            run_level <- function(retention) {
                fits <- lapply(seq_along(ids), function(i) sparq_run(dats[[i]],
                  analysis_function, comparator, model, retention,
                  d$n_iterations, ids[i], x_col, y_col, reference_scale = d$reference_scale,
                  full_result = references[[i]], seed = sparq_seed(d$seed +
                    r, match(ids[i], sort(ids)), 12), min_iterations = d$min_iterations,
                  min_reproducibility = d$min_reproducibility,
                  max_relative_bias = d$max_relative_bias, max_relative_uncertainty = d$max_relative_uncertainty,
                  instability_bootstrap_B = d$instability_bootstrap_B,
                  confidence = d$confidence))
                rows <- sparq_bind_rows(lapply(seq_along(fits),
                  function(i) {
                    z <- fits[[i]]$summary
                    z$benchmark_id <- ids[i]
                    z$replicate <- r
                    z$seed <- d$seed + r
                    z$split <- if (i %in% split$training)
                      "training"
                    else "held_out"
                    z$truth <- if (labels)
                      truth[i]
                    else NA_integer_
                    z
                  }))
                comp <- sparq_bind_rows(lapply(seq_along(fits),
                  function(i) {
                    z <- fits[[i]]$comparisons
                    z$benchmark_id <- ids[i]
                    z$replicate <- r
                    z$stress_model <- model
                    z$retention <- retention
                    z
                  }))
                list(rows = rows, comparisons = comp)
            }
            base <- NULL
            threshold <- NA_real_
            threshold_low <- NA_real_
            threshold_high <- NA_real_
            origin <- "unlabeled; no classification"
            if (labels) {
                if (!is.null(calibration_results) && model ==
                  "uniform_random") {
                  if (d$retention != d$stress_calibration_retention)
                    stop("C–E retention differs from the stress calibration condition.")
                  src <- calibration_results
                  b <- src$benchmark_data
                  if (!identical(src$design, d) || !setequal(b$benchmark_id,
                    ids) || !identical(as.integer(b$truth[match(ids,
                    b$benchmark_id)]), as.integer(truth)))
                    stop("Calibration result does not match the stress design or sample labels.")
                  rows <- src$all_results[src$all_results$replicate ==
                    r, , drop = FALSE]
                  if (!setequal(rows$benchmark_id[rows$split ==
                    "training"], ids[split$training]))
                    stop("Calibration split mismatch.")
                  fingerprints <- vapply(dats, sparq_fingerprint,
                    character(1))
                  if (!"observation_fingerprint" %in% names(rows) ||
                    !identical(unname(rows$observation_fingerprint[match(ids,
                      rows$benchmark_id)]), unname(fingerprints)))
                    stop("Calibration observations do not match this stress replicate.")
                  m <- src$replicate_metrics[src$replicate_metrics$replicate ==
                    r, , drop = FALSE]
                  if (nrow(m) != 1L || !isTRUE(m$informative_training_scores))
                    stop("Invalid baseline calibration.")
                  threshold <- m$threshold
                  threshold_low <- m$threshold_ci_low
                  threshold_high <- m$threshold_ci_high
                  origin <- paste("reused Figure 3C/D training threshold at retention",
                    d$retention)
                }
                else {
                  base <- run_level(d$stress_calibration_retention)
                  tr <- base$rows[split$training, , drop = FALSE]
                  if (any(tr$n_failed > 0L) || any(!is.finite(tr$support_quality)))
                    stop("Calibration has failed training assessments.")
                  cal <- sparq_calibrate_threshold(tr$support_quality,
                    tr$truth, B = d$threshold_bootstrap_B, confidence = d$confidence,
                    seed = sparq_seed(d$seed + r, 0, 13))
                  if (!cal$informative)
                    stop("Training scores are constant; no informative stress threshold can be learned.")
                  threshold <- cal$threshold
                  threshold_low <- cal$threshold_interval$lower
                  threshold_high <- cal$threshold_interval$upper
                  origin <- paste("training only at retention",
                    d$stress_calibration_retention, "within",
                    model)
                }
            }
            thresholds[[length(thresholds) + 1L]] <- data.frame(replicate = r,
                stress_model = model, threshold = threshold,
                threshold_ci_low = threshold_low, threshold_ci_high = threshold_high,
                origin = origin)
            for (ret in d$stress_retentions) {
                fit <- if (!is.null(base) && ret == d$stress_calibration_retention)
                  base
                else run_level(ret)
                fit$rows$threshold <- threshold
                fit$rows$predicted_stable <- if (labels)
                  as.integer(fit$rows$support_quality >= threshold)
                else NA_integer_
                summaries[[length(summaries) + 1L]] <- fit$rows
                comparisons[[length(comparisons) + 1L]] <- fit$comparisons
            }
        }
        z <- list(rows = sparq_bind_rows(summaries), comparisons = sparq_bind_rows(comparisons),
            thresholds = sparq_bind_rows(thresholds))
        if (!is.null(path))
            sparq_atomic_rds(z, path)
        z
    }, n_workers)
    rows <- sparq_bind_rows(lapply(outputs, `[[`, "rows"))
    comps <- sparq_bind_rows(lapply(outputs, `[[`, "comparisons"))
    groups <- split(rows[rows$split == "held_out", , drop = FALSE],
        interaction(rows$replicate[rows$split == "held_out"],
            rows$stress_model[rows$split == "held_out"], rows$retention[rows$split ==
                "held_out"], drop = TRUE))
    met <- sparq_bind_rows(lapply(groups, function(z) {
        complete <- all(z$n_failed == 0L & is.finite(z$support_quality))
        out <- data.frame(replicate = z$replicate[1], stress_model = z$stress_model[1],
            retention = z$retention[1], threshold = z$threshold[1],
            n_expected = nrow(z), n_complete = sum(z$n_failed ==
                0L & is.finite(z$support_quality)), n_failed_iterations = sum(z$n_failed),
            complete = complete, support_quality = if (complete)
                mean(z$support_quality)
            else NA_real_, instability_median = if (complete)
                mean(z$instability_median)
            else NA_real_, uncertainty = if (complete)
                mean(z$uncertainty)
            else NA_real_)
        if (labels && complete)
            out <- cbind(out, sparq_binary_metrics(z$truth, z$support_quality,
                z$threshold[1]))
        out
    }))
    summaries <- list()
    paired <- list()
    keys <- intersect(c("support_quality", "instability_median",
        "uncertainty", "sensitivity", "specificity", "balanced_accuracy"),
        names(met))
    for (model in stress_models) for (k in keys) {
        matrix_values <- vapply(d$stress_retentions, function(ret) {
            z <- met[met$stress_model == model & met$retention ==
                ret, , drop = FALSE]
            z[[k]][match(seq_len(d$n_replicates), z$replicate)]
        }, numeric(d$n_replicates))
        if (d$n_replicates == 1L)
            matrix_values <- matrix(matrix_values, nrow = 1L)
        for (j in seq_along(d$stress_retentions)) {
            x <- matrix_values[, j]
            ci <- if (all(is.finite(x)))
                sparq_bootstrap_ci(x, mean, d$stress_bootstrap_B,
                  d$confidence, sparq_seed(d$seed, 0, 15))
            else c(estimate = NA_real_, lower = NA_real_, upper = NA_real_)
            summaries[[length(summaries) + 1L]] <- data.frame(stress_model = model,
                retention = d$stress_retentions[j], metric = k,
                mean = unname(ci["estimate"]), ci_low = unname(ci["lower"]),
                ci_high = unname(ci["upper"]), n_replicates = length(x),
                n_complete = sum(is.finite(x)))
            baseline <- match(d$stress_calibration_retention,
                d$stress_retentions)
            if (!is.na(baseline)) {
                delta <- x - matrix_values[, baseline]
                dc <- if (all(is.finite(delta)))
                  sparq_bootstrap_ci(delta, mean, d$stress_bootstrap_B,
                    d$confidence, sparq_seed(d$seed, 0, 15))
                else c(estimate = NA_real_, lower = NA_real_,
                  upper = NA_real_)
                paired[[length(paired) + 1L]] <- data.frame(stress_model = model,
                  retention = d$stress_retentions[j], metric = k,
                  reference_retention = d$stress_calibration_retention,
                  mean_difference = unname(dc["estimate"]), ci_low = unname(dc["lower"]),
                  ci_high = unname(dc["upper"]))
            }
        }
    }
    expected <- length(ids) * nr * nm * d$n_replicates
    if (nrow(rows) != expected || nrow(comps) != expected * d$n_iterations)
        stop("Stress output count audit failed.")
    cat("Actual:", length(unique(rows$benchmark_id)), "samples;",
        nr * nm, "conditions;", nrow(rows), "summaries;", sum(rows$n_failed),
        "failed iterations; missing samples:", length(setdiff(ids,
            rows$benchmark_id)), "\n")
    list(sample_results = rows, replicate_metrics = met, summary = sparq_bind_rows(summaries),
        paired_differences = sparq_bind_rows(paired), thresholds = sparq_bind_rows(lapply(outputs,
            `[[`, "thresholds")), comparisons = comps, failures = comps[comps$status !=
            "ok", , drop = FALSE], design = d, input_fingerprint = key,
        implementation_fingerprint = sparq_engine_signature(),
        session_info = utils::sessionInfo())
}
sparq_support_from_comparisons <-
function (comparison_table, n_iterations = nrow(comparison_table),
    reference_scale = 1, min_iterations = 50, min_reproducibility = 0.5,
    max_relative_bias = 0.5, max_relative_uncertainty = 0.25,
    instability_bootstrap_B = 2000, confidence = 0.94999999999999996,
    seed = 1)
{
    sparq_check_number(n_iterations, "n_iterations", 0, integer = TRUE)
    sparq_check_number(min_iterations, "min_iterations", 1, integer = TRUE)
    sparq_check_number(min_reproducibility, "min_reproducibility",
        .Machine$double.eps, 1)
    sparq_check_number(reference_scale, "reference_scale")
    sparq_check_number(max_relative_bias, "max_relative_bias",
        0)
    sparq_check_number(max_relative_uncertainty, "max_relative_uncertainty",
        0)
    sparq_check_number(instability_bootstrap_B, "instability_bootstrap_B",
        1, integer = TRUE)
    if (!is.data.frame(comparison_table) || !"instability" %in%
        names(comparison_table))
        stop("comparison_table must contain instability.")
    if (!is.numeric(comparison_table$instability))
        stop("instability must be numeric.")
    if ("iteration" %in% names(comparison_table) && (anyNA(comparison_table$iteration) ||
        anyDuplicated(comparison_table$iteration)))
        stop("Comparison iteration IDs must be unique and nonmissing.")
    if (nrow(comparison_table) > n_iterations)
        stop("More comparison rows than requested iterations.")
    valid <- is.finite(comparison_table$instability) & comparison_table$instability >=
        0
    if ("delta" %in% names(comparison_table))
        valid <- valid & is.finite(comparison_table$delta)
    if ("status" %in% names(comparison_table))
        valid <- valid & !is.na(comparison_table$status) & comparison_table$status ==
            "ok"
    d <- comparison_table[valid, , drop = FALSE]
    n <- nrow(d)
    signed <- if ("delta" %in% names(d) && n)
        mean(d$delta)
    else NA_real_
    bias <- if (!n)
        NA_real_
    else if ("delta" %in% names(d))
        abs(signed)
    else stats::median(d$instability)
    uncertainty <- if (n >= 2L)
        stats::sd(d$instability)
    else NA_real_
    scale <- max(abs(reference_scale), 9.9999999999999995e-07)
    med <- if (n)
        stats::median(d$instability)
    else NA_real_
    q <- sparq_bound01(1 - med)
    sb <- 1/(1 + bias/scale)
    su <- 1/(1 + uncertainty/scale)
    si <- sparq_bound01(n/min_iterations)
    sr <- sparq_bound01(q/min_reproducibility)
    score <- si * sb * su * sr
    reasons <- c(if (n < min_iterations) "few_iterations", if (n <
        n_iterations) "failed_iterations", if (n < 2L) "insufficient_finite_comparisons",
        if (is.finite(q) && q < min_reproducibility) "low_reproducibility",
        if (is.finite(bias) && bias/scale > max_relative_bias) "high_bias",
        if (is.finite(uncertainty) && uncertainty/scale > max_relative_uncertainty) "high_uncertainty")
    ci <- sparq_bootstrap_ci(d$instability, B = instability_bootstrap_B,
        conf = confidence, seed = seed)
    data.frame(n_iterations = n, n_requested = n_iterations,
        n_failed = n_iterations - n, signed_displacement = signed,
        bias = bias, uncertainty = uncertainty, reference_scale = scale,
        relative_bias = bias/scale, relative_uncertainty = uncertainty/scale,
        reproducibility_score = q, iteration_support = si, bias_support = sb,
        uncertainty_support = su, reproducibility_support = sr,
        support_quality = score, instability_median = med, instability_ci_low = unname(ci["lower"]),
        instability_ci_high = unname(ci["upper"]), insufficient_support = length(reasons) >
            0, support_reason = if (length(reasons))
            paste(reasons, collapse = ";")
        else "adequate_support", stringsAsFactors = FALSE)
}
sparq_threshold_candidates <-
function (score, truth)
{
    sparq_validate_scores(score, truth)
    if (length(unique(truth)) != 2L)
        stop("Both training classes are required.")
    o <- order(score)
    score <- score[o]
    truth <- truth[o]
    ends <- c(which(diff(score) != 0), length(score))
    u <- score[ends]
    n1 <- sum(truth == 1)
    n0 <- sum(truth == 0)
    below1 <- c(0, cumsum(truth == 1)[ends])
    below0 <- c(0, cumsum(truth == 0)[ends])
    mids <- if (length(u) > 1L)
        (head(u, -1L) + tail(u, -1L))/2
    else numeric()
    upper <- max(u) + max(.Machine$double.eps, abs(max(u)) *
        .Machine$double.eps)
    thresholds <- c(0, mids, upper)
    data.frame(threshold = thresholds, sensitivity = (n1 - below1)/n1,
        specificity = below0/n0, balanced_accuracy = ((n1 - below1)/n1 +
            below0/n0)/2)
}
sparq_threshold_summary <-
function (thresholds, B = 2000, confidence = 0.94999999999999996,
    seed = 1)
{
    if (!length(thresholds) || any(!is.finite(thresholds)))
        stop("Thresholds must be finite.")
    ci <- sparq_bootstrap_ci(thresholds, B = B, conf = confidence,
        seed = seed)
    spread <- stats::quantile(thresholds, c((1 - confidence)/2,
        1 - (1 - confidence)/2))
    data.frame(threshold_median = stats::median(thresholds),
        threshold_ci_low = unname(ci["lower"]), threshold_ci_high = unname(ci["upper"]),
        replicate_interval_low = unname(spread[1]), replicate_interval_high = unname(spread[2]),
        n_replicates = length(thresholds), B = B, confidence = confidence,
        ci_target = "median threshold across computational replicates")
}
sparq_validate_benchmark <-
function (benchmark_data, design)
{
    sparq_validate_design(design)
    need <- c("benchmark_id", "difficulty", "true_instability")
    if (!is.data.frame(benchmark_data) || !all(need %in% names(benchmark_data)))
        stop("Missing Figure 3B columns.")
    if (anyDuplicated(names(benchmark_data)) || any(!names(benchmark_data) %in%
        c(need, "truth")))
        stop("Figure 3B input must contain only benchmark_id, difficulty, true_instability and optional truth; precomputed scores are not accepted.")
    b <- benchmark_data
    b$benchmark_id <- as.character(b$benchmark_id)
    b$difficulty <- as.character(b$difficulty)
    if (anyNA(b$benchmark_id) || any(!nzchar(b$benchmark_id)) ||
        anyDuplicated(b$benchmark_id))
        stop("Invalid benchmark IDs.")
    classes <- c("Stable", "Mild instability", "Moderate instability",
        "Severe instability")
    if (anyNA(b$difficulty) || !setequal(unique(b$difficulty),
        classes) || nrow(b) != design$n_samples || any(table(factor(b$difficulty,
        classes)) != design$n_per_class))
        stop("Figure 3B class counts do not match design.")
    if (!is.numeric(b$true_instability) || any(!is.finite(b$true_instability)) ||
        any(b$true_instability < 0))
        stop("true_instability must be finite and nonnegative.")
    expected <- as.integer(b$difficulty == "Stable")
    if ("truth" %in% names(b) && !identical(as.integer(b$truth),
        expected))
        stop("truth contradicts difficulty labels.")
    b$truth <- expected
    b
}
sparq_validate_design <-
function (d)
{
    required <- names(sparq_methods_design())
    if (!is.list(d) || !setequal(names(d), required))
        stop("Design fields do not match the versioned protocol.")
    if (!identical(d$protocol_version, "methods-v1"))
        stop("Unsupported protocol_version.")
    ints <- c("seed", "n_samples", "n_per_class", "n_observations",
        "n_iterations", "n_replicates", "min_iterations", "threshold_bootstrap_B",
        "instability_bootstrap_B", "reliability_bootstrap_B",
        "reliability_permutation_B", "reliability_bins", "min_observations_per_bin",
        "stress_bootstrap_B")
    for (k in ints) sparq_check_number(d[[k]], k, if (k == "seed")
        0
    else 1, .Machine$integer.max, TRUE)
    if (d$seed + d$n_replicates > .Machine$integer.max)
        stop("Replicate seeds exceed R's seed range.")
    if (d$n_samples != 4 * d$n_per_class || d$n_per_class < 4L ||
        d$n_observations < 2L)
        stop("Design requires four equally sized classes and at least two observations per sample.")
    for (k in c("retention", "min_reproducibility", "stress_calibration_retention")) sparq_check_number(d[[k]],
        k, .Machine$double.eps, 1)
    for (k in c("train_fraction", "confidence")) sparq_check_number(d[[k]],
        k, .Machine$double.eps, 1 - .Machine$double.eps)
    for (k in c("reference_scale", "observation_sd_multiplier")) sparq_check_number(d[[k]],
        k, .Machine$double.eps)
    for (k in c("max_relative_bias", "max_relative_uncertainty",
        "observation_sd_floor")) sparq_check_number(d[[k]], k,
        0)
    sparq_check_number(d$reference_mean, "reference_mean")
    if (!is.numeric(d$stress_retentions) || !length(d$stress_retentions) ||
        any(!is.finite(d$stress_retentions)) || any(d$stress_retentions <=
        0 | d$stress_retentions > 1) || anyDuplicated(d$stress_retentions))
        stop("Invalid stress_retentions.")
    if (d$stress_calibration_retention == 1)
        stop("100% retention cannot calibrate a deterministic mean analysis.")
    invisible(d)
}
sparq_validate_scores <-
function (score, truth)
{
    if (!is.numeric(score) || length(score) != length(truth) ||
        !length(score) || any(!is.finite(score)) || any(score <
        0 | score > 1) || anyNA(truth) || !all(truth %in% c(0,
        1)))
        stop("Scores must be finite in [0,1], paired with binary truth.")
    invisible(TRUE)
}
sparq_weighted_median <-
function (x, w)
{
    keep <- w > 0 & is.finite(x)
    if (!any(keep))
        return(NA_real_)
    x <- x[keep]
    w <- w[keep]
    o <- order(x)
    x <- x[o]
    w <- w[o]
    cw <- cumsum(w)
    n <- sum(w)
    mean(c(x[which(cw >= floor((n + 1)/2))[1]], x[which(cw >=
        ceiling((n + 1)/2))[1]]))
}
sparq_with_seed <-
function (seed, code)
{
    sparq_check_number(seed, "seed", 0, .Machine$integer.max,
        TRUE)
    had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
    old_seed <- if (had_seed)
        get(".Random.seed", .GlobalEnv)
    else NULL
    old_kind <- RNGkind()
    on.exit({
        do.call(RNGkind, as.list(old_kind))
        if (had_seed) assign(".Random.seed", old_seed, .GlobalEnv) else if (exists(".Random.seed",
            .GlobalEnv, inherits = FALSE)) rm(".Random.seed",
            envir = .GlobalEnv)
    })
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(as.integer(seed))
    force(code)
}
sparq_write_design <-
function (d, path)
{
    sparq_validate_design(d)
    encode <- function(x) if (is.character(x))
        paste0("\"", x, "\"")
    else if (length(x) > 1L)
        paste0("[", paste(format(x, trim = TRUE, scientific = FALSE),
            collapse = ", "), "]")
    else format(x, trim = TRUE, scientific = FALSE)
    writeLines(paste0(names(d), ": ", vapply(d, encode, character(1))),
        path)
    invisible(path)
}
