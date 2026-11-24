# Presets, conventional result methods, and cohort-facing helpers.

sparq_preset <- function(preset = c("quick", "standard", "thorough")) {
  preset <- match.arg(preset)
  values <- switch(
    preset,
    quick = list(
      retention = 0.75,
      n_iterations = 25L,
      min_iterations = 15L,
      instability_bootstrap_B = 200L,
      support_bootstrap_B = 200L,
      batch_size = 5L,
      progress_every = 5L,
      verbose = TRUE
    ),
    standard = list(
      retention = 0.75,
      n_iterations = 100L,
      min_iterations = 50L,
      instability_bootstrap_B = 2000L,
      support_bootstrap_B = 1000L,
      batch_size = 10L,
      progress_every = 10L,
      verbose = TRUE
    ),
    thorough = list(
      retention = 0.75,
      n_iterations = 500L,
      min_iterations = 100L,
      instability_bootstrap_B = 10000L,
      support_bootstrap_B = 5000L,
      batch_size = 25L,
      progress_every = 25L,
      verbose = TRUE
    )
  )
  c(list(name = preset), values)
}

sparq_resolve_preset <- function(preset = "standard", retention = NULL,
                                 n_iterations = NULL, min_iterations = NULL,
                                 instability_bootstrap_B = NULL,
                                 support_bootstrap_B = NULL, batch_size = NULL,
                                 progress_every = NULL, verbose = NULL) {
  values <- sparq_preset(preset)
  resolved <- list(
    preset = values$name,
    retention = if (is.null(retention)) values$retention else retention,
    n_iterations = if (is.null(n_iterations)) values$n_iterations else n_iterations,
    min_iterations = if (is.null(min_iterations)) values$min_iterations else min_iterations,
    instability_bootstrap_B = if (is.null(instability_bootstrap_B)) values$instability_bootstrap_B else instability_bootstrap_B,
    support_bootstrap_B = if (is.null(support_bootstrap_B)) values$support_bootstrap_B else support_bootstrap_B,
    batch_size = if (is.null(batch_size)) values$batch_size else batch_size,
    progress_every = if (is.null(progress_every)) values$progress_every else progress_every,
    verbose = if (is.null(verbose)) values$verbose else verbose
  )
  sparq_check_number(resolved$retention, "retention", .Machine$double.eps, 1)
  sparq_check_number(resolved$n_iterations, "n_iterations", 1, integer = TRUE)
  sparq_check_number(resolved$min_iterations, "min_iterations", 1, integer = TRUE)
  sparq_check_number(resolved$instability_bootstrap_B, "instability_bootstrap_B", 1, integer = TRUE)
  sparq_check_number(resolved$support_bootstrap_B, "support_bootstrap_B", 1, integer = TRUE)
  sparq_check_number(resolved$batch_size, "batch_size", 1, integer = TRUE)
  sparq_check_number(resolved$progress_every, "progress_every", 1, integer = TRUE)
  if (!is.logical(resolved$verbose) || length(resolved$verbose) != 1L || is.na(resolved$verbose)) {
    stop("verbose must be TRUE or FALSE.", call. = FALSE)
  }
  resolved
}

sparq_inform_preset <- function(values) {
  sparq_inform(
    paste0(
      "Preset `", values$preset, "`: ", values$n_iterations,
      " perturbations; ", round(values$retention * 100), "% retention; minimum ",
      values$min_iterations, " completed; ", format(values$instability_bootstrap_B, big.mark = ","),
      " instability bootstrap draws; progress every ", values$progress_every, "."
    ),
    values$verbose
  )
}

sparq_validate_workflow <- function(data, analysis_function,
                                    output_type = c("scalar", "ranked", "partition", "graph"),
                                    stress_model = c("uniform_random", "contiguous_hole", "none", "custom"),
                                    x_col = NULL, y_col = NULL,
                                    spot_id_column = NULL, custom_function = NULL,
                                    partition_graph = NULL,
                                    preflight_retention = 0.95, seed = 1,
                                    verbose = TRUE) {
  if (!is.data.frame(data) || !nrow(data)) {
    stop("data must be a nonempty data.frame.", call. = FALSE)
  }
  if (!is.function(analysis_function)) {
    stop("analysis_function must be a function.", call. = FALSE)
  }
  output_type <- match.arg(output_type)
  stress_model <- match.arg(stress_model)
  sparq_check_number(preflight_retention, "preflight_retention", .Machine$double.eps, 1 - .Machine$double.eps)
  sparq_check_number(seed, "seed", 0, integer = TRUE)
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE.", call. = FALSE)
  }

  coordinate_status <- "not required"
  if (xor(is.null(x_col), is.null(y_col))) {
    stop("Supply both x_col and y_col, or neither.", call. = FALSE)
  }
  if (!is.null(x_col)) {
    if (!is.character(x_col) || !is.character(y_col) || length(x_col) != 1L || length(y_col) != 1L ||
        is.na(x_col) || is.na(y_col) || !x_col %in% names(data) || !y_col %in% names(data)) {
      stop("x_col and y_col must name columns in data.", call. = FALSE)
    }
    if (!is.numeric(data[[x_col]]) || !is.numeric(data[[y_col]]) ||
        any(!is.finite(data[[x_col]])) || any(!is.finite(data[[y_col]]))) {
      stop("x_col and y_col must contain finite numeric coordinates.", call. = FALSE)
    }
    coordinate_status <- "valid"
  } else if (stress_model == "contiguous_hole") {
    stop("contiguous_hole preflight requires finite numeric x_col and y_col.", call. = FALSE)
  }

  input_ids <- if (!is.null(spot_id_column)) {
    if (!is.character(spot_id_column) || length(spot_id_column) != 1L || !spot_id_column %in% names(data)) {
      stop("spot_id_column must name a column in data.", call. = FALSE)
    }
    as.character(data[[spot_id_column]])
  } else if (!is.null(rownames(data)) && !anyDuplicated(rownames(data)) && all(nzchar(rownames(data)))) {
    rownames(data)
  } else {
    as.character(seq_len(nrow(data)))
  }
  if (anyNA(input_ids) || any(!nzchar(input_ids)) || anyDuplicated(input_ids)) {
    stop("Input spot IDs must be unique and nonmissing.", call. = FALSE)
  }
  if (output_type == "partition") {
    if (is.null(partition_graph)) {
      stop("partition_graph is required for a partition assessment; supply the native spatial-neighbor graph used by the analysis.", call. = FALSE)
    }
    graph <- .sparq_spatial_partition_edges(partition_graph)
    vertices <- unique(c(graph$source, graph$target))
    missing_vertices <- setdiff(input_ids, vertices)
    if (length(missing_vertices)) {
      stop("partition_graph lacks input spot IDs: ", paste(head(missing_vertices, 5L), collapse = ", "), call. = FALSE)
    }
  }

  validate_output <- function(result, expected_ids, stage) {
    fail <- function(message) {
      stop(paste0(stage, " analysis returned an invalid ", output_type, " output: ", message), call. = FALSE)
    }
    if (output_type == "scalar") {
      value <- suppressWarnings(as.numeric(result))
      if (length(value) != 1L || !is.finite(value)) fail("expected one finite numeric value")
      return(list(n_returned = 1L, id_status = "not applicable"))
    }
    values <- as.character(result)
    if (!length(values) || anyNA(values) || any(!nzchar(values))) {
      fail("expected a nonempty vector of nonmissing values")
    }
    if (output_type %in% c("ranked", "graph")) {
      if (anyDuplicated(values)) fail("expected unique feature or edge IDs")
      return(list(n_returned = length(values), id_status = "not applicable"))
    }
    returned_ids <- names(result)
    if (is.null(returned_ids) || anyNA(returned_ids) || any(!nzchar(returned_ids)) || anyDuplicated(returned_ids)) {
      fail("partition labels must be named by unique, nonmissing spot IDs")
    }
    if (!setequal(returned_ids, expected_ids)) {
      missing_ids <- length(setdiff(expected_ids, returned_ids))
      extra_ids <- length(setdiff(returned_ids, expected_ids))
      fail(paste0("spot IDs do not match retained data (", missing_ids, " missing; ", extra_ids, " unexpected)"))
    }
    list(n_returned = length(values), id_status = "exact retained-spot match")
  }

  sparq_inform(paste0("Preflight: validating ", nrow(data), " input observations and a small spot-removal rerun."), verbose)
  full_result <- tryCatch(
    sparq_with_seed(sparq_seed(seed, 0L, 21L), analysis_function(data)),
    error = function(error) stop(paste0("Reference analysis failed during preflight: ", conditionMessage(error)), call. = FALSE)
  )
  full_check <- validate_output(full_result, input_ids, "Reference")

  if (stress_model == "none") {
    perturbed_data <- data
    retained_ids <- input_ids
    perturbation_status <- "not run: stress_model is none"
  } else {
    perturbed_data <- tryCatch(
      sparq_stress_model(
        data = data, retention = preflight_retention, model = "uniform_random",
        x_col = x_col, y_col = y_col, custom_function = custom_function,
        iteration = 1L, seed = sparq_seed(seed, 1L, 21L)
      ),
      error = function(error) stop(paste0("Could not create the preflight perturbation: ", conditionMessage(error)), call. = FALSE)
    )
    if (nrow(perturbed_data) >= nrow(data)) {
      stop("Preflight did not remove any observations. Use more input observations or lower preflight_retention.", call. = FALSE)
    }
    retained_ids <- if (!is.null(spot_id_column)) {
      as.character(perturbed_data[[spot_id_column]])
    } else if (!is.null(rownames(perturbed_data)) && all(nzchar(rownames(perturbed_data)))) {
      rownames(perturbed_data)
    } else {
      as.character(seq_len(nrow(perturbed_data)))
    }
    perturbation_status <- "small uniform-random spot removal"
  }
  perturbed_result <- tryCatch(
    sparq_with_seed(sparq_seed(seed, 1L, 22L), analysis_function(perturbed_data)),
    error = function(error) stop(paste0("Perturbed analysis failed during preflight: ", conditionMessage(error)), call. = FALSE)
  )
  perturbed_check <- validate_output(perturbed_result, retained_ids, "Perturbed")

  checks <- data.frame(
    check = c("input observations", "coordinate columns", "reference output", "perturbed output", "spot-removal rerun"),
    status = c("passed", coordinate_status, "passed", "passed", "passed"),
    detail = c(
      as.character(nrow(data)), coordinate_status,
      paste0(full_check$n_returned, " values; ", full_check$id_status),
      paste0(perturbed_check$n_returned, " values; ", perturbed_check$id_status),
      paste0(nrow(perturbed_data), "/", nrow(data), " retained; ", perturbation_status)
    ),
    stringsAsFactors = FALSE
  )
  out <- list(
    passed = TRUE, input_size = nrow(data), retained_size = nrow(perturbed_data),
    output_type = output_type, stress_model = stress_model,
    preflight_retention = preflight_retention, seed = seed,
    coordinate_status = coordinate_status, checks = checks,
    reference_result = full_result, perturbed_result = perturbed_result
  )
  class(out) <- "sparq_workflow_preflight"
  sparq_inform(paste0("Preflight passed: ", nrow(data), " input observations; ", nrow(perturbed_data), " retained; ", output_type, " output contract confirmed."), verbose)
  out
}

print.sparq_workflow_preflight <- function(x, ...) {
  cat("<SPARQ workflow preflight>\n")
  cat("  Input observations: ", x$input_size, "; retained for dry run: ", x$retained_size, "\n", sep = "")
  cat("  Output type: ", x$output_type, "; coordinate columns: ", x$coordinate_status, "\n", sep = "")
  cat("  Result: passed\n")
  invisible(x)
}

summary.sparq_assessment <- function(object, ...) {
  out <- object$summary
  if (!is.data.frame(out)) {
    stop("This assessment has no summary table.", call. = FALSE)
  }
  attr(out, "result_id") <- object$result_id
  attr(out, "settings") <- object$settings
  class(out) <- unique(c("summary.sparq_assessment", class(out)))
  out
}

summary.sparq_workflow_assessment <- function(object, ...) {
  NextMethod("summary")
}

summary.sparq_script_assessment <- function(object, ...) {
  NextMethod("summary")
}

summary.sparq_spatial_assessment <- function(object, ...) {
  NextMethod("summary")
}

summary.sparq_adaptive_assessment <- function(object, ...) {
  out <- NextMethod("summary")
  if (is.data.frame(object$adaptive_decision) && nrow(object$adaptive_decision) == 1L) {
    for (name in names(object$adaptive_decision)) {
      out[[paste0("adaptive_", name)]] <- object$adaptive_decision[[name]][1]
    }
  }
  out
}

as.data.frame.sparq_assessment <- function(x, row.names = NULL, optional = FALSE, ...) {
  as.data.frame(summary(x), row.names = row.names, optional = optional, ...)
}

as.data.frame.sparq_workflow_assessment <- function(x, row.names = NULL, optional = FALSE, ...) {
  NextMethod("as.data.frame")
}

as.data.frame.sparq_script_assessment <- function(x, row.names = NULL, optional = FALSE, ...) {
  NextMethod("as.data.frame")
}

as.data.frame.sparq_spatial_assessment <- function(x, row.names = NULL, optional = FALSE, ...) {
  NextMethod("as.data.frame")
}

as.data.frame.sparq_adaptive_assessment <- function(x, row.names = NULL, optional = FALSE, ...) {
  NextMethod("as.data.frame")
}

sparq_plot_assessment <- function(x, ...) {
  comparisons <- x$comparisons
  summary_table <- x$summary
  if (!is.data.frame(comparisons) || !is.data.frame(summary_table) || !nrow(summary_table)) {
    stop("This assessment does not contain plottable comparison and summary tables.", call. = FALSE)
  }
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::layout(matrix(c(1, 2), nrow = 1L))
  valid <- comparisons[comparisons$status == "ok" & is.finite(comparisons$instability), , drop = FALSE]
  if (nrow(valid)) {
    graphics::plot(valid$iteration, valid$instability, pch = 16, col = "#2F6B9A",
      xlab = "Perturbation iteration", ylab = "Instability",
      main = "Instability across perturbations", ...)
    graphics::abline(h = stats::median(valid$instability), lty = 2, col = "#B13A3A")
  } else {
    graphics::plot.new()
    graphics::title(main = "Instability across perturbations")
    graphics::text(0.5, 0.5, "No successful perturbation comparisons")
  }
  components <- c("iteration_support", "bias_support", "uncertainty_support", "reproducibility_support")
  components <- components[components %in% names(summary_table)]
  if (length(components)) {
    graphics::barplot(as.numeric(summary_table[1, components, drop = TRUE]), ylim = c(0, 1),
      col = c("#5B8DB8", "#F2A65A", "#73A857", "#9A6FB0")[seq_along(components)],
      names.arg = gsub("_support", "", components), ylab = "Component value",
      main = "Support-quality components")
  } else {
    graphics::plot.new()
    graphics::title(main = "Support-quality components")
  }
  invisible(x)
}

plot.sparq_assessment <- function(x, ...) sparq_plot_assessment(x, ...)
plot.sparq_workflow_assessment <- function(x, ...) sparq_plot_assessment(x, ...)
plot.sparq_script_assessment <- function(x, ...) sparq_plot_assessment(x, ...)
plot.sparq_spatial_assessment <- function(x, ...) sparq_plot_assessment(x, ...)
plot.sparq_adaptive_assessment <- function(x, ...) sparq_plot_assessment(x, ...)

sparq_format_elapsed <- function(seconds) {
  if (!is.finite(seconds) || seconds < 0) return("unknown")
  if (seconds < 60) return(paste0(round(seconds), "s"))
  if (seconds < 3600) return(paste0(round(seconds / 60, 1), "m"))
  paste0(round(seconds / 3600, 1), "h")
}

sparq_assess_cohort_workflow <- function(sample_data, analysis_function,
                                         output_type = c("scalar", "ranked", "partition", "graph"),
                                         preset = "standard", output_dir = NULL,
                                         progress_every = 1L, seed = 1, verbose = TRUE, ...) {
  if (!is.list(sample_data) || !length(sample_data)) {
    stop("sample_data must be a nonempty named list of sample-level data frames.", call. = FALSE)
  }
  if (is.null(names(sample_data)) || anyNA(names(sample_data)) || any(!nzchar(names(sample_data))) || anyDuplicated(names(sample_data))) {
    stop("sample_data must have unique, nonmissing sample IDs.", call. = FALSE)
  }
  if (!is.function(analysis_function)) stop("analysis_function must be a function.", call. = FALSE)
  output_type <- match.arg(output_type)
  sparq_check_number(progress_every, "progress_every", 1, integer = TRUE)
  sparq_check_number(seed, "seed", 0, integer = TRUE)
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) stop("verbose must be TRUE or FALSE.", call. = FALSE)
  if (!is.null(output_dir) && (!is.character(output_dir) || length(output_dir) != 1L || is.na(output_dir) || !nzchar(output_dir))) {
    stop("output_dir must be NULL or one nonempty path.", call. = FALSE)
  }
  if (!is.null(output_dir) && dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
    stop("output_dir already exists and is not empty.", call. = FALSE)
  }
  if (!is.null(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  n_samples <- length(sample_data)
  started <- Sys.time()
  manifest <- vector("list", n_samples)
  assessments <- vector("list", n_samples)
  names(assessments) <- names(sample_data)
  sparq_inform(paste0("Cohort assessment started: 0/", n_samples, " sections complete."), verbose)

  for (i in seq_along(sample_data)) {
    sample_id <- names(sample_data)[i]
    sample_started <- Sys.time()
    fit <- tryCatch(
      sparq_assess_workflow(
        data = sample_data[[i]], analysis_function = analysis_function,
        output_type = output_type, preset = preset, result_id = sample_id,
        seed = sparq_seed(seed, i), verbose = FALSE, ...
      ),
      error = function(error) error
    )
    elapsed <- as.numeric(difftime(Sys.time(), sample_started, units = "secs"))
    if (inherits(fit, "error")) {
      manifest[[i]] <- data.frame(
        sample_id = sample_id, status = "failed", error_message = conditionMessage(fit),
        elapsed_seconds = elapsed, seed = NA_integer_, stringsAsFactors = FALSE
      )
    } else {
      assessments[[i]] <- fit
      manifest[[i]] <- data.frame(
        sample_id = sample_id, status = "completed", error_message = "",
        elapsed_seconds = elapsed, seed = fit$settings$seed,
        support_quality = fit$summary$support_quality[1],
        n_completed = fit$summary$n_iterations[1], n_failed = fit$summary$n_failed[1],
        stringsAsFactors = FALSE
      )
    }
    if (verbose && (i %% progress_every == 0L || i == n_samples)) {
      elapsed_total <- as.numeric(difftime(Sys.time(), started, units = "secs"))
      failed <- sum(vapply(manifest[seq_len(i)], function(x) identical(x$status, "failed"), logical(1)))
      remaining <- if (i) elapsed_total / i * (n_samples - i) else NA_real_
      sparq_inform(
        paste0(i, "/", n_samples, " sections complete; ", failed,
          " failed; estimated remaining ", sparq_format_elapsed(remaining), "."),
        TRUE
      )
    }
  }
  manifest <- sparq_bind_rows(manifest)
  summaries <- sparq_bind_rows(lapply(assessments, function(x) if (is.null(x)) NULL else x$summary))
  failures <- manifest[manifest$status == "failed", c("sample_id", "error_message"), drop = FALSE]
  package_version <- tryCatch(as.character(utils::packageVersion("SPARQ")),
    error = function(error) "source checkout")
  manifest$preset <- preset
  manifest$output_type <- output_type
  manifest$cohort_seed <- seed
  manifest$R_version <- R.version.string
  manifest$SPARQ_version <- package_version
  settings <- list(
    preset = preset, output_type = output_type, cohort_seed = seed, started_at = as.character(started),
    completed_at = as.character(Sys.time()), R_version = R.version.string,
    package_version = package_version
  )
  reliability <- NULL
  completed_assessments <- Filter(Negate(is.null), assessments)
  if (length(completed_assessments) && all(vapply(
      completed_assessments,
      function(assessment) is.data.frame(assessment$reliability_decision) &&
        nrow(assessment$reliability_decision) == 1L,
      logical(1)
    ))) {
    reliability <- sparq_infer_cohort_reliability(completed_assessments)
  }
  files <- list()
  if (!is.null(output_dir)) {
    files$manifest <- file.path(output_dir, "cohort_manifest.tsv")
    files$summaries <- file.path(output_dir, "cohort_support_summary.tsv")
    files$failures <- file.path(output_dir, "cohort_failures.tsv")
    files$settings <- file.path(output_dir, "cohort_settings.rds")
    if (!is.null(reliability)) {
      files$reliability <- file.path(output_dir, "cohort_reliability.tsv")
    }
    utils::write.table(manifest, files$manifest, sep = "\t", row.names = FALSE, quote = FALSE)
    utils::write.table(summaries, files$summaries, sep = "\t", row.names = FALSE, quote = FALSE)
    utils::write.table(failures, files$failures, sep = "\t", row.names = FALSE, quote = FALSE)
    if (!is.null(reliability)) {
      utils::write.table(
        reliability$section_reliability,
        files$reliability,
        sep = "\t", row.names = FALSE, quote = FALSE
      )
    }
    saveRDS(settings, files$settings)
  }
  out <- list(
    summaries = summaries, failures = failures, manifest = manifest,
    settings = settings, assessments = assessments, output_files = files,
    reliability = reliability
  )
  class(out) <- "sparq_cohort_assessment"
  out
}

print.sparq_cohort_assessment <- function(x, ...) {
  cat("<SPARQ cohort assessment>\n")
  cat("  Sections completed: ", sum(x$manifest$status == "completed"), "/", nrow(x$manifest), "\n", sep = "")
  cat("  Sections failed: ", sum(x$manifest$status == "failed"), "\n", sep = "")
  if (length(x$output_files)) cat("  Manifest: ", x$output_files$manifest, "\n", sep = "")
  invisible(x)
}

summary.sparq_cohort_assessment <- function(object, ...) object$summaries
as.data.frame.sparq_cohort_assessment <- function(x, row.names = NULL, optional = FALSE, ...) x$summaries
