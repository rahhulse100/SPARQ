#' Assess several analyses against one shared set of spatial perturbations
#'
#' Runs each supplied analysis on the same reference data and on the same
#' retained spot set at every iteration. This makes cross-analysis comparisons
#' paired, avoids repeated loss generation, and records one shared loss audit.
#'
#' @param data Nonempty spot-level data frame.
#' @param analyses Named list. Each element must contain `analysis_function`
#'   and `output_type`; optional `top_k` and `reference_scale` are supported.
#' @param stress_model,retention,n_iterations,x_col,y_col,loss_model,seed
#'   Standard SPARQ perturbation controls.
#' @param failure_action Whether a failed method iteration is recorded or stops.
#' @param keep_perturbed_results Retain method outputs for every iteration.
#' @param ... Arguments passed to [sparq_support_from_comparisons()].
#'
#' @return A `sparq_multi_assessment` containing one `sparq_assessment` per
#'   named analysis plus shared loss audits and an error ledger.
#' @export
sparq_assess_many <- function(
    data, analyses, stress_model = "uniform_random", retention = NULL,
    n_iterations = NULL, x_col = NULL, y_col = NULL, loss_model = NULL,
    seed = 1, failure_action = c("record", "stop"),
    keep_perturbed_results = FALSE, preset = "standard", verbose = NULL, ...) {
  failure_action <- match.arg(failure_action)
  if (!is.data.frame(data) || !nrow(data)) stop("data must be a nonempty data.frame.", call. = FALSE)
  if (!is.list(analyses) || !length(analyses) || is.null(names(analyses)) ||
      anyNA(names(analyses)) || any(!nzchar(names(analyses))) || anyDuplicated(names(analyses))) {
    stop("analyses must be a nonempty named list with unique names.", call. = FALSE)
  }
  if (!is.null(loss_model) && !inherits(loss_model, "sparq_loss_model")) {
    stop("loss_model must be created by sparq_define_loss_model().", call. = FALSE)
  }
  resolved <- sparq_resolve_preset(preset = preset, retention = retention,
    n_iterations = n_iterations, verbose = verbose)
  retention <- resolved$retention; n_iterations <- resolved$n_iterations
  verbose <- resolved$verbose
  comparator_for <- function(spec) {
    type <- match.arg(spec$output_type, c("scalar", "ranked", "partition", "graph"))
    if (type == "scalar") return(sparq_compare_scalar)
    if (type == "partition") return(sparq_compare_partition)
    if (type == "graph") return(sparq_compare_graph)
    k <- if (is.null(spec$top_k)) 100L else as.integer(spec$top_k)
    function(full, perturbed) sparq_compare_ranked(full, perturbed, k = k)
  }
  specs <- lapply(analyses, function(spec) {
    if (!is.list(spec) || !is.function(spec$analysis_function) || is.null(spec$output_type)) {
      stop("Each analysis must provide analysis_function and output_type.", call. = FALSE)
    }
    spec$comparator <- comparator_for(spec)
    spec$reference_scale <- if (is.null(spec$reference_scale)) 1 else spec$reference_scale
    sparq_check_number(spec$reference_scale, "reference_scale")
    spec
  })
  references <- lapply(seq_along(specs), function(i) {
    sparq_with_seed(sparq_seed(seed, i, 701L), specs[[i]]$analysis_function(data))
  })
  names(references) <- names(specs)
  comparisons <- setNames(vector("list", length(specs)), names(specs))
  outputs <- setNames(vector("list", length(specs)), names(specs))
  for (name in names(specs)) { comparisons[[name]] <- vector("list", n_iterations); outputs[[name]] <- vector("list", n_iterations) }
  audits <- vector("list", n_iterations)
  for (iteration in seq_len(n_iterations)) {
    perturbed <- sparq_stress_model(data, retention = retention, model = stress_model,
      x_col = x_col, y_col = y_col, iteration = iteration, seed = seed, loss_model = loss_model)
    audits[[iteration]] <- attr(perturbed, "sparq_loss_audit")
    for (name in names(specs)) {
      spec <- specs[[name]]
      outcome <- tryCatch({
        result <- spec$analysis_function(perturbed)
        list(comparison = sparq_comparison_row(spec$comparator, references[[name]], result, iteration), result = result)
      }, error = function(error) {
        if (failure_action == "stop") stop(paste0(name, ", iteration ", iteration, ": ", conditionMessage(error)), call. = FALSE)
        list(comparison = data.frame(instability = NA_real_, iteration = iteration, status = "failed",
          error_message = conditionMessage(error), stringsAsFactors = FALSE), result = NULL)
      })
      outcome$comparison$n_retained <- nrow(perturbed)
      outcome$comparison$loss_model <- if (is.null(loss_model)) stress_model else loss_model$type
      comparisons[[name]][[iteration]] <- outcome$comparison
      if (keep_perturbed_results) outputs[[name]][[iteration]] <- outcome$result
    }
    sparq_progress(iteration, n_iterations, verbose, prefix = "Shared perturbations")
  }
  assessments <- lapply(names(specs), function(name) {
    spec <- specs[[name]]; tab <- sparq_bind_rows(comparisons[[name]])
    summary <- sparq_support_from_comparisons(tab, n_iterations = n_iterations,
      reference_scale = spec$reference_scale, seed = sparq_seed(seed, match(name, names(specs)), 702L), ...)
    summary$result_id <- name
    structure(list(result_id = name, output_type = spec$output_type,
      full_result = references[[name]],
      perturbed_results = if (keep_perturbed_results) outputs[[name]] else NULL,
      comparisons = tab, summary = summary, loss_audits = audits,
      failures = tab[tab$status != "ok", , drop = FALSE],
      settings = list(stress_model = stress_model, retention = retention, n_iterations = n_iterations,
        seed = seed, shared_perturbations = TRUE, loss_model = loss_model),
      input_provenance = sparq_input_provenance(data)), class = "sparq_assessment")
  })
  names(assessments) <- names(specs)
  errors <- sparq_bind_rows(lapply(names(assessments), function(name) {
    x <- assessments[[name]]$failures
    if (!nrow(x)) return(NULL)
    x$analysis <- name; x
  }))
  out <- list(assessments = assessments, loss_audits = audits, errors = errors,
    settings = list(retention = retention, n_iterations = n_iterations, stress_model = stress_model,
      loss_model = loss_model, seed = seed, shared_perturbations = TRUE),
    input_provenance = sparq_input_provenance(data))
  class(out) <- "sparq_multi_assessment"
  out
}

#' @export
summary.sparq_multi_assessment <- function(object, ...) {
  out <- sparq_bind_rows(lapply(names(object$assessments), function(name) {
    row <- object$assessments[[name]]$summary; row$analysis <- name; row
  }))
  out[, c("analysis", setdiff(names(out), "analysis")), drop = FALSE]
}

#' @export
print.sparq_multi_assessment <- function(x, ...) {
  cat("<SPARQ multi-analysis assessment>\n")
  cat("  Analyses: ", length(x$assessments), "; shared perturbations: ",
    x$settings$n_iterations, "\n", sep = "")
  cat("  Failed analysis iterations: ", nrow(x$errors), "\n", sep = "")
  invisible(x)
}

#' Write a complete audit bundle for a multi-analysis assessment
#'
#' @param assessment A result from [sparq_assess_many()].
#' @param output_dir New directory for the audit bundle.
#' @param include_plots,verbose Passed to [sparq_report()].
#' @return A list of per-analysis report objects and output paths.
#' @export
sparq_report_many <- function(assessment, output_dir, include_plots = TRUE,
                              verbose = TRUE) {
  if (!inherits(assessment, "sparq_multi_assessment")) {
    stop("assessment must be returned by sparq_assess_many().", call. = FALSE)
  }
  if (!is.character(output_dir) || length(output_dir) != 1L || is.na(output_dir) ||
      !nzchar(output_dir) || dir.exists(output_dir)) {
    stop("output_dir must be one new, nonempty directory path.", call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  summary_table <- summary(assessment)
  utils::write.table(summary_table, file.path(output_dir, "multi_analysis_summary.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.table(assessment$errors, file.path(output_dir, "multi_analysis_errors.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE)
  saveRDS(assessment, file.path(output_dir, "multi_assessment.rds"))
  reports <- lapply(names(assessment$assessments), function(name) {
    sparq_report(assessment$assessments[[name]],
      output_dir = file.path(output_dir, name), include_plots = include_plots,
      verbose = verbose)
  })
  names(reports) <- names(assessment$assessments)
  manifest <- data.frame(
    analysis = names(reports),
    report_directory = file.path(output_dir, names(reports)),
    stringsAsFactors = FALSE
  )
  utils::write.table(manifest, file.path(output_dir, "analysis_report_manifest.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE)
  list(summary = file.path(output_dir, "multi_analysis_summary.tsv"),
    errors = file.path(output_dir, "multi_analysis_errors.tsv"),
    manifest = file.path(output_dir, "analysis_report_manifest.tsv"), reports = reports)
}
