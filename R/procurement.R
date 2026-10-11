#' Assess procurement-aware spatial support sensitivity
#'
#' Evaluates one downstream spatial-transcriptomic result against a declared
#' collection of plausible support-loss scenarios. This workflow is designed to
#' make sensitivity to section coverage, edge loss, local tissue loss, and
#' spot-level quality variation transparent. It does not correct batch effects,
#' impute unobserved tissue, or recover biological information absent from the
#' specimen.
#'
#' The reference analysis is run once and reused for every scenario. Each
#' scenario is otherwise an independent SPARQ assessment with its own
#' perturbation comparisons and loss audit.
#'
#' @param data A nonempty spot-level data frame.
#' @param analysis_function A function rerunning the user's analysis on a
#'   perturbed data frame.
#' @param output_type One of `"scalar"`, `"ranked"`, `"partition"`, or
#'   `"graph"`.
#' @param x_col,y_col Optional names of finite numeric coordinate columns.
#'   They are required by `"contiguous_hole"` and `"edge_weighted"` scenarios.
#' @param quality_col Optional name of a finite numeric spot-level quality
#'   metric, such as library size or detected-gene count. It is required by the
#'   `"quality_weighted"` scenario.
#' @param scenarios Character vector of scenario names. Allowed names are
#'   `"uniform_random"`, `"contiguous_hole"`, `"edge_weighted"`, and
#'   `"quality_weighted"`. If `NULL`, SPARQ uses every scenario supported by
#'   the supplied columns.
#' @param retention,n_iterations,preset Perturbation settings passed to
#'   [sparq_assess_workflow()].
#' @param quality_direction Whether higher or lower values of `quality_col` are
#'   preferentially retained under quality-weighted loss.
#' @param quality_strength,edge_strength Nonnegative loss-model strengths.
#' @param top_k Ranking depth for ranked results.
#' @param result_id Identifier for the assessed section/result.
#' @param reference_scale Comparator reference scale.
#' @param seed Master seed. Each scenario receives a deterministic derived seed.
#' @param preflight Whether to run [sparq_validate_workflow()] before the
#'   scenario assessments.
#' @param ... Additional arguments passed to [sparq_assess_workflow()], such as
#'   `reliability_tolerance`, `required_reliability`, or `cache_dir`.
#'
#' @return A `sparq_procurement_assessment` object containing the shared
#'   reference result, scenario-specific assessments, a compact scenario table,
#'   loss audits, quality profile, settings, and an explicitly non-corrective
#'   interpretation.
#' @export
#'
#' @examples
#' spots <- data.frame(
#'   value = rnorm(60), x = runif(60), y = runif(60),
#'   library_size = rexp(60, rate = 1 / 1000)
#' )
#' assessment <- sparq_assess_procurement(
#'   data = spots,
#'   analysis_function = function(x) mean(x$value),
#'   output_type = "scalar",
#'   x_col = "x", y_col = "y", quality_col = "library_size",
#'   preset = "quick", n_iterations = 10, min_iterations = 5,
#'   instability_bootstrap_B = 20, verbose = FALSE
#' )
#' summary(assessment)
sparq_assess_procurement <- function(
    data, analysis_function,
    output_type = c("scalar", "ranked", "partition", "graph"),
    x_col = NULL, y_col = NULL, quality_col = NULL, scenarios = NULL,
    retention = NULL, n_iterations = NULL, preset = "standard",
    quality_direction = c("higher_quality_retained", "lower_quality_retained"),
    quality_strength = 1, edge_strength = 1, top_k = 100,
    result_id = "result", reference_scale = 1, seed = 1,
    preflight = TRUE, verbose = NULL, ...) {

  output_type <- match.arg(output_type)
  quality_direction <- match.arg(quality_direction)
  sparq_check_number(reference_scale, "reference_scale")
  sparq_check_number(seed, "seed", 0, integer = TRUE)
  sparq_check_number(quality_strength, "quality_strength", 0)
  sparq_check_number(edge_strength, "edge_strength", 0)

  if (!is.data.frame(data) || !nrow(data)) {
    stop("data must be a nonempty data.frame.", call. = FALSE)
  }
  if (!is.function(analysis_function)) {
    stop("analysis_function must be a function.", call. = FALSE)
  }
  if (xor(is.null(x_col), is.null(y_col))) {
    stop("Supply both x_col and y_col, or neither.", call. = FALSE)
  }

  has_coordinates <- !is.null(x_col)
  if (has_coordinates) {
    if (!is.character(x_col) || !is.character(y_col) || length(x_col) != 1L ||
        length(y_col) != 1L || !all(c(x_col, y_col) %in% names(data)) ||
        !is.numeric(data[[x_col]]) || !is.numeric(data[[y_col]]) ||
        any(!is.finite(data[[x_col]])) || any(!is.finite(data[[y_col]]))) {
      stop("x_col and y_col must name finite numeric columns in data.", call. = FALSE)
    }
  }

  has_quality <- !is.null(quality_col)
  if (has_quality && (!is.character(quality_col) || length(quality_col) != 1L ||
      !quality_col %in% names(data) || !is.numeric(data[[quality_col]]) ||
      any(!is.finite(data[[quality_col]])))) {
    stop("quality_col must name one finite numeric column in data.", call. = FALSE)
  }

  supported_scenarios <- "uniform_random"
  if (has_coordinates) {
    supported_scenarios <- c(supported_scenarios, "contiguous_hole", "edge_weighted")
  }
  if (has_quality) {
    supported_scenarios <- c(supported_scenarios, "quality_weighted")
  }
  if (is.null(scenarios)) {
    scenarios <- supported_scenarios
  }
  allowed_scenarios <- c(
    "uniform_random", "contiguous_hole", "edge_weighted", "quality_weighted"
  )
  if (!is.character(scenarios) || !length(scenarios) || anyNA(scenarios) ||
      anyDuplicated(scenarios) || any(!scenarios %in% allowed_scenarios)) {
    stop("scenarios must be unique allowed scenario names.", call. = FALSE)
  }
  unavailable <- setdiff(scenarios, supported_scenarios)
  if (length(unavailable)) {
    stop(
      "Requested scenario(s) require unavailable input columns: ",
      paste(unavailable, collapse = ", "), ".",
      call. = FALSE
    )
  }

  if (!is.logical(preflight) || length(preflight) != 1L || is.na(preflight)) {
    stop("preflight must be TRUE or FALSE.", call. = FALSE)
  }
  resolved <- sparq_resolve_preset(
    preset = preset, retention = retention, n_iterations = n_iterations,
    verbose = verbose
  )
  retention <- resolved$retention
  n_iterations <- resolved$n_iterations
  verbose <- resolved$verbose

  sparq_inform(
    paste0("Procurement assessment: running one reference analysis and ",
      length(scenarios), " declared support-loss scenario(s)."),
    verbose
  )

  preflight_result <- NULL
  if (preflight) {
    preflight_result <- sparq_validate_workflow(
      data = data, analysis_function = analysis_function,
      output_type = output_type, stress_model = "uniform_random",
      x_col = x_col, y_col = y_col, seed = seed, verbose = verbose
    )
  }

  reference_result <- if (!is.null(preflight_result)) {
    preflight_result$reference_result
  } else {
    sparq_with_seed(
      sparq_seed(seed, 0L, 401L),
      analysis_function(data)
    )
  }

  make_loss_model <- function(scenario) {
    switch(
      scenario,
      uniform_random = NULL,
      contiguous_hole = sparq_define_loss_model(type = "contiguous_hole"),
      edge_weighted = sparq_define_loss_model(
        type = "edge_weighted", edge_strength = edge_strength
      ),
      quality_weighted = sparq_define_loss_model(
        type = "quality_weighted", quality_col = quality_col,
        quality_direction = quality_direction,
        quality_strength = quality_strength
      )
    )
  }

  assessments <- vector("list", length(scenarios))
  names(assessments) <- scenarios
  scenario_rows <- vector("list", length(scenarios))

  for (index in seq_along(scenarios)) {
    scenario <- scenarios[[index]]
    model <- make_loss_model(scenario)
    sparq_inform(
      paste0("Procurement scenario ", index, "/", length(scenarios), ": ", scenario, "."),
      verbose
    )
    assessment <- sparq_assess_workflow(
      data = data,
      analysis_function = analysis_function,
      output_type = output_type,
      stress_model = "uniform_random",
      retention = retention,
      n_iterations = n_iterations,
      top_k = top_k,
      x_col = x_col,
      y_col = y_col,
      result_id = paste0(result_id, "__", scenario),
      reference_scale = reference_scale,
      seed = sparq_seed(seed, index, 402L),
      preset = preset,
      verbose = verbose,
      loss_model = model,
      full_result = reference_result,
      keep_perturbed_results = TRUE,
      ...
    )
    assessment$procurement_scenario <- scenario
    assessment$procurement_loss_model <- model
    assessments[[scenario]] <- assessment

    summary_row <- assessment$summary[1, , drop = FALSE]
    summary_row$scenario <- scenario
    summary_row$loss_model <- if (is.null(model)) "uniform_random" else model$type
    summary_row$scenario_seed <- assessment$settings$seed
    if (is.data.frame(assessment$reliability_decision) &&
        nrow(assessment$reliability_decision) == 1L) {
      summary_row$reliability_decision <- assessment$reliability_decision$decision[1]
      summary_row$reliability_lower_confidence_bound <-
        assessment$reliability_decision$lower_confidence_bound[1]
      summary_row$reliability_upper_confidence_bound <-
        assessment$reliability_decision$upper_confidence_bound[1]
    } else {
      summary_row$reliability_decision <- NA_character_
      summary_row$reliability_lower_confidence_bound <- NA_real_
      summary_row$reliability_upper_confidence_bound <- NA_real_
    }
    scenario_rows[[index]] <- summary_row
  }

  scenario_summary <- sparq_bind_rows(scenario_rows)
  scenario_summary <- scenario_summary[match(scenarios, scenario_summary$scenario), , drop = FALSE]
  rownames(scenario_summary) <- NULL

  quality_profile <- if (has_quality) {
    values <- data[[quality_col]]
    data.frame(
      quality_column = quality_col,
      n_spots = length(values),
      minimum = min(values),
      q25 = unname(stats::quantile(values, 0.25)),
      median = stats::median(values),
      q75 = unname(stats::quantile(values, 0.75)),
      maximum = max(values),
      stringsAsFactors = FALSE
    )
  } else {
    NULL
  }

  decisions <- scenario_summary$reliability_decision
  interpretation <- if (all(is.na(decisions))) {
    "Scenario support is descriptive because no tolerance-based reliability decision was requested. SPARQ does not correct batch effects or impute absent tissue."
  } else if (all(decisions == "supported", na.rm = TRUE)) {
    "The result was supported under every evaluated procurement-loss scenario. This is sensitivity evidence conditional on those declared scenarios, not a correction of missing tissue or batch effects."
  } else {
    "At least one evaluated procurement-loss scenario was not supported or remained inconclusive. The result should be interpreted conditional on the available spatial support; SPARQ does not correct batch effects or impute absent tissue."
  }

  out <- list(
    result_id = result_id,
    output_type = output_type,
    reference_result = reference_result,
    assessments = assessments,
    scenario_summary = scenario_summary,
    scenario_loss_audits = lapply(assessments, `[[`, "loss_audits"),
    input_provenance = assessments[[1]]$input_provenance,
    preflight = preflight_result,
    quality_profile = quality_profile,
    settings = list(
      scenarios = scenarios,
      retention = retention,
      n_iterations = n_iterations,
      preset = preset,
      x_col = x_col,
      y_col = y_col,
      quality_col = quality_col,
      quality_direction = quality_direction,
      quality_strength = quality_strength,
      edge_strength = edge_strength,
      reference_scale = reference_scale,
      seed = seed
    ),
    interpretation = interpretation
  )
  class(out) <- "sparq_procurement_assessment"
  sparq_inform("Procurement assessment completed; no correction was applied.", verbose)
  out
}

#' @export
print.sparq_procurement_assessment <- function(x, ...) {
  cat("<SPARQ procurement assessment: ", x$result_id, ">\n", sep = "")
  cat("  Output type: ", x$output_type, "; scenarios: ",
    paste(x$scenario_summary$scenario, collapse = ", "), "\n", sep = "")
  cat("  Reference result was shared across scenarios.\n")
  for (index in seq_len(nrow(x$scenario_summary))) {
    row <- x$scenario_summary[index, , drop = FALSE]
    decision <- if (is.na(row$reliability_decision)) "not requested" else row$reliability_decision
    cat("  ", row$scenario, ": support ",
      format(row$support_quality, digits = 4), "; decision ", decision, "\n", sep = "")
  }
  invisible(x)
}

#' @export
summary.sparq_procurement_assessment <- function(object, ...) {
  object$scenario_summary
}

#' @export
as.data.frame.sparq_procurement_assessment <- function(x, row.names = NULL,
                                                       optional = FALSE, ...) {
  as.data.frame(summary(x), row.names = row.names, optional = optional, ...)
}

#' @export
plot.sparq_procurement_assessment <- function(x, ...) {
  table <- x$scenario_summary
  if (!nrow(table)) {
    stop("No procurement scenarios are available to plot.", call. = FALSE)
  }
  values <- as.numeric(table$support_quality)
  graphics::barplot(
    values,
    names.arg = table$scenario,
    ylim = c(0, max(1, values, na.rm = TRUE)),
    col = "#2F6B9A",
    ylab = "Support quality",
    main = "Support across declared procurement-loss scenarios",
    ...
  )
  invisible(x)
}

#' Export a procurement-aware SPARQ assessment
#'
#' Writes the shared reference result, scenario summary, loss audits, settings,
#' and one complete [sparq_report()] bundle per scenario. No files are written
#' when `output_dir` is `NULL`.
#'
#' @param assessment A `sparq_procurement_assessment` returned by
#'   [sparq_assess_procurement()].
#' @param output_dir Optional new directory for exported files.
#' @param include_plots Whether each scenario report should include a PDF plot.
#' @param verbose Whether to print export progress.
#'
#' @return A list containing the in-memory scenario reports and written paths.
#' @export
sparq_report_procurement <- function(assessment, output_dir = NULL,
                                     include_plots = TRUE, verbose = TRUE) {
  if (!inherits(assessment, "sparq_procurement_assessment")) {
    stop("assessment must be a sparq_procurement_assessment.", call. = FALSE)
  }
  if (!is.logical(include_plots) || length(include_plots) != 1L || is.na(include_plots)) {
    stop("include_plots must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(output_dir) && (!is.character(output_dir) || length(output_dir) != 1L ||
      is.na(output_dir) || !nzchar(output_dir))) {
    stop("output_dir must be one nonempty path or NULL.", call. = FALSE)
  }
  if (!is.null(output_dir) && dir.exists(output_dir) &&
      length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
    stop("output_dir already exists and is not empty.", call. = FALSE)
  }

  reports <- vector("list", length(assessment$assessments))
  names(reports) <- names(assessment$assessments)
  output_files <- list()
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.table(
      assessment$scenario_summary,
      file.path(output_dir, "procurement_scenario_summary.tsv"),
      sep = "\t", row.names = FALSE, quote = FALSE
    )
    output_files$scenario_summary <- file.path(output_dir, "procurement_scenario_summary.tsv")
    if (is.data.frame(assessment$quality_profile)) {
      utils::write.table(
        assessment$quality_profile,
        file.path(output_dir, "quality_profile.tsv"),
        sep = "\t", row.names = FALSE, quote = FALSE
      )
      output_files$quality_profile <- file.path(output_dir, "quality_profile.tsv")
    }
    utils::write.table(
      data.frame(setting = names(assessment$settings), value = vapply(
        assessment$settings, function(value) paste(as.character(value), collapse = ";"), character(1)
      )),
      file.path(output_dir, "procurement_settings.tsv"),
      sep = "\t", row.names = FALSE, quote = FALSE
    )
    output_files$settings <- file.path(output_dir, "procurement_settings.tsv")
    saveRDS(assessment$reference_result, file.path(output_dir, "shared_reference_result.rds"))
    output_files$reference_result <- file.path(output_dir, "shared_reference_result.rds")
    provenance_assessment <- list(
      result_id = assessment$result_id,
      output_type = assessment$output_type,
      full_result = assessment$reference_result,
      summary = assessment$scenario_summary,
      settings = assessment$settings,
      input_provenance = assessment$input_provenance
    )
    output_files <- c(
      output_files,
      sparq_write_reproducibility_bundle(provenance_assessment, output_dir)
    )
  }

  for (scenario in names(assessment$assessments)) {
    scenario_dir <- if (is.null(output_dir)) NULL else file.path(output_dir, scenario)
    reports[[scenario]] <- sparq_report(
      assessment$assessments[[scenario]], output_dir = scenario_dir,
      include_plots = include_plots, verbose = verbose
    )
  }
  if (!is.null(output_dir)) {
    saveRDS(assessment, file.path(output_dir, "procurement_assessment.rds"))
    output_files$assessment <- file.path(output_dir, "procurement_assessment.rds")
  }
  structure(
    list(
      result_id = assessment$result_id,
      scenario_summary = assessment$scenario_summary,
      scenario_reports = reports,
      reproducibility_manifest = sparq_reproducibility_manifest(list(
        result_id = assessment$result_id,
        output_type = assessment$output_type,
        full_result = assessment$reference_result,
        summary = assessment$scenario_summary,
        settings = assessment$settings,
        input_provenance = assessment$input_provenance
      )),
      output_files = output_files,
      interpretation = assessment$interpretation
    ),
    class = "sparq_procurement_report"
  )
}
