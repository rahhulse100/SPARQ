sparq_validate_results <- function(data,
                                   required_cols = c("result_id", "theta_full", "theta_pert", "scenario", "iteration")) {
  missing <- setdiff(required_cols, colnames(data))
  if (length(missing)) {
    stop("Missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  if (!is.numeric(data$theta_full)) {
    stop("theta_full must be numeric.", call. = FALSE)
  }

  if (!is.numeric(data$theta_pert)) {
    stop("theta_pert must be numeric.", call. = FALSE)
  }

  invisible(TRUE)
}

sparq_validate_columns <- function(data, required_cols) {
  missing <- setdiff(required_cols, names(data))
  if (length(missing) > 0) {
    stop(
      "SPARQ input missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

sparq_safe_cor <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  if (length(unique(x[ok])) < 2) return(NA_real_)
  if (length(unique(y[ok])) < 2) return(NA_real_)

  suppressWarnings(
    stats::cor(x[ok], y[ok], method = method)
  )
}

sparq_mode_sign_fraction <- function(x) {
  x <- x[is.finite(x)]
  x <- sign(x[x != 0])
  if (length(x) == 0) return(0)
  max(mean(x > 0), mean(x < 0))
}


sparq_topk_jaccard <- function(full_ids, pert_ids, k = 100) {
  full_ids <- as.character(full_ids)
  pert_ids <- as.character(pert_ids)

  full_top <- unique(full_ids[seq_len(min(k, length(full_ids)))])
  pert_top <- unique(pert_ids[seq_len(min(k, length(pert_ids)))])

  if (length(full_top) == 0 && length(pert_top) == 0) return(NA_real_)

  length(intersect(full_top, pert_top)) / length(union(full_top, pert_top))
}



sparq_topk_replacement <- function(full_ids, pert_ids, k = 100) {
  full_ids <- as.character(full_ids)
  pert_ids <- as.character(pert_ids)

  full_top <- unique(full_ids[seq_len(min(k, length(full_ids)))])
  pert_top <- unique(pert_ids[seq_len(min(k, length(pert_ids)))])

  if (length(full_top) == 0) return(NA_real_)

  1 - (length(intersect(full_top, pert_top)) / length(full_top))
}

# User-facing progress and result-display helpers. These deliberately use base
# R so a standard assessment remains lightweight and usable on a cluster.

sparq_inform <- function(message, verbose = TRUE) {
  if (isTRUE(verbose)) {
    message("SPARQ | ", message)
  }
  invisible(NULL)
}

sparq_progress <- function(completed, total, verbose = TRUE,
                           prefix = "Perturbations") {
  if (!isTRUE(verbose) || total < 1L) {
    return(invisible(NULL))
  }

  width <- 24L
  filled <- min(width, floor(width * completed / total))
  bar <- paste0(
    "[",
    paste(rep("=", filled), collapse = ""),
    paste(rep(".", width - filled), collapse = ""),
    "]"
  )
  message(
    "SPARQ | ", prefix, " ", bar, " ",
    completed, "/", total
  )
  invisible(NULL)
}

# Internal provenance helpers. Fingerprints are hashes rather than copies of
# sensitive input values; row identifiers are retained separately only when an
# assessment was run directly from spot-level input data.

sparq_input_provenance <- function(data) {
  if (!is.data.frame(data) || !nrow(data)) {
    return(NULL)
  }
  row_ids <- rownames(data)
  if (is.null(row_ids) || length(row_ids) != nrow(data) ||
      anyNA(row_ids) || any(!nzchar(row_ids)) || anyDuplicated(row_ids)) {
    row_ids <- as.character(seq_len(nrow(data)))
    row_id_source <- "generated_row_indices"
  } else {
    row_id_source <- "data_frame_row_names"
  }
  list(
    n_observations = nrow(data),
    n_columns = ncol(data),
    column_names = names(data),
    column_name_fingerprint = sparq_fingerprint(names(data)),
    row_id_source = row_id_source,
    row_ids = as.character(row_ids),
    row_id_fingerprint = sparq_fingerprint(as.character(row_ids)),
    data_fingerprint = sparq_fingerprint(data)
  )
}

sparq_git_commit <- function() {
  from_environment <- c(
    Sys.getenv("SPARQ_GIT_COMMIT", unset = ""),
    Sys.getenv("GITHUB_SHA", unset = "")
  )
  from_environment <- from_environment[nzchar(from_environment)]
  if (length(from_environment)) {
    return(from_environment[[1]])
  }

  description <- tryCatch(
    suppressWarnings(utils::packageDescription("SPARQ")),
    error = function(e) NULL
  )
  if (is.list(description)) {
    fields <- c("RemoteSha", "GithubSHA1", "GitHubSHA1", "Commit")
    values <- unlist(description[intersect(fields, names(description))], use.names = FALSE)
    values <- values[!is.na(values) & nzchar(values)]
    if (length(values)) {
      return(as.character(values[[1]]))
    }
  }

  candidates <- c(getwd(), tryCatch(find.package("SPARQ"), error = function(e) ""))
  candidates <- candidates[nzchar(candidates)]
  for (candidate in candidates) {
    path <- normalizePath(candidate, mustWork = FALSE)
    repeat {
      if (dir.exists(file.path(path, ".git"))) {
        commit <- tryCatch(
          suppressWarnings(system2("git", c("-C", path, "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)),
          error = function(e) character()
        )
        if (length(commit) && nzchar(commit[[1]])) {
          return(commit[[1]])
        }
      }
      parent <- dirname(path)
      if (identical(parent, path)) {
        break
      }
      path <- parent
    }
  }
  "unavailable"
}

sparq_loaded_package_versions <- function() {
  packages <- sort(unique(loadedNamespaces()))
  versions <- vapply(packages, function(package) {
    tryCatch(as.character(utils::packageVersion(package)), error = function(e) NA_character_)
  }, character(1))
  data.frame(
    package = packages,
    version = unname(versions),
    stringsAsFactors = FALSE
  )
}

sparq_reproducibility_manifest <- function(assessment) {
  settings <- if (is.list(assessment$settings)) assessment$settings else list()
  setting_rows <- data.frame(
    key = paste0("setting.", names(settings)),
    value = vapply(settings, function(value) {
      if (is.null(value)) return("")
      if (is.atomic(value)) return(paste(as.character(value), collapse = ";"))
      paste(utils::capture.output(str(value, give.attr = FALSE)), collapse = " ")
    }, character(1)),
    stringsAsFactors = FALSE
  )
  input <- assessment$input_provenance
  input_rows <- if (is.list(input)) {
    data.frame(
      key = c(
        "input.n_observations", "input.n_columns", "input.row_id_source",
        "input.row_id_fingerprint", "input.column_name_fingerprint",
        "input.data_fingerprint"
      ),
      value = c(
        input$n_observations, input$n_columns, input$row_id_source,
        input$row_id_fingerprint, input$column_name_fingerprint,
        input$data_fingerprint
      ),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      key = c("input.n_observations", "input.n_columns", "input.row_id_source",
        "input.row_id_fingerprint", "input.column_name_fingerprint", "input.data_fingerprint"),
      value = rep("unavailable", 6L), stringsAsFactors = FALSE
    )
  }
  summary_fingerprint <- if (is.data.frame(assessment$summary)) {
    sparq_fingerprint(assessment$summary)
  } else {
    "unavailable"
  }
  comparison_fingerprint <- if (is.data.frame(assessment$comparisons)) {
    sparq_fingerprint(assessment$comparisons)
  } else {
    "unavailable"
  }
  package_description <- tryCatch(
    suppressWarnings(utils::packageDescription("SPARQ")),
    error = function(e) NULL
  )
  package_version <- if (!is.list(package_description)) {
    "source_or_uninstalled"
  } else {
    as.character(package_description$Version)
  }
  result_id <- if (is.null(assessment$result_id)) "unavailable" else assessment$result_id
  output_type <- if (is.null(assessment$output_type)) "unavailable" else assessment$output_type
  core <- data.frame(
    key = c(
      "manifest_schema", "created_at_utc", "package.name", "package.version",
      "package.git_commit", "r.version", "platform", "result_id", "output_type",
      "reference_result_fingerprint", "summary_fingerprint", "comparisons_fingerprint"
    ),
    value = c(
      "sparq_reproducibility_manifest_v1", format(Sys.time(), tz = "UTC", usetz = TRUE),
      "SPARQ", package_version, sparq_git_commit(), R.version.string,
      R.version$platform, as.character(result_id), as.character(output_type),
      if (is.null(assessment$full_result)) "unavailable" else sparq_fingerprint(assessment$full_result),
      summary_fingerprint, comparison_fingerprint
    ),
    stringsAsFactors = FALSE
  )
  out <- rbind(core, input_rows, setting_rows)
  rownames(out) <- NULL
  out
}

sparq_write_reproducibility_bundle <- function(assessment, output_dir) {
  manifest <- sparq_reproducibility_manifest(assessment)
  manifest_path <- file.path(output_dir, "reproducibility_manifest.tsv")
  package_path <- file.path(output_dir, "package_versions.tsv")
  session_path <- file.path(output_dir, "session_info.txt")
  utils::write.table(manifest, manifest_path, sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.table(sparq_loaded_package_versions(), package_path,
    sep = "\t", row.names = FALSE, quote = FALSE)
  writeLines(utils::capture.output(sessionInfo()), session_path)
  paths <- list(
    reproducibility_manifest = manifest_path,
    package_versions = package_path,
    session_info = session_path
  )
  input <- assessment$input_provenance
  if (is.list(input) && is.character(input$row_ids) && length(input$row_ids)) {
    input_ids_path <- file.path(output_dir, "input_row_ids.tsv")
    utils::write.table(
      data.frame(row_id = input$row_ids, stringsAsFactors = FALSE),
      input_ids_path, sep = "\t", row.names = FALSE, quote = FALSE
    )
    paths$input_row_ids <- input_ids_path
  }
  paths
}

sparq_assessment_status <- function(summary_table) {
  if (!is.data.frame(summary_table) || nrow(summary_table) != 1L) {
    return("unavailable")
  }
  if (isTRUE(summary_table$insufficient_support[1])) {
    return("insufficient support")
  }
  "adequate support"
}

sparq_print_assessment <- function(x, label = "assessment") {
  summary_table <- x$summary
  if (!is.data.frame(summary_table) || nrow(summary_table) != 1L) {
    cat("<SPARQ ", label, ">\n", sep = "")
    return(invisible(x))
  }

  result_id <- if (!is.null(x$result_id)) as.character(x$result_id) else "result"
  output_type <- if (!is.null(x$output_type)) as.character(x$output_type) else "custom"
  cat("<SPARQ ", label, ": ", result_id, ">\n", sep = "")
  cat("  Output type: ", output_type, "\n", sep = "")
  cat(
    "  Support quality: ", format(summary_table$support_quality[1], digits = 4),
    " (", sparq_assessment_status(summary_table), ")\n", sep = ""
  )
  cat(
    "  Completed perturbations: ", summary_table$n_iterations[1],
    "/", summary_table$n_requested[1],
    "; failed: ", summary_table$n_failed[1], "\n", sep = ""
  )
  cat(
    "  Median instability: ", format(summary_table$instability_median[1], digits = 4),
    "; reproducibility: ", format(summary_table$reproducibility_score[1], digits = 4),
    "\n", sep = ""
  )
  if (isTRUE(summary_table$insufficient_support[1])) {
    cat("  Reason: ", summary_table$support_reason[1], "\n", sep = "")
  }
  if (is.data.frame(x$reliability_decision) && nrow(x$reliability_decision) == 1L) {
    cat(
      "  Reliability decision: ", x$reliability_decision$decision[1],
      " (", format(x$reliability_decision$lower_confidence_bound[1], digits = 3),
      "-", format(x$reliability_decision$upper_confidence_bound[1], digits = 3),
      " confidence interval)\n", sep = ""
    )
  }
  invisible(x)
}

print.sparq_assessment <- function(x, ...) {
  sparq_print_assessment(x, label = "assessment")
}

print.sparq_workflow_assessment <- function(x, ...) {
  sparq_print_assessment(x, label = "workflow assessment")
}

print.sparq_script_assessment <- function(x, ...) {
  sparq_print_assessment(x, label = "script assessment")
}

print.sparq_adaptive_assessment <- function(x, ...) {
  sparq_print_assessment(x, label = "adaptive assessment")
  if (is.data.frame(x$adaptive_decision) && nrow(x$adaptive_decision) == 1L) {
    cat(
      "  Adaptive decision: ", x$adaptive_decision$decision[1],
      " (", x$adaptive_decision$stop_reason[1], ")\n", sep = ""
    )
  }
  invisible(x)
}

print.sparq_report <- function(x, ...) {
  cat("<SPARQ report: ", x$result_id, ">\n", sep = "")
  if (is.data.frame(x$support_summary) && nrow(x$support_summary) == 1L) {
    cat(
      "  Support quality: ", format(x$support_summary$support_quality[1], digits = 4),
      "\n", sep = ""
    )
  }
  if (length(x$output_files)) {
    cat("  Exported files: ", length(x$output_files), "\n", sep = "")
  } else {
    cat("  Exported files: none\n")
  }
  invisible(x)
}
