# High-level assessment, fragility localization, and reporting helpers.

sparq_assess_workflow <-
function (data, analysis_function, output_type = c("scalar",
    "ranked", "partition", "graph"), stress_model = c("uniform_random",
    "contiguous_hole", "none", "custom"), retention = NULL, n_iterations = NULL,
    top_k = 100, x_col = NULL, y_col = NULL, custom_function = NULL,
    result_id = "result", reference_scale = 1, seed = 1, cache_dir = NULL,
    resume = TRUE, cache_key = NULL, failure_action = c("record",
        "stop"), verbose = NULL, progress_every = NULL, preset = "standard",
    min_iterations = NULL, instability_bootstrap_B = NULL,
    reliability_tolerance = NULL, required_reliability = 0.9,
    reliability_confidence = 0.95, loss_model = NULL,
    full_result = NULL, keep_perturbed_results = FALSE, ...)
{
    if (!is.data.frame(data) || !nrow(data)) {
        stop("data must be a nonempty data.frame.", call. = FALSE)
    }
    if (!is.function(analysis_function)) {
        stop("analysis_function must be a function.", call. = FALSE)
    }
    output_type <- match.arg(output_type)
    stress_model <- match.arg(stress_model)
    failure_action <- match.arg(failure_action)
    sparq_check_number(reference_scale, "reference_scale")
    sparq_check_number(seed, "seed", 0, integer = TRUE)
    if (output_type == "ranked") {
        sparq_check_number(top_k, "top_k", 1, integer = TRUE)
    }
    comparator <- switch(output_type, scalar = sparq_compare_scalar,
        ranked = function(full_result, perturbed_result) {
            sparq_compare_ranked(full_rank = full_result, perturbed_rank = perturbed_result,
                k = top_k)
        }, partition = sparq_compare_partition, graph = sparq_compare_graph)
    sparq_inform(paste0("Preparing a ", output_type, " workflow assessment."), verbose)
    fit <- sparq_run(data = data, analysis_function = analysis_function,
        comparator = comparator, stress_model = stress_model,
        retention = retention, n_iterations = n_iterations, result_id = result_id,
        x_col = x_col, y_col = y_col, custom_function = custom_function,
        reference_scale = reference_scale, seed = seed, cache_dir = cache_dir,
        resume = resume, cache_key = cache_key, failure_action = failure_action,
        verbose = verbose, progress_every = progress_every, preset = preset,
        min_iterations = min_iterations,
        instability_bootstrap_B = instability_bootstrap_B,
        reliability_tolerance = reliability_tolerance,
        required_reliability = required_reliability,
        reliability_confidence = reliability_confidence,
        loss_model = loss_model,
        full_result = full_result,
        keep_perturbed_results = keep_perturbed_results,
        ...)
    fit$output_type <- output_type
    fit$workflow <- list(output_type = output_type, stress_model = stress_model,
        retention = fit$settings$retention, n_iterations = fit$settings$n_iterations,
        preset = fit$settings$preset, top_k = if (output_type ==
            "ranked") {
            top_k
        } else {
            NA_integer_
        }, x_col = x_col, y_col = y_col, seed = seed,
        loss_model = loss_model)
    fit$fragility <- NULL
    class(fit) <- unique(c("sparq_workflow_assessment", class(fit)))
    fit
}
sparq_localize_fragility <-
function (full_result, perturbed_results, output_type = c("scalar",
    "ranked", "partition", "graph"), top_k = 100, min_agreement = 0.80000000000000004)
{
    output_type <- match.arg(output_type)
    sparq_check_number(min_agreement, "min_agreement", 0, 1)
    if (output_type == "ranked") {
        sparq_check_number(top_k, "top_k", 1, integer = TRUE)
    }
    if (!is.list(perturbed_results)) {
        if (output_type == "scalar") {
            perturbed_results <- as.list(perturbed_results)
        }
        else {
            stop("perturbed_results must be a list for ranked, partition, and graph outputs.",
                call. = FALSE)
        }
    }
    perturbed_results <- Filter(Negate(is.null), perturbed_results)
    if (!length(perturbed_results)) {
        stop("perturbed_results contains no completed perturbation outputs.",
            call. = FALSE)
    }
    mode_value <- function(x) {
        x <- as.character(x)
        x <- x[!is.na(x)]
        if (!length(x)) {
            return(NA_character_)
        }
        counts <- table(x)
        candidates <- names(counts)[counts == max(counts)]
        sort(candidates)[1]
    }
    normalized_entropy <- function(x) {
        x <- as.character(x)
        x <- x[!is.na(x)]
        if (length(unique(x)) < 2L) {
            return(0)
        }
        probabilities <- as.numeric(table(x)/length(x))
        -sum(probabilities * log(probabilities))/log(length(probabilities))
    }
    map_partition_to_reference <- function(reference_labels,
        perturbed_labels) {
        common_ids <- intersect(names(reference_labels), names(perturbed_labels))
        if (length(common_ids) < 2L) {
            return(setNames(rep(NA_character_, length(perturbed_labels)),
                names(perturbed_labels)))
        }
        reference_common <- as.character(reference_labels[common_ids])
        perturbed_common <- as.character(perturbed_labels[common_ids])
        overlap <- table(reference_common, perturbed_common)
        mapping <- vapply(colnames(overlap), function(perturbed_label) {
            counts <- overlap[, perturbed_label, drop = TRUE]
            candidates <- rownames(overlap)[counts == max(counts)]
            sort(candidates)[1]
        }, character(1))
        mapped <- unname(mapping[as.character(perturbed_labels)])
        names(mapped) <- names(perturbed_labels)
        mapped
    }
    if (output_type == "scalar") {
        full_value <- as.numeric(full_result)
        if (length(full_value) != 1L || !is.finite(full_value)) {
            stop("full_result must be one finite numeric value for scalar output.",
                call. = FALSE)
        }
        perturbed_values <- suppressWarnings(as.numeric(unlist(perturbed_results,
            use.names = FALSE)))
        perturbed_values <- perturbed_values[is.finite(perturbed_values)]
        if (!length(perturbed_values)) {
            stop("No finite scalar perturbed results were supplied.",
                call. = FALSE)
        }
        delta <- perturbed_values - full_value
        instability <- abs(delta)
        global <- data.frame(full_result = full_value, n_perturbations = length(perturbed_values),
            signed_displacement_mean = mean(delta), instability_median = stats::median(instability),
            instability_mean = mean(instability), instability_sd = if (length(instability) >
                1L) {
                stats::sd(instability)
            }
            else {
                NA_real_
            }, stringsAsFactors = FALSE)
        out <- list(output_type = "scalar", min_agreement = min_agreement,
            global = global, localization = NULL, interpretation = paste("Scalar outputs have no feature- or spot-level localization.",
                "Use the global perturbation summary."))
        class(out) <- "sparq_fragility"
        return(out)
    }
    if (output_type == "ranked") {
        full_rank <- as.character(full_result)
        if (!length(full_rank) || anyNA(full_rank) || anyDuplicated(full_rank)) {
            stop("full_result must be a nonempty ranked vector of unique feature IDs.",
                call. = FALSE)
        }
        perturbed_ranks <- lapply(perturbed_results, function(x) {
            x <- as.character(x)
            if (!length(x) || anyNA(x) || anyDuplicated(x)) {
                stop("Each perturbed ranked result must contain unique feature IDs.",
                  call. = FALSE)
            }
            x
        })
        feature_ids <- sort(unique(c(full_rank, unlist(perturbed_ranks,
            use.names = FALSE))))
        effective_top_k <- min(top_k, length(full_rank))
        full_rank_position <- match(feature_ids, full_rank)
        rank_matrix <- sapply(perturbed_ranks, function(x) {
            match(feature_ids, x)
        })
        if (is.null(dim(rank_matrix))) {
            rank_matrix <- matrix(rank_matrix, ncol = 1L)
        }
        top_k_matrix <- sapply(perturbed_ranks, function(x) {
            feature_ids %in% head(x, effective_top_k)
        })
        if (is.null(dim(top_k_matrix))) {
            top_k_matrix <- matrix(top_k_matrix, ncol = 1L)
        }
        selection_frequency <- rowMeans(top_k_matrix)
        median_rank <- apply(rank_matrix, 1, function(x) {
            x <- x[is.finite(x)]
            if (!length(x)) {
                NA_real_
            }
            else {
                stats::median(x)
            }
        })
        mean_rank <- apply(rank_matrix, 1, function(x) {
            x <- x[is.finite(x)]
            if (!length(x)) {
                NA_real_
            }
            else {
                mean(x)
            }
        })
        present_frequency <- rowMeans(is.finite(rank_matrix))
        order_index <- order(-selection_frequency, median_rank,
            full_rank_position, feature_ids, na.last = TRUE)
        consensus_rank <- integer(length(feature_ids))
        consensus_rank[order_index] <- seq_along(order_index)
        feature_table <- data.frame(feature_id = feature_ids,
            full_rank = full_rank_position, reference_top_k = feature_ids %in%
                head(full_rank, effective_top_k), perturbation_presence_frequency = present_frequency,
            top_k_selection_frequency = selection_frequency,
            median_perturbed_rank = median_rank, mean_perturbed_rank = mean_rank,
            consensus_rank = consensus_rank, consensus_top_k = consensus_rank <=
                effective_top_k, stable_reference_feature = (feature_ids %in%
                head(full_rank, effective_top_k)) & (selection_frequency >=
                min_agreement), stringsAsFactors = FALSE)
        feature_table <- feature_table[order(feature_table$consensus_rank),
            , drop = FALSE]
        out <- list(output_type = "ranked", top_k = effective_top_k,
            min_agreement = min_agreement, features = feature_table,
            consensus_ranked_features = feature_table$feature_id[feature_table$consensus_top_k],
            stable_reference_features = feature_table$feature_id[feature_table$stable_reference_feature],
            interpretation = paste("top_k_selection_frequency is the fraction of perturbations",
                "in which a feature remained in the top-ranked set."))
        class(out) <- "sparq_fragility"
        return(out)
    }
    if (output_type == "partition") {
        full_labels <- as.character(full_result)
        if (is.null(names(full_result)) || anyNA(names(full_result)) ||
            anyDuplicated(names(full_result)) || anyNA(full_labels)) {
            stop("Partition full_result must be a named vector with unique observation IDs.",
                call. = FALSE)
        }
        names(full_labels) <- names(full_result)
        perturbed_partitions <- lapply(perturbed_results, function(x) {
            labels <- as.character(x)
            if (is.null(names(x)) || anyNA(names(x)) || anyDuplicated(names(x)) ||
                anyNA(labels)) {
                stop("Each perturbed partition must be named by unique observation IDs.",
                  call. = FALSE)
            }
            names(labels) <- names(x)
            map_partition_to_reference(reference_labels = full_labels,
                perturbed_labels = labels)
        })
        spot_table <- lapply(names(full_labels), function(spot_id) {
            mapped_labels <- vapply(perturbed_partitions, function(x) {
                if (spot_id %in% names(x)) {
                  x[[spot_id]]
                }
                else {
                  NA_character_
                }
            }, character(1))
            observed <- mapped_labels[!is.na(mapped_labels)]
            reference_label <- full_labels[[spot_id]]
            consensus_label <- mode_value(observed)
            data.frame(observation_id = spot_id, reference_label = reference_label,
                consensus_label = consensus_label, n_observed = length(observed),
                retention_frequency = (length(observed)/length(perturbed_partitions)),
                reference_agreement = if (length(observed)) {
                  mean(observed == reference_label)
                }
                else {
                  NA_real_
                }, consensus_support = if (length(observed)) {
                  mean(observed == consensus_label)
                }
                else {
                  NA_real_
                }, label_entropy = normalized_entropy(observed),
                fragile = if (length(observed)) {
                  mean(observed == reference_label) < min_agreement
                }
                else {
                  NA
                }, stringsAsFactors = FALSE)
        })
        spot_table <- do.call(rbind, spot_table)
        out <- list(output_type = "partition", min_agreement = min_agreement,
            spots = spot_table, consensus_partition = setNames(spot_table$consensus_label,
                spot_table$observation_id), fragile_spots = spot_table[isTRUE(spot_table$fragile) |
                spot_table$fragile %in% TRUE, , drop = FALSE],
            interpretation = paste("Perturbed labels are mapped to the reference label with",
                "which they overlap most before spot-level agreement is calculated."))
        class(out) <- "sparq_fragility"
        return(out)
    }
    full_edges <- unique(as.character(full_result))
    if (!length(full_edges) || anyNA(full_edges)) {
        stop("full_result must be a nonempty vector of canonical edge IDs.",
            call. = FALSE)
    }
    perturbed_graphs <- lapply(perturbed_results, function(x) {
        x <- unique(as.character(x))
        if (anyNA(x)) {
            stop("Perturbed graph edge IDs cannot contain missing values.",
                call. = FALSE)
        }
        x
    })
    edge_ids <- sort(unique(c(full_edges, unlist(perturbed_graphs,
        use.names = FALSE))))
    edge_matrix <- sapply(perturbed_graphs, function(x) {
        edge_ids %in% x
    })
    if (is.null(dim(edge_matrix))) {
        edge_matrix <- matrix(edge_matrix, ncol = 1L)
    }
    edge_table <- data.frame(edge_id = edge_ids, reference_edge = edge_ids %in%
        full_edges, n_selected = rowSums(edge_matrix), selection_frequency = rowMeans(edge_matrix),
        consensus_edge = rowMeans(edge_matrix) >= min_agreement,
        fragile_reference_edge = (edge_ids %in% full_edges) &
            (rowMeans(edge_matrix) < min_agreement), stringsAsFactors = FALSE)
    edge_table <- edge_table[order(-edge_table$selection_frequency,
        edge_table$edge_id), , drop = FALSE]
    out <- list(output_type = "graph", min_agreement = min_agreement,
        edges = edge_table, consensus_edges = edge_table$edge_id[edge_table$consensus_edge],
        fragile_reference_edges = edge_table$edge_id[edge_table$fragile_reference_edge],
        interpretation = paste("Edge selection frequency is calculated across all perturbations.",
            "It is conservative when an edge cannot be evaluated because",
            "one or both nodes were removed."))
    class(out) <- "sparq_fragility"
    out
}
sparq_report <-
function (assessment, output_dir = NULL, include_plots = TRUE, verbose = TRUE)
{
    if (!is.list(assessment)) {
        stop("assessment must be a SPARQ assessment list.", call. = FALSE)
    }
    required_fields <- c("full_result", "comparisons", "summary",
        "failures", "settings")
    missing_fields <- setdiff(required_fields, names(assessment))
    if (length(missing_fields)) {
        stop(paste("assessment is missing:", paste(missing_fields,
            collapse = ", ")), call. = FALSE)
    }
    if (!is.data.frame(assessment$comparisons) || !is.data.frame(assessment$summary) ||
        !is.data.frame(assessment$failures)) {
        stop(paste("assessment comparisons, summary, and failures",
            "must be data frames."), call. = FALSE)
    }
    if (!is.logical(include_plots) || length(include_plots) !=
        1L || is.na(include_plots)) {
        stop("include_plots must be TRUE or FALSE.", call. = FALSE)
    }
    if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
        stop("verbose must be TRUE or FALSE.", call. = FALSE)
    }
    sparq_inform("Preparing an assessment report.", verbose)
    report_assessment <- assessment
    fragility_note <- NULL
    if (is.null(report_assessment$fragility)) {
        can_localize <- (exists("sparq_localize_fragility", mode = "function") &&
            !is.null(report_assessment$output_type) && !is.null(report_assessment$perturbed_results))
        if (can_localize) {
            workflow_top_k <- 100L
            if (!is.null(report_assessment$workflow) && is.finite(report_assessment$workflow$top_k)) {
                workflow_top_k <- as.integer(report_assessment$workflow$top_k)
            }
            report_assessment$fragility <- tryCatch(sparq_localize_fragility(full_result = report_assessment$full_result,
                perturbed_results = report_assessment$perturbed_results,
                output_type = report_assessment$output_type,
                top_k = workflow_top_k), error = function(error) {
                fragility_note <<- conditionMessage(error)
                NULL
            })
        }
        else {
            fragility_note <- paste("Fragility localization was not calculated because",
                "output_type or perturbed_results is unavailable.")
        }
    }
    summary_table <- report_assessment$summary
    comparison_table <- report_assessment$comparisons
    failure_table <- report_assessment$failures
    loss_audits <- report_assessment$loss_audits
    reliability_decision <- report_assessment$reliability_decision
    reliability_sequence <- report_assessment$reliability_confidence_sequence
    settings_table <- data.frame(setting = names(report_assessment$settings),
        value = vapply(report_assessment$settings, function(x) {
            if (is.null(x)) {
                return("")
            }
            if (is.atomic(x)) {
                return(paste(as.character(x), collapse = ";"))
            }
            paste(utils::capture.output(str(x, give.attr = FALSE)),
                collapse = " ")
        }, character(1)), stringsAsFactors = FALSE)
    output_type <- if (!is.null(report_assessment$output_type)) {
        as.character(report_assessment$output_type)
    }
    else {
        NA_character_
    }
    fragility_tables <- list()
    if (is.list(report_assessment$fragility)) {
        possible_tables <- c("global", "features", "spots", "edges")
        for (table_name in possible_tables) {
            if (is.data.frame(report_assessment$fragility[[table_name]])) {
                fragility_tables[[table_name]] <- report_assessment$fragility[[table_name]]
            }
        }
    }
    draw_summary_plot <- function() {
        component_names <- c("iteration_support", "bias_support",
            "uncertainty_support", "reproducibility_support")
        component_names <- component_names[component_names %in%
            names(summary_table)]
        graphics::layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE))
        old_par <- graphics::par(no.readonly = TRUE)
        on.exit(graphics::par(old_par), add = TRUE)
        valid_comparisons <- comparison_table[comparison_table$status ==
            "ok" & is.finite(comparison_table$instability), ,
            drop = FALSE]
        if (nrow(valid_comparisons)) {
            graphics::plot(valid_comparisons$iteration, valid_comparisons$instability,
                pch = 16, col = "#2F6B9A", xlab = "Perturbation iteration",
                ylab = "Instability", main = "Instability across perturbations")
            graphics::abline(h = stats::median(valid_comparisons$instability),
                lty = 2, col = "#B13A3A")
        }
        else {
            graphics::plot.new()
            graphics::title(main = "Instability across perturbations")
            graphics::text(0.5, 0.5, "No successful perturbation comparisons")
        }
        if (length(component_names)) {
            component_values <- as.numeric(summary_table[1, component_names,
                drop = TRUE])
            graphics::barplot(component_values, ylim = c(0, 1),
                col = c("#5B8DB8", "#F2A65A", "#73A857", "#9A6FB0")[seq_along(component_values)],
                names.arg = gsub("_support", "", component_names),
                ylab = "Component value", main = "Support-quality components")
            graphics::abline(h = 1, lty = 3, col = "grey50")
        }
        else {
            graphics::plot.new()
            graphics::title(main = "Support-quality components")
        }
        if (output_type == "ranked" && !is.null(fragility_tables$features)) {
            graphics::hist(fragility_tables$features$top_k_selection_frequency,
                breaks = seq(0, 1, by = 0.050000000000000003),
                xlim = c(0, 1), col = "#8CB9D6", border = "white",
                xlab = "Top-k selection frequency", main = "Ranked-feature stability")
            graphics::abline(v = 0.80000000000000004, lty = 2,
                col = "#B13A3A")
        }
        else if (output_type == "partition" && !is.null(fragility_tables$spots)) {
            graphics::hist(fragility_tables$spots$reference_agreement,
                breaks = seq(0, 1, by = 0.050000000000000003),
                xlim = c(0, 1), col = "#8CB9D6", border = "white",
                xlab = "Reference-label agreement", main = "Spatial-domain stability")
            graphics::abline(v = 0.80000000000000004, lty = 2,
                col = "#B13A3A")
        }
        else if (output_type == "graph" && !is.null(fragility_tables$edges)) {
            graphics::hist(fragility_tables$edges$selection_frequency,
                breaks = seq(0, 1, by = 0.050000000000000003),
                xlim = c(0, 1), col = "#8CB9D6", border = "white",
                xlab = "Edge selection frequency", main = "Graph-edge stability")
            graphics::abline(v = 0.80000000000000004, lty = 2,
                col = "#B13A3A")
        }
        else {
            graphics::plot.new()
            graphics::title(main = "Assessment summary")
            support_value <- if ("support_quality" %in% names(summary_table)) {
                summary_table$support_quality[1]
            }
            else {
                NA_real_
            }
            graphics::text(0.5, 0.59999999999999998, paste0("Output type: ",
                ifelse(is.na(output_type), "not recorded", output_type)))
            graphics::text(0.5, 0.45000000000000001, paste0("Support quality: ",
                format(support_value, digits = 4)))
            graphics::text(0.5, 0.29999999999999999, paste0("Completed perturbations: ",
                summary_table$n_iterations[1]))
        }
    }
    output_files <- list()
    if (!is.null(output_dir)) {
        if (!is.character(output_dir) || length(output_dir) !=
            1L || is.na(output_dir) || !nzchar(output_dir)) {
            stop("output_dir must be one nonempty path or NULL.",
                call. = FALSE)
        }
        if (dir.exists(output_dir) && length(list.files(output_dir,
            all.files = TRUE, no.. = TRUE))) {
            stop("output_dir already exists and is not empty.",
                call. = FALSE)
        }
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        write_table <- function(table, filename) {
            path <- file.path(output_dir, filename)
            utils::write.table(table, path, sep = "\t", row.names = FALSE,
                quote = FALSE)
            path
        }
        output_files$summary <- write_table(summary_table, "support_summary.tsv")
        output_files$comparisons <- write_table(comparison_table,
            "perturbation_comparisons.tsv")
        output_files$failures <- write_table(failure_table, "failed_iterations.tsv")
        output_files$settings <- write_table(settings_table,
            "settings.tsv")
        if (is.data.frame(reliability_decision)) {
            output_files$reliability_decision <- write_table(reliability_decision,
                "reliability_decision.tsv")
        }
        if (is.data.frame(reliability_sequence)) {
            output_files$reliability_confidence_sequence <- write_table(reliability_sequence,
                "reliability_confidence_sequence.tsv")
        }
        if (is.list(loss_audits) && length(loss_audits)) {
            audit_tables <- lapply(seq_along(loss_audits), function(iteration) {
                audit <- loss_audits[[iteration]]
                if (!is.data.frame(audit)) {
                    return(NULL)
                }
                audit$iteration <- iteration
                audit
            })
            audit_tables <- Filter(Negate(is.null), audit_tables)
            if (length(audit_tables)) {
                output_files$spot_loss_audit <- write_table(sparq_bind_rows(audit_tables),
                    "spot_loss_audit.tsv")
            }
        }
        for (table_name in names(fragility_tables)) {
            output_files[[paste0("fragility_", table_name)]] <- write_table(fragility_tables[[table_name]],
                paste0("fragility_", table_name, ".tsv"))
        }
        output_files$reference_result <- file.path(output_dir,
            "reference_result.rds")
        saveRDS(report_assessment$full_result, output_files$reference_result)
        if (!is.null(report_assessment$perturbed_results)) {
            output_files$perturbed_results <- file.path(output_dir,
                "perturbed_results.rds")
            saveRDS(report_assessment$perturbed_results, output_files$perturbed_results)
        }
        if (include_plots) {
            output_files$plot <- file.path(output_dir, "support_report.pdf")
            grDevices::pdf(output_files$plot, width = 10, height = 7)
            draw_summary_plot()
            grDevices::dev.off()
        }
    }
    if (include_plots && is.null(output_dir)) {
        draw_summary_plot()
    }
    report <- list(result_id = report_assessment$result_id, output_type = output_type,
        reference_result = report_assessment$full_result, support_summary = summary_table,
        comparisons = comparison_table, failures = failure_table,
        settings = report_assessment$settings, fragility = report_assessment$fragility,
        fragility_tables = fragility_tables, fragility_note = fragility_note,
        loss_audits = loss_audits, reliability_decision = reliability_decision,
        reliability_confidence_sequence = reliability_sequence,
        output_files = output_files)
    class(report) <- "sparq_report"
    if (!is.null(output_dir)) {
        report_path <- file.path(output_dir, "sparq_report.rds")
        report$output_files$report <- report_path
        saveRDS(report, report_path)
        sparq_inform(paste0("Report written to ", normalizePath(output_dir)), verbose)
    }
    if (is.null(output_dir)) {
        sparq_inform("Report returned in memory; no files were written.", verbose)
    }
    report
}
