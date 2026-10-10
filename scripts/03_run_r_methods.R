#!/usr/bin/env Rscript
# Figure 4 R-native analyses and SPARQ assessments.
# Pathology annotations are retained for later validation only; none are model input.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop("Usage: 03_run_r_methods.R RUN_DIR R_ENV ITERATIONS RETENTION SEED")
run_dir <- normalizePath(args[[1]], mustWork = TRUE)
r_env <- normalizePath(args[[2]], mustWork = TRUE)
n_iter <- as.integer(args[[3]])
retention <- as.numeric(args[[4]])
seed <- as.integer(args[[5]])
lib <- file.path(run_dir, "software", "R_library")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))

need <- c("Matrix", "data.table", "SPARK", "Banksy", "progeny", "mistyR",
          "SpatialExperiment", "SingleCellExperiment", "SPARQ")
absent <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(absent)) stop("Missing R packages: ", paste(absent, collapse = ", "))

stamp <- function() format(Sys.time(), tz = "UTC", usetz = TRUE)
write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(data.table::as.data.table(x), path, sep = "\t",
    compress = if (grepl("[.]gz$", path)) "gzip" else "none", quote = FALSE)
}
ledger <- function(path, row) data.table::fwrite(data.table::as.data.table(row), path,
  sep = "\t", append = file.exists(path), col.names = !file.exists(path), quote = FALSE)
call_formals <- function(fun, values) do.call(fun, values[names(values) %in% names(formals(fun))])
first_present <- function(x, choices) {
  hit <- choices[choices %in% x]
  if (length(hit)) hit[[1]] else NA_character_
}
canonical <- file.path(run_dir, "inputs", "canonical")
raw <- file.path(run_dir, "raw", "r")
sparq_out <- file.path(run_dir, "sparq")
manifests <- file.path(run_dir, "manifests")
logs <- file.path(run_dir, "logs")
dir.create(raw, recursive = TRUE, showWarnings = FALSE)
dir.create(sparq_out, recursive = TRUE, showWarnings = FALSE)
dir.create(manifests, recursive = TRUE, showWarnings = FALSE)
sections <- sort(dirname(list.files(canonical, pattern = "section.json$", recursive = TRUE, full.names = TRUE)))

read_section <- function(folder) {
  spots <- data.table::fread(file.path(folder, "spots.tsv.gz"), data.table = FALSE)
  spots$spot_id <- as.character(spots$spot_id)
  if (file.exists(file.path(folder, "counts.tsv.gz"))) {
    tab <- data.table::fread(file.path(folder, "counts.tsv.gz"), data.table = FALSE)
    ids <- as.character(tab[[1]])
    exp <- as.matrix(tab[, -1, drop = FALSE]); storage.mode(exp) <- "numeric"
    rownames(exp) <- ids
    counts <- Matrix::Matrix(t(exp), sparse = TRUE)
    colnames(counts) <- ids; rownames(counts) <- colnames(exp)
  } else {
    counts <- Matrix::readMM(file.path(folder, "matrix.mtx.gz"))
    colnames(counts) <- as.character(data.table::fread(file.path(folder, "barcodes.tsv.gz"), header = FALSE)[[1]])
    f <- data.table::fread(file.path(folder, "features.tsv.gz"), header = FALSE)
    rownames(counts) <- make.unique(as.character(if (ncol(f) > 1L) f[[2]] else f[[1]]))
  }
  common <- intersect(colnames(counts), spots$spot_id)
  if (length(common) < 10L) stop("Fewer than 10 expression spots remain after alignment")
  counts <- counts[, common, drop = FALSE]
  spots <- spots[match(common, spots$spot_id), , drop = FALSE]
  xc <- first_present(names(spots), c("x", "pixel_col", "array_col"))
  yc <- first_present(names(spots), c("y", "pixel_row", "array_row"))
  if (is.na(xc) || is.na(yc)) stop("No coordinate columns in spots.tsv.gz")
  xy <- as.matrix(spots[, c(xc, yc), drop = FALSE]); storage.mode(xy) <- "numeric"; rownames(xy) <- common
  list(counts = counts, spots = spots, coords = xy)
}
normalise <- function(counts) {
  x <- t(t(as.matrix(counts)) / pmax(Matrix::colSums(counts), 1) * 1e4)
  log1p(x)
}
partition_banksy <- function(counts, coords, output, tag, local_seed) {
  set.seed(local_seed)
  se <- SpatialExperiment::SpatialExperiment(assays = list(counts = counts), spatialCoords = coords)
  saveRDS(se, file.path(output, paste0(tag, "_input.rds")))
  se <- call_formals(Banksy::computeBanksy, list(se = se, assay_name = "counts", M = 1L, k_geom = 15L, verbose = FALSE))
  se <- call_formals(Banksy::runBanksyPCA, list(se = se, M = 1L, lambda = 0.2, npcs = 20L, seed = local_seed))
  se <- call_formals(Banksy::clusterBanksy, list(se = se, M = 1L, lambda = 0.2, resolution = 1, seed = local_seed))
  saveRDS(se, file.path(output, paste0(tag, "_object.rds")))
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  candidate <- grep("cluster|leiden|banksy", names(cd), value = TRUE, ignore.case = TRUE)
  if (!length(candidate)) stop("BANKSY result has no detectable cluster label")
  ans <- data.frame(spot_id = colnames(se), label = as.character(cd[[tail(candidate, 1L)]]))
  write_tsv(ans, file.path(output, paste0(tag, "_labels.tsv.gz")))
}
rank_sparkx <- function(counts, coords, output, tag) {
  fit <- call_formals(SPARK::sparkx, list(count_in = counts, locus_in = coords, numCores = 1L, option = "mixture", verbose = FALSE))
  saveRDS(fit, file.path(output, paste0(tag, "_object.rds")))
  tab <- as.data.frame(fit$res_mixture)
  tab$feature_id <- rownames(tab)
  pcol <- first_present(names(tab), c("adjusted_pvalue", "adjusted.pvalue", "padj", "pvalue", "pval"))
  if (is.na(pcol)) pcol <- names(tab)[vapply(tab, is.numeric, logical(1))][1L]
  tab$ranking_score <- -log10(pmax(as.numeric(tab[[pcol]]), .Machine$double.xmin))
  tab <- tab[order(tab$ranking_score, decreasing = TRUE), , drop = FALSE]
  tab$rank <- seq_len(nrow(tab)); tab$selected <- tab$rank <= min(100L, nrow(tab))
  write_tsv(tab, file.path(output, paste0(tag, "_full_results.tsv.gz")))
  write_tsv(tab[, c("feature_id", "ranking_score", "rank", "selected")], file.path(output, paste0(tag, "_ranked_features.tsv.gz")))
}
scalar_progeny <- function(counts, output, tag) {
  fit <- call_formals(progeny::progeny, list(expr = normalise(counts), scale = TRUE,
    organism = "Human", top = 100L, perm = 1L, verbose = FALSE))
  saveRDS(fit, file.path(output, paste0(tag, "_object.rds")))
  scores <- as.data.frame(fit)
  if (nrow(scores) == ncol(counts)) rownames(scores) <- colnames(counts)
  write_tsv(cbind(spot_id = rownames(scores), scores), file.path(output, paste0(tag, "_spot_pathway_scores.tsv.gz")))
  summary <- data.frame(result_id = colnames(scores), value = vapply(scores, median, numeric(1), na.rm = TRUE),
    summary = "median_spot_pathway_activity")
  write_tsv(summary, file.path(output, paste0(tag, "_scalar_pathway_activity.tsv.gz")))
}
flatten_number <- function(x) {
  z <- suppressWarnings(as.numeric(unlist(x, recursive = TRUE, use.names = TRUE)))
  z[is.finite(z)]
}
scalar_misty <- function(counts, coords, output, tag) {
  expr <- t(normalise(counts))
  vars <- apply(expr, 2L, var)
  keep <- names(sort(vars, decreasing = TRUE))[seq_len(min(50L, length(vars)))]
  view <- mistyR::create_initial_view(as.data.frame(expr[, keep, drop = FALSE]))
  pos <- data.frame(x = coords[,1], y = coords[,2], row.names = rownames(coords))
  model <- call_formals(mistyR::add_paraview, list(misty = view, positions = pos, l = 10))
  folder <- file.path(output, paste0(tag, "_results")); dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  fit <- call_formals(mistyR::run_misty, list(misty = model, results_folder = folder, verbose = FALSE))
  saveRDS(fit, file.path(output, paste0(tag, "_object.rds")))
  gathered <- tryCatch(mistyR::collect_results(folder), error = function(e) list(error = conditionMessage(e)))
  saveRDS(gathered, file.path(output, paste0(tag, "_collected.rds")))
  components <- flatten_number(gathered)
  write_tsv(data.frame(component = names(components), value = unname(components)), file.path(output, paste0(tag, "_numeric_components.tsv.gz")))
  write_tsv(data.frame(result_id = "misty_mean_numeric_output", value = if (length(components)) mean(components) else NA_real_,
    summary = "mean_of_all_saved_numeric_misty_outputs"), file.path(output, paste0(tag, "_scalar_summary.tsv.gz")))
}
r_ledger <- file.path(manifests, "r_method_ledger.tsv")
for (folder in sections) {
  section <- basename(folder)
  d <- tryCatch(read_section(folder), error = identity)
  if (inherits(d, "error")) {
    ledger(r_ledger, data.frame(section, iteration = "setup", method = "all", status = "failed", started_at = stamp(), finished_at = stamp(), error = conditionMessage(d)))
    next
  }
  dir.create(file.path(raw, "input_objects"), recursive = TRUE, showWarnings = FALSE)
  saveRDS(d, file.path(raw, "input_objects", paste0(section, ".rds")))
  plan <- list(baseline = seq_len(ncol(d$counts)))
  for (i in seq_len(n_iter)) {
    ids <- data.table::fread(file.path(run_dir, "inputs", "perturbations", section, sprintf("iteration_%03d_retained_spot_ids.tsv.gz", i)))[[1]]
    plan[[sprintf("iteration_%03d", i)]] <- match(ids, colnames(d$counts))
  }
  for (tag in names(plan)) for (method in c("banksy", "sparkx", "progeny", "misty")) {
    began <- stamp(); output <- file.path(raw, method, section); dir.create(output, recursive = TRUE, showWarnings = FALSE)
    err <- tryCatch({
      ix <- plan[[tag]]
      if (method == "banksy") partition_banksy(d$counts[,ix,drop=FALSE], d$coords[ix,,drop=FALSE], output, tag, seed + match(tag, names(plan)))
      if (method == "sparkx") rank_sparkx(d$counts[,ix,drop=FALSE], d$coords[ix,,drop=FALSE], output, tag)
      if (method == "progeny") scalar_progeny(d$counts[,ix,drop=FALSE], output, tag)
      if (method == "misty") scalar_misty(d$counts[,ix,drop=FALSE], d$coords[ix,,drop=FALSE], output, tag)
      NULL
    }, error = identity)
    if (inherits(err, "error")) {
      writeLines(conditionMessage(err), file.path(logs, paste0("r_", method, "_", section, "_", tag, ".error.txt")))
      ledger(r_ledger, data.frame(section, iteration = tag, method, status = "failed", started_at = began, finished_at = stamp(), error = conditionMessage(err)))
    } else ledger(r_ledger, data.frame(section, iteration = tag, method, status = "ok", started_at = began, finished_at = stamp(), error = ""))
  }
}

read_out <- function(path) data.table::fread(path, data.table = FALSE)
assess <- function(method, section, output_type, full, perturbed, comparator, top_k = 100L) {
  target <- file.path(sparq_out, method, section); dir.create(target, recursive = TRUE, showWarnings = FALSE)
  began <- stamp()
  result <- tryCatch({
    f <- SPARQ::sparq_assess(full_result = full, perturbed_results = perturbed,
      comparator = comparator, result_id = paste(method, section, sep = "__"),
      reference_scale = 1, seed = seed, min_iterations = n_iter)
    saveRDS(f, file.path(target, "assessment.rds"))
    write_tsv(f$summary, file.path(target, "summary.tsv"))
    write_tsv(f$comparisons, file.path(target, "comparisons.tsv.gz"))
    write_tsv(f$failures, file.path(target, "failures.tsv.gz"))
    loc <- SPARQ::sparq_localize_fragility(full, perturbed, output_type = output_type, top_k = top_k)
    saveRDS(loc, file.path(target, "fragility_localization.rds"))
    tab <- loc$features; if (is.null(tab)) tab <- loc$spots; if (is.null(tab)) tab <- loc$edges; if (is.null(tab)) tab <- loc$global
    if (!is.null(tab)) write_tsv(tab, file.path(target, "fragility_localization.tsv.gz"))
    f$summary
  }, error = identity)
  if (inherits(result, "error")) {
    ledger(file.path(manifests, "sparq_assessment_ledger.tsv"), data.frame(method, section, output_type, status = "failed", started_at = began, finished_at = stamp(), error = conditionMessage(result)))
    return(NULL)
  }
  ledger(file.path(manifests, "sparq_assessment_ledger.tsv"), data.frame(method, section, output_type, status = "ok", started_at = began, finished_at = stamp(), error = ""))
  result
}
all_summary <- list()
for (folder in sections) {
  section <- basename(folder)
  partition_paths <- list(banksy = file.path(raw, "banksy", section),
    pca_leiden = file.path(run_dir, "raw", "python", "pca_leiden", section),
    spagcn = file.path(run_dir, "raw", "python", "spagcn", section))
  for (method in names(partition_paths)) {
    p <- partition_paths[[method]]; files <- sort(Sys.glob(file.path(p, "iteration_*_labels.tsv.gz")))
    if (file.exists(file.path(p, "baseline_labels.tsv.gz")) && length(files)) {
      b <- read_out(file.path(p, "baseline_labels.tsv.gz")); full <- setNames(as.character(b$label), b$spot_id)
      pert <- lapply(files, function(f) { x <- read_out(f); setNames(as.character(x$label), x$spot_id) })
      z <- assess(method, section, "partition", full, pert, SPARQ::sparq_compare_partition)
      if (!is.null(z)) all_summary[[length(all_summary)+1L]] <- cbind(method, section, z)
    }
  }
  ranked_paths <- list(sparkx = list(file.path(raw, "sparkx", section), "ranked_features", "feature_id"),
    spatialdm = list(file.path(run_dir, "raw", "python", "spatialdm", section), "ranked_interactions", "interaction_id"))
  for (method in names(ranked_paths)) {
    item <- ranked_paths[[method]]; p <- item[[1]]; suffix <- item[[2]]; idcol <- item[[3]]
    files <- sort(Sys.glob(file.path(p, paste0("iteration_*_", suffix, ".tsv.gz")))); base <- file.path(p, paste0("baseline_", suffix, ".tsv.gz"))
    if (file.exists(base) && length(files)) {
      full <- as.character(read_out(base)[[idcol]]); pert <- lapply(files, function(f) as.character(read_out(f)[[idcol]]))
      z <- assess(method, section, "ranked", full, pert, function(a,b) SPARQ::sparq_compare_ranked(a,b,k=100L))
      if (!is.null(z)) all_summary[[length(all_summary)+1L]] <- cbind(method, section, z)
    }
  }
  p <- file.path(run_dir, "raw", "python", "squidpy_neighborhood", section)
  files <- sort(Sys.glob(file.path(p, "iteration_*_positive_edges.tsv.gz"))); base <- file.path(p, "baseline_positive_edges.tsv.gz")
  if (file.exists(base) && length(files)) {
    z <- assess("squidpy_neighborhood", section, "graph", read_out(base)$edge_id,
      lapply(files, function(f) read_out(f)$edge_id), SPARQ::sparq_compare_graph)
    if (!is.null(z)) all_summary[[length(all_summary)+1L]] <- cbind(method = "squidpy_neighborhood", section, z)
  }
  scalars <- list(progeny = list(file.path(raw, "progeny", section), "scalar_pathway_activity"),
    misty = list(file.path(raw, "misty", section), "scalar_summary"))
  for (method in names(scalars)) {
    item <- scalars[[method]]; p <- item[[1]]; suffix <- item[[2]]; base <- file.path(p, paste0("baseline_", suffix, ".tsv.gz"))
    files <- sort(Sys.glob(file.path(p, paste0("iteration_*_", suffix, ".tsv.gz"))))
    if (!file.exists(base) || !length(files)) next
    b <- read_out(base)
    for (id in b$result_id) {
      full <- b$value[match(id, b$result_id)]
      pert <- unlist(lapply(files, function(f) { x <- read_out(f); x$value[match(id, x$result_id)] }), use.names = FALSE)
      pert <- pert[is.finite(pert)]
      if (!is.finite(full) || !length(pert)) next
      z <- assess(paste0(method, "__", id), section, "scalar", full, as.list(pert), SPARQ::sparq_compare_scalar)
      if (!is.null(z)) all_summary[[length(all_summary)+1L]] <- cbind(method = paste0(method, "__", id), section, z)
    }
  }
}
if (length(all_summary)) write_tsv(do.call(rbind, all_summary), file.path(sparq_out, "all_assessment_summaries.tsv.gz"))
info <- sessionInfo()
versions <- data.frame(package = c("R", names(info$otherPkgs)),
  version = c(as.character(getRversion()), vapply(info$otherPkgs, function(x) as.character(x$Version), character(1))))
write_tsv(versions, file.path(manifests, "r_package_versions.tsv"))
writeLines(stamp(), file.path(manifests, "r_methods_completed.txt"))

