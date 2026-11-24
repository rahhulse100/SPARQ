# Spatially aware diagnostics for discrete partition results.

.sparq_spatial_partition_edges <- function(edges) {
  if (is.matrix(edges)) edges <- as.data.frame(edges, stringsAsFactors = FALSE)
  if (!is.data.frame(edges) || ncol(edges) < 2L || !nrow(edges)) {
    stop("edges must be a nonempty two-column data frame or matrix.", call. = FALSE)
  }
  out <- data.frame(
    source = as.character(edges[[1L]]),
    target = as.character(edges[[2L]]),
    stringsAsFactors = FALSE
  )
  if (anyNA(out) || any(!nzchar(out$source)) || any(!nzchar(out$target))) {
    stop("edges must contain nonmissing, nonempty spot IDs.", call. = FALSE)
  }
  if (any(out$source == out$target)) {
    stop("Self-loop edges are not valid neighborhood edges.", call. = FALSE)
  }
  left <- pmin(out$source, out$target)
  right <- pmax(out$source, out$target)
  out <- unique(data.frame(source = left, target = right, stringsAsFactors = FALSE))
  rbind(out, data.frame(source = out$target, target = out$source, stringsAsFactors = FALSE))
}

.sparq_align_partitions <- function(full_labels, perturbed_labels) {
  a_names <- names(full_labels)
  b_names <- names(perturbed_labels)
  if (xor(is.null(a_names), is.null(b_names))) {
    stop("Name both partitions by observation ID, or neither.", call. = FALSE)
  }
  if (is.null(a_names)) {
    if (length(full_labels) != length(perturbed_labels)) {
      stop("Unequal partitions need observation-ID names for alignment.", call. = FALSE)
    }
    ids <- as.character(seq_along(full_labels))
  } else {
    if (anyNA(c(a_names, b_names)) || anyDuplicated(a_names) || anyDuplicated(b_names)) {
      stop("Partition IDs must be unique and nonmissing.", call. = FALSE)
    }
    ids <- intersect(a_names, b_names)
    full_labels <- full_labels[match(ids, a_names)]
    perturbed_labels <- perturbed_labels[match(ids, b_names)]
  }
  if (length(full_labels) < 2L || anyNA(full_labels) || anyNA(perturbed_labels)) {
    stop("At least two nonmissing aligned partition labels are required.", call. = FALSE)
  }
  list(
    ids = as.character(ids),
    full = stats::setNames(as.character(full_labels), as.character(ids)),
    perturbed = stats::setNames(as.character(perturbed_labels), as.character(ids))
  )
}

# Solve the maximum-weight one-to-one label matching problem. This is a local
# implementation of the Hungarian algorithm so the package does not require a
# new optimization dependency merely to compare spatial partitions.
.sparq_maximum_overlap_matching <- function(overlap) {
  if (!is.matrix(overlap) || !nrow(overlap) || !ncol(overlap)) {
    stop("overlap must be a nonempty matrix.", call. = FALSE)
  }
  nr <- nrow(overlap)
  nc <- ncol(overlap)
  n <- max(nr, nc)
  padded <- matrix(0, n, n)
  padded[seq_len(nr), seq_len(nc)] <- overlap
  cost <- max(padded) - padded

  # p[j + 1] is the row currently assigned to column j; column zero is the
  # augmenting-path sentinel. Indices are offset because R is one-based.
  u <- numeric(n + 1L)
  v <- numeric(n + 1L)
  p <- integer(n + 1L)
  way <- integer(n + 1L)
  for (row in seq_len(n)) {
    p[[1L]] <- row
    j0 <- 1L
    minv <- rep(Inf, n + 1L)
    used <- rep(FALSE, n + 1L)
    repeat {
      used[[j0]] <- TRUE
      i0 <- p[[j0]]
      delta <- Inf
      j1 <- NA_integer_
      for (j in 2:(n + 1L)) {
        if (!used[[j]]) {
          current <- cost[i0, j - 1L] - u[i0 + 1L] - v[[j]]
          if (current < minv[[j]]) {
            minv[[j]] <- current
            way[[j]] <- j0
          }
          if (minv[[j]] < delta) {
            delta <- minv[[j]]
            j1 <- j
          }
        }
      }
      for (j in seq_len(n + 1L)) {
        if (used[[j]]) {
          u[p[[j]] + 1L] <- u[p[[j]] + 1L] + delta
          v[[j]] <- v[[j]] - delta
        } else {
          minv[[j]] <- minv[[j]] - delta
        }
      }
      j0 <- j1
      if (p[[j0]] == 0L) break
    }
    repeat {
      j1 <- way[[j0]]
      p[[j0]] <- p[[j1]]
      j0 <- j1
      if (j0 == 1L) break
    }
  }
  assignment <- integer(n)
  for (j in 2:(n + 1L)) {
    if (p[[j]] > 0L) assignment[[p[[j]]]] <- j - 1L
  }
  assignment[seq_len(nr)]
}

.sparq_partition_label_map <- function(full_labels, perturbed_labels) {
  overlap <- table(full_labels, perturbed_labels)
  target <- rownames(overlap)
  source <- colnames(overlap)
  assignment <- .sparq_maximum_overlap_matching(unclass(overlap))
  map <- stats::setNames(rep(NA_character_, length(source)), source)
  for (row in seq_along(assignment)) {
    col <- assignment[[row]]
    if (col >= 1L && col <= length(source)) map[[col]] <- target[[row]]
  }
  # A perturbed label can be unmatched when a reference domain splits or a new
  # domain appears. Keep each such label distinct rather than folding it into
  # its largest-overlap reference label; folding would hide precisely the
  # novel adjacency created by a fragmented spatial domain.
  for (col in which(is.na(map))) {
    map[[col]] <- paste0("__sparq_unmatched__", source[[col]])
  }
  map
}

.sparq_partition_boundary_spots <- function(labels, undirected_edges) {
  differs <- labels[undirected_edges$source] != labels[undirected_edges$target]
  unique(c(undirected_edges$source[differs], undirected_edges$target[differs]))
}

.sparq_one_hop_closure <- function(ids, directed_edges) {
  if (!length(ids)) return(character())
  unique(c(ids, directed_edges$target[directed_edges$source %in% ids]))
}

.sparq_partition_adjacencies <- function(labels, undirected_edges) {
  left <- labels[undirected_edges$source]
  right <- labels[undirected_edges$target]
  differs <- left != right
  if (!any(differs)) return(character())
  left <- left[differs]
  right <- right[differs]
  paste(pmin(left, right), pmax(left, right), sep = "\r")
}

#' Compare two spatial partitions using boundaries and domain adjacency
#'
#' This optional comparator supplements global adjusted Rand index (ARI) with
#' two spatial diagnostics. Boundary agreement asks whether spots at a domain
#' interface remain at, or within one neighborhood hop of, an interface after
#' a native rerun. Domain-adjacency agreement asks whether pairs of touching
#' domains remain touching after a one-to-one maximum-overlap matching of
#' rerun labels to reference labels. Unmatched rerun domains remain distinct,
#' so a domain split introduces a detectable new adjacency.
#'
#' The comparator aligns partitions by their shared named spot IDs, then
#' induces the supplied spatial graph on those spots. It therefore works when
#' a perturbation rerun contains only retained spots. No pathology labels or
#' clean-reference result are used.
#'
#' `instability` is one minus a declared weighted mean of ARI, tolerant
#' boundary F1, and adjacency Jaccard. The weights are deliberately explicit:
#' users should pre-specify them for a study rather than tune them against
#' held-out outcomes.
#'
#' @param full_labels Named reference partition labels.
#' @param perturbed_labels Named rerun partition labels.
#' @param edges A two-column spot-neighborhood edge table or matrix. Edges are
#'   treated as undirected and may contain each relationship once or twice.
#' @param weights Named nonnegative weights for `ari`, `boundary_f1`, and
#'   `adjacency_jaccard`. They must sum to a positive value and are normalized
#'   internally. The default gives each component equal weight.
#'
#' @return A one-row data frame containing ARI, tolerant boundary precision,
#' recall and F1, domain-adjacency Jaccard, component weights, and the combined
#' nonnegative `instability`.
#' @export
sparq_compare_spatial_partition <- function(
    full_labels, perturbed_labels, edges,
    weights = c(ari = 1, boundary_f1 = 1, adjacency_jaccard = 1)) {
  aligned <- .sparq_align_partitions(full_labels, perturbed_labels)
  graph <- .sparq_spatial_partition_edges(edges)
  graph <- graph[graph$source %in% aligned$ids & graph$target %in% aligned$ids, , drop = FALSE]
  if (!nrow(graph)) {
    stop("No supplied neighborhood edges remain after partition alignment.", call. = FALSE)
  }
  undirected <- graph[graph$source < graph$target, , drop = FALSE]
  if (!nrow(undirected)) {
    stop("No undirected neighborhood edges remain after partition alignment.", call. = FALSE)
  }
  needed <- c("ari", "boundary_f1", "adjacency_jaccard")
  if (is.null(names(weights)) || !setequal(names(weights), needed) ||
      any(!is.finite(weights)) || any(weights < 0) || sum(weights) <= 0) {
    stop("weights must be named ari, boundary_f1, and adjacency_jaccard, with nonnegative values and positive sum.", call. = FALSE)
  }
  weights <- weights[needed] / sum(weights)

  ari <- .sparq_partition_ari_component(aligned$full, aligned$perturbed)$ari
  full_boundary <- .sparq_partition_boundary_spots(aligned$full, undirected)
  perturbed_boundary <- .sparq_partition_boundary_spots(aligned$perturbed, undirected)
  full_reachable <- .sparq_one_hop_closure(full_boundary, graph)
  perturbed_reachable <- .sparq_one_hop_closure(perturbed_boundary, graph)
  boundary_precision <- if (!length(perturbed_boundary)) {
    if (!length(full_boundary)) 1 else 0
  } else mean(perturbed_boundary %in% full_reachable)
  boundary_recall <- if (!length(full_boundary)) {
    if (!length(perturbed_boundary)) 1 else 0
  } else mean(full_boundary %in% perturbed_reachable)
  boundary_f1 <- if (boundary_precision + boundary_recall == 0) 0 else {
    2 * boundary_precision * boundary_recall / (boundary_precision + boundary_recall)
  }

  label_map <- .sparq_partition_label_map(aligned$full, aligned$perturbed)
  remapped_perturbed <- unname(label_map[aligned$perturbed])
  names(remapped_perturbed) <- names(aligned$perturbed)
  full_adjacency <- unique(.sparq_partition_adjacencies(aligned$full, undirected))
  perturbed_adjacency <- unique(.sparq_partition_adjacencies(remapped_perturbed, undirected))
  adjacency_union <- union(full_adjacency, perturbed_adjacency)
  adjacency_jaccard <- if (!length(adjacency_union)) 1 else {
    length(intersect(full_adjacency, perturbed_adjacency)) / length(adjacency_union)
  }
  agreement <- c(ari = ari, boundary_f1 = boundary_f1, adjacency_jaccard = adjacency_jaccard)
  instability <- 1 - sum(weights * agreement)
  data.frame(
    ari = ari,
    boundary_precision = boundary_precision,
    boundary_recall = boundary_recall,
    boundary_f1 = boundary_f1,
    adjacency_jaccard = adjacency_jaccard,
    ari_weight = unname(weights[["ari"]]),
    boundary_f1_weight = unname(weights[["boundary_f1"]]),
    adjacency_jaccard_weight = unname(weights[["adjacency_jaccard"]]),
    instability = instability,
    n_aligned = length(aligned$ids),
    n_reference_boundary_spots = length(full_boundary),
    n_rerun_boundary_spots = length(perturbed_boundary),
    stringsAsFactors = FALSE
  )
}

#' Compare spatial-domain partitions in their native spatial output space
#'
#' This is SPARQ's public comparator for partition-valued spatial analyses.
#' It requires the analysis's declared native neighborhood graph and combines
#' global label agreement, tolerant boundary agreement, and domain-adjacency
#' agreement. It does not use pathology or an external clean reference.
#'
#' @inheritParams sparq_compare_spatial_partition
#' @return A one-row spatial partition comparison table. See
#'   [sparq_compare_spatial_partition()] for its columns.
#' @export
sparq_compare_partition <- function(
    full_labels, perturbed_labels, edges,
    weights = c(ari = 1, boundary_f1 = 1, adjacency_jaccard = 1)) {
  sparq_compare_spatial_partition(
    full_labels = full_labels,
    perturbed_labels = perturbed_labels,
    edges = edges,
    weights = weights
  )
}
