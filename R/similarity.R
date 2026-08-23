#' Filter and canonicalize a similarity edge table
#'
#' @param edges Data frame containing `Cluster1`, `Cluster2`, and `Similarity`.
#' @param min_similarity Strict similarity cutoff.
#' @return A canonical undirected edge table with one row per unordered pair.
#' @export
prepare_similarity_edges <- function(edges, min_similarity = 0.10) {
  .assert_scalar_number(min_similarity, "min_similarity", lower = 0, upper = 1)
  .assert_columns(edges, c("Cluster1", "Cluster2", "Similarity"), "similarity edges")
  out <- as.data.frame(edges, stringsAsFactors = FALSE)
  out$Cluster1 <- as.character(out$Cluster1)
  out$Cluster2 <- as.character(out$Cluster2)
  out$Similarity <- suppressWarnings(as.numeric(out$Similarity))
  keep <- is.finite(out$Similarity) & out$Similarity > min_similarity &
    !is.na(out$Cluster1) & !is.na(out$Cluster2) &
    nzchar(out$Cluster1) & nzchar(out$Cluster2) & out$Cluster1 != out$Cluster2
  out <- out[keep, , drop = FALSE]
  if (nrow(out) == 0L) {
    .hybs_abort("No similarity edges remained after filtering.")
  }

  swap <- out$Cluster1 > out$Cluster2
  old1 <- out$Cluster1
  out$Cluster1[swap] <- out$Cluster2[swap]
  out$Cluster2[swap] <- old1[swap]
  if (all(c("Nuclei1", "Nuclei2") %in% names(out))) {
    old_n1 <- out$Nuclei1
    out$Nuclei1[swap] <- out$Nuclei2[swap]
    out$Nuclei2[swap] <- old_n1[swap]
  }
  pair_key <- paste(out$Cluster1, out$Cluster2, sep = "\r")
  if (anyDuplicated(pair_key)) {
    duplicates <- unique(pair_key[duplicated(pair_key)])
    .hybs_abort(
      "Similarity edges contain duplicate unordered pairs; examples: ",
      paste(utils::head(duplicates, 3L), collapse = ", "), "."
    )
  }
  out <- out[order(-out$Similarity, out$Cluster1, out$Cluster2), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Calculate pairwise Jaccard similarity between DEG sets
#'
#' The default parameters reproduce the scientific boundary used for Figure 3:
#' DEG records are prepared separately with adjusted `P < 0.05`, and only edges
#' with Jaccard similarity strictly greater than `0.10` are returned.
#'
#' @param deg_data Output from [prepare_deg_sets()].
#' @param min_similarity Strict Jaccard cutoff.
#' @param cross_region_only If `TRUE`, retain only pairs from different regions.
#' @return An edge table including overlap direction categories.
#' @export
calculate_pairwise_jaccard <- function(
    deg_data,
    min_similarity = 0.10,
    cross_region_only = FALSE) {
  .assert_columns(
    deg_data,
    c("gene", "avg_log2FC", "nuclei", "region_cluster"),
    "prepared DEG data"
  )
  .assert_scalar_number(min_similarity, "min_similarity", lower = 0, upper = 1)
  .assert_flag(cross_region_only, "cross_region_only")

  clusters <- sort(unique(as.character(deg_data$region_cluster)))
  if (length(clusters) < 2L) {
    .hybs_abort("At least two subclusters are required for pairwise similarity.")
  }
  split_rows <- split(seq_len(nrow(deg_data)), deg_data$region_cluster)
  genes <- lapply(split_rows, function(i) as.character(deg_data$gene[i]))
  logfc <- lapply(split_rows, function(i) {
    stats::setNames(as.numeric(deg_data$avg_log2FC[i]), as.character(deg_data$gene[i]))
  })

  pairs <- utils::combn(clusters, 2L, simplify = FALSE)
  rows <- vector("list", length(pairs))
  used <- 0L
  for (pair in pairs) {
    c1 <- pair[[1L]]
    c2 <- pair[[2L]]
    n1 <- parse_region(c1)
    n2 <- parse_region(c2)
    if (cross_region_only && identical(n1, n2)) {
      next
    }
    g1 <- genes[[c1]]
    g2 <- genes[[c2]]
    common <- intersect(g1, g2)
    similarity <- length(common) / length(union(g1, g2))
    if (!is.finite(similarity) || similarity <= min_similarity) {
      next
    }
    fc1 <- logfc[[c1]][common]
    fc2 <- logfc[[c2]][common]
    both_up <- common[fc1 > 0 & fc2 > 0]
    both_down <- common[fc1 < 0 & fc2 < 0]
    opposite <- common[fc1 * fc2 < 0]
    used <- used + 1L
    rows[[used]] <- data.frame(
      Cluster1 = c1,
      Cluster2 = c2,
      Similarity = similarity,
      Nuclei1 = n1,
      Nuclei2 = n2,
      Both_Up_Genes = paste(both_up, collapse = "; "),
      Both_Down_Genes = paste(both_down, collapse = "; "),
      Opposite_Genes = paste(opposite, collapse = "; "),
      Both_Up_Count = length(both_up),
      Both_Down_Count = length(both_down),
      Opposite_Count = length(opposite),
      Total_Overlap = length(both_up) + length(both_down) + length(opposite),
      stringsAsFactors = FALSE
    )
  }
  if (used == 0L) {
    .hybs_abort("No Jaccard edges exceeded `min_similarity`.")
  }
  out <- do.call(rbind, rows[seq_len(used)])
  out <- out[order(-out$Similarity, out$Cluster1, out$Cluster2), , drop = FALSE]
  rownames(out) <- NULL
  out
}
