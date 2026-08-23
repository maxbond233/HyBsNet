#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: validate_current_figure3.R PROJECT_ROOT", call. = FALSE)
}
project_root <- normalizePath(args[[1L]], mustWork = TRUE)
suppressPackageStartupMessages(library(HyBsNet))

fail <- function(...) stop(paste0(...), call. = FALSE)
expect_identical_set <- function(x, y, label) {
  if (!identical(sort(unique(as.character(x))), sort(unique(as.character(y))))) {
    fail(label, " differs from the frozen baseline.")
  }
}
edge_key <- function(edges, from = "Cluster1", to = "Cluster2") {
  a <- as.character(edges[[from]])
  b <- as.character(edges[[to]])
  paste(pmin(a, b), pmax(a, b), sep = "\r")
}

deg_raw <- read_hybs_csv(file.path(
  project_root,
  "results", "current", "01_differential_expression", "MAST_deg_summary.csv"
))
baseline_dir <- file.path(
  project_root,
  "Figures", "manuscript", "drafts", "Fig3_refactored", "source_data"
)
expected <- list(
  Obesity = list(panel = "A", nodes = 96L, edges = 632L),
  Diabetes = list(panel = "B", nodes = 99L, edges = 602L)
)
networks <- list()
degs <- list()

for (condition in names(expected)) {
  spec <- expected[[condition]]
  deg <- prepare_deg_sets(deg_raw, condition = condition, padj_cutoff = 0.05)
  frozen_edges <- read_hybs_csv(file.path(
    project_root,
    "results", "current", "02_similarity_network",
    paste0("deg_similarity_all_clusters_jaccard_", tolower(condition), ".csv")
  ))
  network <- build_similarity_network(
    frozen_edges,
    deg_data = deg,
    min_similarity = 0.10,
    condition = condition
  )
  if (nrow(network_nodes(network)) != spec$nodes ||
      nrow(network_edges(network)) != spec$edges) {
    fail("Figure 3", spec$panel, " node/edge counts differ from baseline.")
  }

  baseline_nodes <- read_hybs_csv(file.path(
    baseline_dir, paste0("Fig3", spec$panel, "_nodes.csv")
  ))
  baseline_edges <- read_hybs_csv(file.path(
    baseline_dir, paste0("Fig3", spec$panel, "_edges.csv")
  ))
  expect_identical_set(
    network_nodes(network)$name,
    baseline_nodes$name,
    paste0("Figure 3", spec$panel, " nodes")
  )
  expect_identical_set(
    edge_key(network_edges(network)),
    edge_key(baseline_edges),
    paste0("Figure 3", spec$panel, " edges")
  )

  recalculated <- calculate_pairwise_jaccard(deg, min_similarity = 0.10)
  frozen_clean <- prepare_similarity_edges(frozen_edges, min_similarity = 0.10)
  expect_identical_set(
    edge_key(recalculated),
    edge_key(frozen_clean),
    paste0("Recalculated ", condition, " edges")
  )
  comparison <- merge(
    recalculated[, c("Cluster1", "Cluster2", "Similarity")],
    frozen_clean[, c("Cluster1", "Cluster2", "Similarity")],
    by = c("Cluster1", "Cluster2"),
    suffixes = c("_new", "_frozen")
  )
  if (nrow(comparison) != nrow(frozen_clean) ||
      max(abs(comparison$Similarity_new - comparison$Similarity_frozen)) > 1e-12) {
    fail("Recalculated ", condition, " similarities differ from frozen values.")
  }
  networks[[condition]] <- network
  degs[[condition]] <- deg
  message(
    "Fig3", spec$panel, " OK: ", spec$nodes, " nodes, ", spec$edges, " edges"
  )
}

local <- extract_local_network(
  networks$Diabetes,
  center_node = "SH_IN_10_ESR1",
  order = 1,
  label_top_n = 15
)
if (nrow(network_nodes(local)) != 28L || nrow(network_edges(local)) != 180L) {
  fail("Figure 3G node/edge counts differ from baseline.")
}
baseline_local_nodes <- read_hybs_csv(file.path(
  baseline_dir, "Fig3G_SH_IN_10_ESR1_nodes.csv"
))
baseline_local_edges <- read_hybs_csv(file.path(
  baseline_dir, "Fig3G_SH_IN_10_ESR1_edges.csv"
))
expect_identical_set(network_nodes(local)$name, baseline_local_nodes$name, "Figure 3G nodes")
expect_identical_set(
  edge_key(network_edges(local)),
  edge_key(baseline_local_edges, "from", "to"),
  "Figure 3G edges"
)
if (sum(!is.na(network_nodes(local)$label)) != 16L) {
  fail("Figure 3G labelled-node count differs from baseline.")
}
message("Fig3G OK: 28 nodes, 180 edges, 16 labelled nodes")

core_rows <- 0L
for (condition in names(degs)) {
  centrality <- read_hybs_csv(file.path(
    project_root,
    "results", "current", "03_community",
    paste0("jaccard_centrality_measures_filtered_", tolower(condition), ".csv")
  ))
  prepared <- prepare_core_periphery_data(
    centrality,
    degs[[condition]],
    condition = condition
  )
  if (anyNA(prepared$total_degs)) {
    fail("Figure 3F contains missing total_degs values.")
  }
  core_rows <- core_rows + nrow(prepared)
}
if (core_rows != 172L) {
  fail("Figure 3F row count differs from baseline: ", core_rows)
}
message("Fig3F OK: 172 rows, no missing total_degs")
message("BASELINE_OK")
