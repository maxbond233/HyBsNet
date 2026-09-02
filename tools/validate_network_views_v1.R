#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: validate_network_views_v1.R /path/to/HyBs_Fig3", call. = FALSE)
}

project_root <- normalizePath(args[[1L]], mustWork = TRUE)
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))

deg_path <- file.path(
  project_root,
  "results", "current", "01_differential_expression", "MAST_deg_summary.csv"
)
deg_raw <- read_hybs_csv(
  deg_path,
  c(
    "gene", "avg_log2FC", "p_val_adj", "group", "nuclei", "cluster",
    "pct.1"
  )
)
workflow_root <- file.path(
  project_root,
  "manual_rebuild",
  "Fig3_network_log2fc_1p5"
)
source_dir <- file.path(workflow_root, "outputs", "source_data")
redesign_dir <- file.path(workflow_root, "visual_redesign", "outputs")
expected <- list(
  Obesity = list(panel = "Fig3A", nodes = 96L, edges = 699L),
  Diabetes = list(panel = "Fig3B", nodes = 99L, edges = 485L)
)

qa <- lapply(names(expected), function(condition) {
  settings <- expected[[condition]]
  prepared <- prepare_deg_sets(
    deg_raw[
      is.finite(deg_raw$avg_log2FC) & abs(deg_raw$avg_log2FC) >= log2(1.5),
      ,
      drop = FALSE
    ],
    condition = condition,
    padj_cutoff = 0.05
  )
  edges <- read_hybs_csv(
    file.path(source_dir, paste0(settings$panel, "_log2fc_1p5_edges.csv")),
    c("Cluster1", "Cluster2", "Similarity")
  )
  network <- build_similarity_network(
    edges,
    prepared,
    min_similarity = 0.10,
    condition = condition
  )
  frozen_nodes <- read_hybs_csv(
    file.path(redesign_dir, paste0(tolower(condition), "_redesign_nodes.csv")),
    c("name", "x_redesign", "y_redesign", "graph_signature")
  )
  layout <- as_network_layout(
    frozen_nodes,
    network,
    x_col = "x_redesign",
    y_col = "y_redesign",
    layout_id = "log2fc_1p5_community_first_v1",
    strict_signature = FALSE
  )
  layout_check <- compare_network_layouts(
    data.frame(
      name = frozen_nodes$name,
      x = frozen_nodes$x_redesign,
      y = frozen_nodes$y_redesign
    ),
    layout
  )
  local <- extract_subnetwork(
    network,
    mode = "ego",
    center_node = "SH_IN_10_ESR1",
    order = 1L
  )
  local_layout <- calculate_subnetwork_layout(
    local,
    mode = "frozen",
    parent_layout = layout
  )
  gene_overlay <- prepare_gene_overlay(
    deg_raw,
    gene = "HLA-E",
    condition = condition
  )
  plot <- plot_network_overlay(
    network,
    gene_overlay,
    layout = layout,
    labels = "SH_IN_10_ESR1"
  )
  render_path <- tempfile(fileext = ".png")
  ggplot2::ggsave(render_path, plot, width = 100, height = 85, units = "mm", dpi = 100)
  data.frame(
    condition = condition,
    nodes = nrow(network_nodes(network)),
    expected_nodes = settings$nodes,
    edges = nrow(network_edges(network)),
    expected_edges = settings$edges,
    legacy_signature_present =
      length(unique(frozen_nodes$graph_signature)) == 1L &&
      nzchar(unique(frozen_nodes$graph_signature)),
    current_signature_assigned =
      unique(layout$graph_signature) == network_signature(network),
    frozen_coordinates_identical =
      isTRUE(attr(layout_check, "all_shared_coordinates_identical")),
    local_layout_valid = tryCatch({
      validate_network_layout(local_layout, local)
      TRUE
    }, error = function(e) FALSE),
    render_nonempty = file.exists(render_path) && file.info(render_path)$size > 0,
    stringsAsFactors = FALSE
  )
})
qa <- do.call(rbind, qa)
qa$pass <- with(
  qa,
  nodes == expected_nodes & edges == expected_edges & legacy_signature_present &
    current_signature_assigned &
    frozen_coordinates_identical & local_layout_valid & render_nonempty
)
print(qa, row.names = FALSE)
if (!all(qa$pass)) {
  stop("HyBsNet v1 real-data regression failed.", call. = FALSE)
}
message("HyBsNet v1 real-data regression passed.")
