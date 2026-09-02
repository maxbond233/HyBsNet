#!/usr/bin/env Rscript

# HyBsNet v1 fixed-layout views.
#
# Edit this list and Source the complete file in RStudio. All major objects are
# intentionally retained in the Global Environment for inspection.
network_view_config <- list(
  project_root = Sys.getenv("HYBS_FIG3_ROOT", unset = "/path/to/HyBs_Fig3"),
  condition = "Diabetes",
  panel = "Fig3B",
  gene = "HLA-E",
  pathway = "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  subnetwork = list(
    mode = "ego",
    center_node = "SH_IN_10_ESR1",
    order = 1L,
    edge_mode = "induced"
  ),
  local_layout = "frozen", # or "compact"
  padj_cutoff = 0.05,
  abs_log2fc_cutoff = log2(1.5),
  jaccard_cutoff = 0.10,
  layout_id = "log2fc_1p5_community_first_v1",
  layout_seed = 42L,
  output_dir = {
    value <- Sys.getenv("HYBS_NETWORK_VIEW_OUTPUT", unset = "")
    if (nzchar(value)) value else NULL
  },
  overwrite = FALSE
)

suppressPackageStartupMessages(library(HyBsNet))

project_root <- normalizePath(network_view_config$project_root, mustWork = TRUE)
condition_slug <- tolower(network_view_config$condition)
workflow_root <- file.path(
  project_root,
  "manual_rebuild",
  "Fig3_network_log2fc_1p5"
)
source_data_dir <- file.path(workflow_root, "outputs", "source_data")
redesign_dir <- file.path(workflow_root, "visual_redesign", "outputs")
output_dir <- if (is.null(network_view_config$output_dir)) {
  file.path(workflow_root, "network_views_v1", condition_slug)
} else {
  network_view_config$output_dir
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

deg_path <- file.path(
  project_root,
  "results", "current", "01_differential_expression", "MAST_deg_summary.csv"
)
gsea_path <- file.path(dirname(project_root), "Macaque_HyBs_GSEA", "gsea_all_results.csv")
edge_path <- file.path(
  source_data_dir,
  paste0(network_view_config$panel, "_log2fc_1p5_edges.csv")
)
frozen_node_path <- file.path(
  redesign_dir,
  paste0(condition_slug, "_redesign_nodes.csv")
)

deg_raw <- read_hybs_csv(
  deg_path,
  c(
    "gene", "avg_log2FC", "p_val_adj", "group", "nuclei", "cluster",
    "pct.1"
  )
)
deg_for_network <- deg_raw[
  is.finite(deg_raw$avg_log2FC) &
    abs(deg_raw$avg_log2FC) >= network_view_config$abs_log2fc_cutoff,
  ,
  drop = FALSE
]
prepared_deg <- prepare_deg_sets(
  deg_for_network,
  condition = network_view_config$condition,
  padj_cutoff = network_view_config$padj_cutoff
)
edge_data <- read_hybs_csv(
  edge_path,
  c("Cluster1", "Cluster2", "Similarity")
)
network <- build_similarity_network(
  edge_data,
  deg_data = prepared_deg,
  min_similarity = network_view_config$jaccard_cutoff,
  condition = network_view_config$condition
)

frozen_nodes <- read_hybs_csv(
  frozen_node_path,
  c("name", "x_redesign", "y_redesign", "graph_signature")
)
frozen_layout <- as_network_layout(
  frozen_nodes,
  network,
  x_col = "x_redesign",
  y_col = "y_redesign",
  layout_id = network_view_config$layout_id,
  layout_method = "community_first_frozen",
  layout_seed = network_view_config$layout_seed,
  # One-time migration from legacy CSV signatures generated before export.
  # The newly saved layout below receives the current reproducible signature.
  strict_signature = FALSE
)

gene_overlay <- prepare_gene_overlay(
  deg_raw,
  gene = network_view_config$gene,
  condition = network_view_config$condition,
  padj_cutoff = network_view_config$padj_cutoff,
  min_abs_score = network_view_config$abs_log2fc_cutoff
)
gene_plot <- plot_network_overlay(
  network,
  gene_overlay,
  layout = frozen_layout,
  labels = network_view_config$subnetwork$center_node,
  size_by = "overlay",
  title = paste(network_view_config$condition, network_view_config$gene)
)

subnetwork <- do.call(
  extract_subnetwork,
  c(list(network = network), network_view_config$subnetwork)
)
subnetwork_layout <- calculate_subnetwork_layout(
  subnetwork,
  mode = network_view_config$local_layout,
  parent_layout = if (network_view_config$local_layout == "frozen") {
    frozen_layout
  } else {
    NULL
  },
  seed = network_view_config$layout_seed
)
subnetwork_gene_plot <- plot_subnetwork(
  subnetwork,
  subnetwork_layout,
  overlay = gene_overlay,
  title = paste(network_view_config$gene, "subnetwork")
)
subnetwork_context_plot <- plot_subnetwork_context(
  network,
  frozen_layout,
  subnetwork,
  overlay = gene_overlay,
  title = paste(network_view_config$gene, "network context")
)

pathway_overlay <- NULL
pathway_plot <- NULL
if (!is.null(network_view_config$pathway) && nzchar(network_view_config$pathway)) {
  gsea_data <- read_hybs_csv(
    gsea_path,
    c("group", "nuclei", "cluster", "pathway", "NES", "padj")
  )
  pathway_overlay <- prepare_pathway_overlay(
    gsea_data,
    pathway = network_view_config$pathway,
    condition = network_view_config$condition,
    method = "preranked GSEA; verify ranked-universe provenance"
  )
  pathway_plot <- plot_network_overlay(
    network,
    pathway_overlay,
    layout = frozen_layout,
    size_by = "deg_count",
    title = paste(network_view_config$condition, network_view_config$pathway)
  )
}

export_network_panel(
  gene_plot,
  output_dir,
  paste0(condition_slug, "_gene_", gsub("[^A-Za-z0-9]+", "_", network_view_config$gene)),
  width_mm = 180,
  height_mm = 150,
  formats = c("pdf", "png", "svg"),
  overwrite = network_view_config$overwrite
)
export_network_panel(
  subnetwork_gene_plot,
  output_dir,
  paste0(condition_slug, "_gene_subnetwork"),
  width_mm = 125,
  height_mm = 110,
  formats = c("pdf", "png", "svg"),
  overwrite = network_view_config$overwrite
)
export_network_panel(
  subnetwork_context_plot,
  output_dir,
  paste0(condition_slug, "_subnetwork_context"),
  width_mm = 150,
  height_mm = 130,
  formats = c("pdf", "png", "svg"),
  overwrite = network_view_config$overwrite
)
if (!is.null(pathway_plot)) {
  export_network_panel(
    pathway_plot,
    output_dir,
    paste0(condition_slug, "_pathway_", gsub("[^A-Za-z0-9]+", "_", network_view_config$pathway)),
    width_mm = 180,
    height_mm = 150,
    formats = c("pdf", "png", "svg"),
    overwrite = network_view_config$overwrite
  )
}

source_output_dir <- file.path(output_dir, "source_data")
dir.create(source_output_dir, recursive = TRUE, showWarnings = FALSE)
write_overlay <- function(data, path, overwrite = FALSE) {
  if (file.exists(path) && !overwrite) {
    stop("Refusing to overwrite existing file: ", path, call. = FALSE)
  }
  utils::write.csv(as.data.frame(data), path, row.names = FALSE, na = "")
  invisible(path)
}
save_network_layout(
  frozen_layout,
  file.path(source_output_dir, paste0(condition_slug, "_frozen_layout.csv")),
  overwrite = network_view_config$overwrite
)
export_network_source_data(
  subnetwork,
  subnetwork_layout,
  source_output_dir,
  paste0(condition_slug, "_subnetwork"),
  overwrite = network_view_config$overwrite
)
write_overlay(
  gene_overlay,
  file.path(source_output_dir, paste0(condition_slug, "_gene_overlay.csv")),
  overwrite = network_view_config$overwrite
)
if (!is.null(pathway_overlay)) {
  write_overlay(
    pathway_overlay,
    file.path(source_output_dir, paste0(condition_slug, "_pathway_overlay.csv")),
    overwrite = network_view_config$overwrite
  )
}

message("HyBsNet fixed-layout network views completed: ", output_dir)
