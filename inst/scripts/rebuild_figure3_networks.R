#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop(
    "Usage: rebuild_figure3_networks.R PROJECT_ROOT NEW_OUTPUT_DIR [STYLE_YAML]",
    call. = FALSE
  )
}

project_root <- normalizePath(args[[1L]], mustWork = TRUE)
output_dir <- normalizePath(args[[2L]], mustWork = FALSE)
style_path <- if (length(args) >= 3L) {
  normalizePath(args[[3L]], mustWork = TRUE)
} else {
  system.file("config", "figure3_style.yml", package = "HyBsNet")
}
analysis_path <- system.file("config", "figure3_analysis.yml", package = "HyBsNet")
if (!nzchar(analysis_path) || !nzchar(style_path)) {
  stop("HyBsNet configuration files were not found.", call. = FALSE)
}

suppressPackageStartupMessages(library(HyBsNet))
analysis_config <- read_network_config(analysis_path)
style <- style_from_config(read_network_config(style_path))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
source_data_dir <- file.path(output_dir, "source_data")
dir.create(source_data_dir, recursive = TRUE, showWarnings = FALSE)

deg_path <- file.path(
  project_root,
  "results", "current", "01_differential_expression", "MAST_deg_summary.csv"
)
deg_raw <- read_hybs_csv(
  deg_path,
  c("gene", "avg_log2FC", "p_val_adj", "group", "nuclei", "cluster")
)

networks <- list()
prepared_degs <- list()
run_summary <- list()
baseline_source_dir <- file.path(
  project_root,
  "Figures", "manuscript", "drafts", "Fig3_refactored", "source_data"
)
for (condition in names(analysis_config$panels)) {
  panel <- analysis_config$panels[[condition]]
  deg <- prepare_deg_sets(
    deg_raw,
    condition = condition,
    padj_cutoff = analysis_config$analysis$deg_padj_cutoff
  )
  edge_path <- file.path(
    project_root,
    "results", "current", "02_similarity_network",
    paste0("deg_similarity_all_clusters_jaccard_", tolower(condition), ".csv")
  )
  edges <- read_hybs_csv(
    edge_path,
    c("Cluster1", "Cluster2", "Similarity", "Total_Overlap")
  )
  network <- build_similarity_network(
    edges,
    deg_data = deg,
    min_similarity = analysis_config$analysis$jaccard_cutoff,
    condition = condition
  )
  baseline_layout_path <- file.path(
    baseline_source_dir, paste0("Fig3", panel$panel, "_nodes.csv")
  )
  if (file.exists(baseline_layout_path)) {
    layout <- read_hybs_csv(baseline_layout_path, c("name", "x", "y"))
    validate_network_layout(layout, network)
    layout$graph_signature <- network_signature(network)
    layout$layout_method <- "cached_current"
    layout$layout_seed <- analysis_config$layout$seed
  } else {
    layout <- calculate_network_layout(
      network,
      method = analysis_config$layout$method,
      seed = analysis_config$layout$seed,
      iterations = analysis_config$layout$iterations
    )
  }
  labels <- unlist(panel$labels, use.names = FALSE)
  plot <- plot_similarity_network(
    network,
    layout = layout,
    labels = labels,
    style = style,
    title = paste0("DEG similarity network (Jaccard) - ", condition),
    seed = analysis_config$layout$seed
  )
  stem <- paste0(
    "Fig3", panel$panel, "_similarity_network_", tolower(condition)
  )
  export_network_panel(
    plot,
    output_dir = output_dir,
    stem = stem,
    width_mm = analysis_config$output$full_network_width_mm,
    height_mm = analysis_config$output$full_network_height_mm,
    dpi = analysis_config$output$dpi,
    formats = unlist(analysis_config$output$formats, use.names = FALSE)
  )
  export_network_source_data(
    network,
    layout,
    output_dir = source_data_dir,
    stem = paste0("Fig3", panel$panel)
  )
  networks[[condition]] <- network
  prepared_degs[[condition]] <- deg
  run_summary[[condition]] <- list(
    panel = panel$panel,
    nodes = nrow(network_nodes(network)),
    edges = nrow(network_edges(network)),
    graph_signature = network_signature(network)
  )
}

local_config <- analysis_config$local_network
local_network <- extract_local_network(
  networks[[local_config$condition]],
  center_node = local_config$center_node,
  order = local_config$neighborhood_order,
  label_top_n = local_config$label_top_n
)
local_baseline_path <- file.path(
  baseline_source_dir,
  paste0("Fig3", local_config$panel, "_", local_config$center_node, "_nodes.csv")
)
if (file.exists(local_baseline_path)) {
  local_layout <- read_hybs_csv(local_baseline_path, c("name", "x", "y"))
  validate_network_layout(local_layout, local_network)
  local_layout$graph_signature <- network_signature(local_network)
  local_layout$layout_method <- "cached_current"
  local_layout$layout_seed <- analysis_config$layout$seed
} else {
  local_layout <- calculate_network_layout(
    local_network,
    method = analysis_config$layout$method,
    seed = analysis_config$layout$seed,
    weights = NULL,
    iterations = analysis_config$layout$iterations
  )
}
local_plot <- plot_local_network(
  local_network,
  layout = local_layout,
  style = style,
  title = bquote(
    italic(.(local_config$center_node)) ~ "local network -" ~ .(local_config$condition)
  ),
  seed = analysis_config$layout$seed
)
local_stem <- paste0(
  "Fig3", local_config$panel, "_", local_config$center_node,
  "_local_network_", tolower(local_config$condition)
)
export_network_panel(
  local_plot,
  output_dir = output_dir,
  stem = local_stem,
  width_mm = analysis_config$output$local_network_width_mm,
  height_mm = analysis_config$output$local_network_height_mm,
  dpi = analysis_config$output$dpi,
  formats = unlist(analysis_config$output$formats, use.names = FALSE)
)
export_network_source_data(
  local_network,
  local_layout,
  output_dir = source_data_dir,
  stem = paste0("Fig3", local_config$panel, "_", local_config$center_node)
)
run_summary[[local_config$condition]]$local_network <- list(
  panel = local_config$panel,
  center_node = local_config$center_node,
  nodes = nrow(network_nodes(local_network)),
  edges = nrow(network_edges(local_network)),
  graph_signature = network_signature(local_network)
)

core_periphery <- list()
for (condition in names(prepared_degs)) {
  centrality_path <- file.path(
    project_root,
    "results", "current", "03_community",
    paste0("jaccard_centrality_measures_filtered_", tolower(condition), ".csv")
  )
  centrality <- read_hybs_csv(
    centrality_path,
    c("Cluster", "Nuclei", "Centrality_Score")
  )
  core_periphery[[condition]] <- prepare_core_periphery_data(
    centrality,
    prepared_degs[[condition]],
    condition = condition,
    cutoff = analysis_config$analysis$core_cutoff
  )
}
core_periphery <- do.call(rbind, core_periphery)
rownames(core_periphery) <- NULL
core_path <- file.path(source_data_dir, "Fig3F_core_periphery_data.csv")
if (file.exists(core_path)) {
  stop("Refusing to overwrite existing file: ", core_path, call. = FALSE)
}
utils::write.csv(core_periphery, core_path, row.names = FALSE)

manifest <- list(
  package_version = as.character(utils::packageVersion("HyBsNet")),
  project_root = project_root,
  output_dir = output_dir,
  analysis_config = analysis_config,
  style_config = normalizePath(style_path, mustWork = TRUE),
  networks = run_summary,
  fig3f_rows = nrow(core_periphery)
)
yaml::write_yaml(manifest, file.path(output_dir, "HyBsNet_run_manifest.yml"))
message("HyBsNet Figure 3 network workflow completed successfully: ", output_dir)
