#' Export an editable and reviewable network panel
#'
#' @param plot A ggplot object.
#' @param output_dir Output directory.
#' @param stem File stem without extension.
#' @param width_mm,height_mm Figure size in millimetres.
#' @param dpi PNG resolution.
#' @param formats Any of `"pdf"`, `"png"`, and `"svg"`.
#' @param overwrite Whether existing files may be replaced.
#' @return A named character vector of output paths, invisibly.
#' @export
export_network_panel <- function(
    plot,
    output_dir,
    stem,
    width_mm,
    height_mm,
    dpi = 300,
    formats = c("pdf", "png"),
    overwrite = FALSE) {
  if (!inherits(plot, "ggplot")) {
    .hybs_abort("`plot` must be a ggplot object.")
  }
  .assert_nonempty_string(output_dir, "output_dir")
  .assert_nonempty_string(stem, "stem")
  .assert_scalar_number(width_mm, "width_mm", lower = 1)
  .assert_scalar_number(height_mm, "height_mm", lower = 1)
  .assert_scalar_number(dpi, "dpi", lower = 72)
  formats <- unique(tolower(as.character(formats)))
  invalid <- setdiff(formats, c("pdf", "png", "svg"))
  if (length(invalid) > 0L) {
    .hybs_abort("Unsupported output format(s): ", paste(invalid, collapse = ", "))
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- stats::setNames(file.path(output_dir, paste0(stem, ".", formats)), formats)
  existing <- paths[file.exists(paths)]
  if (length(existing) > 0L && !overwrite) {
    .hybs_abort(
      "Refusing to overwrite existing output(s): ",
      paste(existing, collapse = ", ")
    )
  }
  for (format in formats) {
    device <- switch(
      format,
      pdf = grDevices::cairo_pdf,
      png = "png",
      svg = grDevices::svg
    )
    ggplot2::ggsave(
      filename = paths[[format]],
      plot = plot,
      width = width_mm,
      height = height_mm,
      units = "mm",
      dpi = dpi,
      device = device,
      bg = "white"
    )
  }
  invisible(paths)
}

#' Export node, edge, layout, and parameter source data
#'
#' @param network A `hybs_network` object.
#' @param layout Validated layout data.
#' @param output_dir Output directory.
#' @param stem File stem.
#' @param overwrite Whether existing files may be replaced.
#' @return Paths to the three source-data files, invisibly.
#' @export
export_network_source_data <- function(
    network,
    layout,
    output_dir,
    stem,
    overwrite = FALSE) {
  .assert_hybs_network(network)
  validate_network_layout(layout, network)
  .assert_nonempty_string(output_dir, "output_dir")
  .assert_nonempty_string(stem, "stem")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  nodes <- network$nodes
  index <- match(nodes$name, layout$name)
  nodes$x <- layout$x[index]
  nodes$y <- layout$y[index]
  nodes$graph_signature <- network_signature(network)
  paths <- c(
    nodes = file.path(output_dir, paste0(stem, "_nodes.csv")),
    edges = file.path(output_dir, paste0(stem, "_edges.csv")),
    parameters = file.path(output_dir, paste0(stem, "_parameters.yml"))
  )
  existing <- paths[file.exists(paths)]
  if (length(existing) > 0L && !overwrite) {
    .hybs_abort(
      "Refusing to overwrite existing source-data output(s): ",
      paste(existing, collapse = ", ")
    )
  }
  .atomic_write_csv(nodes, paths[["nodes"]], overwrite = overwrite)
  .atomic_write_csv(network$edges, paths[["edges"]], overwrite = overwrite)
  metadata <- utils::modifyList(
    network$parameters,
    list(
      graph_signature = network_signature(network),
      n_nodes = nrow(network$nodes),
      n_edges = nrow(network$edges),
      layout_method = unique(as.character(layout$layout_method %||% NA_character_)),
      layout_seed = unique(as.integer(layout$layout_seed %||% NA_integer_))
    )
  )
  yaml::write_yaml(metadata, paths[["parameters"]])
  invisible(paths)
}
