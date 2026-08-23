#' Create a reusable network figure style
#'
#' @param region_colors Named colours for BS, SH, TH, or other regions.
#' @param edge_colors Named colours for `Within region` and `Cross-region`.
#' @param local_edge_color Edge colour for focal induced networks.
#' @param center_color Colour for a focal local-network node.
#' @param font_family Font family.
#' @param node_size_range Numeric minimum and maximum node sizes.
#' @param edge_width_range Numeric minimum and maximum edge widths.
#' @param node_alpha,edge_alpha Point and edge opacity.
#' @param label_size,center_label_size Label sizes in mm.
#' @param title_size,legend_title_size,legend_text_size Text sizes in points.
#' @param legend_position A ggplot2 legend position.
#' @return An object of class `hybs_network_style`.
#' @export
hybs_network_style <- function(
    region_colors = c(BS = "#316B9D", SH = "#DC7E39", TH = "#4A9267"),
    edge_colors = c(
      `Within region` = "#B7B7B7",
      `Cross-region` = "#F39C12"
    ),
    local_edge_color = "#A8A8A8",
    center_color = "#C51B29",
    font_family = "Arial",
    node_size_range = c(1.6, 6.0),
    edge_width_range = c(0.2, 1.25),
    node_alpha = 0.98,
    edge_alpha = 0.82,
    label_size = 2.2,
    center_label_size = 2.6,
    title_size = 10,
    legend_title_size = 7.5,
    legend_text_size = 7,
    legend_position = "right") {
  if (is.null(names(region_colors)) || any(!nzchar(names(region_colors)))) {
    .hybs_abort("`region_colors` must be a named character vector.")
  }
  needed_edges <- c("Within region", "Cross-region")
  if (is.null(names(edge_colors)) || !all(needed_edges %in% names(edge_colors))) {
    .hybs_abort(
      "`edge_colors` must define: ", paste(needed_edges, collapse = ", "), "."
    )
  }
  for (item in list(
    node_size_range = node_size_range,
    edge_width_range = edge_width_range
  )) {
    if (length(item) != 2L || any(!is.finite(item)) || any(item <= 0)) {
      .hybs_abort("Size and width ranges must contain two positive finite values.")
    }
  }
  .assert_scalar_number(node_alpha, "node_alpha", lower = 0, upper = 1)
  .assert_scalar_number(edge_alpha, "edge_alpha", lower = 0, upper = 1)
  out <- list(
    region_colors = region_colors,
    edge_colors = edge_colors,
    local_edge_color = local_edge_color,
    center_color = center_color,
    font_family = font_family,
    node_size_range = node_size_range,
    edge_width_range = edge_width_range,
    node_alpha = node_alpha,
    edge_alpha = edge_alpha,
    label_size = label_size,
    center_label_size = center_label_size,
    title_size = title_size,
    legend_title_size = legend_title_size,
    legend_text_size = legend_text_size,
    legend_position = legend_position
  )
  class(out) <- "hybs_network_style"
  out
}

#' Current Figure 3 network style
#'
#' @return A `hybs_network_style` object reproducing the current palette and
#'   main sizing defaults.
#' @export
theme_fig3_current <- function() {
  hybs_network_style()
}

#' Colour-blind-friendly Figure 3 network style
#'
#' @return A `hybs_network_style` object.
#' @export
theme_fig3_colorblind <- function() {
  hybs_network_style(
    region_colors = c(BS = "#0072B2", SH = "#D55E00", TH = "#009E73"),
    edge_colors = c(
      `Within region` = "#B3B3B3",
      `Cross-region` = "#E69F00"
    ),
    local_edge_color = "#A6A6A6",
    center_color = "#CC79A7"
  )
}

.assert_style <- function(style) {
  if (!inherits(style, "hybs_network_style")) {
    .hybs_abort("`style` must be created by `hybs_network_style()`.")
  }
  invisible(style)
}

#' Read an editable YAML configuration
#'
#' @param path YAML path.
#' @return A named list.
#' @export
read_network_config <- function(path) {
  .assert_nonempty_string(path, "path")
  if (!file.exists(path)) {
    .hybs_abort("Missing configuration file: ", path)
  }
  config <- yaml::read_yaml(path)
  if (!is.list(config)) {
    .hybs_abort("Network configuration must decode to a named list.")
  }
  config
}

#' Convert a style configuration list into a figure style
#'
#' @param config A list, typically read from `figure3_style.yml`.
#' @return A `hybs_network_style` object.
#' @export
style_from_config <- function(config) {
  if (!is.list(config)) {
    .hybs_abort("`config` must be a list.")
  }
  args <- config$style %||% config
  if (!is.list(args)) {
    .hybs_abort("Style configuration must be a named list.")
  }
  vector_fields <- c("region_colors", "edge_colors", "node_size_range", "edge_width_range")
  for (field in intersect(vector_fields, names(args))) {
    args[[field]] <- unlist(args[[field]], use.names = TRUE)
    if (field %in% c("node_size_range", "edge_width_range")) {
      args[[field]] <- as.numeric(args[[field]])
    }
  }
  do.call(hybs_network_style, args)
}
