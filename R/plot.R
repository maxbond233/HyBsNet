.plot_data <- function(network, layout) {
  validate_network_layout(layout, network)
  nodes <- network$nodes
  index <- match(nodes$name, layout$name)
  nodes$x <- layout$x[index]
  nodes$y <- layout$y[index]
  edges <- network$edges
  from <- match(edges$Cluster1, nodes$name)
  to <- match(edges$Cluster2, nodes$name)
  edges$x <- nodes$x[from]
  edges$y <- nodes$y[from]
  edges$xend <- nodes$x[to]
  edges$yend <- nodes$y[to]
  edges$connection <- ifelse(
    parse_region(edges$Cluster1) == parse_region(edges$Cluster2),
    "Within region",
    "Cross-region"
  )
  list(nodes = nodes, edges = edges)
}

.validate_observed_regions <- function(nodes, style) {
  missing <- setdiff(unique(as.character(nodes$region)), names(style$region_colors))
  if (length(missing) > 0L) {
    .hybs_abort(
      "No colour was supplied for observed region(s): ",
      paste(missing, collapse = ", "), "."
    )
  }
  invisible(nodes)
}

.network_plot_theme <- function(style) {
  ggplot2::theme_void(base_family = style$font_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = style$title_size, hjust = 0),
      legend.title = ggplot2::element_text(size = style$legend_title_size),
      legend.text = ggplot2::element_text(size = style$legend_text_size),
      legend.position = style$legend_position,
      plot.margin = ggplot2::margin(3, 3, 3, 3, unit = "mm")
    )
}

#' Plot a full DEG similarity network
#'
#' @param network A `hybs_network` object.
#' @param layout Cached or newly calculated layout. If `NULL`, an FR layout is
#'   calculated with `seed`.
#' @param labels Node names to label.
#' @param style A `hybs_network_style` object.
#' @param title Optional plot title.
#' @param seed Seed for layout and label placement.
#' @return A ggplot object.
#' @export
plot_similarity_network <- function(
    network,
    layout = NULL,
    labels = character(),
    style = theme_fig3_current(),
    title = NULL,
    seed = 42L) {
  .assert_hybs_network(network)
  .assert_style(style)
  layout <- layout %||% calculate_network_layout(network, seed = seed)
  data <- .plot_data(network, layout)
  .validate_observed_regions(data$nodes, style)
  data$nodes$label <- ifelse(data$nodes$name %in% labels, data$nodes$name, NA_character_)
  within <- data$edges[data$edges$connection == "Within region", , drop = FALSE]
  cross <- data$edges[data$edges$connection == "Cross-region", , drop = FALSE]

  plot <- ggplot2::ggplot()
  if (nrow(within) > 0L) {
    plot <- plot + ggplot2::geom_segment(
      data = within,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend, linewidth = Similarity),
      colour = unname(style$edge_colors[["Within region"]]),
      alpha = style$edge_alpha,
      lineend = "round"
    )
  }
  if (nrow(cross) > 0L) {
    plot <- plot + ggplot2::geom_segment(
      data = cross,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend, linewidth = Similarity),
      colour = unname(style$edge_colors[["Cross-region"]]),
      alpha = style$edge_alpha,
      lineend = "round"
    )
  }
  plot <- plot +
    ggplot2::geom_point(
      data = data$nodes,
      ggplot2::aes(x = x, y = y, size = deg_count, colour = region),
      alpha = style$node_alpha
    ) +
    ggrepel::geom_text_repel(
      data = data$nodes[!is.na(data$nodes$label), , drop = FALSE],
      ggplot2::aes(x = x, y = y, label = label),
      family = style$font_family,
      size = style$label_size,
      max.overlaps = Inf,
      seed = as.integer(seed),
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(values = style$region_colors, name = "Region") +
    ggplot2::scale_size_continuous(range = style$node_size_range, name = "DEG count") +
    ggplot2::scale_linewidth_continuous(
      range = style$edge_width_range,
      name = "Jaccard\nsimilarity"
    ) +
    ggplot2::guides(
      linewidth = ggplot2::guide_legend(
        override.aes = list(colour = "#404040", alpha = 1)
      )
    ) +
    ggplot2::labs(title = title) +
    ggplot2::coord_equal(clip = "off") +
    .network_plot_theme(style)
  plot
}

#' Plot a focal induced network
#'
#' @param network A local `hybs_network` returned by [extract_local_network()].
#' @param layout Cached or newly calculated layout.
#' @param style A `hybs_network_style` object.
#' @param title Optional title.
#' @param seed Seed for layout and labels.
#' @return A ggplot object.
#' @export
plot_local_network <- function(
    network,
    layout = NULL,
    style = theme_fig3_current(),
    title = NULL,
    seed = 42L) {
  .assert_hybs_network(network)
  .assert_style(style)
  center_node <- network$parameters$center_node
  if (is.null(center_node)) {
    .hybs_abort("Local network metadata do not contain `center_node`.")
  }
  layout <- layout %||% calculate_network_layout(network, seed = seed, weights = NULL)
  data <- .plot_data(network, layout)
  .validate_observed_regions(data$nodes, style)
  if (!"label" %in% names(data$nodes)) {
    data$nodes$label <- NA_character_
  }
  center <- data$nodes[data$nodes$name == center_node, , drop = FALSE]
  neighbors <- data$nodes[data$nodes$name != center_node, , drop = FALSE]
  edge_colour <- style$local_edge_color

  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = data$edges,
      ggplot2::aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        linewidth = Similarity,
        alpha = Similarity
      ),
      colour = edge_colour,
      lineend = "round"
    ) +
    ggplot2::geom_point(
      data = neighbors,
      ggplot2::aes(x = x, y = y, colour = region),
      size = style$node_size_range[[2L]] * 0.63,
      alpha = style$node_alpha
    ) +
    ggplot2::geom_point(
      data = center,
      ggplot2::aes(x = x, y = y),
      colour = style$center_color,
      size = style$node_size_range[[2L]] * 1.17,
      alpha = style$node_alpha
    ) +
    ggrepel::geom_text_repel(
      data = neighbors[!is.na(neighbors$label), , drop = FALSE],
      ggplot2::aes(x = x, y = y, label = label),
      family = style$font_family,
      size = style$label_size,
      max.overlaps = Inf,
      force = 3,
      box.padding = 0.35,
      point.padding = 0.15,
      segment.colour = "#777777",
      segment.alpha = 0.7,
      seed = as.integer(seed),
      show.legend = FALSE
    ) +
    ggrepel::geom_text_repel(
      data = center,
      ggplot2::aes(x = x, y = y, label = label),
      family = style$font_family,
      size = style$center_label_size,
      colour = style$center_color,
      fontface = "bold",
      max.overlaps = Inf,
      segment.colour = style$center_color,
      segment.alpha = 0.7,
      seed = as.integer(seed),
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(values = style$region_colors, name = "Region") +
    ggplot2::scale_linewidth_continuous(range = style$edge_width_range, guide = "none") +
    ggplot2::scale_alpha_continuous(range = c(0.18, 0.78), guide = "none") +
    ggplot2::labs(title = title) +
    ggplot2::coord_equal(clip = "off") +
    .network_plot_theme(style)
}
