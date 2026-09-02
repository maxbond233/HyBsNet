# Overlay and subnetwork plotting ------------------------------------------

.overlay_is_diverging <- function(score_type) {
  any(grepl(
    "log2fc|nes|signed|direction|contrast|difference",
    tolower(as.character(score_type))
  ))
}

.overlay_limits <- function(score, diverging, limits = NULL) {
  if (!is.null(limits)) {
    if (length(limits) != 2L || any(!is.finite(limits)) || limits[[1L]] >= limits[[2L]]) {
      .hybs_abort("`limits` must contain two increasing finite numbers.")
    }
    return(as.numeric(limits))
  }
  finite <- score[is.finite(score)]
  if (length(finite) == 0L) return(if (diverging) c(-1, 1) else c(0, 1))
  if (diverging) {
    bound <- max(abs(finite))
    if (!is.finite(bound) || bound == 0) bound <- 1
    return(c(-bound, bound))
  }
  observed <- range(finite)
  if (diff(observed) == 0) {
    padding <- max(abs(observed[[1L]]) * 0.05, 0.5)
    observed <- observed + c(-padding, padding)
  }
  observed
}

.add_overlay_fill_scale <- function(
    plot,
    nodes,
    score_label,
    diverging,
    limits,
    missing_color,
    low,
    mid,
    high) {
  if (diverging) {
    plot + ggplot2::scale_fill_gradient2(
      low = low,
      mid = mid,
      high = high,
      midpoint = 0,
      limits = limits,
      oob = scales::squish,
      na.value = missing_color,
      name = score_label
    )
  } else {
    plot + ggplot2::scale_fill_gradient(
      low = mid,
      high = high,
      limits = limits,
      oob = scales::squish,
      na.value = missing_color,
      name = score_label
    )
  }
}

#' Plot a fixed-layout network with a node overlay
#'
#' Overlay values alter node aesthetics only. They never change topology,
#' communities, or coordinates.
#'
#' @param network A `hybs_network` object.
#' @param overlay A `hybs_node_overlay` object.
#' @param layout Cached or newly calculated layout.
#' @param condition Optional overlay condition.
#' @param strict_overlay_nodes If `TRUE`, overlay records outside the plotted
#'   network are rejected. The default records and ignores them because a full
#'   DEG or pathway table commonly includes nodes with no retained network edge.
#' @param labels Node names to label.
#' @param style A `hybs_network_style` object.
#' @param title Optional plot title.
#' @param size_by One of `"overlay"`, `"deg_count"`, or `"constant"`.
#' @param diverging Whether to use a zero-centred fill scale. When `NULL`, this
#'   is inferred from `score_type`.
#' @param limits Optional two-number fill limits.
#' @param missing_color Fill used for missing or untested nodes.
#' @param low,mid,high Continuous fill colours.
#' @param significant_outline,nonsignificant_outline Outline colours.
#' @param seed Seed for layout and label placement.
#' @return A ggplot object.
#' @export
plot_network_overlay <- function(
    network,
    overlay,
    layout = NULL,
    condition = NULL,
    strict_overlay_nodes = FALSE,
    labels = character(),
    style = theme_fig3_current(),
    title = NULL,
    size_by = c("overlay", "deg_count", "constant"),
    diverging = NULL,
    limits = NULL,
    missing_color = "#D9DDE0",
    low = "#2F6F9F",
    mid = "#F7F7F7",
    high = "#C94B40",
    significant_outline = "#1F2933",
    nonsignificant_outline = "#FFFFFF",
    seed = 42L) {
  .assert_hybs_network(network)
  .assert_style(style)
  .assert_flag(strict_overlay_nodes, "strict_overlay_nodes")
  size_by <- match.arg(size_by)
  attached <- attach_node_overlay(
    network,
    overlay,
    condition = condition,
    strict = strict_overlay_nodes
  )
  layout <- layout %||% calculate_network_layout(attached, seed = seed)
  data <- .plot_data(attached, layout)
  data$nodes$label <- ifelse(
    data$nodes$name %in% labels,
    data$nodes$name,
    NA_character_
  )
  data$nodes$overlay_score[!data$nodes$overlay_available] <- NA_real_
  if (is.null(diverging)) {
    diverging <- .overlay_is_diverging(data$nodes$overlay_score_type)
  }
  .assert_flag(diverging, "diverging")
  fill_limits <- .overlay_limits(data$nodes$overlay_score, diverging, limits)
  if (size_by == "overlay") {
    finite_size <- is.finite(data$nodes$overlay_size)
    if (any(data$nodes$overlay_size[finite_size] < 0)) {
      .hybs_abort("Overlay size values must be non-negative.")
    }
    data$nodes$plot_size <- ifelse(finite_size, data$nodes$overlay_size, 0)
    size_label <- attached$parameters$overlay$size_label %||% "Overlay size"
  } else if (size_by == "deg_count") {
    data$nodes$plot_size <- data$nodes$deg_count
    size_label <- "DEG count"
  } else {
    data$nodes$plot_size <- 1
    size_label <- NULL
  }
  data$nodes$overlay_outline <- ifelse(
    data$nodes$overlay_significant,
    significant_outline,
    nonsignificant_outline
  )
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
      ggplot2::aes(
        x = x,
        y = y,
        size = plot_size,
        fill = overlay_score,
        colour = overlay_outline
      ),
      shape = 21,
      stroke = 0.55,
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
    ggplot2::scale_colour_identity() +
    ggplot2::scale_linewidth_continuous(
      range = style$edge_width_range,
      name = "Jaccard\nsimilarity"
    ) +
    ggplot2::labs(title = title) +
    ggplot2::coord_equal(clip = "off") +
    .network_plot_theme(style)
  plot <- .add_overlay_fill_scale(
    plot,
    data$nodes,
    attached$parameters$overlay$score_label %||% "Score",
    diverging,
    fill_limits,
    missing_color,
    low,
    mid,
    high
  )
  plot + ggplot2::scale_size_continuous(
    range = style$node_size_range,
    name = size_label,
    guide = if (size_by == "constant") "none" else "legend"
  )
}

#' Plot an extracted subnetwork
#'
#' @param subnetwork A `hybs_subnetwork` object.
#' @param layout A validated subnetwork layout.
#' @param overlay Optional node overlay.
#' @param labels Optional node labels. Sensible seed-focused defaults are used
#'   when omitted.
#' @param style,title,seed Plot settings.
#' @param ... Additional arguments passed to [plot_network_overlay()] when an
#'   overlay is supplied.
#' @return A ggplot object.
#' @export
plot_subnetwork <- function(
    subnetwork,
    layout,
    overlay = NULL,
    labels = NULL,
    style = theme_fig3_current(),
    title = NULL,
    seed = 42L,
    ...) {
  if (!inherits(subnetwork, "hybs_subnetwork")) {
    .hybs_abort("`subnetwork` must be returned by `extract_subnetwork()`.")
  }
  validate_network_layout(layout, subnetwork)
  if (is.null(labels)) {
    labels <- if (nrow(subnetwork$nodes) <= 20L) {
      as.character(subnetwork$nodes$name)
    } else {
      as.character(subnetwork$nodes$name[subnetwork$nodes$is_seed])
    }
  }
  plot <- if (is.null(overlay)) {
    plot_similarity_network(
      subnetwork,
      layout = layout,
      labels = labels,
      style = style,
      title = title,
      seed = seed
    )
  } else {
    plot_network_overlay(
      subnetwork,
      overlay,
      layout = layout,
      labels = labels,
      style = style,
      title = title,
      seed = seed,
      ...
    )
  }
  node_data <- .plot_data(subnetwork, layout)$nodes
  seeds <- node_data[node_data$is_seed %in% TRUE, , drop = FALSE]
  if (nrow(seeds) > 0L) {
    plot <- plot + ggplot2::geom_point(
      data = seeds,
      ggplot2::aes(x = x, y = y),
      shape = 21,
      fill = NA,
      colour = style$center_color,
      size = style$node_size_range[[2L]] + 0.7,
      stroke = 0.7,
      inherit.aes = FALSE,
      show.legend = FALSE
    )
  }
  plot
}

#' Plot a subnetwork in its full-network context
#'
#' @param network Parent `hybs_network` object.
#' @param layout Parent frozen layout.
#' @param subnetwork Extracted `hybs_subnetwork`.
#' @param overlay Optional overlay used only for highlighted subnetwork nodes.
#' @param style,title,seed Plot settings.
#' @return A ggplot object.
#' @export
plot_subnetwork_context <- function(
    network,
    layout,
    subnetwork,
    overlay = NULL,
    style = theme_fig3_current(),
    title = NULL,
    seed = 42L) {
  .assert_hybs_network(network)
  if (!inherits(subnetwork, "hybs_subnetwork")) {
    .hybs_abort("`subnetwork` must be returned by `extract_subnetwork()`.")
  }
  validate_network_layout(layout, network)
  parent_signature <- subnetwork$parameters$parent_graph_signature
  if (!identical(parent_signature, network_signature(network))) {
    .hybs_abort("Subnetwork was not extracted from the supplied parent network.")
  }
  data <- .plot_data(network, layout)
  selected <- data$nodes[
    data$nodes$name %in% subnetwork$nodes$name,
    , drop = FALSE
  ]
  seed_names <- subnetwork$nodes$name[subnetwork$nodes$is_seed %in% TRUE]
  selected$is_seed <- selected$name %in% seed_names
  plot <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = data$edges,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      colour = "#C5CBD3",
      alpha = 0.32,
      linewidth = 0.25,
      lineend = "round"
    ) +
    ggplot2::geom_point(
      data = data$nodes,
      ggplot2::aes(x = x, y = y),
      colour = "#D7DCE2",
      size = 1.8
    )
  if (is.null(overlay)) {
    .validate_observed_regions(selected, style)
    plot <- plot +
      ggplot2::geom_point(
        data = selected,
        ggplot2::aes(x = x, y = y, fill = region),
        shape = 21,
        colour = "white",
        stroke = 0.45,
        size = 3.4
      ) +
      ggplot2::scale_fill_manual(values = style$region_colors, name = "Region")
  } else {
    attached <- attach_node_overlay(network, overlay, strict = FALSE)
    attached_data <- .plot_data(attached, layout)$nodes
    selected <- attached_data[
      attached_data$name %in% subnetwork$nodes$name,
      , drop = FALSE
    ]
    selected$overlay_score[!selected$overlay_available] <- NA_real_
    diverging <- .overlay_is_diverging(selected$overlay_score_type)
    fill_limits <- .overlay_limits(selected$overlay_score, diverging)
    plot <- plot + ggplot2::geom_point(
      data = selected,
      ggplot2::aes(x = x, y = y, fill = overlay_score),
      shape = 21,
      colour = "white",
      stroke = 0.45,
      size = 3.4
    )
    plot <- .add_overlay_fill_scale(
      plot,
      selected,
      attached$parameters$overlay$score_label %||% "Score",
      diverging,
      fill_limits,
      "#D9DDE0",
      "#2F6F9F",
      "#F7F7F7",
      "#C94B40"
    )
  }
  seeds <- selected[selected$name %in% seed_names, , drop = FALSE]
  plot +
    ggplot2::geom_point(
      data = seeds,
      ggplot2::aes(x = x, y = y),
      shape = 21,
      fill = NA,
      colour = style$center_color,
      size = 4.4,
      stroke = 0.8,
      inherit.aes = FALSE
    ) +
    ggrepel::geom_text_repel(
      data = seeds,
      ggplot2::aes(x = x, y = y, label = name),
      family = style$font_family,
      size = style$label_size,
      seed = as.integer(seed),
      show.legend = FALSE
    ) +
    ggplot2::labs(title = title) +
    ggplot2::coord_equal(clip = "off") +
    .network_plot_theme(style)
}
