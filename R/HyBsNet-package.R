#' HyBsNet: reusable DEG similarity networks
#'
#' `HyBsNet` separates four concerns that were previously interleaved in panel
#' scripts: scientific filtering, graph analysis, deterministic coordinates,
#' and visual styling. The current Figure 3 defaults remain adjusted `P < 0.05`,
#' Jaccard `> 0.10`, and random seed `42`.
#'
#' @keywords internal
"_PACKAGE"

utils::globalVariables(c(
  "x", "y", "xend", "yend", "Similarity", "deg_count", "region", "label",
  "plot_size", "overlay_score", "overlay_outline", "name"
))
