#' Calculate a deterministic network layout
#'
#' @param network A `hybs_network` object.
#' @param method One of `"fr"`, `"kk"`, `"circle"`, or `"nicely"`.
#' @param seed Random seed used without changing the caller's RNG state.
#' @param weights Edge attribute used as layout weights, or `NULL`.
#' @param iterations Iterations for the Fruchterman-Reingold layout.
#' @return A data frame with node coordinates and a topology signature.
#' @export
calculate_network_layout <- function(
    network,
    method = c("fr", "kk", "circle", "nicely"),
    seed = 42L,
    weights = "Similarity",
    iterations = 500L) {
  .assert_hybs_network(network)
  method <- match.arg(method)
  .assert_scalar_number(iterations, "iterations", lower = 1)
  graph <- network$graph
  weight_values <- NULL
  if (!is.null(weights)) {
    .assert_nonempty_string(weights, "weights")
    weight_values <- igraph::edge_attr(graph, weights)
    if (is.null(weight_values) || !is.numeric(weight_values) || any(!is.finite(weight_values))) {
      .hybs_abort("Requested layout weights are absent or non-finite: ", weights)
    }
  }
  coords <- .with_seed(seed, {
    switch(
      method,
      fr = igraph::layout_with_fr(
        graph,
        weights = weight_values,
        niter = as.integer(iterations)
      ),
      kk = igraph::layout_with_kk(graph, weights = weight_values),
      circle = igraph::layout_in_circle(graph),
      nicely = igraph::layout_nicely(graph)
    )
  })
  out <- data.frame(
    name = igraph::V(graph)$name,
    x = coords[, 1L],
    y = coords[, 2L],
    stringsAsFactors = FALSE
  )
  node_index <- match(out$name, network$nodes$name)
  extra <- setdiff(names(network$nodes), "name")
  for (column in extra) {
    out[[column]] <- network$nodes[[column]][node_index]
  }
  out$graph_signature <- network_signature(network)
  out$layout_method <- method
  out$layout_seed <- as.integer(seed)
  class(out) <- c("hybs_network_layout", "data.frame")
  out
}

#' Validate layout coordinates against network topology
#'
#' @param layout Layout data frame.
#' @param network A `hybs_network` object.
#' @return The layout invisibly.
#' @export
validate_network_layout <- function(layout, network) {
  .assert_hybs_network(network)
  .assert_columns(layout, c("name", "x", "y"), "network layout")
  if (anyDuplicated(layout$name)) {
    .hybs_abort("Network layout contains duplicate node names.")
  }
  if (any(!is.finite(layout$x)) || any(!is.finite(layout$y))) {
    .hybs_abort("Network layout coordinates must be finite.")
  }
  expected <- sort(igraph::V(network$graph)$name)
  observed <- sort(as.character(layout$name))
  if (!identical(expected, observed)) {
    .hybs_abort("Network layout nodes do not match network nodes.")
  }
  if ("graph_signature" %in% names(layout)) {
    signatures <- unique(as.character(layout$graph_signature))
    if (length(signatures) != 1L || signatures != network_signature(network)) {
      .hybs_abort("Cached layout topology signature does not match the network.")
    }
  }
  invisible(layout)
}

#' Save reusable network coordinates
#'
#' @param layout Layout returned by [calculate_network_layout()].
#' @param path Output CSV path.
#' @param overwrite Whether to replace an existing file.
#' @return `path`, invisibly.
#' @export
save_network_layout <- function(layout, path, overwrite = FALSE) {
  .assert_columns(layout, c("name", "x", "y", "graph_signature"), "network layout")
  .atomic_write_csv(as.data.frame(layout), path, overwrite = overwrite)
}

#' Load and validate cached network coordinates
#'
#' @param path Layout CSV path.
#' @param network A `hybs_network` object.
#' @return A `hybs_network_layout` data frame.
#' @export
load_network_layout <- function(path, network) {
  layout <- read_hybs_csv(path, c("name", "x", "y", "graph_signature"))
  validate_network_layout(layout, network)
  class(layout) <- c("hybs_network_layout", "data.frame")
  layout
}
