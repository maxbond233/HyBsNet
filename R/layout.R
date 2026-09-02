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

#' Standardize external coordinates as a network layout
#'
#' This is the preferred entry point for manually adjusted or publication
#' layouts. The coordinates are copied into a small, topology-aware object;
#' node annotations in the input table are intentionally not required.
#'
#' @param data A data frame containing node identifiers and coordinates.
#' @param network A `hybs_network` object.
#' @param node_col,x_col,y_col Columns containing node names and coordinates.
#' @param layout_id Stable identifier for this frozen layout version.
#' @param layout_method Description of how the coordinates were produced.
#' @param layout_seed Optional layout seed. Use `NA_integer_` when not known.
#' @param strict_signature If `TRUE`, an existing `graph_signature` column must
#'   match `network`.
#' @return A `hybs_network_layout` data frame.
#' @export
as_network_layout <- function(
    data,
    network,
    node_col = "name",
    x_col = "x",
    y_col = "y",
    layout_id = "external",
    layout_method = "frozen_external",
    layout_seed = NA_integer_,
    strict_signature = TRUE) {
  .assert_hybs_network(network)
  .assert_flag(strict_signature, "strict_signature")
  .assert_nonempty_string(node_col, "node_col")
  .assert_nonempty_string(x_col, "x_col")
  .assert_nonempty_string(y_col, "y_col")
  .assert_nonempty_string(layout_id, "layout_id")
  .assert_nonempty_string(layout_method, "layout_method")
  .assert_columns(data, c(node_col, x_col, y_col), "external layout")

  if (strict_signature && "graph_signature" %in% names(data)) {
    signatures <- unique(as.character(data$graph_signature))
    signatures <- signatures[!is.na(signatures) & nzchar(signatures)]
    if (length(signatures) > 0L &&
        (length(signatures) != 1L || signatures != network_signature(network))) {
      .hybs_abort("External layout topology signature does not match the network.")
    }
  }
  out <- data.frame(
    name = as.character(data[[node_col]]),
    x = suppressWarnings(as.numeric(data[[x_col]])),
    y = suppressWarnings(as.numeric(data[[y_col]])),
    stringsAsFactors = FALSE
  )
  out$graph_signature <- network_signature(network)
  out$layout_id <- layout_id
  out$layout_method <- layout_method
  out$layout_seed <- as.integer(layout_seed)
  validate_network_layout(out, network)
  class(out) <- c("hybs_network_layout", "data.frame")
  out
}

#' Freeze and version network coordinates
#'
#' @param layout A layout data frame.
#' @param network A `hybs_network` object.
#' @param layout_id Stable identifier for this layout version.
#' @param layout_method Optional method label. The value already stored in
#'   `layout` is used when omitted.
#' @param layout_seed Optional seed. The value already stored in `layout` is
#'   used when omitted.
#' @return A minimal, validated `hybs_network_layout` data frame.
#' @export
freeze_network_layout <- function(
    layout,
    network,
    layout_id,
    layout_method = NULL,
    layout_seed = NULL) {
  .assert_nonempty_string(layout_id, "layout_id")
  method <- layout_method %||%
    if ("layout_method" %in% names(layout)) {
      unique(as.character(layout$layout_method))[[1L]]
    } else {
      "frozen"
    }
  seed <- layout_seed %||%
    if ("layout_seed" %in% names(layout)) {
      unique(as.integer(layout$layout_seed))[[1L]]
    } else {
      NA_integer_
    }
  as_network_layout(
    layout,
    network,
    layout_id = layout_id,
    layout_method = method,
    layout_seed = seed,
    strict_signature = TRUE
  )
}

#' Compare coordinates shared by two layouts
#'
#' @param reference,candidate Layout data frames.
#' @param tolerance Maximum absolute coordinate difference considered equal.
#' @return A data frame with coordinate deltas for shared nodes.
#' @export
compare_network_layouts <- function(reference, candidate, tolerance = 0) {
  .assert_columns(reference, c("name", "x", "y"), "reference layout")
  .assert_columns(candidate, c("name", "x", "y"), "candidate layout")
  .assert_scalar_number(tolerance, "tolerance", lower = 0)
  if (anyDuplicated(reference$name) || anyDuplicated(candidate$name)) {
    .hybs_abort("Layouts must contain unique node names.")
  }
  shared <- sort(intersect(as.character(reference$name), as.character(candidate$name)))
  ref_index <- match(shared, reference$name)
  candidate_index <- match(shared, candidate$name)
  out <- data.frame(
    name = shared,
    reference_x = as.numeric(reference$x[ref_index]),
    reference_y = as.numeric(reference$y[ref_index]),
    candidate_x = as.numeric(candidate$x[candidate_index]),
    candidate_y = as.numeric(candidate$y[candidate_index]),
    stringsAsFactors = FALSE
  )
  out$delta_x <- out$candidate_x - out$reference_x
  out$delta_y <- out$candidate_y - out$reference_y
  out$coordinates_identical <-
    abs(out$delta_x) <= tolerance & abs(out$delta_y) <= tolerance
  attr(out, "reference_only") <- setdiff(as.character(reference$name), shared)
  attr(out, "candidate_only") <- setdiff(as.character(candidate$name), shared)
  attr(out, "all_shared_coordinates_identical") <-
    nrow(out) > 0L && all(out$coordinates_identical)
  out
}

#' Reconcile a frozen layout with a changed node set
#'
#' Shared nodes retain their exact coordinates. New nodes can either trigger an
#' error or be placed near the similarity-weighted centre of positioned
#' neighbours. Removed nodes are dropped and recorded as an attribute.
#'
#' @param reference_layout Frozen coordinates for a related network.
#' @param network Target `hybs_network` object.
#' @param new_nodes Either `"error"` or `"anchor_neighbors"`.
#' @param anchor_offset_fraction Deterministic offset relative to the coordinate
#' span, used to keep multiple new nodes from overlapping exactly.
#' @param layout_id Identifier for the reconciled layout.
#' @return A validated `hybs_network_layout` data frame.
#' @export
reconcile_network_layout <- function(
    reference_layout,
    network,
    new_nodes = c("error", "anchor_neighbors"),
    anchor_offset_fraction = 0.035,
    layout_id = "reconciled") {
  .assert_hybs_network(network)
  new_nodes <- match.arg(new_nodes)
  .assert_scalar_number(
    anchor_offset_fraction, "anchor_offset_fraction", lower = 0
  )
  .assert_nonempty_string(layout_id, "layout_id")
  .assert_columns(reference_layout, c("name", "x", "y"), "reference layout")
  if (anyDuplicated(reference_layout$name)) {
    .hybs_abort("Reference layout contains duplicate node names.")
  }
  if (any(!is.finite(reference_layout$x)) || any(!is.finite(reference_layout$y))) {
    .hybs_abort("Reference layout coordinates must be finite.")
  }

  wanted <- as.character(network$nodes$name)
  reference_names <- as.character(reference_layout$name)
  shared <- intersect(wanted, reference_names)
  added <- sort(setdiff(wanted, reference_names))
  removed <- sort(setdiff(reference_names, wanted))
  if (length(shared) == 0L) {
    .hybs_abort("Reference layout and target network have no shared nodes.")
  }
  if (length(added) > 0L && new_nodes == "error") {
    .hybs_abort(
      "Target network contains node(s) absent from the frozen layout: ",
      paste(utils::head(added, 5L), collapse = ", "), "."
    )
  }

  coords <- reference_layout[
    match(shared, reference_names), c("name", "x", "y"), drop = FALSE
  ]
  if (length(added) > 0L) {
    edges <- network_edges(network)
    span <- max(diff(range(coords$x)), diff(range(coords$y)))
    if (!is.finite(span) || span == 0) span <- 1
    offset <- anchor_offset_fraction * span
    remaining <- added
    for (pass in seq_len(max(1L, length(added)))) {
      placed <- character()
      for (node in remaining) {
        incident <- edges$Cluster1 == node | edges$Cluster2 == node
        node_edges <- edges[incident, , drop = FALSE]
        neighbours <- ifelse(
          node_edges$Cluster1 == node,
          node_edges$Cluster2,
          node_edges$Cluster1
        )
        positioned <- match(neighbours, coords$name)
        usable <- !is.na(positioned)
        if (!any(usable)) next
        weights <- suppressWarnings(as.numeric(node_edges$Similarity[usable]))
        weights[!is.finite(weights) | weights <= 0] <- 1
        neighbour_coords <- coords[positioned[usable], , drop = FALSE]
        angle_index <- match(node, added) - 1L
        angle <- 2 * pi * angle_index / max(1L, length(added))
        coords <- rbind(
          coords,
          data.frame(
            name = node,
            x = stats::weighted.mean(neighbour_coords$x, weights) +
              offset * cos(angle),
            y = stats::weighted.mean(neighbour_coords$y, weights) +
              offset * sin(angle),
            stringsAsFactors = FALSE
          )
        )
        placed <- c(placed, node)
      }
      remaining <- setdiff(remaining, placed)
      if (length(placed) == 0L) break
    }
    if (length(remaining) > 0L) {
      .hybs_abort(
        "Could not anchor new node(s) to positioned neighbours: ",
        paste(remaining, collapse = ", "), "."
      )
    }
  }
  coords <- coords[match(wanted, coords$name), , drop = FALSE]
  coords$graph_signature <- network_signature(network)
  coords$layout_id <- layout_id
  coords$layout_method <- if (length(added) > 0L) {
    "frozen_shared_plus_anchored_new"
  } else {
    "frozen_shared_subset"
  }
  coords$layout_seed <- NA_integer_
  validate_network_layout(coords, network)
  class(coords) <- c("hybs_network_layout", "data.frame")
  attr(coords, "shared_nodes") <- sort(shared)
  attr(coords, "added_nodes") <- added
  attr(coords, "removed_nodes") <- removed
  coords
}

#' Extract the node coordinates from a layout
#'
#' @param layout A network layout.
#' @return A data frame with `name`, `x`, and `y`.
#' @export
layout_nodes <- function(layout) {
  .assert_columns(layout, c("name", "x", "y"), "network layout")
  as.data.frame(layout[, c("name", "x", "y"), drop = FALSE])
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
