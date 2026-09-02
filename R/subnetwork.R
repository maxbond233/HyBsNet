# Subnetwork extraction -----------------------------------------------------

.hybs_cell_class <- function(node) {
  parts <- strsplit(as.character(node), "_", fixed = TRUE)
  vapply(parts, function(x) {
    if (length(x) < 2L) NA_character_ else x[[2L]]
  }, character(1))
}

.select_subnetwork_seeds <- function(
    network,
    mode,
    center_node,
    nodes,
    pattern,
    community,
    community_col,
    group,
    group_col,
    neuronal_classes) {
  node_data <- network$nodes
  available <- as.character(node_data$name)
  selected <- switch(
    mode,
    ego = {
      .assert_nonempty_string(center_node, "center_node")
      center_node
    },
    nodes = {
      if (length(nodes) == 0L) {
        .hybs_abort("`nodes` must contain at least one node in nodes mode.")
      }
      unique(as.character(nodes))
    },
    pattern = {
      .assert_nonempty_string(pattern, "pattern")
      available[grepl(pattern, available)]
    },
    community = {
      .assert_nonempty_string(community_col, "community_col")
      .assert_columns(node_data, community_col, "network nodes")
      available[as.character(node_data[[community_col]]) %in% as.character(community)]
    },
    group = {
      .assert_nonempty_string(group, "group")
      if (!is.null(group_col)) {
        .assert_nonempty_string(group_col, "group_col")
        .assert_columns(node_data, group_col, "network nodes")
        available[as.character(node_data[[group_col]]) == group]
      } else {
        alias <- tolower(trimws(group))
        cell_class <- .hybs_cell_class(available)
        neuronal <- cell_class %in% neuronal_classes
        if (alias %in% c("non-neuron", "nonneuron", "non neuron")) {
          available[!neuronal]
        } else if (alias == "neuron") {
          available[neuronal]
        } else if (grepl("^(bs|sh|th)\\s+neuron$", alias)) {
          region <- toupper(sub("\\s+neuron$", "", alias))
          available[node_data$region == region & neuronal]
        } else if (toupper(alias) %in% unique(as.character(node_data$region))) {
          available[as.character(node_data$region) == toupper(alias)]
        } else {
          .hybs_abort(
            "Unknown group alias: ", group,
            ". Supply `group_col` for a custom node grouping."
          )
        }
      }
    }
  )
  selected <- unique(as.character(selected))
  missing <- setdiff(selected, available)
  if (length(missing) > 0L) {
    .hybs_abort(
      "Subnetwork seed node(s) are absent from the network: ",
      paste(utils::head(missing, 5L), collapse = ", "), "."
    )
  }
  if (length(selected) == 0L) {
    .hybs_abort("Subnetwork query selected no nodes.")
  }
  selected
}

.network_graph_with_nodes <- function(network) {
  graph <- network$graph
  graph_names <- igraph::V(graph)$name
  index <- match(graph_names, network$nodes$name)
  for (column in setdiff(names(network$nodes), "name")) {
    graph <- igraph::set_vertex_attr(
      graph,
      column,
      value = network$nodes[[column]][index]
    )
  }
  graph
}

#' Extract a reusable subnetwork view
#'
#' @param network A `hybs_network` object.
#' @param mode One of `"ego"`, `"nodes"`, `"pattern"`, `"community"`, or
#'   `"group"`.
#' @param center_node Focal node for ego mode.
#' @param order Maximum graph distance used for neighbourhood expansion.
#' @param nodes Explicit seed nodes for nodes mode.
#' @param pattern Regular expression for pattern mode.
#' @param community Community value(s) for community mode.
#' @param community_col Community column in `network$nodes`.
#' @param group Group alias for group mode. Built-in aliases include regions,
#'   `"neuron"`, `"non-neuron"`, and values such as `"BS neuron"`.
#' @param group_col Optional node column for custom group matching.
#' @param neuronal_classes Node-name cell classes treated as neuronal.
#' @param include_neighbors Whether to expand non-ego seeds by `order` steps.
#' @param edge_mode `"induced"` retains all edges among selected nodes;
#'   `"incident"` retains only edges touching an original seed.
#' @param min_similarity Optional additional inclusive edge cutoff.
#' @param keep_isolates Whether selected nodes without retained edges remain.
#' @return A `hybs_subnetwork`, which also inherits from `hybs_network`.
#' @export
extract_subnetwork <- function(
    network,
    mode = c("ego", "nodes", "pattern", "community", "group"),
    center_node = NULL,
    order = 1L,
    nodes = character(),
    pattern = NULL,
    community = NULL,
    community_col = "Community",
    group = NULL,
    group_col = NULL,
    neuronal_classes = c("CHAT", "IN", "EX", "AVP", "OXT", "SST", "HDC"),
    include_neighbors = FALSE,
    edge_mode = c("induced", "incident"),
    min_similarity = NULL,
    keep_isolates = TRUE) {
  .assert_hybs_network(network)
  mode <- match.arg(mode)
  edge_mode <- match.arg(edge_mode)
  .assert_scalar_number(order, "order", lower = 0)
  .assert_flag(include_neighbors, "include_neighbors")
  .assert_flag(keep_isolates, "keep_isolates")
  if (!is.null(min_similarity)) {
    .assert_scalar_number(min_similarity, "min_similarity", lower = 0, upper = 1)
  }
  seeds <- .select_subnetwork_seeds(
    network, mode, center_node, nodes, pattern, community, community_col,
    group, group_col, neuronal_classes
  )
  graph <- .network_graph_with_nodes(network)
  if (!is.null(min_similarity) && igraph::ecount(graph) > 0L) {
    similarity <- igraph::edge_attr(graph, "Similarity")
    if (is.null(similarity)) {
      .hybs_abort("The network has no `Similarity` edge attribute.")
    }
    graph <- igraph::delete_edges(graph, which(similarity < min_similarity))
  }

  expand <- mode == "ego" || include_neighbors
  selected <- seeds
  if (expand && order > 0L) {
    neighbourhoods <- igraph::ego(
      graph,
      order = as.integer(order),
      nodes = seeds
    )
    selected <- unique(unlist(lapply(neighbourhoods, igraph::as_ids)))
  }
  local_graph <- igraph::induced_subgraph(graph, vids = selected)
  if (edge_mode == "incident" && igraph::ecount(local_graph) > 0L) {
    ends <- igraph::ends(local_graph, igraph::E(local_graph), names = TRUE)
    keep <- ends[, 1L] %in% seeds | ends[, 2L] %in% seeds
    local_graph <- igraph::delete_edges(local_graph, which(!keep))
  }
  if (!keep_isolates && igraph::vcount(local_graph) > 0L) {
    isolates <- which(igraph::degree(local_graph) == 0L)
    if (length(isolates) > 0L) {
      local_graph <- igraph::delete_vertices(local_graph, isolates)
    }
  }
  if (igraph::vcount(local_graph) == 0L) {
    .hybs_abort("No nodes remained after applying the subnetwork query.")
  }

  local_nodes <- igraph::as_data_frame(local_graph, what = "vertices")
  local_nodes$is_seed <- local_nodes$name %in% seeds
  local_nodes$is_center <- if (is.null(center_node)) {
    rep(FALSE, nrow(local_nodes))
  } else {
    local_nodes$name == center_node
  }
  local_nodes$role <- ifelse(
    local_nodes$is_center,
    "Center",
    ifelse(local_nodes$is_seed, "Seed", "Neighbour")
  )
  local_edges <- igraph::as_data_frame(local_graph, what = "edges")
  names(local_edges)[1:2] <- c("Cluster1", "Cluster2")
  if (!"connection" %in% names(local_edges) && nrow(local_edges) > 0L) {
    local_edges$connection <- ifelse(
      parse_region(local_edges$Cluster1) == parse_region(local_edges$Cluster2),
      "Within region",
      "Cross-region"
    )
  }
  parameters <- utils::modifyList(
    network$parameters,
    list(
      parent_graph_signature = network_signature(network),
      subnetwork_mode = mode,
      center_node = center_node,
      seed_nodes = seeds,
      neighbourhood_order = as.integer(order),
      include_neighbors = expand,
      edge_mode = edge_mode,
      subnetwork_min_similarity = min_similarity,
      keep_isolates = keep_isolates
    )
  )
  out <- list(
    graph = local_graph,
    nodes = local_nodes,
    edges = local_edges,
    parameters = parameters
  )
  class(out) <- c("hybs_subnetwork", "hybs_network")
  out
}

#' Calculate frozen or compact coordinates for a subnetwork
#'
#' @param subnetwork A `hybs_subnetwork` object.
#' @param mode `"frozen"` subsets the parent layout; `"compact"` calculates a
#'   deterministic layout for the extracted graph.
#' @param parent_layout Required for frozen mode.
#' @param seed,weights,iterations Compact-layout settings.
#' @return A `hybs_network_layout` object.
#' @export
calculate_subnetwork_layout <- function(
    subnetwork,
    mode = c("frozen", "compact"),
    parent_layout = NULL,
    seed = 42L,
    weights = "Similarity",
    iterations = 2000L) {
  if (!inherits(subnetwork, "hybs_subnetwork")) {
    .hybs_abort("`subnetwork` must be returned by `extract_subnetwork()`.")
  }
  mode <- match.arg(mode)
  if (mode == "frozen") {
    if (is.null(parent_layout)) {
      .hybs_abort("`parent_layout` is required for frozen subnetwork layout.")
    }
    .assert_columns(parent_layout, c("name", "x", "y"), "parent layout")
    wanted <- as.character(subnetwork$nodes$name)
    missing <- setdiff(wanted, as.character(parent_layout$name))
    if (length(missing) > 0L) {
      .hybs_abort(
        "Parent layout is missing subnetwork node(s): ",
        paste(utils::head(missing, 5L), collapse = ", "), "."
      )
    }
    subset <- parent_layout[
      match(wanted, parent_layout$name), c("name", "x", "y"), drop = FALSE
    ]
    return(as_network_layout(
      subset,
      subnetwork,
      layout_id = "parent_frozen_subset",
      layout_method = "frozen_parent_subset",
      layout_seed = NA_integer_,
      strict_signature = FALSE
    ))
  }
  out <- calculate_network_layout(
    subnetwork,
    method = "fr",
    seed = seed,
    weights = weights,
    iterations = iterations
  )
  out$layout_id <- "compact_subnetwork"
  out$layout_method <- "compact_fr"
  out
}

#' @export
print.hybs_subnetwork <- function(x, ...) {
  cat("<hybs_subnetwork>\n")
  cat("  mode:", x$parameters$subnetwork_mode, "\n")
  cat("  seeds:", length(x$parameters$seed_nodes), "\n")
  cat("  nodes:", igraph::vcount(x$graph), "\n")
  cat("  edges:", igraph::ecount(x$graph), "\n")
  invisible(x)
}
