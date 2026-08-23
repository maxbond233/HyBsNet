#' Build a DEG similarity network
#'
#' @param edges Similarity edge table. Additional columns are retained as edge
#'   attributes.
#' @param deg_data Optional prepared DEG table used to calculate node DEG counts.
#' @param min_similarity Strict edge cutoff.
#' @param condition Optional condition stored in network metadata.
#' @return An object of class `hybs_network`.
#' @export
build_similarity_network <- function(
    edges,
    deg_data = NULL,
    min_similarity = 0.10,
    condition = NULL) {
  clean_edges <- prepare_similarity_edges(edges, min_similarity = min_similarity)
  if (!"Total_Overlap" %in% names(clean_edges)) {
    clean_edges$Total_Overlap <- NA_integer_
  }
  node_names <- sort(unique(c(clean_edges$Cluster1, clean_edges$Cluster2)))
  nodes <- data.frame(
    name = node_names,
    region = parse_region(node_names),
    deg_count = 0L,
    stringsAsFactors = FALSE
  )
  if (!is.null(deg_data)) {
    .assert_columns(deg_data, c("gene", "region_cluster"), "prepared DEG data")
    counts <- table(as.character(deg_data$region_cluster))
    matched <- unname(counts[nodes$name])
    matched[is.na(matched)] <- 0L
    nodes$deg_count <- as.integer(matched)
  }
  graph <- igraph::graph_from_data_frame(
    clean_edges,
    directed = FALSE,
    vertices = nodes
  )
  endpoint_regions <- parse_region(igraph::ends(graph, igraph::E(graph)))
  endpoint_regions <- matrix(endpoint_regions, ncol = 2L)
  igraph::E(graph)$connection <- ifelse(
    endpoint_regions[, 1L] == endpoint_regions[, 2L],
    "Within region",
    "Cross-region"
  )
  if (is.null(condition) && !is.null(deg_data)) {
    observed <- unique(as.character(deg_data$group))
    if (length(observed) == 1L) {
      condition <- observed
    }
  }
  out <- list(
    graph = graph,
    nodes = nodes,
    edges = clean_edges,
    parameters = list(
      condition = condition,
      min_similarity = min_similarity
    )
  )
  class(out) <- "hybs_network"
  out
}

#' @export
print.hybs_network <- function(x, ...) {
  cat("<hybs_network>\n")
  cat("  nodes:", igraph::vcount(x$graph), "\n")
  cat("  edges:", igraph::ecount(x$graph), "\n")
  condition <- x$parameters$condition %||% NA_character_
  if (!is.na(condition)) {
    cat("  condition:", condition, "\n")
  }
  cat("  similarity cutoff: >", x$parameters$min_similarity, "\n")
  invisible(x)
}

.assert_hybs_network <- function(x) {
  if (!inherits(x, "hybs_network") || !inherits(x$graph, "igraph")) {
    .hybs_abort("`network` must be a `hybs_network` object.")
  }
  invisible(x)
}

#' Extract network node data
#'
#' @param network A `hybs_network` object.
#' @return A node data frame.
#' @export
network_nodes <- function(network) {
  .assert_hybs_network(network)
  as.data.frame(network$nodes, stringsAsFactors = FALSE)
}

#' Extract network edge data
#'
#' @param network A `hybs_network` object.
#' @return An edge data frame.
#' @export
network_edges <- function(network) {
  .assert_hybs_network(network)
  as.data.frame(network$edges, stringsAsFactors = FALSE)
}

#' Create a stable topology signature
#'
#' The signature changes when nodes, endpoints, or similarity values change. It
#' intentionally ignores colours and other figure styling.
#'
#' @param network A `hybs_network` object.
#' @return A SHA-256 string.
#' @export
network_signature <- function(network) {
  .assert_hybs_network(network)
  edges <- network$edges[, c("Cluster1", "Cluster2", "Similarity"), drop = FALSE]
  edges <- edges[order(edges$Cluster1, edges$Cluster2), , drop = FALSE]
  payload <- list(
    nodes = sort(as.character(network$nodes$name)),
    edges = data.frame(
      Cluster1 = as.character(edges$Cluster1),
      Cluster2 = as.character(edges$Cluster2),
      Similarity = sprintf("%.17g", as.numeric(edges$Similarity)),
      stringsAsFactors = FALSE
    )
  )
  digest::digest(payload, algo = "sha256", serialize = TRUE)
}

#' Extract an induced neighbourhood network
#'
#' The current Figure 3G design uses `order = 1`, then retains all edges among
#' the centre and its first-order neighbours.
#'
#' @param network A `hybs_network` object.
#' @param center_node Name of the focal node.
#' @param order Maximum graph distance from the focal node.
#' @param label_top_n Number of strongest direct neighbours to label.
#' @param label_nodes Optional explicit node labels. When supplied, it overrides
#'   `label_top_n`; the centre node is always labelled.
#' @return A `hybs_network` object containing the induced subgraph.
#' @export
extract_local_network <- function(
    network,
    center_node,
    order = 1L,
    label_top_n = 15L,
    label_nodes = NULL) {
  .assert_hybs_network(network)
  .assert_nonempty_string(center_node, "center_node")
  .assert_scalar_number(order, "order", lower = 1)
  .assert_scalar_number(label_top_n, "label_top_n", lower = 0)
  graph <- network$graph
  if (!(center_node %in% igraph::V(graph)$name)) {
    .hybs_abort("Center node is absent from the network: ", center_node)
  }
  vertex_ids <- igraph::ego(graph, order = as.integer(order), nodes = center_node)[[1L]]
  local_graph <- igraph::induced_subgraph(graph, vids = vertex_ids)
  distances <- as.numeric(igraph::distances(local_graph, v = center_node)[1L, ])
  names(distances) <- igraph::V(local_graph)$name

  local_edges <- igraph::as_data_frame(local_graph, what = "edges")
  names(local_edges)[1:2] <- c("Cluster1", "Cluster2")
  direct <- local_edges[
    local_edges$Cluster1 == center_node | local_edges$Cluster2 == center_node,
    , drop = FALSE
  ]
  direct$neighbor <- ifelse(
    direct$Cluster1 == center_node,
    direct$Cluster2,
    direct$Cluster1
  )
  similarity <- if ("Similarity" %in% names(direct)) direct$Similarity else rep(1, nrow(direct))
  direct <- direct[order(-similarity, direct$neighbor), , drop = FALSE]
  if (is.null(label_nodes)) {
    selected <- utils::head(direct$neighbor, as.integer(label_top_n))
  } else {
    selected <- as.character(label_nodes)
  }
  selected <- unique(c(center_node, selected))

  nodes <- igraph::as_data_frame(local_graph, what = "vertices")
  if (!"region" %in% names(nodes)) {
    nodes$region <- parse_region(nodes$name)
  }
  if (!"deg_count" %in% names(nodes)) {
    nodes$deg_count <- 0L
  }
  nodes$distance <- unname(distances[nodes$name])
  nodes$node_type <- ifelse(
    nodes$name == center_node,
    "Center",
    paste0("Neighbor order ", nodes$distance)
  )
  nodes$label <- ifelse(nodes$name %in% selected, nodes$name, NA_character_)
  if ("connection" %in% names(local_edges)) {
    local_edges$connection <- as.character(local_edges$connection)
  }
  out <- list(
    graph = local_graph,
    nodes = nodes,
    edges = local_edges,
    parameters = utils::modifyList(
      network$parameters,
      list(
        center_node = center_node,
        neighborhood_order = as.integer(order),
        label_top_n = as.integer(label_top_n)
      )
    )
  )
  class(out) <- "hybs_network"
  out
}
