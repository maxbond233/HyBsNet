#' Detect network communities
#'
#' @param network A `hybs_network` object.
#' @param method Community algorithm: `"louvain"` or `"fast_greedy"`.
#' @param weights `NULL` for the current unweighted Figure 3 behaviour, or the
#'   name of a numeric edge attribute such as `"Similarity"`.
#' @param seed Random seed used without changing the caller's global RNG state.
#' @return A data frame with node, region, and community membership.
#' @export
detect_network_communities <- function(
    network,
    method = c("louvain", "fast_greedy"),
    weights = NULL,
    seed = 42L) {
  .assert_hybs_network(network)
  method <- match.arg(method)
  graph <- network$graph
  weight_values <- NULL
  if (!is.null(weights)) {
    .assert_nonempty_string(weights, "weights")
    weight_values <- igraph::edge_attr(graph, weights)
    if (is.null(weight_values) || !is.numeric(weight_values) || any(!is.finite(weight_values))) {
      .hybs_abort("Requested community weights are absent or non-finite: ", weights)
    }
  }
  community <- .with_seed(seed, {
    if (method == "louvain") {
      igraph::cluster_louvain(graph, weights = weight_values)
    } else {
      igraph::cluster_fast_greedy(graph, weights = weight_values)
    }
  })
  data.frame(
    Cluster = igraph::V(graph)$name,
    Nuclei = parse_region(igraph::V(graph)$name),
    Community = as.integer(igraph::membership(community)),
    stringsAsFactors = FALSE
  )
}

#' Calculate the Figure 3 composite centrality score
#'
#' The defaults preserve the existing logic: detect unweighted Louvain
#' communities, retain communities with at least three nodes and nodes with at
#' least two links, calculate four unweighted centralities in the induced
#' subgraph, z-standardize them, and average the standardized values.
#'
#' @param network A `hybs_network` object.
#' @param min_community_size Minimum community size.
#' @param min_degree Minimum degree in the full graph.
#' @param community_weights Optional edge attribute for community detection.
#' @param centrality_weights Optional edge attribute for centrality. `NULL`
#'   preserves the current unweighted workflow.
#' @param seed Random seed.
#' @return A data frame of raw, standardized, and composite centrality values.
#' @export
calculate_network_centrality <- function(
    network,
    min_community_size = 3L,
    min_degree = 2L,
    community_weights = NULL,
    centrality_weights = NULL,
    seed = 42L) {
  .assert_hybs_network(network)
  .assert_scalar_number(min_community_size, "min_community_size", lower = 1)
  .assert_scalar_number(min_degree, "min_degree", lower = 0)
  graph <- network$graph
  membership <- detect_network_communities(
    network,
    method = "louvain",
    weights = community_weights,
    seed = seed
  )
  sizes <- table(membership$Community)
  main_communities <- as.integer(names(sizes)[sizes >= as.integer(min_community_size)])
  full_degree <- igraph::degree(graph)
  keep <- membership$Cluster[
    membership$Community %in% main_communities &
      full_degree[membership$Cluster] >= as.integer(min_degree)
  ]
  if (length(keep) < 2L) {
    .hybs_abort("Fewer than two nodes remained for centrality analysis.")
  }
  subgraph <- igraph::induced_subgraph(graph, vids = keep)
  weight_values <- NULL
  if (!is.null(centrality_weights)) {
    .assert_nonempty_string(centrality_weights, "centrality_weights")
    weight_values <- igraph::edge_attr(subgraph, centrality_weights)
    if (is.null(weight_values) || any(!is.finite(weight_values)) || any(weight_values <= 0)) {
      .hybs_abort("Centrality weights must be present, finite, and positive.")
    }
  }
  degree_raw <- as.numeric(igraph::degree(subgraph))
  betweenness_raw <- as.numeric(igraph::betweenness(subgraph, weights = weight_values))
  closeness_raw <- as.numeric(igraph::closeness(subgraph, weights = weight_values))
  eigenvector_raw <- as.numeric(igraph::eigen_centrality(subgraph, weights = weight_values)$vector)
  degree_z <- .scaled_vector(degree_raw)
  betweenness_z <- .scaled_vector(betweenness_raw)
  closeness_z <- .scaled_vector(closeness_raw)
  eigenvector_z <- .scaled_vector(eigenvector_raw)
  score_matrix <- cbind(degree_z, betweenness_z, closeness_z, eigenvector_z)
  score <- rowMeans(score_matrix, na.rm = TRUE)
  score[rowSums(is.finite(score_matrix)) == 0L] <- NA_real_
  community_map <- stats::setNames(membership$Community, membership$Cluster)
  out <- data.frame(
    Cluster = igraph::V(subgraph)$name,
    Nuclei = parse_region(igraph::V(subgraph)$name),
    Community = unname(community_map[igraph::V(subgraph)$name]),
    Degree_raw = degree_raw,
    Betweenness_raw = betweenness_raw,
    Closeness_raw = closeness_raw,
    Eigenvector_raw = eigenvector_raw,
    Degree = degree_z,
    Betweenness = betweenness_z,
    Closeness = closeness_z,
    Eigenvector = eigenvector_z,
    Centrality_Score = score,
    stringsAsFactors = FALSE
  )
  out <- out[order(-out$Centrality_Score, out$Cluster), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Prepare network centrality versus DEG-burden data
#'
#' @param centrality_data A table containing `Cluster`, `Nuclei`, and
#'   `Centrality_Score`.
#' @param deg_data Prepared DEG records for one condition.
#' @param condition Optional condition label.
#' @param cutoff `"median"` or a numeric centrality threshold.
#' @return A panel-ready data frame with `total_degs` and `network_position`.
#' @export
prepare_core_periphery_data <- function(
    centrality_data,
    deg_data,
    condition = NULL,
    cutoff = "median") {
  .assert_columns(
    centrality_data,
    c("Cluster", "Nuclei", "Centrality_Score"),
    "centrality data"
  )
  .assert_columns(deg_data, c("region_cluster", "gene"), "prepared DEG data")
  counts <- table(as.character(deg_data$region_cluster))
  out <- as.data.frame(centrality_data, stringsAsFactors = FALSE)
  out$total_degs <- as.integer(unname(counts[as.character(out$Cluster)]))
  out$total_degs[is.na(out$total_degs)] <- 0L
  if (identical(cutoff, "median")) {
    threshold <- stats::median(out$Centrality_Score, na.rm = TRUE)
  } else {
    .assert_scalar_number(cutoff, "cutoff")
    threshold <- cutoff
  }
  out$network_position <- ifelse(
    out$Centrality_Score > threshold,
    "Core",
    "Periphery"
  )
  condition <- condition %||% attr(deg_data, "condition")
  out$condition <- condition %||% NA_character_
  out
}
