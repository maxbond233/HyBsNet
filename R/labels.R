#' Select network labels with a reusable rule
#'
#' @param network A `hybs_network` object.
#' @param mode One of `"explicit"`, `"top_degree"`, `"top_strength"`,
#'   `"existing"`, or `"none"`.
#' @param nodes Explicit node names.
#' @param n Number of automatically selected nodes.
#' @return A character vector of node names.
#' @export
select_network_labels <- function(
    network,
    mode = c("explicit", "top_degree", "top_strength", "existing", "none"),
    nodes = character(),
    n = 10L) {
  .assert_hybs_network(network)
  mode <- match.arg(mode)
  .assert_scalar_number(n, "n", lower = 0)
  graph <- network$graph
  available <- igraph::V(graph)$name
  selected <- switch(
    mode,
    explicit = as.character(nodes),
    top_degree = {
      score <- igraph::degree(graph)
      names(sort(score, decreasing = TRUE))[seq_len(min(as.integer(n), length(score)))]
    },
    top_strength = {
      weight <- igraph::edge_attr(graph, "Similarity")
      if (is.null(weight)) {
        .hybs_abort("The network has no `Similarity` edge attribute.")
      }
      score <- igraph::strength(graph, weights = weight)
      names(sort(score, decreasing = TRUE))[seq_len(min(as.integer(n), length(score)))]
    },
    existing = {
      if (!"label" %in% names(network$nodes)) character() else {
        as.character(network$nodes$name[!is.na(network$nodes$label) & nzchar(network$nodes$label)])
      }
    },
    none = character()
  )
  unique(intersect(selected, available))
}
