annotated_triangle_network <- function() {
  deg <- prepare_deg_sets(synthetic_degs(), condition = "Disease")
  net <- build_similarity_network(triangle_edges(), deg, condition = "Disease")
  net$nodes$Community <- c(1L, 1L, 2L)
  net
}

test_that("five subnetwork query modes select expected seeds", {
  net <- annotated_triangle_network()
  ego <- extract_subnetwork(net, mode = "ego", center_node = "BS_A", order = 1)
  explicit <- extract_subnetwork(net, mode = "nodes", nodes = c("BS_A", "SH_B"))
  pattern <- extract_subnetwork(net, mode = "pattern", pattern = "^(BS|SH)_")
  community <- extract_subnetwork(net, mode = "community", community = 1)
  group <- extract_subnetwork(net, mode = "group", group = "BS")
  expect_s3_class(ego, "hybs_subnetwork")
  expect_equal(nrow(network_nodes(ego)), 3L)
  expect_setequal(explicit$parameters$seed_nodes, c("BS_A", "SH_B"))
  expect_setequal(pattern$parameters$seed_nodes, c("BS_A", "SH_B"))
  expect_setequal(community$parameters$seed_nodes, c("BS_A", "SH_B"))
  expect_identical(group$parameters$seed_nodes, "BS_A")
})

test_that("incident and induced edge modes remain explicit", {
  net <- annotated_triangle_network()
  induced <- extract_subnetwork(
    net,
    mode = "ego",
    center_node = "BS_A",
    edge_mode = "induced"
  )
  incident <- extract_subnetwork(
    net,
    mode = "ego",
    center_node = "BS_A",
    edge_mode = "incident"
  )
  expect_equal(nrow(network_edges(induced)), 3L)
  expect_equal(nrow(network_edges(incident)), 2L)
})

test_that("subnetwork layouts separate frozen and compact coordinates", {
  net <- annotated_triangle_network()
  parent_layout <- calculate_network_layout(net, method = "circle")
  sub <- extract_subnetwork(net, mode = "nodes", nodes = c("BS_A", "SH_B"))
  frozen <- calculate_subnetwork_layout(
    sub,
    mode = "frozen",
    parent_layout = parent_layout
  )
  compact <- calculate_subnetwork_layout(sub, mode = "compact", seed = 42)
  expected <- parent_layout[match(frozen$name, parent_layout$name), c("x", "y")]
  expect_equal(frozen[, c("x", "y")], expected)
  expect_identical(unique(frozen$layout_method), "frozen_parent_subset")
  expect_identical(unique(compact$layout_method), "compact_fr")
})

test_that("subnetwork and context plots return ggplot objects", {
  net <- annotated_triangle_network()
  parent_layout <- calculate_network_layout(net, method = "circle")
  sub <- extract_subnetwork(net, mode = "ego", center_node = "BS_A")
  local_layout <- calculate_subnetwork_layout(
    sub,
    mode = "frozen",
    parent_layout = parent_layout
  )
  expect_s3_class(plot_subnetwork(sub, local_layout), "ggplot")
  expect_s3_class(
    plot_subnetwork_context(net, parent_layout, sub),
    "ggplot"
  )
})
