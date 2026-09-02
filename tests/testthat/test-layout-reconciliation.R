test_that("external and frozen layouts retain exact coordinates", {
  deg <- prepare_deg_sets(synthetic_degs(), condition = "Disease")
  net <- build_similarity_network(triangle_edges(), deg, condition = "Disease")
  external <- data.frame(
    node = network_nodes(net)$name,
    x_redesign = c(-1, 0, 1),
    y_redesign = c(0.5, -0.5, 0.5),
    stringsAsFactors = FALSE
  )
  layout <- as_network_layout(
    external,
    net,
    node_col = "node",
    x_col = "x_redesign",
    y_col = "y_redesign",
    layout_id = "community_first_v1"
  )
  frozen <- freeze_network_layout(layout, net, "community_first_v1")
  expect_s3_class(frozen, "hybs_network_layout")
  expect_equal(layout_nodes(frozen)$x, external$x_redesign)
  expect_equal(unique(frozen$layout_id), "community_first_v1")
  expect_invisible(validate_network_layout(frozen, net))
})

test_that("layout reconciliation freezes shared nodes and anchors new nodes", {
  deg <- prepare_deg_sets(synthetic_degs(), condition = "Disease")
  baseline <- build_similarity_network(triangle_edges(), deg)
  baseline_layout <- calculate_network_layout(baseline, method = "circle")
  changed_edges <- rbind(
    triangle_edges(),
    data.frame(
      Cluster1 = "BS_A",
      Cluster2 = "SH_D",
      Similarity = 0.3,
      Total_Overlap = 1L
    )
  )
  changed <- build_similarity_network(changed_edges)
  expect_error(
    reconcile_network_layout(baseline_layout, changed),
    "absent from the frozen layout"
  )
  reconciled <- reconcile_network_layout(
    baseline_layout,
    changed,
    new_nodes = "anchor_neighbors",
    layout_id = "anchored_v1"
  )
  comparison <- compare_network_layouts(baseline_layout, reconciled)
  expect_true(all(comparison$coordinates_identical))
  expect_identical(attr(reconciled, "added_nodes"), "SH_D")
  expect_true(all(is.finite(reconciled$x)))
  expect_invisible(validate_network_layout(reconciled, changed))
})

test_that("layout comparison records changed node sets", {
  reference <- data.frame(name = c("A", "B"), x = c(0, 1), y = c(0, 1))
  candidate <- data.frame(name = c("B", "C"), x = c(1, 2), y = c(1, 2))
  comparison <- compare_network_layouts(reference, candidate)
  expect_identical(comparison$name, "B")
  expect_true(comparison$coordinates_identical)
  expect_identical(attr(comparison, "reference_only"), "A")
  expect_identical(attr(comparison, "candidate_only"), "C")
})
