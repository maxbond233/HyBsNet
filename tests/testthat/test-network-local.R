test_that("network construction and signatures are deterministic", {
  deg <- prepare_deg_sets(synthetic_degs(), condition = "Disease")
  net1 <- build_similarity_network(triangle_edges(), deg, condition = "Disease")
  net2 <- build_similarity_network(triangle_edges()[3:1, ], deg, condition = "Disease")
  expect_s3_class(net1, "hybs_network")
  expect_equal(nrow(network_nodes(net1)), 3L)
  expect_equal(nrow(network_edges(net1)), 3L)
  expect_identical(network_signature(net1), network_signature(net2))
})

test_that("local network is induced and labels strongest neighbours", {
  deg <- prepare_deg_sets(synthetic_degs(), condition = "Disease")
  net <- build_similarity_network(triangle_edges(), deg)
  local <- extract_local_network(net, "BS_A", label_top_n = 1L)
  expect_equal(nrow(network_nodes(local)), 3L)
  expect_equal(nrow(network_edges(local)), 3L)
  expect_setequal(
    network_nodes(local)$name[!is.na(network_nodes(local)$label)],
    c("BS_A", "SH_B")
  )
  expect_error(extract_local_network(net, "missing"), "absent")
})
