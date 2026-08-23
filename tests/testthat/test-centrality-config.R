test_that("centrality and core-periphery helpers return explicit contracts", {
  deg <- prepare_deg_sets(synthetic_degs(), condition = "Disease")
  net <- build_similarity_network(triangle_edges(), deg)
  centrality <- calculate_network_centrality(
    net,
    min_community_size = 1,
    min_degree = 0
  )
  expect_true(all(c(
    "Cluster", "Nuclei", "Degree_raw", "Degree", "Centrality_Score"
  ) %in% names(centrality)))
  prepared <- prepare_core_periphery_data(
    centrality,
    deg,
    condition = "Disease"
  )
  expect_false(anyNA(prepared$total_degs))
  expect_true(all(prepared$network_position %in% c("Core", "Periphery")))
  expect_false(anyNA(prepared$network_position))
})

test_that("YAML style configuration is reusable", {
  path <- system.file("config", "figure3_style.yml", package = "HyBsNet")
  config <- read_network_config(path)
  style <- style_from_config(config)
  expect_s3_class(style, "hybs_network_style")
  expect_identical(unname(style$region_colors[["BS"]]), "#316B9D")
})
