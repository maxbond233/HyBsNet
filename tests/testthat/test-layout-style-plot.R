test_that("layout cache is tied to topology and preserves caller RNG", {
  deg <- prepare_deg_sets(synthetic_degs(), condition = "Disease")
  net <- build_similarity_network(triangle_edges(), deg)
  set.seed(100)
  expected_next <- runif(1)
  set.seed(100)
  layout <- calculate_network_layout(net, method = "fr", seed = 42)
  observed_next <- runif(1)
  expect_equal(observed_next, expected_next)
  expect_invisible(validate_network_layout(layout, net))
  path <- tempfile(fileext = ".csv")
  save_network_layout(layout, path)
  loaded <- load_network_layout(path, net)
  expect_equal(loaded[, c("name", "x", "y")], layout[, c("name", "x", "y")])
})

test_that("styles change appearance without changing topology or coordinates", {
  deg <- prepare_deg_sets(synthetic_degs(), condition = "Disease")
  net <- build_similarity_network(triangle_edges(), deg)
  layout <- calculate_network_layout(net, method = "circle")
  signature <- network_signature(net)
  p1 <- plot_similarity_network(net, layout, style = theme_fig3_current())
  p2 <- plot_similarity_network(net, layout, style = theme_fig3_colorblind())
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_identical(network_signature(net), signature)
  expect_equal(layout$x, layout$x)
})

test_that("local plots and exports are non-overwriting", {
  deg <- prepare_deg_sets(synthetic_degs(), condition = "Disease")
  net <- build_similarity_network(triangle_edges(), deg)
  local <- extract_local_network(net, "BS_A", label_top_n = 1L)
  layout <- calculate_network_layout(local, method = "circle", weights = NULL)
  plot <- plot_local_network(local, layout)
  expect_s3_class(plot, "ggplot")
  out <- tempfile(pattern = "hybsnet-")
  dir.create(out)
  paths <- export_network_panel(
    plot, out, "panel", 50, 50, formats = "png"
  )
  expect_true(file.exists(paths[["png"]]))
  expect_error(
    export_network_panel(plot, out, "panel", 50, 50, formats = "png"),
    "Refusing to overwrite"
  )
})
