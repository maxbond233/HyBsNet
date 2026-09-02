gene_overlay_data <- function() {
  data.frame(
    gene = rep("SST", 3),
    group = rep("Disease", 3),
    nuclei = c("BS", "SH", "TH"),
    cluster = c("A", "B", "C"),
    avg_log2FC = c(1.2, -0.8, NA),
    pct.1 = c(0.7, 0.4, NA),
    p_val_adj = c(0.001, 0.2, NA),
    stringsAsFactors = FALSE
  )
}

test_that("gene overlays distinguish unavailable from nonsignificant nodes", {
  overlay <- prepare_gene_overlay(
    gene_overlay_data(),
    gene = "SST",
    condition = "Disease"
  )
  expect_s3_class(overlay, "hybs_node_overlay")
  expect_identical(overlay$node, c("BS_A", "SH_B", "TH_C"))
  expect_identical(overlay$available, c(TRUE, TRUE, FALSE))
  expect_identical(overlay$significant, c(TRUE, FALSE, FALSE))
  expect_true(is.na(overlay$score[overlay$node == "TH_C"]))
})

test_that("attaching an overlay does not change topology", {
  deg <- prepare_deg_sets(synthetic_degs(), condition = "Disease")
  net <- build_similarity_network(triangle_edges(), deg, condition = "Disease")
  signature <- network_signature(net)
  overlay <- prepare_gene_overlay(
    gene_overlay_data(),
    gene = "SST",
    condition = "Disease"
  )
  attached <- attach_node_overlay(net, overlay)
  expect_identical(network_signature(attached), signature)
  expect_true(all(c("overlay_score", "overlay_tested") %in% names(attached$nodes)))
  expect_equal(attached$nodes$overlay_score[attached$nodes$name == "BS_A"], 1.2)
})

test_that("partial overlays retain missing network nodes as untested NA", {
  deg <- prepare_deg_sets(synthetic_degs(), condition = "Disease")
  net <- build_similarity_network(triangle_edges(), deg, condition = "Disease")
  partial <- prepare_node_overlay(
    data.frame(node = "BS_A", value = 2),
    node_col = "node",
    value_col = "value",
    feature = "custom",
    score_type = "custom"
  )
  attached <- attach_node_overlay(net, partial)
  missing <- attached$nodes$name != "BS_A"
  expect_true(all(is.na(attached$nodes$overlay_score[missing])))
  expect_true(all(!attached$nodes$overlay_tested[missing]))
})

test_that("duplicate overlay rows require explicit aggregation", {
  duplicate <- data.frame(node = c("BS_A", "BS_A"), value = c(1, 3))
  expect_error(
    prepare_node_overlay(
      duplicate,
      node_col = "node",
      value_col = "value",
      feature = "custom",
      score_type = "custom"
    ),
    "duplicate"
  )
  overlay <- prepare_node_overlay(
    duplicate,
    node_col = "node",
    value_col = "value",
    feature = "custom",
    score_type = "custom",
    aggregate = "mean"
  )
  expect_equal(overlay$score, 2)
})

test_that("pathway overlays retain method and score semantics", {
  pathway_data <- data.frame(
    group = "Disease",
    nuclei = c("BS", "SH"),
    cluster = c("A", "B"),
    pathway = "HALLMARK_OXPHOS",
    NES = c(2.1, -1.4),
    padj = c(0.01, 0.2)
  )
  overlay <- prepare_pathway_overlay(
    pathway_data,
    pathway = "HALLMARK_OXPHOS",
    condition = "Disease"
  )
  expect_identical(unique(overlay$score_type), "gsea_nes")
  expect_identical(attr(overlay, "method"), "GSEA")
})

test_that("expression overlays are not labelled as differential effects", {
  expression_data <- data.frame(
    gene = "SST",
    group = "Disease",
    nuclei = c("BS", "SH"),
    cluster = c("A", "B"),
    avg_expression = c(2.5, 1.1),
    pct_expressing = c(0.8, 0.3)
  )
  overlay <- prepare_expression_overlay(
    expression_data,
    gene = "SST",
    condition = "Disease"
  )
  expect_identical(unique(overlay$score_type), "expression")
  expect_true(all(!overlay$significant))
})

test_that("overlay plots use an existing frozen layout", {
  deg <- prepare_deg_sets(synthetic_degs(), condition = "Disease")
  net <- build_similarity_network(triangle_edges(), deg, condition = "Disease")
  layout <- calculate_network_layout(net, method = "circle")
  overlay <- prepare_gene_overlay(
    gene_overlay_data(),
    gene = "SST",
    condition = "Disease"
  )
  plot <- plot_network_overlay(
    net,
    overlay,
    layout = layout,
    size_by = "overlay"
  )
  expect_s3_class(plot, "ggplot")
  expect_invisible(validate_network_layout(layout, net))
})
