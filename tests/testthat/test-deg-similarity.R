test_that("DEG preparation uses nuclei plus cluster and strict padj cutoff", {
  raw <- synthetic_degs()
  raw <- rbind(
    raw,
    data.frame(
      gene = "excluded",
      avg_log2FC = 1,
      p_val_adj = 0.05,
      group = "Disease",
      nuclei = "BS",
      cluster = "A"
    )
  )
  deg <- prepare_deg_sets(raw, condition = "Disease", padj_cutoff = 0.05)
  expect_s3_class(deg, "hybs_deg_data")
  expect_setequal(unique(deg$region_cluster), c("BS_A", "SH_B", "TH_C"))
  expect_false("excluded" %in% deg$gene)
})

test_that("pairwise Jaccard and direction categories are exact", {
  deg <- prepare_deg_sets(synthetic_degs(), condition = "Disease")
  edges <- calculate_pairwise_jaccard(deg, min_similarity = 0.10)
  expect_equal(nrow(edges), 2L)
  row <- edges[edges$Cluster1 == "BS_A" & edges$Cluster2 == "SH_B", ]
  expect_equal(row$Similarity, 0.5)
  expect_equal(row$Both_Up_Count, 1L)
  expect_equal(row$Opposite_Count, 1L)
  expect_equal(row$Total_Overlap, 2L)
})

test_that("duplicate DEG keys and duplicate unordered edges fail loudly", {
  raw <- synthetic_degs()
  expect_error(prepare_deg_sets(rbind(raw, raw[1, ])), "duplicate")
  edges <- rbind(triangle_edges(), triangle_edges()[1, c(2, 1, 3, 4)])
  names(edges) <- c("Cluster1", "Cluster2", "Similarity", "Total_Overlap")
  expect_error(prepare_similarity_edges(edges), "duplicate unordered")
})
