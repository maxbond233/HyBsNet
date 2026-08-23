synthetic_degs <- function() {
  data.frame(
    gene = c("a", "b", "c", "b", "c", "d", "a", "e"),
    avg_log2FC = c(1, 1, -1, 1, 2, -1, -1, 1),
    p_val_adj = rep(0.01, 8),
    group = rep("Disease", 8),
    nuclei = c(rep("BS", 3), rep("SH", 3), rep("TH", 2)),
    cluster = c(rep("A", 3), rep("B", 3), rep("C", 2)),
    stringsAsFactors = FALSE
  )
}

triangle_edges <- function() {
  data.frame(
    Cluster1 = c("BS_A", "BS_A", "SH_B"),
    Cluster2 = c("SH_B", "TH_C", "TH_C"),
    Similarity = c(0.5, 0.25, 0.2),
    Total_Overlap = c(2L, 1L, 1L),
    stringsAsFactors = FALSE
  )
}
