#' Parse the anatomical region from a subcluster name
#'
#' The formal HyBs identifier is constructed from `nuclei + cluster`, for
#' example `SH_IN_10_ESR1`. The region is the prefix before the first underscore.
#'
#' @param cluster_name Character vector of subcluster identifiers.
#' @return A character vector containing region prefixes.
#' @export
parse_region <- function(cluster_name) {
  cluster_name <- as.character(cluster_name)
  sub("_.*$", "", cluster_name)
}

#' Read and validate a CSV table
#'
#' @param path Path to a CSV file.
#' @param required Required column names.
#' @return A data frame.
#' @export
read_hybs_csv <- function(path, required = character()) {
  .assert_nonempty_string(path, "path")
  if (!file.exists(path)) {
    .hybs_abort("Missing input file: ", path)
  }
  data <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (length(required) > 0L) {
    .assert_columns(data, required, path)
  }
  data
}

#' Prepare significant DEG records for network analysis
#'
#' This function creates the formal `region_cluster` identifier using both
#' `nuclei` and `cluster`, filters by BH-adjusted P value, and optionally limits
#' records to one disease condition. It does not run differential expression.
#'
#' @param data DEG result table.
#' @param condition Optional condition value such as `"Obesity"`.
#' @param padj_cutoff Strict adjusted-P cutoff. The current Figure 3 default is
#'   `0.05`.
#' @param gene_col,logfc_col,padj_col,condition_col,nuclei_col,cluster_col Column
#'   names in `data`.
#' @return A standardized data frame with class `hybs_deg_data`.
#' @export
prepare_deg_sets <- function(
    data,
    condition = NULL,
    padj_cutoff = 0.05,
    gene_col = "gene",
    logfc_col = "avg_log2FC",
    padj_col = "p_val_adj",
    condition_col = "group",
    nuclei_col = "nuclei",
    cluster_col = "cluster") {
  .assert_scalar_number(padj_cutoff, "padj_cutoff", lower = 0, upper = 1)
  requested <- c(gene_col, logfc_col, padj_col, condition_col, nuclei_col, cluster_col)
  .assert_columns(data, requested, "DEG data")

  out <- data.frame(
    gene = as.character(data[[gene_col]]),
    avg_log2FC = suppressWarnings(as.numeric(data[[logfc_col]])),
    p_val_adj = suppressWarnings(as.numeric(data[[padj_col]])),
    group = as.character(data[[condition_col]]),
    nuclei = as.character(data[[nuclei_col]]),
    cluster = as.character(data[[cluster_col]]),
    stringsAsFactors = FALSE
  )
  keep <- is.finite(out$p_val_adj) & out$p_val_adj < padj_cutoff
  if (!is.null(condition)) {
    .assert_nonempty_string(condition, "condition")
    keep <- keep & out$group == condition
  }
  out <- out[keep, , drop = FALSE]
  if (nrow(out) == 0L) {
    .hybs_abort("No DEG records remained after filtering.")
  }
  essential <- c("gene", "group", "nuclei", "cluster")
  bad <- Reduce(`|`, lapply(out[essential], function(x) is.na(x) | !nzchar(x)))
  if (any(bad)) {
    .hybs_abort("Filtered DEG data contain missing or empty identifiers.")
  }
  out$region_cluster <- paste(out$nuclei, out$cluster, sep = "_")
  key <- paste(out$region_cluster, out$gene, sep = "\r")
  if (anyDuplicated(key)) {
    examples <- unique(key[duplicated(key)])
    .hybs_abort(
      "DEG data contain duplicate gene/subcluster records; examples: ",
      paste(utils::head(examples, 3L), collapse = ", "), "."
    )
  }
  rownames(out) <- NULL
  class(out) <- c("hybs_deg_data", "data.frame")
  attr(out, "padj_cutoff") <- padj_cutoff
  attr(out, "condition") <- condition
  out
}
