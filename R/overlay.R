# Node overlays -------------------------------------------------------------

.overlay_numeric_summary <- function(x, method) {
  x <- suppressWarnings(as.numeric(x))
  finite <- is.finite(x)
  if (!any(finite)) return(NA_real_)
  x <- x[finite]
  switch(
    method,
    mean = mean(x),
    median = stats::median(x),
    first = x[[1L]],
    .hybs_abort("Unsupported overlay aggregation method: ", method)
  )
}

.collapse_overlay_rows <- function(data, method) {
  key <- paste(data$condition, data$node, sep = "\r")
  if (!anyDuplicated(key)) return(data)
  if (method == "error") {
    duplicates <- unique(key[duplicated(key)])
    .hybs_abort(
      "Overlay contains duplicate condition/node records; examples: ",
      paste(utils::head(duplicates, 3L), collapse = ", "),
      ". Select an explicit aggregation method."
    )
  }
  pieces <- lapply(split(data, key), function(x) {
    data.frame(
      condition = x$condition[[1L]],
      node = x$node[[1L]],
      feature = x$feature[[1L]],
      score = .overlay_numeric_summary(x$score, method),
      score_type = x$score_type[[1L]],
      padj = if (any(is.finite(x$padj))) min(x$padj, na.rm = TRUE) else NA_real_,
      size_value = .overlay_numeric_summary(x$size_value, method),
      tested = any(x$tested %in% TRUE),
      available = any(x$available %in% TRUE),
      significant = any(x$significant %in% TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}

#' Prepare a generic node overlay
#'
#' The returned object separates a numeric score from testing status. Missing
#' or untested nodes remain `NA`; they are never converted to zero.
#'
#' @param data Input data frame.
#' @param node_col,value_col Node identifier and score columns.
#' @param feature Feature label, such as a gene or pathway name.
#' @param score_type Machine-readable score type, such as
#'   `"differential_log2fc"`, `"gsea_nes"`, or `"ucell_activity"`.
#' @param condition Optional condition to retain or assign.
#' @param condition_col Optional condition column.
#' @param padj_col,size_col,tested_col,significant_col Optional columns.
#' @param padj_cutoff Strict adjusted-P cutoff used when `significant_col` is
#'   omitted.
#' @param min_abs_score Optional minimum absolute score required for
#'   significance.
#' @param aggregate How duplicate condition/node records are handled.
#' @param score_label,size_label Human-readable legend labels.
#' @param method Optional method description retained as metadata.
#' @return A `hybs_node_overlay` data frame.
#' @export
prepare_node_overlay <- function(
    data,
    node_col,
    value_col,
    feature,
    score_type,
    condition = NULL,
    condition_col = NULL,
    padj_col = NULL,
    size_col = NULL,
    tested_col = NULL,
    significant_col = NULL,
    padj_cutoff = 0.05,
    min_abs_score = NULL,
    aggregate = c("error", "mean", "median", "first"),
    score_label = value_col,
    size_label = size_col,
    method = NULL) {
  .assert_data_frame(data, "data")
  .assert_nonempty_string(node_col, "node_col")
  .assert_nonempty_string(value_col, "value_col")
  .assert_nonempty_string(feature, "feature")
  .assert_nonempty_string(score_type, "score_type")
  .assert_scalar_number(padj_cutoff, "padj_cutoff", lower = 0, upper = 1)
  aggregate <- match.arg(aggregate)
  optional_columns <- c(
    condition_col, padj_col, size_col, tested_col, significant_col
  )
  optional_columns <- optional_columns[
    !vapply(optional_columns, is.null, logical(1))
  ]
  .assert_columns(data, c(node_col, value_col, optional_columns), "overlay data")

  input <- data
  if (!is.null(condition)) {
    .assert_nonempty_string(condition, "condition")
    if (!is.null(condition_col)) {
      input <- input[as.character(input[[condition_col]]) == condition, , drop = FALSE]
      if (nrow(input) == 0L) {
        .hybs_abort("No overlay records matched condition: ", condition, ".")
      }
    }
  }
  condition_values <- if (!is.null(condition_col)) {
    as.character(input[[condition_col]])
  } else {
    rep(condition %||% NA_character_, nrow(input))
  }
  node_values <- as.character(input[[node_col]])
  bad_node <- is.na(node_values) | !nzchar(node_values)
  if (any(bad_node)) {
    .hybs_abort("Overlay node identifiers must be non-missing and non-empty.")
  }
  score <- suppressWarnings(as.numeric(input[[value_col]]))
  padj <- if (is.null(padj_col)) {
    rep(NA_real_, nrow(input))
  } else {
    suppressWarnings(as.numeric(input[[padj_col]]))
  }
  invalid_padj <- !is.na(padj) & (!is.finite(padj) | padj < 0 | padj > 1)
  if (any(invalid_padj)) {
    .hybs_abort("Overlay adjusted P values must be missing or finite in [0, 1].")
  }
  size_value <- if (is.null(size_col)) {
    rep(NA_real_, nrow(input))
  } else {
    suppressWarnings(as.numeric(input[[size_col]]))
  }
  tested <- if (is.null(tested_col)) {
    rep(TRUE, nrow(input))
  } else {
    as.logical(input[[tested_col]])
  }
  tested[is.na(tested)] <- FALSE
  available <- is.finite(score)
  significant <- if (!is.null(significant_col)) {
    as.logical(input[[significant_col]])
  } else if (!is.null(padj_col)) {
    tested & available & is.finite(padj) & padj < padj_cutoff
  } else {
    rep(FALSE, nrow(input))
  }
  significant[is.na(significant)] <- FALSE
  if (!is.null(min_abs_score)) {
    .assert_scalar_number(min_abs_score, "min_abs_score", lower = 0)
    significant <- significant & available & abs(score) >= min_abs_score
  }

  out <- data.frame(
    condition = condition_values,
    node = node_values,
    feature = rep(feature, nrow(input)),
    score = score,
    score_type = rep(score_type, nrow(input)),
    padj = padj,
    size_value = size_value,
    tested = tested,
    available = available,
    significant = significant,
    stringsAsFactors = FALSE
  )
  out <- .collapse_overlay_rows(out, aggregate)
  out$available <- out$tested & is.finite(out$score)
  if (is.null(significant_col)) {
    out$significant <- if (!is.null(padj_col)) {
      out$tested & out$available & is.finite(out$padj) & out$padj < padj_cutoff
    } else {
      FALSE
    }
    if (!is.null(min_abs_score)) {
      out$significant <- out$significant & abs(out$score) >= min_abs_score
    }
  }
  class(out) <- c("hybs_node_overlay", "data.frame")
  attr(out, "score_label") <- score_label %||% value_col
  attr(out, "size_label") <- size_label
  attr(out, "method") <- method
  attr(out, "padj_cutoff") <- padj_cutoff
  attr(out, "min_abs_score") <- min_abs_score
  out
}

#' Prepare a gene-level network overlay
#'
#' By default node fill represents differential `avg_log2FC`, not absolute
#' expression. To plot an expression summary, supply its column through
#' `value_col` and set `score_type = "expression"`.
#'
#' @param data Gene-level table.
#' @param gene Gene to retain.
#' @param condition Optional condition to retain.
#' @param gene_col,condition_col,nuclei_col,cluster_col Input columns used to
#'   identify the feature, condition, and formal node key.
#' @param value_col,padj_col,size_col Score, adjusted-P, and optional size
#'   columns.
#' @param score_type Score interpretation.
#' @param padj_cutoff,min_abs_score Significance thresholds.
#' @param aggregate Duplicate handling passed to [prepare_node_overlay()].
#' @return A `hybs_node_overlay` object.
#' @export
prepare_gene_overlay <- function(
    data,
    gene,
    condition = NULL,
    gene_col = "gene",
    condition_col = "group",
    nuclei_col = "nuclei",
    cluster_col = "cluster",
    value_col = "avg_log2FC",
    padj_col = "p_val_adj",
    size_col = "pct.1",
    score_type = "differential_log2fc",
    padj_cutoff = 0.05,
    min_abs_score = log2(1.5),
    aggregate = c("error", "mean", "median", "first")) {
  .assert_nonempty_string(gene, "gene")
  aggregate <- match.arg(aggregate)
  required <- c(
    gene_col, condition_col, nuclei_col, cluster_col, value_col,
    padj_col, size_col
  )
  required <- required[!vapply(required, is.null, logical(1))]
  .assert_columns(data, required, "gene overlay data")
  input <- data[as.character(data[[gene_col]]) == gene, , drop = FALSE]
  if (nrow(input) == 0L) {
    .hybs_abort("Gene is absent from the overlay table: ", gene, ".")
  }
  input$.hybs_node <- paste(
    as.character(input[[nuclei_col]]),
    as.character(input[[cluster_col]]),
    sep = "_"
  )
  prepare_node_overlay(
    input,
    node_col = ".hybs_node",
    value_col = value_col,
    feature = gene,
    score_type = score_type,
    condition = condition,
    condition_col = condition_col,
    padj_col = padj_col,
    size_col = size_col,
    padj_cutoff = padj_cutoff,
    min_abs_score = min_abs_score,
    aggregate = aggregate,
    score_label = value_col,
    size_label = size_col,
    method = "gene"
  )
}

#' Prepare an average-expression network overlay
#'
#' This convenience wrapper is for an actual expression summary (for example
#' pseudobulk mean expression), rather than a differential-expression effect.
#'
#' @param data Gene-level expression summary.
#' @param gene Gene to retain.
#' @param condition Optional condition to retain.
#' @param gene_col,condition_col,nuclei_col,cluster_col Input identifier columns.
#' @param value_col Expression summary column.
#' @param size_col Optional detection-fraction or other non-negative size column.
#' @param tested_col Optional logical column distinguishing tested from untested
#'   records.
#' @param score_type Score interpretation stored in the overlay.
#' @param aggregate Duplicate handling passed to [prepare_node_overlay()].
#' @return A `hybs_node_overlay` object.
#' @export
prepare_expression_overlay <- function(
    data,
    gene,
    condition = NULL,
    gene_col = "gene",
    condition_col = "group",
    nuclei_col = "nuclei",
    cluster_col = "cluster",
    value_col = "avg_expression",
    size_col = "pct_expressing",
    tested_col = NULL,
    score_type = "expression",
    aggregate = c("error", "mean", "median", "first")) {
  .assert_nonempty_string(gene, "gene")
  aggregate <- match.arg(aggregate)
  required <- c(
    gene_col, condition_col, nuclei_col, cluster_col, value_col, size_col,
    tested_col
  )
  required <- required[!vapply(required, is.null, logical(1))]
  .assert_columns(data, required, "expression overlay data")
  input <- data[as.character(data[[gene_col]]) == gene, , drop = FALSE]
  if (nrow(input) == 0L) {
    .hybs_abort("Gene is absent from the expression table: ", gene, ".")
  }
  input$.hybs_node <- paste(
    as.character(input[[nuclei_col]]),
    as.character(input[[cluster_col]]),
    sep = "_"
  )
  prepare_node_overlay(
    input,
    node_col = ".hybs_node",
    value_col = value_col,
    feature = gene,
    score_type = score_type,
    condition = condition,
    condition_col = condition_col,
    size_col = size_col,
    tested_col = tested_col,
    aggregate = aggregate,
    score_label = value_col,
    size_label = size_col,
    method = "expression summary"
  )
}

#' Prepare a pathway-level network overlay
#'
#' The default score is GSEA NES. Scores from UCell, AUCell, GSVA, or another
#' method can be supplied by changing `value_col`, `score_type`, and `method`.
#'
#' @param data Pathway result table.
#' @param pathway Pathway to retain.
#' @param condition Optional condition to retain.
#' @param pathway_col,condition_col,nuclei_col,cluster_col Input columns.
#' @param value_col,padj_col Score and adjusted-P columns.
#' @param score_type Score interpretation.
#' @param method Method description retained as metadata.
#' @param padj_cutoff Adjusted-P threshold.
#' @param aggregate Duplicate handling passed to [prepare_node_overlay()].
#' @return A `hybs_node_overlay` object.
#' @export
prepare_pathway_overlay <- function(
    data,
    pathway,
    condition = NULL,
    pathway_col = "pathway",
    condition_col = "group",
    nuclei_col = "nuclei",
    cluster_col = "cluster",
    value_col = "NES",
    padj_col = "padj",
    score_type = "gsea_nes",
    method = "GSEA",
    padj_cutoff = 0.05,
    aggregate = c("error", "mean", "median", "first")) {
  .assert_nonempty_string(pathway, "pathway")
  aggregate <- match.arg(aggregate)
  required <- c(
    pathway_col, condition_col, nuclei_col, cluster_col, value_col, padj_col
  )
  required <- required[!vapply(required, is.null, logical(1))]
  .assert_columns(data, required, "pathway overlay data")
  input <- data[as.character(data[[pathway_col]]) == pathway, , drop = FALSE]
  if (nrow(input) == 0L) {
    .hybs_abort("Pathway is absent from the overlay table: ", pathway, ".")
  }
  input$.hybs_node <- paste(
    as.character(input[[nuclei_col]]),
    as.character(input[[cluster_col]]),
    sep = "_"
  )
  prepare_node_overlay(
    input,
    node_col = ".hybs_node",
    value_col = value_col,
    feature = pathway,
    score_type = score_type,
    condition = condition,
    condition_col = condition_col,
    padj_col = padj_col,
    padj_cutoff = padj_cutoff,
    aggregate = aggregate,
    score_label = value_col,
    method = method
  )
}

#' Validate a node overlay
#'
#' @param overlay A `hybs_node_overlay` object.
#' @param network Optional network used to check node identifiers.
#' @param condition Optional condition to validate.
#' @param allow_partial Whether the overlay may omit network nodes.
#' @return `overlay`, invisibly.
#' @export
validate_node_overlay <- function(
    overlay,
    network = NULL,
    condition = NULL,
    allow_partial = TRUE) {
  .assert_flag(allow_partial, "allow_partial")
  if (!inherits(overlay, "hybs_node_overlay")) {
    .hybs_abort("`overlay` must be returned by an overlay preparation function.")
  }
  required <- c(
    "condition", "node", "feature", "score", "score_type", "padj",
    "size_value", "tested", "available", "significant"
  )
  .assert_columns(overlay, required, "node overlay")
  if (anyDuplicated(paste(overlay$condition, overlay$node, sep = "\r"))) {
    .hybs_abort("Node overlay must contain one record per condition/node.")
  }
  if (!is.null(condition)) {
    .assert_nonempty_string(condition, "condition")
  }
  if (!is.null(network)) {
    .assert_hybs_network(network)
    subset <- overlay
    if (!is.null(condition)) {
      subset <- subset[is.na(subset$condition) | subset$condition == condition, , drop = FALSE]
    }
    extra <- setdiff(as.character(subset$node), as.character(network$nodes$name))
    if (length(extra) > 0L) {
      .hybs_abort(
        "Overlay contains node(s) absent from the network: ",
        paste(utils::head(extra, 5L), collapse = ", "), "."
      )
    }
    if (!allow_partial) {
      missing <- setdiff(as.character(network$nodes$name), as.character(subset$node))
      if (length(missing) > 0L) {
        .hybs_abort(
          "Overlay is missing network node(s): ",
          paste(utils::head(missing, 5L), collapse = ", "), "."
        )
      }
    }
  }
  invisible(overlay)
}

#' Attach one overlay to network nodes
#'
#' @param network A `hybs_network` object.
#' @param overlay A `hybs_node_overlay` object.
#' @param condition Condition to select. By default the condition stored in the
#'   network is used.
#' @param strict If `TRUE`, overlay-only nodes are rejected.
#' @return A copy of `network` with standardized `overlay_*` node attributes.
#' @export
attach_node_overlay <- function(
    network,
    overlay,
    condition = NULL,
    strict = TRUE) {
  .assert_hybs_network(network)
  .assert_flag(strict, "strict")
  validate_node_overlay(overlay)
  target_condition <- condition %||% network$parameters$condition
  observed_conditions <- unique(as.character(overlay$condition))
  observed_conditions <- observed_conditions[
    !is.na(observed_conditions) & nzchar(observed_conditions)
  ]
  if (is.null(target_condition)) {
    if (length(observed_conditions) == 1L) {
      target_condition <- observed_conditions
    } else if (length(observed_conditions) > 1L) {
      .hybs_abort("Specify `condition` when an overlay contains multiple conditions.")
    }
  }
  selected <- overlay
  if (!is.null(target_condition)) {
    .assert_nonempty_string(target_condition, "condition")
    selected <- selected[
      is.na(selected$condition) | selected$condition == target_condition,
      , drop = FALSE
    ]
  }
  if (nrow(selected) == 0L) {
    .hybs_abort("No overlay records matched the target network condition.")
  }
  extra <- setdiff(as.character(selected$node), as.character(network$nodes$name))
  if (length(extra) > 0L && strict) {
    .hybs_abort(
      "Overlay contains node(s) absent from the network: ",
      paste(utils::head(extra, 5L), collapse = ", "), "."
    )
  }
  n_input <- nrow(selected)
  selected <- selected[selected$node %in% network$nodes$name, , drop = FALSE]
  index <- match(network$nodes$name, selected$node)
  out <- network
  assignments <- list(
    overlay_feature = selected$feature[index],
    overlay_score = selected$score[index],
    overlay_score_type = selected$score_type[index],
    overlay_padj = selected$padj[index],
    overlay_size = selected$size_value[index],
    overlay_tested = selected$tested[index],
    overlay_available = selected$available[index],
    overlay_significant = selected$significant[index]
  )
  assignments$overlay_tested[is.na(assignments$overlay_tested)] <- FALSE
  assignments$overlay_available[is.na(assignments$overlay_available)] <- FALSE
  assignments$overlay_significant[is.na(assignments$overlay_significant)] <- FALSE
  for (column in names(assignments)) {
    out$nodes[[column]] <- assignments[[column]]
    out$graph <- igraph::set_vertex_attr(
      out$graph,
      column,
      value = assignments[[column]]
    )
  }
  out$parameters$overlay <- list(
    condition = target_condition,
    feature = unique(as.character(selected$feature)),
    score_type = unique(as.character(selected$score_type)),
    score_label = attr(overlay, "score_label"),
    size_label = attr(overlay, "size_label"),
    method = attr(overlay, "method"),
    n_input = n_input,
    n_network_nodes = nrow(network$nodes),
    n_ignored = length(extra),
    ignored_nodes = sort(extra),
    n_available = sum(assignments$overlay_available),
    n_significant = sum(assignments$overlay_significant)
  )
  out
}

#' @export
print.hybs_node_overlay <- function(x, ...) {
  cat("<hybs_node_overlay>\n")
  cat("  feature:", unique(as.character(x$feature)), "\n")
  cat("  score type:", unique(as.character(x$score_type)), "\n")
  cat("  records:", nrow(x), "\n")
  cat("  available:", sum(x$available), "\n")
  cat("  significant:", sum(x$significant), "\n")
  invisible(x)
}
