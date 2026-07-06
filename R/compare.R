#' Extract edge list from a copula fit result
#'
#' @param result Output of [FitStratumCopula()].
#' @param label Stratum label.
#' @return Data frame with columns `from`, `to`, `pcor`, `stratum`, `edge`.
#' @export
ExtractEdges <- function(result, label) {
  if (is.null(result)) {
    return(NULL)
  }
  pcor <- result$pcor
  adj <- result$adjacency
  vars <- colnames(pcor)

  edges <- which(adj > 0 & upper.tri(adj), arr.ind = TRUE)
  if (nrow(edges) == 0) {
    return(NULL)
  }

  data.frame(
    from = vars[edges[, 1]],
    to = vars[edges[, 2]],
    pcor = pcor[edges],
    stratum = label,
    edge = paste(
      pmin(vars[edges[, 1]], vars[edges[, 2]]),
      pmax(vars[edges[, 1]], vars[edges[, 2]]),
      sep = " -- "
    ),
    stringsAsFactors = FALSE
  )
}

#' Extract all pairwise values from a symmetric matrix
#' @keywords internal
extract_pairwise_values <- function(mat, vars = colnames(mat), edge_only = FALSE, adj = NULL) {
  pairs <- combn(vars, 2, simplify = FALSE)
  if (length(pairs) == 0) {
    return(data.frame(
      varI = character(0),
      varJ = character(0),
      edge = character(0),
      value = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  out <- lapply(pairs, function(variable_pair) {
    i <- variable_pair[1]
    j <- variable_pair[2]
    val <- mat[i, j]
    if (edge_only) {
      if (is.null(adj) || adj[i, j] == 0) {
        val <- 0
      }
    }
    data.frame(
      varI = i,
      varJ = j,
      edge = paste(pmin(i, j), pmax(i, j), sep = " -- "),
      value = val,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}

#' Compare pairwise matrix values across two strata
#'
#' @param resA First stratum fit result.
#' @param resB Second stratum fit result.
#' @param labelA Label for stratum A.
#' @param labelB Label for stratum B.
#' @param matrixName One of `"pcor"` or `"copulaCor"`.
#' @param edgeOnly For `"pcor"`, if `TRUE` only adjacency edges are non-zero.
#'   For `"copulaCor"`, defaults to `FALSE` (all pairs).
#' @param deltaThreshold Threshold for direction classification.
#' @return Data frame with comparison columns.
#' @export
ComparePairwiseMatrices <- function(resA,
                                    resB,
                                    labelA,
                                    labelB,
                                    matrixName = c("pcor", "copulaCor"),
                                    edgeOnly = NULL,
                                    deltaThreshold = 0.05) {
  matrixName <- match.arg(matrixName)

  if (is.null(resA) || is.null(resB)) {
    stop("Both fit results must be non-NULL.", call. = FALSE)
  }

  mat_a <- resA[[matrixName]]
  mat_b <- resB[[matrixName]]
  shared <- intersect(colnames(mat_a), colnames(mat_b))

  if (length(shared) < 2) {
    stop("Fewer than 2 shared variables between strata.", call. = FALSE)
  }

  only_in_a <- setdiff(colnames(mat_a), colnames(mat_b))
  only_in_b <- setdiff(colnames(mat_b), colnames(mat_a))
  if (length(only_in_a) > 0 || length(only_in_b) > 0) {
    warning(
      "Variable mismatch - using intersection (n = ", length(shared), "). ",
      "Only in A: ", paste(only_in_a, collapse = ", "),
      "; only in B: ", paste(only_in_b, collapse = ", "),
      call. = FALSE
    )
  }

  mat_a <- mat_a[shared, shared, drop = FALSE]
  mat_b <- mat_b[shared, shared, drop = FALSE]

  if (is.null(edgeOnly)) {
    edgeOnly <- matrixName == "pcor"
  }

  adj_a <- if (edgeOnly && matrixName == "pcor") resA$adjacency[shared, shared, drop = FALSE] else NULL
  adj_b <- if (edgeOnly && matrixName == "pcor") resB$adjacency[shared, shared, drop = FALSE] else NULL

  pairs_a <- extract_pairwise_values(mat_a, shared, edge_only = edgeOnly, adj = adj_a)
  pairs_b <- extract_pairwise_values(mat_b, shared, edge_only = edgeOnly, adj = adj_b)

  comparison_result <- dplyr::full_join(
    pairs_a %>% dplyr::select(edge, valueA = value),
    pairs_b %>% dplyr::select(edge, valueB = value),
    by = "edge"
  ) %>%
    tidyr::replace_na(list(valueA = 0, valueB = 0)) %>%
    dplyr::mutate(
      delta = valueA - valueB,
      absDelta = abs(delta),
      direction = dplyr::case_when(
        abs(valueA) > 0 & abs(valueB) == 0 ~ paste("Only in", labelA),
        abs(valueA) == 0 & abs(valueB) > 0 ~ paste("Only in", labelB),
        delta > deltaThreshold ~ paste("Stronger in", labelA),
        delta < -deltaThreshold ~ paste("Stronger in", labelB),
        TRUE ~ "Similar"
      ),
      matrix = matrixName,
      labelA = labelA,
      labelB = labelB
    ) %>%
    dplyr::arrange(dplyr::desc(absDelta))

  comparison_result
}

#' Compare two fitted strata across matrix types
#'
#' @param resA First stratum fit result.
#' @param resB Second stratum fit result.
#' @param labelA Label for stratum A.
#' @param labelB Label for stratum B.
#' @param matrices Character vector of matrices to compare (`"pcor"`, `"copulaCor"`).
#' @param deltaThreshold Threshold for direction labels and bar charts.
#' @param edgeOnlyPcor If `TRUE`, pcor comparison uses adjacency edges only.
#' @param edgeOnlyCor If `FALSE` (default), copula cor compares all pairs.
#' @return List with comparison data frames keyed by matrix name.
#' @export
CompareTwoStrata <- function(resA,
                             resB,
                             labelA = "A",
                             labelB = "B",
                             matrices = c("pcor", "copulaCor"),
                             deltaThreshold = 0.05,
                             edgeOnlyPcor = TRUE,
                             edgeOnlyCor = FALSE) {
  out <- list()
  for (mat in matrices) {
    edge_only <- if (mat == "pcor") edgeOnlyPcor else edgeOnlyCor
    out[[mat]] <- ComparePairwiseMatrices(
      resA, resB, labelA, labelB,
      matrixName = mat,
      edgeOnly = edge_only,
      deltaThreshold = deltaThreshold
    )
  }
  structure(out, class = "CopulaStratumComparison")
}
