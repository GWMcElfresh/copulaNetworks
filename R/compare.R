#' Extract edge list from a copula fit result
#'
#' @param result Output of [fit_stratum_copula()].
#' @param label Stratum label.
#' @return Data frame with columns `from`, `to`, `pcor`, `stratum`, `edge`.
#' @export
extract_edges <- function(result, label) {
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
      var_i = character(0),
      var_j = character(0),
      edge = character(0),
      value = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  out <- lapply(pairs, function(pair) {
    i <- pair[1]
    j <- pair[2]
    val <- mat[i, j]
    if (edge_only) {
      if (is.null(adj) || adj[i, j] == 0) {
        val <- 0
      }
    }
    data.frame(
      var_i = i,
      var_j = j,
      edge = paste(pmin(i, j), pmax(i, j), sep = " -- "),
      value = val,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}

#' Compare pairwise matrix values across two strata
#'
#' @param res_a First stratum fit result.
#' @param res_b Second stratum fit result.
#' @param label_a Label for stratum A.
#' @param label_b Label for stratum B.
#' @param matrix_name One of `"pcor"` or `"copula_cor"`.
#' @param edge_only For `"pcor"`, if `TRUE` only adjacency edges are non-zero.
#'   For `"copula_cor"`, defaults to `FALSE` (all pairs).
#' @param delta_threshold Threshold for direction classification.
#' @return Data frame with comparison columns.
#' @export
compare_pairwise_matrices <- function(res_a,
                                      res_b,
                                      label_a,
                                      label_b,
                                      matrix_name = c("pcor", "copula_cor"),
                                      edge_only = NULL,
                                      delta_threshold = 0.05) {
  matrix_name <- match.arg(matrix_name)

  if (is.null(res_a) || is.null(res_b)) {
    stop("Both fit results must be non-NULL.", call. = FALSE)
  }

  mat_a <- res_a[[matrix_name]]
  mat_b <- res_b[[matrix_name]]
  shared <- intersect(colnames(mat_a), colnames(mat_b))

  if (length(shared) < 2) {
    stop("Fewer than 2 shared variables between strata.", call. = FALSE)
  }

  only_a <- setdiff(colnames(mat_a), colnames(mat_b))
  only_b <- setdiff(colnames(mat_b), colnames(mat_a))
  if (length(only_a) > 0 || length(only_b) > 0) {
    warning(
      "Variable mismatch — using intersection (n = ", length(shared), "). ",
      "Only in A: ", paste(only_a, collapse = ", "),
      "; only in B: ", paste(only_b, collapse = ", "),
      call. = FALSE
    )
  }

  mat_a <- mat_a[shared, shared, drop = FALSE]
  mat_b <- mat_b[shared, shared, drop = FALSE]

  if (is.null(edge_only)) {
    edge_only <- matrix_name == "pcor"
  }

  adj_a <- if (edge_only && matrix_name == "pcor") res_a$adjacency[shared, shared, drop = FALSE] else NULL
  adj_b <- if (edge_only && matrix_name == "pcor") res_b$adjacency[shared, shared, drop = FALSE] else NULL

  pairs_a <- extract_pairwise_values(mat_a, shared, edge_only = edge_only, adj = adj_a)
  pairs_b <- extract_pairwise_values(mat_b, shared, edge_only = edge_only, adj = adj_b)

  cmp <- dplyr::full_join(
    pairs_a %>% dplyr::select(edge, value_a = value),
    pairs_b %>% dplyr::select(edge, value_b = value),
    by = "edge"
  ) %>%
    tidyr::replace_na(list(value_a = 0, value_b = 0)) %>%
    dplyr::mutate(
      delta = value_a - value_b,
      abs_delta = abs(delta),
      direction = dplyr::case_when(
        abs(value_a) > 0 & abs(value_b) == 0 ~ paste("Only in", label_a),
        abs(value_a) == 0 & abs(value_b) > 0 ~ paste("Only in", label_b),
        delta > delta_threshold ~ paste("Stronger in", label_a),
        delta < -delta_threshold ~ paste("Stronger in", label_b),
        TRUE ~ "Similar"
      ),
      matrix = matrix_name,
      label_a = label_a,
      label_b = label_b
    ) %>%
    dplyr::arrange(dplyr::desc(abs_delta))

  cmp
}

#' Compare two fitted strata across matrix types
#'
#' @param res_a First stratum fit result.
#' @param res_b Second stratum fit result.
#' @param label_a Label for stratum A.
#' @param label_b Label for stratum B.
#' @param matrices Character vector of matrices to compare (`"pcor"`, `"copula_cor"`).
#' @param delta_threshold Threshold for direction labels and bar charts.
#' @param edge_only_pcor If `TRUE`, pcor comparison uses adjacency edges only.
#' @param edge_only_cor If `FALSE` (default), copula cor compares all pairs.
#' @return List with comparison data frames keyed by matrix name.
#' @export
compare_two_strata <- function(res_a,
                               res_b,
                               label_a = "A",
                               label_b = "B",
                               matrices = c("pcor", "copula_cor"),
                               delta_threshold = 0.05,
                               edge_only_pcor = TRUE,
                               edge_only_cor = FALSE) {
  out <- list()
  for (mat in matrices) {
    edge_only <- if (mat == "pcor") edge_only_pcor else edge_only_cor
    out[[mat]] <- compare_pairwise_matrices(
      res_a, res_b, label_a, label_b,
      matrix_name = mat,
      edge_only = edge_only,
      delta_threshold = delta_threshold
    )
  }
  structure(out, class = "copula_stratum_comparison")
}
