#' Resolve node columns from data and role specifications
#'
#' @param data Input data frame.
#' @param nodeCols Explicit node columns, or NULL to auto-detect numeric columns.
#' @param idCols ID columns to exclude.
#' @param strataCols Stratification columns to exclude.
#' @param excludeCols Additional columns to exclude.
#' @return Character vector of node column names.
#' @keywords internal
resolve_node_cols <- function(data, nodeCols, idCols, strataCols, excludeCols) {
  reserved <- unique(c(idCols, strataCols, excludeCols))
  missing_reserved <- setdiff(reserved, colnames(data))
  if (length(missing_reserved) > 0) {
    stop("Columns not found in data: ", paste(missing_reserved, collapse = ", "), call. = FALSE)
  }

  if (is.null(nodeCols)) {
    nodeCols <- colnames(data)[vapply(data, is.numeric, logical(1))]
    nodeCols <- setdiff(nodeCols, reserved)
  } else {
    missing_nodes <- setdiff(nodeCols, colnames(data))
    if (length(missing_nodes) > 0) {
      stop("nodeCols not found in data: ", paste(missing_nodes, collapse = ", "), call. = FALSE)
    }
  }

  if (length(nodeCols) < 3) {
    stop("At least 3 node columns are required; found ", length(nodeCols), ".", call. = FALSE)
  }

  all_na <- vapply(nodeCols, function(col) all(is.na(data[[col]])), logical(1))
  if (any(all_na)) {
    stop("nodeCols with all NA values: ", paste(nodeCols[all_na], collapse = ", "), call. = FALSE)
  }

  nodeCols
}

#' Prepare data and strata for copula modeling (Step 0)
#'
#' Validates column roles, resolves node variables, and builds named strata from
#' declarative recipes. User supplies a clean data frame; no imputation is performed.
#'
#' @param data Clean input data frame.
#' @param idCols Character vector of identifier columns (excluded from model).
#' @param strataCols Character vector of exogenous stratification columns (excluded from model).
#' @param nodeCols Character vector of endogenous node columns. If `NULL`, all numeric
#'   columns not in `idCols`, `strataCols`, or `excludeCols` are used.
#' @param excludeCols Additional columns to exclude from the model.
#' @param strataSpecs Named list of stratum recipes passed to [BuildStrata()]. Each
#'   recipe may include `mutate`, `filter`, `group_by`, `stratumCol`, `nameSep`, `minN`.
#'   Keys become prefixes: stratum names are `"<spec_key>::<stratum_label>"`.
#' @param outDir Optional directory to save `prep.rds`.
#' @return List with elements `data`, `idCols`, `strataCols`, `nodeCols`, `strata`,
#'   and `meta` (per-stratum row counts).
#' @export
PrepareCopulaData <- function(data,
                              idCols = character(0),
                              strataCols = character(0),
                              nodeCols = NULL,
                              excludeCols = character(0),
                              strataSpecs = list(all = list()),
                              outDir = NULL) {
  if (!is.data.frame(data)) {
    stop("data must be a data.frame.", call. = FALSE)
  }

  nodeCols <- resolve_node_cols(data, nodeCols, idCols, strataCols, excludeCols)

  all_strata <- list()
  meta <- list()

  if (length(strataSpecs) == 0) {
    strataSpecs <- list(all = list())
  }

  for (spec_name in names(strataSpecs)) {
    spec <- strataSpecs[[spec_name]]
    if (is.null(spec$stratumCol) && is.null(spec$group_by) && spec_name == "all" && length(spec) == 0) {
      all_strata[["all"]] <- data
      meta[["all"]] <- list(n = nrow(data), strataSpec = "all")
      next
    }

    built_strata <- BuildStrata(data, spec, specName = spec_name)
    for (stratum_name in names(built_strata)) {
      key <- if (spec_name == stratum_name || grepl(paste0("^", spec_name, "::"), stratum_name)) {
        stratum_name
      } else {
        paste(spec_name, stratum_name, sep = "::")
      }
      all_strata[[key]] <- built_strata[[stratum_name]]
      meta[[key]] <- list(n = nrow(built_strata[[stratum_name]]), strataSpec = spec_name)
    }
  }

  if (length(all_strata) == 0) {
    stop("No strata were built. Check strataSpecs and minN thresholds.", call. = FALSE)
  }

  prepared_data <- list(
    data = data,
    idCols = idCols,
    strataCols = strataCols,
    nodeCols = nodeCols,
    strata = all_strata,
    meta = meta
  )

  if (!is.null(outDir)) {
    SaveCheckpoint(prepared_data, outDir, filename = "prep.rds")
  }

  prepared_data
}
