#' Resolve node columns from data and role specifications
#'
#' @param data Input data frame.
#' @param node_cols Explicit node columns, or NULL to auto-detect numeric columns.
#' @param id_cols ID columns to exclude.
#' @param strata_cols Stratification columns to exclude.
#' @param exclude_cols Additional columns to exclude.
#' @return Character vector of node column names.
#' @keywords internal
resolve_node_cols <- function(data, node_cols, id_cols, strata_cols, exclude_cols) {
  reserved <- unique(c(id_cols, strata_cols, exclude_cols))
  missing_reserved <- setdiff(reserved, colnames(data))
  if (length(missing_reserved) > 0) {
    stop("Columns not found in data: ", paste(missing_reserved, collapse = ", "), call. = FALSE)
  }

  if (is.null(node_cols)) {
    node_cols <- colnames(data)[vapply(data, is.numeric, logical(1))]
    node_cols <- setdiff(node_cols, reserved)
  } else {
    missing_nodes <- setdiff(node_cols, colnames(data))
    if (length(missing_nodes) > 0) {
      stop("node_cols not found in data: ", paste(missing_nodes, collapse = ", "), call. = FALSE)
    }
  }

  if (length(node_cols) < 3) {
    stop("At least 3 node columns are required; found ", length(node_cols), ".", call. = FALSE)
  }

  all_na <- vapply(node_cols, function(col) all(is.na(data[[col]])), logical(1))
  if (any(all_na)) {
    stop("node_cols with all NA values: ", paste(node_cols[all_na], collapse = ", "), call. = FALSE)
  }

  node_cols
}

#' Prepare data and strata for copula modeling (Step 0)
#'
#' Validates column roles, resolves node variables, and builds named strata from
#' declarative recipes. User supplies a clean data frame; no imputation is performed.
#'
#' @param data Clean input data frame.
#' @param id_cols Character vector of identifier columns (excluded from model).
#' @param strata_cols Character vector of exogenous stratification columns (excluded from model).
#' @param node_cols Character vector of endogenous node columns. If `NULL`, all numeric
#'   columns not in `id_cols`, `strata_cols`, or `exclude_cols` are used.
#' @param exclude_cols Additional columns to exclude from the model.
#' @param strata_specs Named list of stratum recipes passed to [build_strata()]. Each
#'   recipe may include `mutate`, `filter`, `group_by`, `stratum_col`, `name_sep`, `min_n`.
#'   Keys become prefixes: stratum names are `"<spec_key>::<stratum_label>"`.
#' @param out_dir Optional directory to save `prep.rds`.
#' @return List with elements `data`, `id_cols`, `strata_cols`, `node_cols`, `strata`,
#'   and `meta` (per-stratum row counts).
#' @export
prepare_copula_data <- function(data,
                                id_cols = character(0),
                                strata_cols = character(0),
                                node_cols = NULL,
                                exclude_cols = character(0),
                                strata_specs = list(all = list()),
                                out_dir = NULL) {
  if (!is.data.frame(data)) {
    stop("data must be a data.frame.", call. = FALSE)
  }

  node_cols <- resolve_node_cols(data, node_cols, id_cols, strata_cols, exclude_cols)

  all_strata <- list()
  meta <- list()

  if (length(strata_specs) == 0) {
    strata_specs <- list(all = list())
  }

  for (spec_name in names(strata_specs)) {
    spec <- strata_specs[[spec_name]]
    if (is.null(spec$stratum_col) && is.null(spec$group_by) && spec_name == "all" && length(spec) == 0) {
      all_strata[["all"]] <- data
      meta[["all"]] <- list(n = nrow(data), spec = "all")
      next
    }

    built <- build_strata(data, spec, spec_name = spec_name)
    for (nm in names(built)) {
      key <- if (spec_name == nm || grepl(paste0("^", spec_name, "::"), nm)) {
        nm
      } else {
        paste(spec_name, nm, sep = "::")
      }
      all_strata[[key]] <- built[[nm]]
      meta[[key]] <- list(n = nrow(built[[nm]]), spec = spec_name)
    }
  }

  if (length(all_strata) == 0) {
    stop("No strata were built. Check strata_specs and min_n thresholds.", call. = FALSE)
  }

  prep <- list(
    data = data,
    id_cols = id_cols,
    strata_cols = strata_cols,
    node_cols = node_cols,
    strata = all_strata,
    meta = meta
  )

  if (!is.null(out_dir)) {
    save_checkpoint(prep, out_dir, filename = "prep.rds")
  }

  prep
}
