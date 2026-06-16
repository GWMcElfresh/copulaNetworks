#' Build named strata from a declarative recipe
#'
#' @param data Input data frame.
#' @param spec Named list with optional elements:
#'   \describe{
#'     \item{`mutate`}{Quoted expression evaluated with [rlang::eval_tidy()].}
#'     \item{`filter`}{Quoted expression for row filtering.}
#'     \item{`group_by`}{Character vector of grouping columns, or a single column name.}
#'     \item{`stratum_col`}{Pre-built stratum column (skips `group_by` split).}
#'     \item{`name_sep`}{Separator when joining multiple `group_by` columns (default `" | "`).}
#'     \item{`min_n`}{Minimum rows per stratum (default 1).}
#'   }
#' @param spec_name Name prefix for strata (used in messages).
#' @return Named list of stratum data frames.
#' @export
build_strata <- function(data, spec, spec_name = "strata") {
  df <- data

  if (!is.null(spec$mutate)) {
    df <- rlang::eval_tidy(
      rlang::expr(dplyr::mutate(df, !!spec$mutate)),
      data = list(df = df)
    )
  }

  if (!is.null(spec$filter)) {
    keep <- rlang::eval_tidy(spec$filter, data = df)
    df <- df[keep, , drop = FALSE]
  }

  min_n <- if (is.null(spec$min_n)) 1L else as.integer(spec$min_n)
  name_sep <- if (is.null(spec$name_sep)) " | " else spec$name_sep

  if (!is.null(spec$stratum_col)) {
    stratum_col <- spec$stratum_col
    if (!stratum_col %in% colnames(df)) {
      stop("stratum_col '", stratum_col, "' not found in data.", call. = FALSE)
    }
    levels <- sort(unique(as.character(df[[stratum_col]])))
    strata <- stats::setNames(vector("list", length(levels)), levels)
    for (lv in levels) {
      sub <- df[df[[stratum_col]] == lv, , drop = FALSE]
      if (nrow(sub) >= min_n) {
        strata[[lv]] <- sub
      } else {
        warning("Dropping stratum '", lv, "' in ", spec_name, " — n = ", nrow(sub),
                " < min_n = ", min_n, ".", call. = FALSE)
      }
    }
    strata <- strata[!vapply(strata, is.null, logical(1))]
    return(strata)
  }

  if (is.null(spec$group_by)) {
    stop("strata spec must include either 'stratum_col' or 'group_by'.", call. = FALSE)
  }

  group_cols <- spec$group_by
  missing <- setdiff(group_cols, colnames(df))
  if (length(missing) > 0) {
    stop("group_by columns not found: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  if (length(group_cols) == 1) {
    stratum_names <- as.character(df[[group_cols[1]]])
  } else {
    stratum_names <- apply(df[, group_cols, drop = FALSE], 1, function(row) {
      paste(row, collapse = name_sep)
    })
  }

  splits <- split(df, stratum_names)
  strata <- splits

  for (nm in names(strata)) {
    strata[[nm]] <- strata[[nm]][, setdiff(colnames(strata[[nm]]), ".stratum_name"), drop = FALSE]
    if (nrow(strata[[nm]]) < min_n) {
      warning("Dropping stratum '", nm, "' in ", spec_name, " — n = ", nrow(strata[[nm]]),
              " < min_n = ", min_n, ".", call. = FALSE)
      strata[[nm]] <- NULL
    }
  }

  strata <- strata[!vapply(strata, is.null, logical(1))]
  strata
}
