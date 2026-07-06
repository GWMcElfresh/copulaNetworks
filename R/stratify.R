#' Build named strata from a declarative recipe
#'
#' @param data Input data frame.
#' @param spec Named list with optional elements:
#'   \describe{
#'     \item{`mutate`}{Quoted expression evaluated with [rlang::eval_tidy()].}
#'     \item{`filter`}{Quoted expression for row filtering.}
#'     \item{`group_by`}{Character vector of grouping columns, or a single column name.}
#'     \item{`stratumCol`}{Pre-built stratum column (skips `group_by` split).}
#'     \item{`nameSep`}{Separator when joining multiple `group_by` columns (default `" | "`).}
#'     \item{`minN`}{Minimum rows per stratum (default 1).}
#'   }
#' @param specName Name prefix for strata (used in messages).
#' @return Named list of stratum data frames.
#' @export
BuildStrata <- function(data, spec, specName = "strata") {
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

  min_n <- if (is.null(spec$minN)) 1L else as.integer(spec$minN)
  name_sep <- if (is.null(spec$nameSep)) " | " else spec$nameSep

  if (!is.null(spec$stratumCol)) {
    stratum_col <- spec$stratumCol
    if (!stratum_col %in% colnames(df)) {
      stop("stratumCol '", stratum_col, "' not found in data.", call. = FALSE)
    }
    levels <- sort(unique(as.character(df[[stratum_col]])))
    strata <- stats::setNames(vector("list", length(levels)), levels)
    for (lv in levels) {
      sub <- df[df[[stratum_col]] == lv, , drop = FALSE]
      if (nrow(sub) >= min_n) {
        strata[[lv]] <- sub
      } else {
        warning("Dropping stratum '", lv, "' in ", specName, " - n = ", nrow(sub),
                " < minN = ", min_n, ".", call. = FALSE)
      }
    }
    strata <- strata[!vapply(strata, is.null, logical(1))]
    return(strata)
  }

  if (is.null(spec$group_by)) {
    stop("strata spec must include either 'stratumCol' or 'group_by'.", call. = FALSE)
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

  for (stratum_name in names(strata)) {
    strata[[stratum_name]] <- strata[[stratum_name]][, setdiff(colnames(strata[[stratum_name]]), ".stratum_name"), drop = FALSE]
    if (nrow(strata[[stratum_name]]) < min_n) {
      warning("Dropping stratum '", stratum_name, "' in ", specName, " - n = ", nrow(strata[[stratum_name]]),
              " < minN = ", min_n, ".", call. = FALSE)
      strata[[stratum_name]] <- NULL
    }
  }

  strata <- strata[!vapply(strata, is.null, logical(1))]
  strata
}
