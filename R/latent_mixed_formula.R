#' Validate and parse a RHS-only latent mixed-model formula
#'
#' Supports a single grouping factor with common random-intercept / random-slope
#' terms: `(1 | g)`, `(0 + x | g)`, `(1 + x | g)`, `(x | g)`.
#'
#' @param latentFormula RHS-only formula for the latent Gaussian layer.
#' @param data Data frame used to resolve variable names.
#' @return List with `fixedFormula`, `groupVar`, `reTerms`, `hasIntercept`,
#'   `slopeVars`, and `nRe`.
#' @keywords internal
#' @noRd
parse_latent_formula <- function(latentFormula, data) {
  if (!inherits(latentFormula, "formula")) {
    stop("latentFormula must be a formula, e.g. ~ time + (1 | subject).", call. = FALSE)
  }
  if (length(latentFormula) != 2L) {
    stop(
      "latentFormula must be RHS-only (no response), e.g. ~ time + (1 | subject). ",
      "Each node in nodeCols is the latent response.",
      call. = FALSE
    )
  }

  formula_chr <- paste(deparse(latentFormula), collapse = "")
  bar_locs <- gregexpr("|", formula_chr, fixed = TRUE)[[1]]
  bar_count <- if (identical(bar_locs, -1L)) 0L else length(bar_locs)
  if (bar_count < 1L) {
    stop(
      "latentFormula must include a random-effect term with '|', e.g. (1 | subject).",
      call. = FALSE
    )
  }
  if (bar_count > 1L) {
    stop(
      "v1 supports a single grouping factor only (one '|' term). ",
      "Nested or crossed random effects are not supported.",
      call. = FALSE
    )
  }

  # Extract ( ... | group ) via terms attributes when available; fallback to regex
  re_match <- regexec("\\(([^|()]+)\\|([^)]+)\\)", formula_chr, perl = TRUE)
  re_cap <- regmatches(formula_chr, re_match)[[1]]
  if (length(re_cap) < 3L) {
    stop("Could not parse random-effect term in latentFormula.", call. = FALSE)
  }
  re_lhs <- trimws(re_cap[2])
  group_var <- trimws(re_cap[3])
  if (!group_var %in% names(data)) {
    stop("Grouping variable '", group_var, "' not found in data.", call. = FALSE)
  }

  # Fixed formula: drop the RE parenthetical
  fixed_chr <- gsub("\\([^|()]+\\|[^)]+\\)", "", formula_chr)
  fixed_chr <- gsub("~\\s*\\+\\s*", "~ ", fixed_chr)
  fixed_chr <- gsub("\\+\\s*$", "", fixed_chr)
  fixed_chr <- gsub("~\\s*$", "~ 1", fixed_chr)
  fixed_chr <- trimws(fixed_chr)
  if (!grepl("^~", fixed_chr)) {
    fixed_chr <- paste0("~ ", fixed_chr)
  }
  fixed_formula <- stats::as.formula(fixed_chr)

  # Parse RE LHS: intercept / slopes
  has_intercept <- TRUE
  slope_vars <- character(0)
  re_lhs_clean <- gsub("\\s+", " ", re_lhs)
  if (grepl("^0\\s*\\+", re_lhs_clean) || grepl("^0\\s*$", re_lhs_clean)) {
    has_intercept <- FALSE
    re_lhs_clean <- sub("^0\\s*\\+\\s*", "", re_lhs_clean)
    re_lhs_clean <- sub("^0\\s*$", "", re_lhs_clean)
  }
  if (nzchar(re_lhs_clean) && re_lhs_clean != "1") {
    parts <- strsplit(re_lhs_clean, "\\+", perl = TRUE)[[1]]
    parts <- trimws(parts)
    parts <- parts[nzchar(parts) & parts != "1"]
    slope_vars <- parts
  }
  for (sv in slope_vars) {
    if (!sv %in% names(data)) {
      stop("Random-slope variable '", sv, "' not found in data.", call. = FALSE)
    }
  }

  n_re <- as.integer(has_intercept) + length(slope_vars)
  if (n_re < 1L) {
    stop("Random-effect term must include an intercept and/or at least one slope.", call. = FALSE)
  }

  list(
    fixedFormula = fixed_formula,
    groupVar = group_var,
    reLhs = re_lhs,
    hasIntercept = has_intercept,
    slopeVars = slope_vars,
    nRe = n_re,
    latentFormula = latentFormula
  )
}

#' Build fixed and random design matrices for the latent mixed model
#'
#' @param parsed Output of [parse_latent_formula()].
#' @param data Data frame.
#' @return List with `X`, `Z` (n x nRe), `groupId` (integer 1..G), `groupLevels`,
#'   `groupSizes`, `fixedNames`, `reNames`.
#' @keywords internal
#' @noRd
build_latent_design_matrices <- function(parsed, data) {
  n_obs <- nrow(data)
  X <- stats::model.matrix(parsed$fixedFormula, data = data)
  if (nrow(X) != n_obs) {
    stop("model.matrix dropped rows (check for NA in covariates).", call. = FALSE)
  }

  group_factor <- factor(data[[parsed$groupVar]])
  group_id <- as.integer(group_factor)
  group_levels <- levels(group_factor)
  n_groups <- length(group_levels)
  group_sizes <- as.integer(table(group_factor))

  Z <- matrix(0, nrow = n_obs, ncol = parsed$nRe)
  re_names <- character(parsed$nRe)
  col_idx <- 1L
  if (isTRUE(parsed$hasIntercept)) {
    Z[, col_idx] <- 1
    re_names[col_idx] <- "(Intercept)"
    col_idx <- col_idx + 1L
  }
  for (sv in parsed$slopeVars) {
    Z[, col_idx] <- as.numeric(data[[sv]])
    re_names[col_idx] <- sv
    col_idx <- col_idx + 1L
  }
  colnames(Z) <- re_names

  list(
    X = X,
    Z = Z,
    groupId = group_id,
    groupLevels = group_levels,
    nGroups = n_groups,
    groupSizes = group_sizes,
    fixedNames = colnames(X),
    reNames = re_names,
    nRe = parsed$nRe,
    groupVar = parsed$groupVar
  )
}

#' Expand group-level random effects to observation level
#'
#' @param b_mat G x q matrix of group random effects.
#' @param group_id Integer group index per observation.
#' @param Z n x q random design matrix.
#' @return Length-n vector of Z_i' b_{g(i)}.
#' @keywords internal
#' @noRd
expand_random_effects <- function(b_mat, group_id, Z) {
  n_obs <- nrow(Z)
  out <- numeric(n_obs)
  for (i in seq_len(n_obs)) {
    out[i] <- sum(Z[i, ] * b_mat[group_id[i], ])
  }
  out
}
