#' Run the full stratified copula pipeline
#'
#' Convenience wrapper that runs prepare, fit, plot, and optional comparisons.
#' Each step writes RDS checkpoints when `out_dir` is provided.
#'
#' @param data Clean input data frame.
#' @param id_cols Identifier columns.
#' @param strata_cols Exogenous stratification columns.
#' @param node_cols Node columns (NULL = auto-detect numeric).
#' @param exclude_cols Columns to exclude from modeling.
#' @param strata_specs Named list of stratum recipes for [prepare_copula_data()].
#' @param compare_pairs Optional list of length-2 character vectors naming strata to compare.
#' @param out_dir Output directory for checkpoints and figures.
#' @param method Lambda selection method.
#' @param nlambda Number of lambda values.
#' @param stars_thresh StARS threshold.
#' @param min_n Minimum observations per stratum.
#' @param include_full Fit unstratified baseline.
#' @param plot_diagnostics If `TRUE`, save per-stratum diagnostic plots.
#' @param node_groups Optional node group mapping for plots.
#' @param width Save width in inches for stratum diagnostic plots.
#' @param height Save height in inches for stratum diagnostic plots.
#' @param dpi Save resolution for stratum diagnostic PNG files.
#' @param comparison_width Save width in inches for comparison plots.
#' @param comparison_height Save height in inches for comparison plots.
#' @param comparison_dpi Save resolution for comparison PNG files.
#' @param ... Additional arguments passed to [plot_stratum_diagnostics()].
#' @return List with `prep`, `fits`, `plots`, and `comparisons`.
#' @export
run_copula_pipeline <- function(data,
                                id_cols = character(0),
                                strata_cols = character(0),
                                node_cols = NULL,
                                exclude_cols = character(0),
                                strata_specs,
                                compare_pairs = list(),
                                out_dir = "checkpoints/copula_run",
                                method = c("stars", "ebic"),
                                nlambda = 40,
                                stars_thresh = 0.1,
                                min_n = 10,
                                include_full = FALSE,
                                plot_diagnostics = TRUE,
                                node_groups = NULL,
                                width = 10,
                                height = 10,
                                dpi = 150,
                                comparison_width = 16,
                                comparison_height = 11,
                                comparison_dpi = 150,
                                ...) {
  method <- match.arg(method)

  prep <- prepare_copula_data(
    data = data,
    id_cols = id_cols,
    strata_cols = strata_cols,
    node_cols = node_cols,
    exclude_cols = exclude_cols,
    strata_specs = strata_specs,
    out_dir = out_dir
  )

  fits <- fit_copula_strata(
    prep = prep,
    method = method,
    nlambda = nlambda,
    stars_thresh = stars_thresh,
    min_n = min_n,
    include_full = include_full,
    out_dir = out_dir
  )

  plots <- NULL
  if (isTRUE(plot_diagnostics)) {
    fig_dir <- file.path(out_dir, "figures")
    plots <- plot_all_strata(
      fits,
      out_dir = fig_dir,
      node_groups = node_groups,
      width = width,
      height = height,
      dpi = dpi,
      ...
    )
  }

  comparisons <- list()
  if (length(compare_pairs) > 0) {
    cmp_dir <- file.path(out_dir, "comparisons")
    for (pair in compare_pairs) {
      if (length(pair) != 2) {
        warning("compare_pairs entries must be length-2 vectors; skipping.", call. = FALSE)
        next
      }
      nm_a <- pair[1]
      nm_b <- pair[2]
      if (!all(c(nm_a, nm_b) %in% names(fits$fits))) {
        warning(
          "Comparison pair not found in fits: ", nm_a, " vs ", nm_b, " — skipping.",
          call. = FALSE
        )
        next
      }
      cmp <- compare_two_strata(
        fits$fits[[nm_a]],
        fits$fits[[nm_b]],
        label_a = nm_a,
        label_b = nm_b
      )
      pair_dir <- file.path(cmp_dir, paste(nm_a, "vs", nm_b, sep = "_"))
      plot_stratum_comparison(
        cmp,
        out_dir = pair_dir,
        width = comparison_width,
        height = comparison_height,
        dpi = comparison_dpi
      )
      comparisons[[paste(nm_a, "vs", nm_b, sep = "_")]] <- cmp
    }
    if (length(comparisons) > 0 && !is.null(out_dir)) {
      save_checkpoint(comparisons, out_dir, filename = "comparisons.rds")
    }
  }

  result <- list(
    prep = prep,
    fits = fits,
    plots = plots,
    comparisons = comparisons
  )

  if (!is.null(out_dir)) {
    save_checkpoint(result, out_dir, filename = "pipeline_result.rds")
  }

  result
}
