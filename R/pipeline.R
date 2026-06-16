#' Run the full stratified copula pipeline
#'
#' Convenience wrapper that runs prepare, fit, plot, and optional comparisons.
#' Each step writes RDS checkpoints when `outDir` is provided.
#'
#' @param data Clean input data frame.
#' @param idCols Identifier columns.
#' @param strataCols Exogenous stratification columns.
#' @param nodeCols Node columns (NULL = auto-detect numeric).
#' @param excludeCols Columns to exclude from modeling.
#' @param strataSpecs Named list of stratum recipes for [PrepareCopulaData()].
#' @param comparePairs Optional list of length-2 character vectors naming strata to compare.
#' @param outDir Output directory for checkpoints and figures.
#' @param method Lambda selection method.
#' @param nlambda Number of lambda values.
#' @param starsThresh StARS threshold.
#' @param minN Minimum observations per stratum.
#' @param includeFull Fit unstratified baseline.
#' @param plotDiagnostics If `TRUE`, save per-stratum diagnostic plots.
#' @param nodeGroups Optional node group mapping for plots.
#' @param width Save width in inches for stratum diagnostic plots.
#' @param height Save height in inches for stratum diagnostic plots.
#' @param dpi Save resolution for stratum diagnostic PNG files.
#' @param comparisonWidth Save width in inches for comparison plots.
#' @param comparisonHeight Save height in inches for comparison plots.
#' @param comparisonDpi Save resolution for comparison PNG files.
#' @param ... Additional arguments passed to [PlotStratumDiagnostics()].
#' @return List with `preparedData`, `fitResults`, `plotArtifacts`, and `comparisons`.
#' @export
RunCopulaPipeline <- function(data,
                              idCols = character(0),
                              strataCols = character(0),
                              nodeCols = NULL,
                              excludeCols = character(0),
                              strataSpecs,
                              comparePairs = list(),
                              outDir = "checkpoints/copula_run",
                              method = c("stars", "ebic"),
                              nlambda = 40,
                              starsThresh = 0.1,
                              minN = 10,
                              includeFull = FALSE,
                              plotDiagnostics = TRUE,
                              nodeGroups = NULL,
                              width = 10,
                              height = 10,
                              dpi = 150,
                              comparisonWidth = 16,
                              comparisonHeight = 11,
                              comparisonDpi = 150,
                              ...) {
  method <- match.arg(method)

  prepared_data <- PrepareCopulaData(
    data = data,
    idCols = idCols,
    strataCols = strataCols,
    nodeCols = nodeCols,
    excludeCols = excludeCols,
    strataSpecs = strataSpecs,
    outDir = outDir
  )

  fit_results <- FitCopulaStrata(
    prep = prepared_data,
    method = method,
    nlambda = nlambda,
    starsThresh = starsThresh,
    minN = minN,
    includeFull = includeFull,
    outDir = outDir
  )

  plot_artifacts <- NULL
  if (isTRUE(plotDiagnostics)) {
    fig_dir <- file.path(outDir, "figures")
    plot_artifacts <- PlotAllStrata(
      fit_results,
      outDir = fig_dir,
      nodeGroups = nodeGroups,
      width = width,
      height = height,
      dpi = dpi,
      ...
    )
  }

  comparisons <- list()
  if (length(comparePairs) > 0) {
    cmp_dir <- file.path(outDir, "comparisons")
    for (pair in comparePairs) {
      if (length(pair) != 2) {
        warning("comparePairs entries must be length-2 vectors; skipping.", call. = FALSE)
        next
      }
      nm_a <- pair[1]
      nm_b <- pair[2]
      if (!all(c(nm_a, nm_b) %in% names(fit_results$fits))) {
        warning(
          "Comparison pair not found in fits: ", nm_a, " vs ", nm_b, " — skipping.",
          call. = FALSE
        )
        next
      }
      cmp <- CompareTwoStrata(
        fit_results$fits[[nm_a]],
        fit_results$fits[[nm_b]],
        labelA = nm_a,
        labelB = nm_b
      )
      pair_dir <- file.path(cmp_dir, paste(nm_a, "vs", nm_b, sep = "_"))
      PlotStratumComparison(
        cmp,
        outDir = pair_dir,
        width = comparisonWidth,
        height = comparisonHeight,
        dpi = comparisonDpi
      )
      comparisons[[paste(nm_a, "vs", nm_b, sep = "_")]] <- cmp
    }
    if (length(comparisons) > 0 && !is.null(outDir)) {
      SaveCheckpoint(comparisons, outDir, filename = "comparisons.rds")
    }
  }

  result <- list(
    preparedData = prepared_data,
    fitResults = fit_results,
    plotArtifacts = plot_artifacts,
    comparisons = comparisons
  )

  if (!is.null(outDir)) {
    SaveCheckpoint(result, outDir, filename = "pipeline_result.rds")
  }

  result
}
