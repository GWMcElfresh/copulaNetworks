#' Build a graphical-fit-like object from a latent mixed-model residual matrix
#'
#' Convenience wrapper around [PseudoGraphicalFitFromCor()] / optional glasso
#' so vignettes can reuse [PlotCopulaCorHeatmap()], [PlotPcorHeatmap()], and
#' [PlotCopulaNetwork()] exactly as in the factor-vine and meta-analysis articles.
#'
#' @param fit Output of [FitCopulaLatentMixedModel()].
#' @param useGlasso If `TRUE` and `fit$residualGraphical` has edges, return
#'   that glasso fit; otherwise correlation-based pseudo-fit.
#' @return List compatible with existing `Plot*` helpers.
#' @export
PseudoGraphicalFitFromLatentMixed <- function(fit, useGlasso = TRUE) {
  if (isTRUE(useGlasso) && !is.null(fit$residualGraphical)) {
    adj <- fit$residualGraphical$adjacency
    if (!is.null(adj) && sum(adj) > 0) {
      return(fit$residualGraphical)
    }
  }
  if (!is.null(fit$impliedCor)) {
    return(PseudoGraphicalFitFromCor(fit$impliedCor, nObs = nrow(fit$residualMatrix)))
  }
  cor_mat <- stats::cor(fit$residualMatrix, use = "pairwise.complete.obs")
  PseudoGraphicalFitFromCor(cor_mat, nObs = nrow(fit$residualMatrix))
}

#' Plot residual copula correlation heatmap from a latent mixed fit
#'
#' Thin wrapper on [PlotCopulaCorHeatmap()] for vignette parity with other articles.
#'
#' @param fit Output of [FitCopulaLatentMixedModel()].
#' @param title Plot title.
#' @param ... Passed to [PlotCopulaCorHeatmap()].
#' @return pheatmap object (draw with grid for knitr capture).
#' @export
PlotLatentMixedResidualCorHeatmap <- function(fit,
                                              title = "Residual copula correlation",
                                              ...) {
  graphical <- PseudoGraphicalFitFromLatentMixed(fit, useGlasso = TRUE)
  PlotCopulaCorHeatmap(graphical, title = title, ...)
}

#' Plot residual partial-correlation heatmap from a latent mixed fit
#'
#' @inheritParams PlotLatentMixedResidualCorHeatmap
#' @return pheatmap object.
#' @export
PlotLatentMixedResidualPcorHeatmap <- function(fit,
                                               title = "Residual partial correlation",
                                               ...) {
  graphical <- PseudoGraphicalFitFromLatentMixed(fit, useGlasso = TRUE)
  PlotPcorHeatmap(graphical, title = title, ...)
}

#' Plot residual graphical network from a latent mixed fit
#'
#' @inheritParams PlotLatentMixedResidualCorHeatmap
#' @param printPlot Passed to [PlotCopulaNetwork()].
#' @param ... Passed to [PlotCopulaNetwork()].
#' @return ggraph / ggplot object.
#' @export
PlotLatentMixedResidualNetwork <- function(fit,
                                           title = "Residual copula network",
                                           printPlot = TRUE,
                                           ...) {
  graphical <- PseudoGraphicalFitFromLatentMixed(fit, useGlasso = TRUE)
  PlotCopulaNetwork(graphical, title = title, printPlot = printPlot, ...)
}
