#' Plot copula correlation heatmap
#'
#' @param result Output of [FitStratumCopula()].
#' @param title Plot title.
#' @param vars Optional subset of variables to display.
#' @return pheatmap grob.
#' @export
PlotCopulaCorHeatmap <- function(result, title = "Copula Correlation Matrix", vars = NULL) {
  if (is.null(result)) {
    stop("result is NULL.", call. = FALSE)
  }
  mat <- result$copulaCor
  if (!is.null(vars)) {
    vars <- intersect(vars, colnames(mat))
    mat <- mat[vars, vars, drop = FALSE]
  }
  mat[!is.finite(mat)] <- 0
  diag(mat) <- 1
  can_cluster <- nrow(mat) > 1L && all(is.finite(mat))
  pheatmap::pheatmap(
    mat,
    cluster_rows = can_cluster,
    cluster_cols = can_cluster,
    color = grDevices::colorRampPalette(c("blue", "white", "red"))(100),
    breaks = seq(-1, 1, length.out = 101),
    main = title,
    silent = TRUE
  )
}

#' Symmetric color breaks centered at 0, bounded by max |value| in mat
#' @keywords internal
symmetric_matrix_breaks <- function(mat, n = 101) {
  lim <- max(abs(mat), na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) {
    lim <- 1e-6
  }
  seq(-lim, lim, length.out = n)
}

#' Plot partial correlation heatmap
#'
#' @param result Output of [FitStratumCopula()].
#' @param title Plot title.
#' @param vars Optional subset of variables to display.
#' @param zeroDiag If `TRUE`, zero the diagonal for display.
#' @return pheatmap grob.
#' @export
PlotPcorHeatmap <- function(result,
                             title = "Partial Correlation Matrix (Glasso)",
                             vars = NULL,
                             zeroDiag = TRUE) {
  if (is.null(result)) {
    stop("result is NULL.", call. = FALSE)
  }
  mat <- result$pcor
  if (isTRUE(zeroDiag)) {
    diag(mat) <- 0
  }
  if (!is.null(vars)) {
    vars <- intersect(vars, colnames(mat))
    mat <- mat[vars, vars, drop = FALSE]
  }
  pheatmap::pheatmap(
    mat,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = grDevices::colorRampPalette(c("blue", "white", "red"))(100),
    breaks = symmetric_matrix_breaks(mat),
    main = title,
    silent = TRUE
  )
}

#' Resolve save dimension with optional override
#' @keywords internal
resolve_dim <- function(default, override = NULL) {
  if (!is.null(override)) override else default
}

#' Save a ggplot as PNG and PDF with white background
#' @keywords internal
copula_ggsave <- function(path, plot, width, height, dpi = 150) {
  base <- sub("\\.[^.]*$", "", path)
  ggplot2::ggsave(
    paste0(base, ".png"),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
  ggplot2::ggsave(
    paste0(base, ".pdf"),
    plot = plot,
    width = width,
    height = height,
    bg = "white"
  )
  invisible(base)
}

#' Save a grob to PNG and PDF via grid
#' @keywords internal
save_grob <- function(grob, path, width = 10, height = 10, dpi = 150) {
  base <- sub("\\.[^.]*$", "", path)
  grDevices::png(
    paste0(base, ".png"),
    width = width,
    height = height,
    units = "in",
    res = dpi,
    bg = "white"
  )
  grid::grid.draw(grob)
  grDevices::dev.off()
  grDevices::pdf(paste0(base, ".pdf"), width = width, height = height, bg = "white")
  grid::grid.draw(grob)
  grDevices::dev.off()
  invisible(base)
}

#' Plot single-stratum diagnostics (network + heatmaps)
#'
#' @param fitResult Output of [FitStratumCopula()].
#' @param stratumLabel Label used in plot titles and file names.
#' @param outDir Directory for saved PNG/PDF files. If `NULL`, plots are not saved.
#' @param minPcor Minimum |partial correlation| for network edges.
#' @param nodeGroups Optional node group mapping.
#' @param seed Layout seed.
#' @param width Default save width in inches (network and heatmaps).
#' @param height Default save height in inches (network and heatmaps).
#' @param networkWidth Optional network plot width override.
#' @param networkHeight Optional network plot height override.
#' @param heatmapWidth Optional heatmap width override.
#' @param heatmapHeight Optional heatmap height override.
#' @param dpi Resolution for saved PNG files.
#' @return List with ggplot/grob objects: `network`, `copulaCorHeatmap`, `pcorHeatmap`.
#' @export
PlotStratumDiagnostics <- function(fitResult,
                                   stratumLabel = "stratum",
                                   outDir = NULL,
                                   minPcor = 0.01,
                                   nodeGroups = NULL,
                                   seed = 42,
                                   width = 10,
                                   height = 10,
                                   networkWidth = NULL,
                                   networkHeight = NULL,
                                   heatmapWidth = NULL,
                                   heatmapHeight = NULL,
                                   dpi = 150) {
  if (is.null(fitResult)) {
    stop("fitResult is NULL.", call. = FALSE)
  }

  title_base <- stratumLabel
  network_plot <- PlotCopulaNetwork(
    fitResult,
    title = paste(title_base, "\u2014 Nonparanormal Copula Network"),
    seed = seed,
    minPcor = minPcor,
    nodeGroups = nodeGroups,
    printPlot = FALSE
  )

  copula_heatmap <- PlotCopulaCorHeatmap(
    fitResult,
    title = paste("Copula Correlation Matrix (", title_base, ")", sep = "")
  )
  pcor_heatmap <- PlotPcorHeatmap(
    fitResult,
    title = paste("Partial Correlation Matrix (Glasso,", title_base, ")", sep = " ")
  )

  if (!is.null(outDir)) {
    if (!dir.exists(outDir)) {
      dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
    }
    network_width <- resolve_dim(width, networkWidth)
    network_height <- resolve_dim(height, networkHeight)
    heatmap_width <- resolve_dim(width, heatmapWidth)
    heatmap_height <- resolve_dim(height, heatmapHeight)
    if (!is.null(network_plot)) {
      copula_ggsave(
        file.path(outDir, "network.png"),
        plot = network_plot,
        width = network_width,
        height = network_height,
        dpi = dpi
      )
    }
    save_grob(copula_heatmap, file.path(outDir, "copula_cor_heatmap.png"), heatmap_width, heatmap_height, dpi)
    save_grob(pcor_heatmap, file.path(outDir, "pcor_heatmap.png"), heatmap_width, heatmap_height, dpi)
    message("Saved diagnostics to: ", outDir)
  }

  invisible(list(
    network = network_plot,
    copulaCorHeatmap = copula_heatmap,
    pcorHeatmap = pcor_heatmap
  ))
}

#' Plot diagnostics for all fitted strata
#'
#' @param fits Output of [FitCopulaStrata()] or a named list of fit results.
#' @param outDir Base output directory. Each stratum gets a subfolder.
#' @param width Save width in inches (passed to [PlotStratumDiagnostics()]).
#' @param height Save height in inches (passed to [PlotStratumDiagnostics()]).
#' @param dpi Save resolution for PNG files (passed to [PlotStratumDiagnostics()]).
#' @param ... Additional arguments passed to [PlotStratumDiagnostics()].
#' @return Named list of per-stratum plot objects.
#' @export
PlotAllStrata <- function(fits,
                          outDir = "figures",
                          width = 10,
                          height = 10,
                          dpi = 150,
                          ...) {
  fit_results_list <- if (!is.null(fits$fits)) fits$fits else fits
  plots <- list()
  extra <- list(...)
  if (!is.null(extra$width)) width <- extra$width
  if (!is.null(extra$height)) height <- extra$height
  if (!is.null(extra$dpi)) dpi <- extra$dpi
  extra$width <- NULL
  extra$height <- NULL
  extra$dpi <- NULL

  for (nm in names(fit_results_list)) {
    stratum_output_dir <- if (!is.null(outDir)) {
      file.path(outDir, gsub("[^A-Za-z0-9._-]+", "_", nm))
    } else {
      NULL
    }
    plots[[nm]] <- do.call(
      PlotStratumDiagnostics,
      c(
        list(
          fitResult = fit_results_list[[nm]],
          stratumLabel = nm,
          outDir = stratum_output_dir,
          width = width,
          height = height,
          dpi = dpi
        ),
        extra
      )
    )
  }

  invisible(plots)
}
