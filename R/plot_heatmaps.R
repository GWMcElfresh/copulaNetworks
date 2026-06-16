#' Plot copula correlation heatmap
#'
#' @param result Output of [fit_stratum_copula()].
#' @param title Plot title.
#' @param vars Optional subset of variables to display.
#' @return pheatmap grob.
#' @export
plot_copula_cor_heatmap <- function(result, title = "Copula Correlation Matrix", vars = NULL) {
  if (is.null(result)) {
    stop("result is NULL.", call. = FALSE)
  }
  mat <- result$copula_cor
  if (!is.null(vars)) {
    vars <- intersect(vars, colnames(mat))
    mat <- mat[vars, vars, drop = FALSE]
  }
  pheatmap::pheatmap(
    mat,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
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
#' @param result Output of [fit_stratum_copula()].
#' @param title Plot title.
#' @param vars Optional subset of variables to display.
#' @param zero_diag If `TRUE`, zero the diagonal for display.
#' @return pheatmap grob.
#' @export
plot_pcor_heatmap <- function(result,
                              title = "Partial Correlation Matrix (Glasso)",
                              vars = NULL,
                              zero_diag = TRUE) {
  if (is.null(result)) {
    stop("result is NULL.", call. = FALSE)
  }
  mat <- result$pcor
  if (isTRUE(zero_diag)) {
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
#' @param fit_result Output of [fit_stratum_copula()].
#' @param stratum_label Label used in plot titles and file names.
#' @param out_dir Directory for saved PNG/PDF files. If `NULL`, plots are not saved.
#' @param min_pcor Minimum |partial correlation| for network edges.
#' @param node_groups Optional node group mapping.
#' @param seed Layout seed.
#' @param width Plot width in inches.
#' @param height Plot height in inches.
#' @param dpi Resolution for saved files.
#' @return List with ggplot/grob objects: `network`, `copula_cor_heatmap`, `pcor_heatmap`.
#' @export
plot_stratum_diagnostics <- function(fit_result,
                                     stratum_label = "stratum",
                                     out_dir = NULL,
                                     min_pcor = 0.01,
                                     node_groups = NULL,
                                     seed = 42,
                                     width = 10,
                                     height = 10,
                                     dpi = 150) {
  if (is.null(fit_result)) {
    stop("fit_result is NULL.", call. = FALSE)
  }

  title_base <- stratum_label
  p_net <- plot_copula_network(
    fit_result,
    title = paste(title_base, "\u2014 Nonparanormal Copula Network"),
    seed = seed,
    min_pcor = min_pcor,
    node_groups = node_groups,
    print_plot = FALSE
  )

  ph_cor <- plot_copula_cor_heatmap(
    fit_result,
    title = paste("Copula Correlation Matrix (", title_base, ")", sep = "")
  )
  ph_pcor <- plot_pcor_heatmap(
    fit_result,
    title = paste("Partial Correlation Matrix (Glasso,", title_base, ")", sep = " ")
  )

  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    if (!is.null(p_net)) {
      copula_ggsave(
        file.path(out_dir, "network.png"),
        plot = p_net,
        width = width,
        height = height,
        dpi = dpi
      )
    }
    save_grob(ph_cor, file.path(out_dir, "copula_cor_heatmap.png"), width, height, dpi)
    save_grob(ph_pcor, file.path(out_dir, "pcor_heatmap.png"), width, height, dpi)
    message("Saved diagnostics to: ", out_dir)
  }

  invisible(list(
    network = p_net,
    copula_cor_heatmap = ph_cor,
    pcor_heatmap = ph_pcor
  ))
}

#' Plot diagnostics for all fitted strata
#'
#' @param fits Output of [fit_copula_strata()] or a named list of fit results.
#' @param out_dir Base output directory. Each stratum gets a subfolder.
#' @param ... Passed to [plot_stratum_diagnostics()].
#' @return Named list of per-stratum plot objects.
#' @export
plot_all_strata <- function(fits, out_dir = "figures", ...) {
  fit_list <- if (!is.null(fits$fits)) fits$fits else fits
  plots <- list()

  for (nm in names(fit_list)) {
    stratum_dir <- if (!is.null(out_dir)) {
      file.path(out_dir, gsub("[^A-Za-z0-9._-]+", "_", nm))
    } else {
      NULL
    }
    plots[[nm]] <- plot_stratum_diagnostics(
      fit_list[[nm]],
      stratum_label = nm,
      out_dir = stratum_dir,
      ...
    )
  }

  invisible(plots)
}
