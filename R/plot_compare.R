#' ggplot theme: egg::theme_article() when available, else theme_bw()
#' @keywords internal
copula_plot_theme <- function(base_size = 12) {
  base <- if (requireNamespace("egg", quietly = TRUE)) {
    egg::theme_article(base_size = base_size)
  } else {
    ggplot2::theme_bw(base_size = base_size)
  }
  base +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.background = ggplot2::element_rect(fill = "white", colour = NA),
      strip.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

#' Shorten edge names for comparison plot labels
#' @keywords internal
format_comparison_edge_label <- function(edges) {
  short <- function(v) {
    v <- format_node_labels(v)
    if (grepl("__", v, fixed = TRUE)) {
      sub("^.*__", "", v)
    } else {
      v
    }
  }
  vapply(strsplit(edges, " -- ", fixed = TRUE), function(parts) {
    parts <- vapply(parts, short, character(1))
    paste(parts, collapse = " -- ")
  }, character(1))
}

#' Select the most extreme points to label on comparison scatter plots
#' @keywords internal
select_comparison_labels <- function(comparison_df, label_threshold = 0.05, max_labels = 25) {
  if (nrow(comparison_df) == 0) {
    return(comparison_df)
  }

  comparison_df <- comparison_df %>%
    dplyr::mutate(
      extreme_score = absDelta + 0.15 * (abs(valueA) + abs(valueB))
    ) %>%
    dplyr::arrange(dplyr::desc(extreme_score))

  label_df <- comparison_df %>% dplyr::filter(absDelta >= label_threshold)
  if (nrow(label_df) == 0) {
    label_df <- comparison_df
  }
  label_df %>%
    dplyr::slice_head(n = max_labels) %>%
    dplyr::mutate(label = format_comparison_edge_label(edge)) %>%
    add_quadrant_label_nudges()
}

#' Bias repelled labels toward Q2/Q4 (sparse for signed correlations)
#' @keywords internal
add_quadrant_label_nudges <- function(label_df) {
  xr <- range(label_df$valueB, na.rm = TRUE)
  yr <- range(label_df$valueA, na.rm = TRUE)
  dx <- max(diff(xr), 0.05)
  dy <- max(diff(yr), 0.05)

  label_df %>%
    dplyr::mutate(
      .q = (dplyr::row_number() - 1L) %% 2L,
      nudge_x = dplyr::if_else(.q == 0L, -0.5 * dx, 0.5 * dx),
      nudge_y = dplyr::if_else(.q == 0L, 0.4 * dy, -0.4 * dy)
    ) %>%
    dplyr::select(-.q)
}

#' Build comparison scatter plot
#' @keywords internal
comparison_scatter_plot <- function(comparison_df,
                                    labelA,
                                    labelB,
                                    title,
                                    label_threshold = 0.05,
                                    max_labels = 25) {
  comparison_df <- comparison_df %>%
    dplyr::mutate(
      direction = factor(direction, levels = unique(direction))
    )

  xr <- range(comparison_df$valueB, na.rm = TRUE)
  yr <- range(comparison_df$valueA, na.rm = TRUE)
  dx <- max(diff(xr), 0.05)
  dy <- max(diff(yr), 0.05)

  scatter_plot <- ggplot2::ggplot(comparison_df, ggplot2::aes(x = valueB, y = valueA, colour = direction)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
    ggplot2::geom_hline(yintercept = 0, colour = "grey80") +
    ggplot2::geom_vline(xintercept = 0, colour = "grey80") +
    ggplot2::geom_point(ggplot2::aes(size = absDelta), alpha = 0.7) +
    ggplot2::scale_size_continuous(range = c(1, 5), name = "|\u0394|") +
    ggplot2::labs(
      title = title,
      x = paste(labelB),
      y = paste(labelA),
      colour = NULL
    ) +
    ggplot2::coord_cartesian(
      xlim = c(xr[1] - 0.5 * dx, xr[2] + 0.5 * dx),
      ylim = c(yr[1] - 0.5 * dy, yr[2] + 0.5 * dy),
      clip = "off"
    ) +
    copula_plot_theme(base_size = 12) +
    ggplot2::theme(plot.margin = ggplot2::margin(20, 20, 20, 20))

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    label_df <- select_comparison_labels(
      comparison_df,
      label_threshold = label_threshold,
      max_labels = max_labels
    )
    if (nrow(label_df) > 0) {
      scatter_plot <- scatter_plot + ggrepel::geom_label_repel(
        data = label_df,
        ggplot2::aes(label = label),
        nudge_x = label_df$nudge_x,
        nudge_y = label_df$nudge_y,
        size = 3.1,
        fill = "white",
        alpha = 0.93,
        label.size = 0.2,
        label.padding = ggplot2::unit(0.2, "lines"),
        label.r = ggplot2::unit(0.12, "lines"),
        max.overlaps = Inf,
        min.segment.length = 0,
        box.padding = 1,
        point.padding = 0.75,
        segment.color = "grey40",
        segment.size = 0.3,
        force = 5,
        force_pull = 0.1,
        max.iter = 20000,
        seed = 42,
        show.legend = FALSE
      )
    }
  }

  scatter_plot
}

#' Build differential bar chart
#' @keywords internal
comparison_bar_plot <- function(comparison_df, labelA, labelB, title, delta_threshold = 0.05) {
  bar_data <- comparison_df %>%
    dplyr::filter(absDelta > delta_threshold) %>%
    tidyr::pivot_longer(
      cols = c(valueA, valueB),
      names_to = "stratum",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      # !! injects function args; bare labelA/labelB collide with comparison_df columns
      stratum = dplyr::recode(stratum, valueA = !!labelA, valueB = !!labelB),
      edge = reorder(edge, absDelta)
    )

  if (nrow(bar_data) == 0) {
    return(NULL)
  }

  ggplot2::ggplot(bar_data, ggplot2::aes(x = value, y = edge, fill = stratum)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7), width = 0.6) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey40") +
    ggplot2::labs(
      title = title,
      x = "Value",
      y = NULL,
      fill = "Stratum"
    ) +
    copula_plot_theme(base_size = 12)
}

#' Plot stratum comparison diagnostics
#'
#' @param cmp Output of [CompareTwoStrata()].
#' @param outDir Optional directory to save PNG files.
#' @param labelThreshold Minimum |delta| for edge labels in scatter plots.
#' @param maxLabels Maximum number of extreme pairs to label in scatter plots.
#' @param deltaThreshold Minimum |delta| for bar charts.
#' @param width Default save width in inches (scatter and bar plots).
#' @param height Default save height in inches (scatter and bar plots).
#' @param scatterWidth Optional scatter plot width override.
#' @param scatterHeight Optional scatter plot height override.
#' @param barWidth Optional bar plot width override.
#' @param barHeight Optional bar plot height override.
#' @param dpi Resolution for saved PNG files.
#' @return List with `plots` keyed by matrix type (each containing `scatter` and `bar`).
#' @export
PlotStratumComparison <- function(cmp,
                                  outDir = NULL,
                                  labelThreshold = 0.05,
                                  maxLabels = 25,
                                  deltaThreshold = 0.05,
                                  width = 16,
                                  height = 11,
                                  scatterWidth = NULL,
                                  scatterHeight = NULL,
                                  barWidth = NULL,
                                  barHeight = NULL,
                                  dpi = 150) {
  if (!inherits(cmp, "CopulaStratumComparison") && !is.list(cmp)) {
    stop("cmp must be output of CompareTwoStrata().", call. = FALSE)
  }

  plots <- list()

  for (mat_name in names(cmp)) {
    comparison_df <- cmp[[mat_name]]
    if (nrow(comparison_df) == 0) {
      next
    }

    label_a <- comparison_df$labelA[1]
    label_b <- comparison_df$labelB[1]
    mat_label <- if (mat_name == "pcor") "Partial Correlation" else "Copula Correlation"

    scatter_plot <- comparison_scatter_plot(
      comparison_df,
      labelA = label_a,
      labelB = label_b,
      title = paste(mat_label, "Comparison:", label_a, "vs.", label_b),
      label_threshold = labelThreshold,
      max_labels = maxLabels
    )

    bar_plot <- comparison_bar_plot(
      comparison_df,
      labelA = label_a,
      labelB = label_b,
      title = paste("Differential Pairs (|", "\u0394", "| > ", deltaThreshold, "): ",
                    label_a, " vs. ", label_b, sep = ""),
      delta_threshold = deltaThreshold
    )

    plots[[mat_name]] <- list(scatter = scatter_plot, bar = bar_plot)

    if (!is.null(outDir)) {
      if (!dir.exists(outDir)) {
        dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
      }
      scatter_width <- resolve_dim(width, scatterWidth)
      scatter_height <- resolve_dim(height, scatterHeight)
      bar_width <- resolve_dim(width, barWidth)
      bar_height <- resolve_dim(height, barHeight)
      copula_ggsave(
        file.path(outDir, paste0(mat_name, "_scatter.png")),
        plot = scatter_plot,
        width = scatter_width,
        height = scatter_height,
        dpi = dpi
      )
      if (!is.null(bar_plot)) {
        copula_ggsave(
          file.path(outDir, paste0(mat_name, "_bar.png")),
          plot = bar_plot,
          width = bar_width,
          height = bar_height,
          dpi = dpi
        )
      }
      message("Saved comparison plots to: ", outDir)
    }
  }

  invisible(list(plots = plots, cmp = cmp))
}
