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
select_comparison_labels <- function(cmp_df, label_threshold = 0.05, max_labels = 25) {
  if (nrow(cmp_df) == 0) {
    return(cmp_df)
  }

  cmp_df <- cmp_df %>%
    dplyr::mutate(
      extreme_score = abs_delta + 0.15 * (abs(value_a) + abs(value_b))
    ) %>%
    dplyr::arrange(dplyr::desc(extreme_score))

  label_df <- cmp_df %>% dplyr::filter(abs_delta >= label_threshold)
  if (nrow(label_df) == 0) {
    label_df <- cmp_df
  }
  label_df %>%
    dplyr::slice_head(n = max_labels) %>%
    dplyr::mutate(label = format_comparison_edge_label(edge)) %>%
    add_quadrant_label_nudges()
}

#' Bias repelled labels toward Q2/Q4 (sparse for signed correlations)
#' @keywords internal
add_quadrant_label_nudges <- function(label_df) {
  xr <- range(label_df$value_b, na.rm = TRUE)
  yr <- range(label_df$value_a, na.rm = TRUE)
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
comparison_scatter_plot <- function(cmp_df,
                                    label_a,
                                    label_b,
                                    title,
                                    label_threshold = 0.05,
                                    max_labels = 25) {
  cmp_df <- cmp_df %>%
    dplyr::mutate(
      direction = factor(direction, levels = unique(direction))
    )

  xr <- range(cmp_df$value_b, na.rm = TRUE)
  yr <- range(cmp_df$value_a, na.rm = TRUE)
  dx <- max(diff(xr), 0.05)
  dy <- max(diff(yr), 0.05)

  p <- ggplot2::ggplot(cmp_df, ggplot2::aes(x = value_b, y = value_a, colour = direction)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
    ggplot2::geom_hline(yintercept = 0, colour = "grey80") +
    ggplot2::geom_vline(xintercept = 0, colour = "grey80") +
    ggplot2::geom_point(ggplot2::aes(size = abs_delta), alpha = 0.7) +
    ggplot2::scale_size_continuous(range = c(1, 5), name = "|\u0394|") +
    ggplot2::labs(
      title = title,
      x = paste(label_b),
      y = paste(label_a),
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
      cmp_df,
      label_threshold = label_threshold,
      max_labels = max_labels
    )
    if (nrow(label_df) > 0) {
      p <- p + ggrepel::geom_label_repel(
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

  p
}

#' Build differential bar chart
#' @keywords internal
comparison_bar_plot <- function(cmp_df, label_a, label_b, title, delta_threshold = 0.05) {
  bar_data <- cmp_df %>%
    dplyr::filter(abs_delta > delta_threshold) %>%
    tidyr::pivot_longer(
      cols = c(value_a, value_b),
      names_to = "stratum",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      stratum = dplyr::recode(stratum, value_a = label_a, value_b = label_b),
      edge = reorder(edge, abs_delta)
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
#' @param cmp Output of [compare_two_strata()].
#' @param out_dir Optional directory to save PNG files.
#' @param label_threshold Minimum |delta| for edge labels in scatter plots.
#' @param max_labels Maximum number of extreme pairs to label in scatter plots.
#' @param delta_threshold Minimum |delta| for bar charts.
#' @param width Plot width in inches.
#' @param height Plot height in inches.
#' @param dpi Resolution.
#' @return List with `plots` keyed by matrix type (each containing `scatter` and `bar`).
#' @export
plot_stratum_comparison <- function(cmp,
                                     out_dir = NULL,
                                     label_threshold = 0.05,
                                     max_labels = 25,
                                     delta_threshold = 0.05,
                                     width = 16,
                                     height = 11,
                                     dpi = 150) {
  if (!inherits(cmp, "copula_stratum_comparison") && !is.list(cmp)) {
    stop("cmp must be output of compare_two_strata().", call. = FALSE)
  }

  plots <- list()

  for (mat_name in names(cmp)) {
    cmp_df <- cmp[[mat_name]]
    if (nrow(cmp_df) == 0) {
      next
    }

    label_a <- cmp_df$label_a[1]
    label_b <- cmp_df$label_b[1]
    mat_label <- if (mat_name == "pcor") "Partial Correlation" else "Copula Correlation"

    p_scatter <- comparison_scatter_plot(
      cmp_df,
      label_a = label_a,
      label_b = label_b,
      title = paste(mat_label, "Comparison:", label_a, "vs.", label_b),
      label_threshold = label_threshold,
      max_labels = max_labels
    )

    p_bar <- comparison_bar_plot(
      cmp_df,
      label_a = label_a,
      label_b = label_b,
      title = paste("Differential Pairs (|", "\u0394", "| > ", delta_threshold, "): ",
                    label_a, " vs. ", label_b, sep = ""),
      delta_threshold = delta_threshold
    )

    plots[[mat_name]] <- list(scatter = p_scatter, bar = p_bar)

    if (!is.null(out_dir)) {
      if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      }
      copula_ggsave(
        file.path(out_dir, paste0(mat_name, "_scatter.png")),
        plot = p_scatter,
        width = width,
        height = height,
        dpi = dpi
      )
      if (!is.null(p_bar)) {
        copula_ggsave(
          file.path(out_dir, paste0(mat_name, "_bar.png")),
          plot = p_bar,
          width = width,
          height = height,
          dpi = dpi
        )
      }
      message("Saved comparison plots to: ", out_dir)
    }
  }

  invisible(list(plots = plots, cmp = cmp))
}
