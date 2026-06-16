#' Resolve node group labels for plotting
#'
#' @param vars Character vector of variable names.
#' @param node_groups Named character vector mapping variable to group, or a
#'   function accepting variable names and returning group labels. Default: all `"Node"`.
#' @return Character vector of group labels.
#' @keywords internal
resolve_node_groups <- function(vars, node_groups = NULL) {
  if (is.null(node_groups)) {
    return(rep("Node", length(vars)))
  }
  if (is.function(node_groups)) {
    return(node_groups(vars))
  }
  if (is.character(node_groups)) {
    out <- node_groups[vars]
    out[is.na(out)] <- "Other"
    return(unname(out))
  }
  stop("node_groups must be NULL, a character named vector, or a function.", call. = FALSE)
}

#' Default palette for node groups
#' @keywords internal
default_group_palette <- function(groups) {
  base <- c(
    "Node" = "grey60",
    "Disease Outcome" = "#E41A1C",
    "T Cell" = "#377EB8",
    "Myeloid" = "#FF7F00",
    "B Cell" = "#984EA3",
    "Spatial" = "#4DAF4A",
    "Other" = "grey60"
  )
  unique_groups <- unique(groups)
  missing <- setdiff(unique_groups, names(base))
  if (length(missing) > 0) {
    extra <- grDevices::colorRampPalette(c("#999999", "#E69F00", "#56B4E9", "#009E73"))(length(missing))
    names(extra) <- missing
    base <- c(base, extra)
  }
  base[unique_groups]
}

#' Format variable names for plot labels
#' @keywords internal
format_node_labels <- function(vars) {
  out <- gsub("Percentage", "%", vars)
  out <- gsub("_Freq$", "", out)
  gsub("_", " ", out)
}

#' Plot a copula network using ggraph
#'
#' Edge width and alpha encode |partial correlation|; edge colour encodes sign.
#'
#' @param result Output of [fit_stratum_copula()].
#' @param title Plot title.
#' @param seed Random seed for FR layout.
#' @param min_pcor Minimum |partial correlation| to display an edge.
#' @param node_groups Optional group mapping (named vector or function).
#' @param print_plot If `TRUE`, print the plot before returning.
#' @return ggplot object, or `NULL` if nothing to plot.
#' @export
plot_copula_network <- function(result,
                                title = "",
                                seed = 42,
                                min_pcor = 0.01,
                                node_groups = NULL,
                                print_plot = TRUE) {
  if (is.null(result)) {
    message("No result to plot.")
    return(invisible(NULL))
  }

  pcor_mat <- result$pcor
  adj <- result$adjacency
  vars <- colnames(pcor_mat)

  edges_idx <- which(adj > 0 & upper.tri(adj), arr.ind = TRUE)
  if (nrow(edges_idx) == 0) {
    message("No edges to plot for: ", title)
    return(invisible(NULL))
  }

  edge_df <- data.frame(
    from = vars[edges_idx[, 1]],
    to = vars[edges_idx[, 2]],
    pcor = pcor_mat[edges_idx],
    stringsAsFactors = FALSE
  )
  edge_df <- edge_df[abs(edge_df$pcor) >= min_pcor, , drop = FALSE]
  edge_df$abs_pcor <- abs(edge_df$pcor)
  edge_df$direction <- ifelse(edge_df$pcor > 0, "Positive", "Negative")

  if (nrow(edge_df) == 0) {
    message("No edges above min_pcor for: ", title)
    return(invisible(NULL))
  }

  groups <- resolve_node_groups(vars, node_groups)
  node_df <- data.frame(
    name = vars,
    group = groups,
    label = format_node_labels(vars),
    stringsAsFactors = FALSE
  )

  g <- tidygraph::tbl_graph(nodes = node_df, edges = edge_df, directed = FALSE)
  palette <- default_group_palette(node_df$group)

  set.seed(seed)
  p <- ggraph(g, layout = "fr") +
    geom_edge_link(
      aes(width = abs_pcor, alpha = abs_pcor, colour = direction),
      lineend = "round"
    ) +
    scale_edge_width_continuous(
      range = c(0.4, 3.5),
      name = "|Partial correlation|",
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    scale_edge_alpha_continuous(range = c(0.25, 1), guide = "none") +
    scale_edge_colour_manual(
      values = c("Positive" = "#C0392B", "Negative" = "#2980B9"),
      name = "Direction"
    ) +
    geom_node_point(aes(fill = group), shape = 21, size = 6, colour = "white", stroke = 0.8) +
    geom_node_text(aes(label = label), repel = TRUE, size = 3, fontface = "bold",
                   bg.colour = "white", bg.r = 0.15) +
    scale_fill_manual(values = palette, name = "Variable group") +
    labs(
      title = title,
      subtitle = bquote(
        n == .(result$n) ~ "|" ~ p == .(length(result$kept_cols)) ~
          "|" ~ edges == .(nrow(edge_df)) ~ "|" ~ lambda == .(round(result$lambda_opt, 4))
      )
    ) +
    theme_graph(base_family = "sans", base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, colour = "grey40"),
      legend.position = "right"
    )

  if (isTRUE(print_plot)) {
    print(p)
  }
  invisible(p)
}
