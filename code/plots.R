plot_top_enrichment_results <- function(data, class, fontsize = 9) {
  max_nes <- data |>
    slice_max(order_by = abs(NES), n = 1) |>
    pull(NES) |>
    abs()

  up_plot <- data |>
    filter(NES >= 0) |>
    ggplot(aes(
      x = NES, y = fct_reorder(geneset, NES),
      fill = as_factor(sign(NES))
    )) +
    geom_col() +
    labs(x = "NES", y = NULL) +
    scale_x_continuous(
      position = "top", limits = c(0, max_nes),
      breaks = ~ unique(floor(pretty(seq(0, (max(.x) + 1) * 1.1))))
    ) +
    theme(legend.position = "none") +
    scale_fill_manual(values = AachenColorPalette::aachen_color("red")) +
    my_theme(grid = "no", fsize = fontsize) +
    geom_text(aes(label = stars),
      vjust = +0.75, hjust = 1.5,
      color = "white", size = fontsize / (14 / 5)
    ) +
    theme(plot.margin = margin(0, -0.15, 0, 0, "cm"))

  down_plot <- data |>
    filter(NES < 0) |>
    ggplot(aes(
      x = -NES, y = fct_reorder(geneset, NES),
      fill = as_factor(sign(NES))
    )) +
    geom_col() +
    labs(x = "NES", y = NULL) +
    scale_x_reverse(
      limits = c(max_nes, 0),
      breaks = ~ unique(floor(pretty(seq(0, (max(.x) + 1) * 1.1)))),
      labels = ~ -.x
    ) +
    scale_y_discrete(position = "right") +
    theme(legend.position = "none") +
    scale_fill_manual(values = AachenColorPalette::aachen_color("blue")) +
    my_theme(grid = "no", fsize = fontsize) +
    geom_text(aes(label = stars),
      vjust = +0.75, hjust = -0.5,
      color = "white", size = fontsize / (14 / 5)
    ) +
    theme(plot.margin = margin(0, 0, 0, -0.15, "cm"))

  p <- up_plot + down_plot + ggtitle(str_glue("Enriched {class}"))

  return(p)
}

plot_top_deg_genes <- function(data, geneset_label, fontsize = 9, ...) {
  data <- data |>
    mutate(tick_color = "black")

  max_logfc <- data |>
    slice_max(order_by = abs(logFC), n = 1) |>
    pull(logFC) |>
    abs()

  n_up <- data |>
    filter(logFC >= 0) |>
    nrow()

  n_down <- data |>
    filter(logFC < 0) |>
    nrow()

  if (n_up != n_down) {
    diff <- abs(n_up - n_down)
    # case 1: less downregulated genes
    if (n_down < n_up) {
      data <- data |>
        add_row(
          gene = unlist(map(seq(diff), ~ str_c(rep(" ", .x), collapse = ""))),
          logFC = 0,
          regulation = factor("down"),
          `sign(logFC)` = -1,
          tick_color = "transparent"
        )
    } else if (n_up < n_down) {
      data <- data |>
        add_row(
          gene = unlist(map(seq(diff), ~ str_c(rep(" ", .x), collapse = ""))),
          logFC = 0,
          regulation = factor("up"),
          `sign(logFC)` = 1,
          tick_color = "transparent"
        )
    }
  }

  up_plot <- data |>
    filter(logFC >= 0 & regulation == "up") |>
    ggplot(aes(
      x = logFC, y = fct_reorder(gene, logFC),
      fill = as_factor(`sign(logFC)`)
    )) +
    geom_col() +
    labs(x = "Upregulated logFC", y = NULL) +
    scale_x_continuous(
      position = "top", limits = c(0, max_logfc),
      breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))
    ) +
    theme(legend.position = "none") +
    scale_fill_manual(values = AachenColorPalette::aachen_color("red")) +
    my_theme(grid = "no", fsize = fz) +
    geom_text(aes(label = stars),
      vjust = +0.75, hjust = 1.15,
      color = "white", size = (fontsize - 1) / (14 / 5)
    ) +
    theme(
      plot.margin = margin(0, -0.5, 0, 0, "cm"),
      axis.ticks.y = element_line(
        colour = rev(pull(
          filter(data, logFC >= 0 & regulation == "up"),
          tick_color
        ))
      )
    )

  # check if the up plot has only transparent labels, i.e no bar
  if (!("black" %in% unique(pull(
    filter(data, logFC >= 0 & regulation == "up"),
    tick_color
  )))) {
    up_plot <- up_plot +
      annotate(
        geom = "text",
        x = max_logfc / 2,
        y = (nrow(filter(data, logFC >= 0 & regulation == "up")) / 2) + 0.5,
        label = "No differentially\nexpressed gene",
        color = "gray",
        size = (fontsize - 1) / (14 / 5)
      )
  }


  down_plot <- data |>
    filter(logFC <= 0 & regulation == "down") |>
    ggplot(
      aes(
        x = -logFC,
        y = fct_reorder(gene, logFC),
        fill = as_factor(`sign(logFC)`)
      )
    ) +
    geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
    labs(x = "Downregulated logFC", y = NULL) +
    scale_y_discrete(position = "right") +
    scale_x_reverse(
      limits = c(max_logfc, 0),
      breaks = ~ unique(floor(pretty(seq(0, (max(.x) + 1) * 1.1)))),
      labels = ~ -.x
    ) +
    theme(legend.position = "none") +
    scale_fill_manual(values = AachenColorPalette::aachen_color("blue")) +
    my_theme(grid = "no", fsize = fz) +
    geom_text(aes(label = stars),
      vjust = +0.75, hjust = -0.15,
      color = "white", size = (fontsize - 1) / (14 / 5)
    ) +
    theme(
      plot.margin = margin(0, 0, 0, -0.5, "cm"),
      axis.ticks.y = element_line(
        colour = pull(
          filter(data, logFC <= 0 & regulation == "down"),
          tick_color
        )
      )
    )

  # check if the down plot has only transparent labels, i.e no bar
  if (!("black" %in% unique(
    pull(filter(data, logFC <= 0 & regulation == "down"), tick_color)
  ))) {
    down_plot <- down_plot +
      annotate(
        geom = "text",
        x = max_logfc / 2,
        y = (nrow(filter(data, logFC <= 0 & regulation == "down")) / 2) + 0.5,
        label = "No differentially\nexpressed gene",
        color = "gray", size = (fontsize - 1) / (14 / 5)
      )
  }


  up_plot + down_plot +
    plot_annotation(
      subtitle = geneset_label,
      theme = theme(
        plot.subtitle = element_text(size = fontsize + 1)
      )
    )
}

plot_correlation <- function(data, fontsize = 9, ...) {
  data |>
    ggplot(aes(x = hs_logFC, y = mm_logFC, label = gene)) +
    geom_point(size = 1.5) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_smooth(
      method = "lm", se = FALSE, size = 0.75, formula = y ~ x,
      color = aachen_color("green")
    ) +
    my_theme(fsize = fontsize) +
    labs(x = "Human logFC", y = "Mouse logFC") +
    geom_text(
      data = distinct(data, anno),
      aes(label = anno),
      x = -Inf, y = Inf,
      hjust = -0.05, vjust = 1.5,
      size = fontsize / (14 / 5)
    ) +
    theme(
      axis.line = element_blank(), axis.ticks = element_blank(),
      legend.position = "none"
    )
}
