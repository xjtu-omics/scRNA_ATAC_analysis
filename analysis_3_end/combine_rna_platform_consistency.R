suppressPackageStartupMessages({
  required_packages <- c("ggplot2", "dplyr", "tidyr", "readr", "stringr", "purrr", "scales")
  missing_packages <- setdiff(required_packages, rownames(installed.packages()))
  if (length(missing_packages) > 0) {
    stop(
      "Missing required R packages: ",
      paste(missing_packages, collapse = ", "),
      ". Install them before running this script.",
      call. = FALSE
    )
  }

  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)

env_or_default <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

base_dir <- if (length(args) >= 1) {
  args[[1]]
} else {
  "/path/to/snv_effect_by_cell_type"
}

output_dir <- if (length(args) >= 2) {
  args[[2]]
} else {
  file.path(base_dir, "combined_rna_platform_consistency")
}

snv_table_dir <- file.path(base_dir, "rna_platform_consistency", "tables")
cell_table_dir <- file.path(base_dir, "cell_perturbation_rna_platform_consistency", "tables")
table_dir <- file.path(output_dir, "tables")
plot_dir <- file.path(output_dir, "plots")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

snv_profile_path <- file.path(snv_table_dir, "cellclass_score_profile_wide.csv")
snv_profile_cor_path <- file.path(snv_table_dir, "cellclass_profile_correlation_summary.csv")
cell_profile_path <- file.path(cell_table_dir, "cell_perturbation_score_profile_wide.csv")
cell_profile_cor_path <- file.path(cell_table_dir, "cell_perturbation_profile_correlation_summary.csv")
cell_distribution_path <- file.path(cell_table_dir, "cell_perturbation_distribution_distance_summary.csv")

primary_snv_metric <- "perturb_euclidean_distance"
primary_cell_metric <- "perturbation_score_euclidean"
snv_metric_labels <- c(
  "perturb_euclidean_distance" = "Perturbation Euclidean distance"
)
cell_metric_labels <- c(
  "perturbation_score_euclidean" = "Euclidean score"
)
primary_snv_label <- unname(snv_metric_labels[[primary_snv_metric]])
primary_cell_label <- unname(cell_metric_labels[[primary_cell_metric]])
primary_stats <- c("mean_score", "median_score")
rank_concordance_cutoff <- 2

drop_cell_classes <- c("Unassigned", "Low information")
cell_class_order <- c(
  "Fibroblast",
  "Keratinocyte",
  "Endothelial",
  "Perivascular smooth muscle",
  "T cell",
  "B cell",
  "Plasma cell",
  "Myeloid/dendritic",
  "NK cell",
  "Mast cell",
  "Other immune",
  "Melanocyte",
  "Adnexal epithelial",
  "Schwann",
  "Sebocyte"
)
nature_palette <- c(
  "Fibroblast" = "#4E79A7",
  "Keratinocyte" = "#E15759",
  "Endothelial" = "#59A14F",
  "Perivascular smooth muscle" = "#F28E2B",
  "T cell" = "#B07AA1",
  "B cell" = "#8CD17D",
  "Plasma cell" = "#FF9DA7",
  "Myeloid/dendritic" = "#A0CBE8",
  "NK cell" = "#D37295",
  "Mast cell" = "#FABFD2",
  "Other immune" = "#BAB0AC",
  "Melanocyte" = "#9C755F",
  "Adnexal epithelial" = "#76B7B2",
  "Schwann" = "#EDC948",
  "Sebocyte" = "#AF7AA1"
)
platform_palette <- c("3-end RNA-seq" = "#0072B2", "Full-length RNA-seq" = "#D55E00")

check_file_exists <- function(path) {
  if (!file.exists(path)) {
    stop("Required input file does not exist: ", path, call. = FALSE)
  }
}

walk(
  c(
    snv_profile_path,
    snv_profile_cor_path,
    cell_profile_path,
    cell_profile_cor_path,
    cell_distribution_path
  ),
  check_file_exists
)

safe_cor <- function(x, y, method, min_n = 3L) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < min_n || length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(NA_real_)
  }
  suppressWarnings(cor(x, y, method = method))
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  mean(x)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  median(x)
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  max(x)
}

safe_ratio <- function(full_length, three_end, eps = 1e-12) {
  log2((full_length + eps) / (three_end + eps))
}

direction_label <- function(full_length, three_end, rel_tol = 0.05) {
  ratio <- (full_length + 1e-12) / (three_end + 1e-12)
  case_when(
    !is.finite(ratio) ~ NA_character_,
    ratio > 1 + rel_tol ~ "full_higher",
    ratio < 1 / (1 + rel_tol) ~ "three_end_higher",
    TRUE ~ "similar"
  )
}

theme_nature <- function(base_size = 10) {
  axis_col <- "#2B2B2B"
  theme_bw(base_size = base_size) %+replace%
    theme(
      text = element_text(colour = axis_col, family = "sans"),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.border = element_rect(fill = NA, colour = axis_col, linewidth = 0.4),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = axis_col, linewidth = 0.35),
      axis.ticks = element_line(colour = axis_col, linewidth = 0.35),
      axis.text = element_text(colour = axis_col, size = base_size * 0.82),
      axis.title = element_text(colour = axis_col, face = "bold", size = base_size * 0.95),
      plot.title = element_text(colour = axis_col, face = "bold", hjust = 0, size = base_size * 1.05),
      plot.subtitle = element_text(colour = axis_col, hjust = 0, size = base_size * 0.86),
      plot.caption = element_text(colour = axis_col, hjust = 0, size = base_size * 0.72),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.key = element_rect(fill = "white", colour = NA),
      legend.title = element_text(colour = axis_col, face = "bold", size = base_size * 0.82),
      legend.text = element_text(colour = axis_col, size = base_size * 0.74),
      strip.background = element_rect(fill = "#F4F4F4", colour = axis_col, linewidth = 0.3),
      strip.text = element_text(colour = axis_col, face = "bold", size = base_size * 0.78),
      plot.margin = margin(5, 6, 5, 5)
    )
}

save_plot_png_pdf <- function(filename_png, plot, width, height, dpi = 300, limitsize = FALSE) {
  ggsave(filename = filename_png, plot = plot, width = width, height = height, dpi = dpi, bg = "white", limitsize = limitsize)
  ggsave(filename = sub("\\.png$", ".pdf", filename_png, ignore.case = TRUE), plot = plot, width = width, height = height, device = grDevices::pdf, bg = "white", limitsize = limitsize)
}

format_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
}

profile_rank_summary <- function(profile_df, metric_name, profile_stats) {
  metric_df <- profile_df %>% filter(.data$metric == metric_name)
  map_dfr(
    profile_stats,
    function(stat) {
      three_col <- paste0(stat, "_three_end")
      full_col <- paste0(stat, "_full_length")
      if (!all(c(three_col, full_col) %in% names(metric_df))) {
        return(tibble(profile_stat = stat, median_abs_rank_delta = NA_real_))
      }
      rank_df <- metric_df %>%
        mutate(
          rank_three_end = min_rank(desc(.data[[three_col]])),
          rank_full_length = min_rank(desc(.data[[full_col]])),
          abs_rank_delta = abs(.data$rank_full_length - .data$rank_three_end),
          normalized_abs_rank_delta = .data$abs_rank_delta / pmax(n() - 1, 1)
        )
      tibble(profile_stat = stat, median_abs_rank_delta = median(rank_df$normalized_abs_rank_delta, na.rm = TRUE))
    }
  )
}

snv_profile_wide <- read_csv(snv_profile_path, show_col_types = FALSE)
snv_profile_cor <- read_csv(snv_profile_cor_path, show_col_types = FALSE)
cell_profile_wide <- read_csv(cell_profile_path, show_col_types = FALSE)
cell_profile_cor <- read_csv(cell_profile_cor_path, show_col_types = FALSE)
cell_distribution <- read_csv(cell_distribution_path, show_col_types = FALSE)

snv_primary <- snv_profile_wide %>%
  filter(.data$metric == primary_snv_metric) %>%
  transmute(
    cell_class = as.character(.data$cell_class),
    snv_n_snv_rows_three_end = .data$n_snv_rows_three_end,
    snv_n_unique_snvs_three_end = .data$n_unique_snvs_three_end,
    snv_mean_three_end = .data$mean_score_three_end,
    snv_median_three_end = .data$median_score_three_end,
    snv_p90_three_end = .data$p90_score_three_end,
    snv_p95_three_end = .data$p95_score_three_end,
    snv_top50_mean_three_end = .data$top50_mean_score_three_end,
    snv_n_snv_rows_full_length = .data$n_snv_rows_full_length,
    snv_n_unique_snvs_full_length = .data$n_unique_snvs_full_length,
    snv_mean_full_length = .data$mean_score_full_length,
    snv_median_full_length = .data$median_score_full_length,
    snv_p90_full_length = .data$p90_score_full_length,
    snv_p95_full_length = .data$p95_score_full_length,
    snv_top50_mean_full_length = .data$top50_mean_score_full_length
  )

cell_primary <- cell_profile_wide %>%
  filter(.data$metric == primary_cell_metric) %>%
  transmute(
    cell_class = as.character(.data$cell_class),
    n_cells_three_end = .data$n_cells_three_end,
    n_cells_full_length = .data$n_cells_full_length,
    cell_mean_three_end = .data$mean_score_three_end,
    cell_median_three_end = .data$median_score_three_end,
    cell_p90_three_end = .data$p90_score_three_end,
    cell_p95_three_end = .data$p95_score_three_end,
    cell_top_n_mean_three_end = .data$top_n_mean_score_three_end,
    cell_mean_full_length = .data$mean_score_full_length,
    cell_median_full_length = .data$median_score_full_length,
    cell_p90_full_length = .data$p90_score_full_length,
    cell_p95_full_length = .data$p95_score_full_length,
    cell_top_n_mean_full_length = .data$top_n_mean_score_full_length
  )

joint_profile <- snv_primary %>%
  inner_join(cell_primary, by = "cell_class") %>%
  mutate(cell_class = factor(.data$cell_class, levels = cell_class_order)) %>%
  arrange(.data$cell_class) %>%
  mutate(
    cell_class = as.character(.data$cell_class),
    snv_mean_rank_three_end = min_rank(desc(.data$snv_mean_three_end)),
    snv_mean_rank_full_length = min_rank(desc(.data$snv_mean_full_length)),
    snv_mean_rank_delta = .data$snv_mean_rank_full_length - .data$snv_mean_rank_three_end,
    snv_median_rank_three_end = min_rank(desc(.data$snv_median_three_end)),
    snv_median_rank_full_length = min_rank(desc(.data$snv_median_full_length)),
    snv_median_rank_delta = .data$snv_median_rank_full_length - .data$snv_median_rank_three_end,
    cell_mean_rank_three_end = min_rank(desc(.data$cell_mean_three_end)),
    cell_mean_rank_full_length = min_rank(desc(.data$cell_mean_full_length)),
    cell_mean_rank_delta = .data$cell_mean_rank_full_length - .data$cell_mean_rank_three_end,
    cell_median_rank_three_end = min_rank(desc(.data$cell_median_three_end)),
    cell_median_rank_full_length = min_rank(desc(.data$cell_median_full_length)),
    cell_median_rank_delta = .data$cell_median_rank_full_length - .data$cell_median_rank_three_end,
    snv_mean_log2_ratio_full_vs_three = safe_ratio(.data$snv_mean_full_length, .data$snv_mean_three_end),
    cell_mean_log2_ratio_full_vs_three = safe_ratio(.data$cell_mean_full_length, .data$cell_mean_three_end),
    joint_mean_abs_rank_delta = (abs(.data$snv_mean_rank_delta) + abs(.data$cell_mean_rank_delta)) / 2,
    joint_median_abs_rank_delta = (abs(.data$snv_median_rank_delta) + abs(.data$cell_median_rank_delta)) / 2,
    snv_mean_rank_concordant = abs(.data$snv_mean_rank_delta) <= rank_concordance_cutoff,
    cell_mean_rank_concordant = abs(.data$cell_mean_rank_delta) <= rank_concordance_cutoff,
    joint_mean_rank_concordant = .data$snv_mean_rank_concordant & .data$cell_mean_rank_concordant,
    snv_direction = direction_label(.data$snv_mean_full_length, .data$snv_mean_three_end),
    cell_direction = direction_label(.data$cell_mean_full_length, .data$cell_mean_three_end),
    direction_concordant = .data$snv_direction == .data$cell_direction | .data$snv_direction == "similar" | .data$cell_direction == "similar"
  )

write_csv(joint_profile, file.path(table_dir, "joint_cellclass_consistency_profile.csv"))

snv_rank_summary <- profile_rank_summary(snv_profile_wide, primary_snv_metric, unique(snv_profile_cor$profile_stat))

joint_profile_cor <- bind_rows(
  snv_profile_cor %>%
    left_join(snv_rank_summary, by = "profile_stat") %>%
    mutate(analysis_layer = "SNV effect profile"),
  cell_profile_cor %>%
    mutate(analysis_layer = "Cell perturbation profile")
) %>%
  mutate(
    evidence_tier = if_else(
      (.data$analysis_layer == "SNV effect profile" & .data$metric == primary_snv_metric & .data$profile_stat %in% primary_stats) |
        (.data$analysis_layer == "Cell perturbation profile" & .data$metric == primary_cell_metric & .data$profile_stat %in% primary_stats),
      "primary",
      "secondary"
    )
  ) %>%
  select("analysis_layer", "metric", "metric_label", "profile_stat", "profile_stat_label", "n_cell_classes", "pearson", "spearman", "median_abs_rank_delta", "evidence_tier")

write_csv(joint_profile_cor, file.path(table_dir, "joint_profile_correlation_summary.csv"))

get_cor_row <- function(cor_df, layer, metric_name, stat_name) {
  cor_df %>% filter(.data$analysis_layer == layer, .data$metric == metric_name, .data$profile_stat == stat_name) %>% slice(1)
}

snv_mean_cor <- get_cor_row(joint_profile_cor, "SNV effect profile", primary_snv_metric, "mean_score")
snv_median_cor <- get_cor_row(joint_profile_cor, "SNV effect profile", primary_snv_metric, "median_score")
cell_mean_cor <- get_cor_row(joint_profile_cor, "Cell perturbation profile", primary_cell_metric, "mean_score")
cell_median_cor <- get_cor_row(joint_profile_cor, "Cell perturbation profile", primary_cell_metric, "median_score")

joint_consistency_summary <- tibble(
  analysis_layer = c("SNV effect profile", "Cell perturbation profile", "Joint descriptive summary"),
  primary_metric = c(primary_snv_metric, primary_cell_metric, paste(primary_snv_metric, "+", primary_cell_metric)),
  n_cell_classes = c(snv_mean_cor$n_cell_classes, cell_mean_cor$n_cell_classes, nrow(joint_profile)),
  mean_profile_pearson = c(snv_mean_cor$pearson, cell_mean_cor$pearson, mean(c(snv_mean_cor$pearson, cell_mean_cor$pearson), na.rm = TRUE)),
  mean_profile_spearman = c(snv_mean_cor$spearman, cell_mean_cor$spearman, mean(c(snv_mean_cor$spearman, cell_mean_cor$spearman), na.rm = TRUE)),
  median_profile_pearson = c(snv_median_cor$pearson, cell_median_cor$pearson, mean(c(snv_median_cor$pearson, cell_median_cor$pearson), na.rm = TRUE)),
  median_profile_spearman = c(snv_median_cor$spearman, cell_median_cor$spearman, mean(c(snv_median_cor$spearman, cell_median_cor$spearman), na.rm = TRUE)),
  mean_of_mean_and_median_spearman = c(
    mean(c(snv_mean_cor$spearman, snv_median_cor$spearman), na.rm = TRUE),
    mean(c(cell_mean_cor$spearman, cell_median_cor$spearman), na.rm = TRUE),
    mean(c(snv_mean_cor$spearman, snv_median_cor$spearman, cell_mean_cor$spearman, cell_median_cor$spearman), na.rm = TRUE)
  ),
  interpretation = c(
    sprintf("SNV-level %s profile correlations are summarized across marker-harmonized coarse cell classes.", primary_snv_label),
    sprintf("Cell-level %s profile correlations and distribution differences are summarized after marker-based cell-class harmonization.", primary_cell_label),
    "Agreement can differ by layer and statistic; this is descriptive integration, not a formal meta-analysis."
  )
)
write_csv(joint_consistency_summary, file.path(table_dir, "joint_consistency_score_summary.csv"))

rank_concordance <- bind_rows(
  joint_profile %>%
    transmute(cell_class, analysis_layer = paste("SNV", primary_snv_label, "mean"), profile_stat = "mean_score", `3-end RNA-seq` = snv_mean_rank_three_end, `Full-length RNA-seq` = snv_mean_rank_full_length, score_three_end = snv_mean_three_end, score_full_length = snv_mean_full_length),
  joint_profile %>%
    transmute(cell_class, analysis_layer = paste("Cell perturbation", primary_cell_label, "mean"), profile_stat = "mean_score", `3-end RNA-seq` = cell_mean_rank_three_end, `Full-length RNA-seq` = cell_mean_rank_full_length, score_three_end = cell_mean_three_end, score_full_length = cell_mean_full_length)
) %>%
  pivot_longer(cols = c("3-end RNA-seq", "Full-length RNA-seq"), names_to = "platform", values_to = "rank") %>%
  mutate(
    score = if_else(.data$platform == "3-end RNA-seq", .data$score_three_end, .data$score_full_length),
    cell_class = factor(.data$cell_class, levels = rev(intersect(cell_class_order, joint_profile$cell_class))),
    platform = factor(.data$platform, levels = names(platform_palette))
  ) %>%
  select("cell_class", "analysis_layer", "profile_stat", "platform", "score", "rank")
write_csv(rank_concordance, file.path(table_dir, "joint_cellclass_rank_concordance.csv"))

primary_cell_distribution <- cell_distribution %>% filter(.data$metric == primary_cell_metric)
largest_distribution_diffs <- primary_cell_distribution %>% arrange(desc(.data$ks_statistic)) %>% slice_head(n = 5)

supplementary_evidence <- bind_rows(
  tibble(
    evidence_group = "Cell perturbation distribution distance",
    statistic = c("median_ks_statistic", "max_ks_statistic", "median_standardized_mean_difference"),
    value = as.character(c(safe_median(primary_cell_distribution$ks_statistic), safe_max(primary_cell_distribution$ks_statistic), safe_median(primary_cell_distribution$standardized_mean_difference))),
    note = "Distribution-level differences are supplementary because cells are unpaired across platforms."
  ),
  largest_distribution_diffs %>%
    transmute(
      evidence_group = "Largest cell distribution differences",
      statistic = paste0(.data$cell_class, "_ks_statistic"),
      value = as.character(.data$ks_statistic),
      note = sprintf("Median delta full-minus-3' = %.3f", .data$median_delta_full_minus_three)
    )
)
write_csv(supplementary_evidence, file.path(table_dir, "joint_supplementary_evidence_summary.csv"))

scatter_df <- bind_rows(
  joint_profile %>%
    transmute(
      cell_class,
      analysis_layer = paste("SNV", primary_snv_label, "mean"),
      score_three_end = snv_mean_three_end,
      score_full_length = snv_mean_full_length,
      pearson = snv_mean_cor$pearson,
      spearman = snv_mean_cor$spearman
    ),
  joint_profile %>%
    transmute(
      cell_class,
      analysis_layer = paste("Cell perturbation", primary_cell_label, "mean"),
      score_three_end = cell_mean_three_end,
      score_full_length = cell_mean_full_length,
      pearson = cell_mean_cor$pearson,
      spearman = cell_mean_cor$spearman
    )
) %>%
  mutate(
    cell_class = factor(.data$cell_class, levels = cell_class_order),
    layer_label = sprintf("%s\nPearson %.2f; Spearman %.2f", .data$analysis_layer, .data$pearson, .data$spearman)
  )

p_overview <- ggplot(scatter_df, aes(x = .data$score_three_end, y = .data$score_full_length, colour = .data$cell_class)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#666666", linewidth = 0.35) +
  geom_point(size = 3.2, alpha = 0.9) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, colour = "#2B2B2B", linewidth = 0.42, linetype = "solid") +
  geom_text(aes(label = .data$cell_class), size = 2.25, nudge_x = 0.01, nudge_y = 0.01, colour = "#2B2B2B", show.legend = FALSE, check_overlap = TRUE) +
  facet_wrap(vars(.data$layer_label), scales = "free", nrow = 1) +
  scale_x_continuous(expand = expansion(mult = c(0.16, 0.22))) +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.08))) +
  scale_colour_manual(values = nature_palette, drop = TRUE) +
  labs(
    x = "3' RNA-seq mean score",
    y = "Full-length RNA-seq mean score",
    colour = "Cell class",
    title = "Orthogonal analyses compare coarse cell-class platform profiles",
    subtitle = sprintf("SNV %s and cell-level %s profiles are compared on the same marker-derived cell classes", primary_snv_label, primary_cell_label)
  ) +
  theme_nature(base_size = 9) +
  theme(legend.position = "right", panel.grid.major = element_line(colour = "#E8E8E8", linewidth = 0.22))
save_plot_png_pdf(file.path(plot_dir, "01_joint_consistency_evidence_overview.png"), p_overview, width = 10.8, height = 4.8)

rank_segment_df <- rank_concordance %>%
  select("cell_class", "analysis_layer", "platform", "rank") %>%
  pivot_wider(names_from = "platform", values_from = "rank")

p_rank <- ggplot(rank_concordance, aes(x = .data$rank, y = .data$cell_class)) +
  geom_segment(
    data = rank_segment_df,
    aes(x = .data$`3-end RNA-seq`, xend = .data$`Full-length RNA-seq`, y = .data$cell_class, yend = .data$cell_class),
    inherit.aes = FALSE,
    colour = "#999999",
    linewidth = 0.48
  ) +
  geom_point(aes(fill = .data$platform), shape = 21, size = 3.2, colour = "#2B2B2B", stroke = 0.25) +
  facet_wrap(vars(.data$analysis_layer), nrow = 1) +
  scale_fill_manual(values = platform_palette) +
  scale_x_reverse(breaks = seq_len(nrow(joint_profile)), limits = c(nrow(joint_profile) + 0.5, 0.5)) +
  labs(
    x = "Rank (1 = highest mean score)",
    y = NULL,
    fill = NULL,
    title = "Cell-class rank concordance across RNA platforms",
    subtitle = "Both evidence layers are ranked independently within each platform"
  ) +
  theme_nature(base_size = 9) +
  theme(panel.grid.major.y = element_line(colour = "#E8E8E8", linewidth = 0.22), legend.position = "top")
save_plot_png_pdf(file.path(plot_dir, "02_joint_cellclass_rank_concordance.png"), p_rank, width = 9.2, height = 5.3)

rank_delta_df <- joint_profile %>%
  transmute(
    cell_class = factor(.data$cell_class, levels = rev(intersect(cell_class_order, joint_profile$cell_class))),
    `SNV mean` = abs(.data$snv_mean_rank_delta),
    `SNV median` = abs(.data$snv_median_rank_delta),
    `Cell mean` = abs(.data$cell_mean_rank_delta),
    `Cell median` = abs(.data$cell_median_rank_delta)
  ) %>%
  pivot_longer(cols = -"cell_class", names_to = "evidence", values_to = "abs_rank_delta") %>%
  mutate(label = sprintf("%.0f", .data$abs_rank_delta))

p_rank_heatmap <- ggplot(rank_delta_df, aes(x = .data$evidence, y = .data$cell_class, fill = .data$abs_rank_delta)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = .data$label), size = 2.8, colour = "#2B2B2B", fontface = "bold") +
  scale_fill_gradient(low = "#F7F7F7", high = "#D95F02", limits = c(0, max(rank_delta_df$abs_rank_delta, na.rm = TRUE)), oob = scales::squish) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Absolute\nrank delta",
    title = "Rank disagreement across evidence layers",
    subtitle = "Lower absolute rank deltas indicate stronger platform consistency"
  ) +
  theme_nature(base_size = 9) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid = element_blank(), legend.position = "right")
save_plot_png_pdf(file.path(plot_dir, "03_joint_rank_delta_heatmap.png"), p_rank_heatmap, width = 6.8, height = 4.8)

cor_plot_df <- joint_profile_cor %>%
  filter(.data$evidence_tier == "primary") %>%
  mutate(
    layer_stat = factor(
      paste(.data$analysis_layer, .data$profile_stat_label, sep = "\n"),
      levels = unique(paste(.data$analysis_layer, .data$profile_stat_label, sep = "\n"))
    )
  ) %>%
  select("layer_stat", "pearson", "spearman") %>%
  pivot_longer(cols = c("pearson", "spearman"), names_to = "correlation_type", values_to = "correlation") %>%
  mutate(correlation_type = recode(.data$correlation_type, "pearson" = "Pearson", "spearman" = "Spearman"))

p_cor <- ggplot(cor_plot_df, aes(x = .data$layer_stat, y = .data$correlation, fill = .data$correlation_type)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62, colour = "#2B2B2B", linewidth = 0.18) +
  geom_hline(yintercept = 0, colour = "#2B2B2B", linewidth = 0.3) +
  scale_fill_manual(values = c("Pearson" = "#0072B2", "Spearman" = "#D55E00")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = NULL,
    y = "Cross-platform correlation across cell classes",
    fill = NULL,
    title = "Primary profile correlations across both analysis layers",
    subtitle = sprintf("%s and %s mean and median profiles are summarized together", primary_snv_label, primary_cell_label)
  ) +
  theme_nature(base_size = 9) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "top", panel.grid.major.y = element_line(colour = "#E8E8E8", linewidth = 0.22))
save_plot_png_pdf(file.path(plot_dir, "04_joint_correlation_summary_barplot.png"), p_cor, width = 7.4, height = 4.6)

shared_cell_classes <- joint_profile$cell_class
joint_rank_concordant_count <- sum(joint_profile$joint_mean_rank_concordant, na.rm = TRUE)
direction_concordant_count <- sum(joint_profile$direction_concordant, na.rm = TRUE)

summary_lines <- c(
  "Combined RNA-platform consistency summary",
  sprintf("Base directory: %s", base_dir),
  sprintf("SNV consistency tables: %s", snv_table_dir),
  sprintf("Cell perturbation consistency tables: %s", cell_table_dir),
  sprintf("Output directory: %s", output_dir),
  sprintf("Shared marker-derived coarse cell classes: %d (%s)", length(shared_cell_classes), paste(shared_cell_classes, collapse = ", ")),
  sprintf("Primary metric mapping: SNV layer = %s; cell perturbation layer = %s", primary_snv_metric, primary_cell_metric),
  "Primary cross-platform profile evidence:",
  sprintf("  - SNV %s mean profile: Pearson %.3f; Spearman %.3f", primary_snv_label, snv_mean_cor$pearson, snv_mean_cor$spearman),
  sprintf("  - SNV %s median profile: Pearson %.3f; Spearman %.3f", primary_snv_label, snv_median_cor$pearson, snv_median_cor$spearman),
  sprintf("  - Cell perturbation %s mean profile: Pearson %.3f; Spearman %.3f", primary_cell_label, cell_mean_cor$pearson, cell_mean_cor$spearman),
  sprintf("  - Cell perturbation %s median profile: Pearson %.3f; Spearman %.3f", primary_cell_label, cell_median_cor$pearson, cell_median_cor$spearman),
  sprintf("Joint mean-rank concordant cell classes with abs rank delta <= %d in both layers: %d/%d", rank_concordance_cutoff, joint_rank_concordant_count, nrow(joint_profile)),
  sprintf("Direction-concordant or similar cell classes by mean score: %d/%d", direction_concordant_count, nrow(joint_profile)),
  "Supplementary evidence:",
  sprintf("  - Cell distribution median KS statistic for the primary metric: %.3f", safe_median(primary_cell_distribution$ks_statistic)),
  "Interpretation:",
  "  - Agreement differs by analysis layer and correlation statistic; the reported Pearson, Spearman, rank-concordance, and distribution-distance results should be interpreted together rather than as uniform cross-platform reproducibility.",
  "  - Because cells are not directly paired, the result should be framed as a marker-harmonized cell-class comparison rather than one-to-one cell reproducibility.",
  sprintf("Tables: %s", table_dir),
  sprintf("Plots: %s", plot_dir)
)
writeLines(summary_lines, file.path(output_dir, "combined_rna_platform_consistency_summary.txt"), useBytes = TRUE)

message("Finished combined RNA-platform consistency analysis.")
message("Tables: ", table_dir)
message("Plots: ", plot_dir)
message("Summary: ", file.path(output_dir, "combined_rna_platform_consistency_summary.txt"))
