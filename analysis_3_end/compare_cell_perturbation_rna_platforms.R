suppressPackageStartupMessages({
  required_packages <- c("ggplot2", "dplyr", "tidyr", "readr", "stringr", "purrr", "scales", "ragg")
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
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)

env_or_default <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

input_dir <- if (length(args) >= 1) {
  args[[1]]
} else {
  "I:/极空间/文献插图/最终版/补充材料/snv_effect_by_cell_type"
}

output_dir <- if (length(args) >= 2) {
  args[[2]]
} else {
  file.path(input_dir, "cell_perturbation_rna_platform_consistency")
}

top_n <- if (length(args) >= 3) {
  as.integer(args[[3]])
} else {
  200L
}
if (is.na(top_n) || top_n <= 0) {
  stop("top_n must be a positive integer.", call. = FALSE)
}

three_end_score_path <- env_or_default(
  "PRISMSNV_THREE_END_CELL_PATH",
  file.path(input_dir, "3-end_cell_perturbation_scores.csv")
)
full_length_score_path <- env_or_default(
  "PRISMSNV_FULL_LENGTH_CELL_PATH",
  file.path(input_dir, "full-length_cell_perturbation_scores.csv")
)
three_end_marker_path <- env_or_default(
  "PRISMSNV_THREE_END_MARKER_PATH",
  file.path(input_dir, "3-end_cell_cluster_marker_genes_top10.csv")
)
full_length_marker_path <- env_or_default(
  "PRISMSNV_FULL_LENGTH_MARKER_PATH",
  file.path(input_dir, "full-length_cell_cluster_marker_genes_top10.csv")
)

table_dir <- file.path(output_dir, "tables")
plot_dir <- file.path(output_dir, "plots")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

score_columns <- "perturbation_score_euclidean"
primary_score_column <- "perturbation_score_euclidean"
min_pairs_for_correlation <- 3L
min_cells_for_distribution_test <- 20L
drop_cell_classes <- c("Unassigned", "Low information")

metric_labels <- c(
  "perturbation_score_euclidean" = "Euclidean score"
)
primary_metric_label <- unname(metric_labels[[primary_score_column]])
profile_stat_labels <- c(
  "mean_score" = "Mean",
  "median_score" = "Median",
  "p25_score" = "P25",
  "p75_score" = "P75",
  "p90_score" = "P90",
  "p95_score" = "P95",
  "p99_score" = "P99",
  "positive_fraction" = "Positive fraction",
  "top_n_mean_score" = sprintf("Top %d mean", top_n)
)

marker_panels <- list(
  "Fibroblast" = c(
    "COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL6A1", "COL6A2", "COL6A3",
    "DCN", "LUM", "CCDC80", "FBN1", "FBLN2", "MMP2", "C3", "CFH", "ABI3BP", "LAMA2",
    "PRRX1", "THBS2", "MEG3", "TWIST2", "ADAMTS2", "ROR2", "DCLK1", "SLIT3", "LTBP2"
  ),
  "Keratinocyte" = c(
    "KRT1", "KRT5", "KRT10", "KRT14", "KRT15", "KRT16", "KRT17", "COL17A1", "SFN",
    "TP63", "DSG1", "DSG3", "DSC3", "DSP", "PKP1", "PKP3", "JUP", "PERP", "IRF6",
    "S100A2", "S100A8", "S100A9", "S100A14", "AQP3", "CSTA", "LAMB3", "SERPINB5"
  ),
  "Endothelial" = c(
    "PECAM1", "VWF", "KDR", "FLT1", "CDH5", "EMCN", "RAMP2", "ADGRL4", "ERG", "ENG",
    "LDB2", "PTPRB", "PTPRM", "RHOJ", "MAGI1", "PLVAP", "ACKR1", "SELE", "SELP",
    "CALCRL", "CLDN5", "ESAM", "PROX1", "CCL21", "MMRN1", "FLT4"
  ),
  "Perivascular smooth muscle" = c(
    "ACTA2", "TAGLN", "MYH11", "MYL9", "MYLK", "CNN1", "CALD1", "RGS5", "MCAM",
    "NOTCH3", "PDGFRB", "PRKG1", "TPM2", "SYNPO2", "COL4A2", "ITGA7", "GPC6",
    "ZFHX3", "EBF1", "NR2F2"
  ),
  "T cell" = c(
    "CD3D", "CD3E", "CD2", "TRAC", "TRBC1", "TRBC2", "IL32", "CD247", "FYN",
    "SKAP1", "ITK", "CD96", "CTLA4", "FOXP3", "IL2RA", "IKZF2", "SPOCK2",
    "MALT1", "CXCR4", "STK17B", "CD52"
  ),
  "B cell" = c(
    "MS4A1", "BANK1", "CD79A", "CD79B", "CD37", "PAX5", "BLK", "MEF2C",
    "FCRL5", "BCL2", "LYN", "SWAP70", "RIPOR2", "CD74", "CD83",
    "HLA-DRA", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1"
  ),
  "Plasma cell" = c(
    "MZB1", "XBP1", "JCHAIN", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGKC",
    "TXNDC5", "SSR4", "DERL3", "FKBP11", "HERPUD1", "SEC11C", "ST6GAL1",
    "TENT5C", "IRF4", "POU2AF1"
  ),
  "Myeloid/dendritic" = c(
    "LYZ", "LCP1", "TYROBP", "FCER1G", "IRF8", "CST3", "HLA-DRA", "HLA-DRB1",
    "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1", "CD74", "CD83", "CSF2RA",
    "FCGR2A", "ITGAX", "PLEK", "WDFY4", "CPVL", "CTSZ", "CTSH", "AIF1"
  ),
  "NK cell" = c("GNLY", "NKG7", "KLRD1", "KLRK1", "CTSW", "GZMB", "IL2RB", "STAT4", "SYTL3"),
  "Mast cell" = c("CPA3", "TPSB2", "TPSAB1", "MS4A2", "KIT", "HDC", "HPGDS", "HPGD", "SLC18A2", "SLC24A3", "IL18R1", "IL1RL1", "GATA2", "VWA5A"),
  "Other immune" = c("PTPRC", "CD52", "SRGN", "SAMSN1", "ARHGAP15", "IRF7", "IL3RA", "RHEX", "SEL1L3", "SLC15A4", "PAG1", "STK4", "RUBCN", "LPXN"),
  "Melanocyte" = c("MLANA", "PMEL", "DCT", "TYR", "TYRP1", "MITF", "TRPM1", "KIT", "FMN1", "LINC01320", "GCNT2", "INPP4B", "PLEKHA5", "SOX5", "CAPN3"),
  "Adnexal epithelial" = c("KRT7", "KRT8", "KRT18", "KRT19", "CLDN4", "ELF3", "AZGP1", "PIP", "DCD", "SCGB1B2P", "MUCL1", "SLC12A2", "TSPAN8", "ESRRG", "PLA2R1", "ZG16B", "NDRG2"),
  "Schwann" = c("SOX10", "ERBB3", "PMP22", "MPZ", "PLP1", "CDH19", "GPM6B", "NTM", "L1CAM", "NRXN1", "NRXN3", "CADM1", "GAS7", "FRMD5"),
  "Sebocyte" = c("FASN", "FADS2", "APOE", "APOC1", "ACSBG1", "ACSL1", "ALOX15B", "HMGCS1", "FDPS", "THRSP", "CYB5A", "MGST1"),
  "Low information" = c(
    "RPS3", "RPS6", "RPS7", "RPS18", "RPS21", "RPS23", "RPS27", "RPL10A", "RPL12",
    "RPL17", "RPL23A", "RPL31", "RPL32", "RPL34", "RPL35A", "RPL36", "RPL36A",
    "RPL39", "MT-CO1", "MT-CO2", "MT-CO3", "MT-ND1", "MT-ND2", "HSPA1A", "HSPA1B",
    "UBC", "GAS5"
  )
)

cell_class_order <- setdiff(names(marker_panels), drop_cell_classes)
nature_palette <- c(
  "Fibroblast" = "#4E79A7", "Keratinocyte" = "#E15759", "Endothelial" = "#59A14F",
  "Perivascular smooth muscle" = "#F28E2B", "T cell" = "#B07AA1", "B cell" = "#8CD17D",
  "Plasma cell" = "#FF9DA7", "Myeloid/dendritic" = "#A0CBE8", "NK cell" = "#D37295",
  "Mast cell" = "#FABFD2", "Other immune" = "#BAB0AC", "Melanocyte" = "#9C755F",
  "Adnexal epithelial" = "#76B7B2", "Schwann" = "#EDC948", "Sebocyte" = "#AF7AA1"
)
platform_palette <- c("3-end RNA-seq" = "#0072B2", "Full-length RNA-seq" = "#D55E00")

check_file_exists <- function(path) {
  if (!file.exists(path)) {
    stop("Required input file does not exist: ", path, call. = FALSE)
  }
}

walk(c(three_end_score_path, full_length_score_path, three_end_marker_path, full_length_marker_path), check_file_exists)

obsolete_plot_files <- list.files(
  plot_dir,
  pattern = "^04_.*\\.(png|pdf)$",
  full.names = TRUE
)
if (length(obsolete_plot_files) > 0) {
  unlink(obsolete_plot_files, force = TRUE)
}

obsolete_table_files <- list.files(
  table_dir,
  pattern = "^(04_.*|cellclass_mean_score_scatter.*)\\.csv$",
  full.names = TRUE
)
if (length(obsolete_table_files) > 0) {
  unlink(obsolete_table_files, force = TRUE)
}

normalize_gene <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    str_replace("\\..*$", "") %>%
    str_to_upper()
}

safe_cor <- function(x, y, method, min_n = min_pairs_for_correlation) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < min_n || length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(NA_real_)
  }
  suppressWarnings(cor(x, y, method = method))
}

safe_numeric_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

safe_numeric_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  sd(x)
}

safe_quantile <- function(x, p) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  as.numeric(quantile(x, p, names = FALSE))
}

top_numeric_mean <- function(x, n = top_n) {
  x <- sort(x[is.finite(x)], decreasing = TRUE)
  if (length(x) == 0) return(NA_real_)
  mean(head(x, min(n, length(x))))
}

safe_ks <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) < min_cells_for_distribution_test || length(y) < min_cells_for_distribution_test) {
    return(tibble(ks_statistic = NA_real_, ks_p_value = NA_real_, test_status = "skipped_too_few_cells"))
  }
  result <- tryCatch(suppressWarnings(ks.test(x, y)), error = function(exc) NULL)
  if (is.null(result)) {
    return(tibble(ks_statistic = NA_real_, ks_p_value = NA_real_, test_status = "ks_failed"))
  }
  tibble(ks_statistic = as.numeric(result$statistic), ks_p_value = result$p.value, test_status = "ok")
}

safe_wilcox_p <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) < min_cells_for_distribution_test || length(y) < min_cells_for_distribution_test) {
    return(NA_real_)
  }
  result <- tryCatch(suppressWarnings(wilcox.test(x, y)), error = function(exc) NULL)
  if (is.null(result)) return(NA_real_)
  result$p.value
}

pooled_sd <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  sqrt(((length(x) - 1) * var(x) + (length(y) - 1) * var(y)) / (length(x) + length(y) - 2))
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

read_marker_table <- function(path, platform_label) {
  marker_df <- read_csv(path, col_types = cols(.default = col_character()), show_col_types = FALSE)
  names(marker_df) <- str_to_lower(names(marker_df))
  required_cols <- c("group", "names")
  missing_cols <- setdiff(required_cols, names(marker_df))
  if (length(missing_cols) > 0) {
    stop("Marker file is missing required columns: ", path, " | ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  marker_df %>%
    transmute(
      platform = platform_label,
      celltype = as.character(.data$group),
      gene = normalize_gene(.data$names),
      marker_score = suppressWarnings(as.numeric(.data$scores)),
      marker_logfc = suppressWarnings(as.numeric(.data$logfoldchanges))
    ) %>%
    filter(!is.na(.data$gene), .data$gene != "")
}

assign_one_celltype <- function(genes) {
  genes <- unique(genes)
  score_df <- imap_dfr(
    marker_panels,
    function(panel_genes, cell_class) {
      panel_genes <- normalize_gene(panel_genes)
      matched <- intersect(genes, panel_genes)
      tibble(cell_class = cell_class, marker_hits = length(matched), matched_marker_genes = paste(matched, collapse = ";"))
    }
  ) %>%
    mutate(panel_priority = match(.data$cell_class, names(marker_panels))) %>%
    arrange(desc(.data$marker_hits), .data$panel_priority)

  top <- score_df %>% slice(1)
  second_hits <- if (nrow(score_df) >= 2) score_df$marker_hits[[2]] else 0L
  if (is.na(top$marker_hits) || top$marker_hits <= 0) {
    return(tibble(cell_class = "Unassigned", marker_hits = 0L, second_marker_hits = second_hits, assignment_margin = 0L, matched_marker_genes = NA_character_))
  }
  tibble(
    cell_class = top$cell_class,
    marker_hits = top$marker_hits,
    second_marker_hits = second_hits,
    assignment_margin = top$marker_hits - second_hits,
    matched_marker_genes = top$matched_marker_genes
  )
}

assign_celltypes <- function(marker_df) {
  marker_df %>%
    group_by(.data$platform, .data$celltype) %>%
    summarise(n_marker_genes = n_distinct(.data$gene), marker_genes = paste(unique(.data$gene), collapse = ";"), .groups = "drop") %>%
    mutate(assignment = map(str_split(.data$marker_genes, ";"), assign_one_celltype)) %>%
    unnest("assignment")
}

read_cell_score_table <- function(path, platform_label, assignment_df) {
  score_df <- read_csv(path, col_types = cols(.default = col_character()), show_col_types = FALSE)
  required_cols <- c("cell_id", "cell_cluster", score_columns)
  missing_cols <- setdiff(required_cols, names(score_df))
  if (length(missing_cols) > 0) {
    stop("Cell perturbation score file is missing required columns: ", path, " | ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  score_df %>%
    mutate(
      platform = platform_label,
      cell_id = as.character(.data$cell_id),
      source_cluster = as.character(.data$cell_cluster),
      across(all_of(score_columns), ~ suppressWarnings(as.numeric(.x)))
    ) %>%
    select("cell_id", "source_cluster", all_of(score_columns), "platform") %>%
    left_join(
      assignment_df %>% select("platform", celltype = "celltype", "cell_class", "marker_hits", "assignment_margin", "matched_marker_genes"),
      by = c("platform", "source_cluster" = "celltype")
    ) %>%
    mutate(cell_class = if_else(is.na(.data$cell_class), "Unassigned", .data$cell_class))
}

make_coverage_summary <- function(cell_scores_all) {
  cell_scores_all %>%
    mutate(retained = !.data$cell_class %in% drop_cell_classes) %>%
    group_by(.data$platform, .data$cell_class, .data$retained) %>%
    summarise(n_cells = n(), n_source_clusters = n_distinct(.data$source_cluster), source_clusters = paste(sort(unique(.data$source_cluster)), collapse = ";"), .groups = "drop") %>%
    arrange(.data$platform, desc(.data$retained), factor(.data$cell_class, levels = c(cell_class_order, drop_cell_classes)))
}

make_cluster_count_summary <- function(cell_scores_all) {
  cell_scores_all %>%
    count(.data$platform, .data$source_cluster, .data$cell_class, .data$marker_hits, .data$assignment_margin, .data$matched_marker_genes, name = "n_cells") %>%
    mutate(retained = !.data$cell_class %in% drop_cell_classes) %>%
    arrange(.data$platform, suppressWarnings(as.integer(.data$source_cluster)))
}

make_cellclass_score_profile <- function(cell_scores) {
  cell_scores %>%
    pivot_longer(cols = all_of(score_columns), names_to = "metric", values_to = "score") %>%
    filter(is.finite(.data$score)) %>%
    group_by(.data$platform, .data$cell_class, .data$metric) %>%
    summarise(
      n_cells = n(),
      n_source_clusters = n_distinct(.data$source_cluster),
      source_clusters = paste(sort(unique(.data$source_cluster)), collapse = ";"),
      mean_score = mean(.data$score),
      median_score = median(.data$score),
      sd_score = safe_numeric_sd(.data$score),
      mad_score = mad(.data$score, na.rm = TRUE),
      p05_score = safe_quantile(.data$score, 0.05),
      p25_score = safe_quantile(.data$score, 0.25),
      p75_score = safe_quantile(.data$score, 0.75),
      p90_score = safe_quantile(.data$score, 0.90),
      p95_score = safe_quantile(.data$score, 0.95),
      p99_score = safe_quantile(.data$score, 0.99),
      min_score = min(.data$score),
      max_score = max(.data$score),
      positive_fraction = mean(.data$score > 0),
      top_n_mean_score = top_numeric_mean(.data$score, top_n),
      .groups = "drop"
    ) %>%
    mutate(
      metric_label = metric_labels[.data$metric],
      cell_class = factor(.data$cell_class, levels = cell_class_order)
    ) %>%
    arrange(.data$metric, .data$cell_class, .data$platform)
}

make_profile_wide <- function(profile_df) {
  inner_join(
    profile_df %>% filter(.data$platform == "3-end RNA-seq"),
    profile_df %>% filter(.data$platform == "Full-length RNA-seq"),
    by = c("cell_class", "metric", "metric_label"),
    suffix = c("_three_end", "_full_length")
  )
}

summarise_profile_correlations <- function(profile_wide_df) {
  profile_stats <- names(profile_stat_labels)
  map_dfr(
    score_columns,
    function(metric_name) {
      metric_df <- profile_wide_df %>% filter(.data$metric == metric_name)
      map_dfr(
        profile_stats,
        function(stat_name) {
          three_col <- paste0(stat_name, "_three_end")
          full_col <- paste0(stat_name, "_full_length")
          rank_three <- percent_rank(metric_df[[three_col]])
          rank_full <- percent_rank(metric_df[[full_col]])
          tibble(
            metric = metric_name,
            metric_label = metric_labels[[metric_name]],
            profile_stat = stat_name,
            profile_stat_label = profile_stat_labels[[stat_name]],
            n_cell_classes = nrow(metric_df),
            pearson = safe_cor(metric_df[[three_col]], metric_df[[full_col]], "pearson"),
            spearman = safe_cor(metric_df[[three_col]], metric_df[[full_col]], "spearman"),
            median_abs_rank_delta = median(abs(rank_three - rank_full), na.rm = TRUE)
          )
        }
      )
    }
  )
}

summarise_distribution_distances <- function(cell_scores) {
  long_df <- cell_scores %>%
    pivot_longer(cols = all_of(score_columns), names_to = "metric", values_to = "score") %>%
    filter(is.finite(.data$score))

  shared_groups <- long_df %>%
    distinct(.data$platform, .data$cell_class, .data$metric) %>%
    count(.data$cell_class, .data$metric, name = "n_platforms") %>%
    filter(.data$n_platforms == 2)

  pmap_dfr(
    shared_groups,
    function(cell_class, metric, n_platforms) {
      x <- long_df %>% filter(.data$platform == "3-end RNA-seq", .data$cell_class == .env$cell_class, .data$metric == .env$metric) %>% pull(.data$score)
      y <- long_df %>% filter(.data$platform == "Full-length RNA-seq", .data$cell_class == .env$cell_class, .data$metric == .env$metric) %>% pull(.data$score)
      ks <- safe_ks(x, y)
      psd <- pooled_sd(x, y)
      tibble(
        cell_class = cell_class,
        metric = metric,
        metric_label = metric_labels[[metric]],
        n_cells_three_end = length(x),
        n_cells_full_length = length(y),
        mean_three_end = safe_numeric_mean(x),
        mean_full_length = safe_numeric_mean(y),
        median_three_end = median(x, na.rm = TRUE),
        median_full_length = median(y, na.rm = TRUE),
        positive_fraction_three_end = mean(x > 0),
        positive_fraction_full_length = mean(y > 0),
        median_delta_full_minus_three = median(y, na.rm = TRUE) - median(x, na.rm = TRUE),
        positive_fraction_delta_full_minus_three = mean(y > 0) - mean(x > 0),
        standardized_mean_difference = if_else(is.na(psd) || psd == 0, NA_real_, (mean(y, na.rm = TRUE) - mean(x, na.rm = TRUE)) / psd),
        wilcox_p_value = safe_wilcox_p(x, y)
      ) %>%
        bind_cols(ks)
    }
  ) %>%
    mutate(cell_class = factor(.data$cell_class, levels = cell_class_order)) %>%
    arrange(.data$metric, .data$cell_class)
}

make_top_cells <- function(cell_scores) {
  cell_scores %>%
    pivot_longer(cols = all_of(score_columns), names_to = "metric", values_to = "score") %>%
    filter(is.finite(.data$score)) %>%
    group_by(.data$platform, .data$cell_class, .data$metric) %>%
    arrange(desc(.data$score), .by_group = TRUE) %>%
    mutate(score_rank = row_number()) %>%
    slice_head(n = top_n) %>%
    ungroup() %>%
    mutate(metric_label = metric_labels[.data$metric]) %>%
    select("platform", "cell_class", "source_cluster", "cell_id", "metric", "metric_label", "score", "score_rank")
}

three_end_markers <- read_marker_table(three_end_marker_path, "3-end RNA-seq")
full_length_markers <- read_marker_table(full_length_marker_path, "Full-length RNA-seq")
marker_df <- bind_rows(three_end_markers, full_length_markers)
celltype_assignment <- assign_celltypes(marker_df)
write_csv(celltype_assignment, file.path(table_dir, "celltype_marker_assignments.csv"))

three_end_cells_all <- read_cell_score_table(three_end_score_path, "3-end RNA-seq", celltype_assignment)
full_length_cells_all <- read_cell_score_table(full_length_score_path, "Full-length RNA-seq", celltype_assignment)
cell_scores_all <- bind_rows(three_end_cells_all, full_length_cells_all)
cell_scores <- cell_scores_all %>%
  filter(!.data$cell_class %in% drop_cell_classes) %>%
  filter(is.finite(.data[[primary_score_column]]), .data[[primary_score_column]] != 0) %>%
  mutate(cell_class = factor(.data$cell_class, levels = cell_class_order))

if (nrow(cell_scores) == 0) {
  stop("No retained cells remained after marker-based cell-class filtering.", call. = FALSE)
}

shared_cell_ids <- intersect(three_end_cells_all$cell_id, full_length_cells_all$cell_id)
coverage_summary <- make_coverage_summary(cell_scores_all)
cluster_count_summary <- make_cluster_count_summary(cell_scores_all)
cellclass_profile <- make_cellclass_score_profile(cell_scores)
cellclass_profile_wide <- make_profile_wide(cellclass_profile)
profile_correlation_summary <- summarise_profile_correlations(cellclass_profile_wide)
distribution_distance_summary <- summarise_distribution_distances(cell_scores)
top_cells <- make_top_cells(cell_scores)

write_csv(cell_scores, file.path(table_dir, "cell_perturbation_scores_with_cellclass.csv"))
write_csv(coverage_summary, file.path(table_dir, "cell_perturbation_coverage_summary.csv"))
write_csv(cluster_count_summary, file.path(table_dir, "cluster_assignment_cell_counts.csv"))
write_csv(cellclass_profile, file.path(table_dir, "cell_perturbation_score_profile.csv"))
write_csv(cellclass_profile_wide, file.path(table_dir, "cell_perturbation_score_profile_wide.csv"))
write_csv(profile_correlation_summary, file.path(table_dir, "cell_perturbation_profile_correlation_summary.csv"))
write_csv(distribution_distance_summary, file.path(table_dir, "cell_perturbation_distribution_distance_summary.csv"))
write_csv(top_cells, file.path(table_dir, sprintf("top%d_cells_by_metric.csv", top_n)))

assignment_plot_df <- cluster_count_summary %>%
  mutate(
    source_cluster_numeric = suppressWarnings(as.integer(.data$source_cluster)),
    source_cluster_plot = factor(.data$source_cluster, levels = unique(.data$source_cluster[order(.data$platform, .data$source_cluster_numeric)])),
    cell_class = factor(.data$cell_class, levels = c(cell_class_order, drop_cell_classes))
  )

p_assignment <- ggplot(
  assignment_plot_df,
  aes(x = .data$source_cluster_plot, y = .data$cell_class, size = .data$marker_hits, fill = .data$platform)
) +
  geom_point(shape = 21, colour = "#2B2B2B", alpha = 0.86, stroke = 0.25) +
  facet_wrap(vars(.data$platform), scales = "free_x", nrow = 1) +
  scale_size_continuous(range = c(1.5, 5.6), breaks = c(1, 3, 5, 10), limits = c(0, NA)) +
  scale_fill_manual(values = platform_palette) +
  labs(x = "Original cluster", y = "Coarse cell class", size = "Marker hits", fill = NULL, title = "Marker-based coarse cell-class assignment") +
  theme_nature(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major.y = element_line(colour = "#E8E8E8", linewidth = 0.22), legend.position = "right")
save_plot_png_pdf(file.path(plot_dir, "01_marker_based_cellclass_assignment.png"), p_assignment, width = 9.0, height = 4.8)

primary_profile <- cellclass_profile_wide %>%
  filter(.data$metric == primary_score_column) %>%
  transmute(
    cell_class = factor(.data$cell_class, levels = cell_class_order),
    score_three_end = .data$median_score_three_end,
    score_full_length = .data$median_score_full_length
  ) %>%
  mutate(
    rank_three_end = min_rank(desc(.data$score_three_end)),
    rank_full_length = min_rank(desc(.data$score_full_length)),
    cell_class_label = recode(
      as.character(.data$cell_class),
      "Adnexal epithelial" = "Adnexal epi.",
      "Perivascular smooth muscle" = "Perivascular",
      "Myeloid/dendritic" = "Myeloid/DC",
      "Other immune" = "Other imm.",
      .default = as.character(.data$cell_class)
    )
  )
primary_rank_max <- max(primary_profile$rank_three_end, primary_profile$rank_full_length, na.rm = TRUE)
primary_rank_spearman <- safe_cor(primary_profile$rank_three_end, primary_profile$rank_full_length, "spearman")

p_rank <- ggplot(primary_profile, aes(x = .data$rank_three_end, y = .data$rank_full_length, colour = .data$cell_class)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#666666", linewidth = 0.4) +
  geom_point(size = 9.9, alpha = 0.92) +
  geom_text(aes(label = .data$cell_class_label), size = 4.9, nudge_x = 0.25, nudge_y = 0.15, colour = "#2B2B2B", show.legend = FALSE) +
  scale_colour_manual(values = nature_palette, drop = FALSE) +
  scale_x_reverse(breaks = seq_len(primary_rank_max), limits = c(primary_rank_max + 0.8, 0.2), expand = expansion(mult = c(0, 0))) +
  scale_y_reverse(breaks = seq_len(primary_rank_max), limits = c(primary_rank_max + 0.8, 0.2), expand = expansion(mult = c(0, 0))) +
  coord_equal(clip = "off") +
  labs(
    x = "3' RNA-seq rank (1 = highest median score)",
    y = "Full-length RNA-seq rank (1 = highest median score)",
    colour = "Cell class",
    title = "Cell-level perturbation profiles show cell-class rank agreement",
    subtitle = sprintf("Primary metric: %s; Spearman rank correlation = %.2f", metric_labels[[primary_score_column]], primary_rank_spearman)
  ) +
  theme_nature(base_size = 10) +
  theme(legend.position = "none", panel.grid.major = element_line(colour = "#E8E8E8", linewidth = 0.22), plot.margin = margin(5, 18, 5, 8))
save_plot_png_pdf(file.path(plot_dir, "02_cellclass_mean_score_rank_consistency.png"), p_rank, width = 7.2, height = 5.4)

profile_cor_plot_df <- profile_correlation_summary %>%
  mutate(
    profile_stat_label = factor(.data$profile_stat_label, levels = profile_stat_labels),
    metric_label = factor(.data$metric_label, levels = metric_labels[score_columns]),
    label = if_else(is.na(.data$spearman), "NA", sprintf("%.2f", .data$spearman))
  )

p_cor_heatmap <- ggplot(profile_cor_plot_df, aes(x = .data$profile_stat_label, y = .data$metric_label, fill = .data$spearman)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = .data$label), colour = "#2B2B2B", size = 2.8, fontface = "bold") +
  scale_fill_gradient2(low = "#3B75AF", mid = "#F7F7F7", high = "#D95F02", midpoint = 0, limits = c(-1, 1), oob = scales::squish, na.value = "#E5E5E5") +
  labs(x = "Cell-class profile statistic", y = NULL, fill = "Spearman", title = "Cell perturbation profile correlations across RNA-seq strategies") +
  theme_nature(base_size = 9) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank())
save_plot_png_pdf(file.path(plot_dir, "03_cellclass_profile_correlation_heatmap.png"), p_cor_heatmap, width = 8.4, height = 4.2)

violin_df <- cell_scores %>%
  select("platform", "cell_class", all_of(primary_score_column)) %>%
  rename(score = all_of(primary_score_column)) %>%
  filter(is.finite(.data$score)) %>%
  mutate(cell_class = factor(.data$cell_class, levels = rev(cell_class_order)))

p_violin <- ggplot(violin_df, aes(x = .data$cell_class, y = .data$score, fill = .data$platform)) +
  geom_violin(position = position_dodge(width = 0.78), scale = "width", trim = TRUE, linewidth = 0.18, colour = "#333333", alpha = 0.82) +
  coord_flip() +
  scale_fill_manual(values = platform_palette) +
  labs(x = NULL, y = metric_labels[[primary_score_column]], fill = NULL, title = "Cell-level perturbation score distributions by coarse cell class") +
  theme_nature(base_size = 9) +
  theme(panel.grid.major.x = element_line(colour = "#E8E8E8", linewidth = 0.22), legend.position = "top")
save_plot_png_pdf(file.path(plot_dir, "05_cell_level_distribution_violin_by_cellclass.png"), p_violin, width = 8.2, height = max(5.0, 0.33 * length(unique(violin_df$cell_class)) + 1.8))

shared_boxplot_cell_classes <- cell_scores %>%
  distinct(.data$platform, .data$cell_class) %>%
  count(.data$cell_class, name = "n_platforms") %>%
  filter(.data$n_platforms == 2) %>%
  pull(.data$cell_class) %>%
  as.character()

top_boxplot_cell_classes <- cell_scores %>%
  filter(as.character(.data$cell_class) %in% shared_boxplot_cell_classes) %>%
  count(.data$cell_class, name = "n_cells") %>%
  arrange(desc(.data$n_cells)) %>%
  slice_head(n = 6) %>%
  pull(.data$cell_class) %>%
  as.character()

boxplot_df <- cell_scores %>%
  filter(as.character(.data$cell_class) %in% top_boxplot_cell_classes) %>%
  select("platform", "cell_class", all_of(primary_score_column)) %>%
  rename(score = all_of(primary_score_column)) %>%
  filter(is.finite(.data$score)) %>%
  mutate(cell_class = factor(.data$cell_class, levels = top_boxplot_cell_classes))

boxplot_score_range <- range(boxplot_df$score, na.rm = TRUE)
boxplot_score_pad <- diff(boxplot_score_range) * 0.05
if (!is.finite(boxplot_score_pad) || boxplot_score_pad == 0) {
  boxplot_score_pad <- 0.1
}
boxplot_y_limits <- boxplot_score_range + c(-boxplot_score_pad, boxplot_score_pad)
boxplot_x_limits <- c(0.4, length(top_boxplot_cell_classes) + 0.6)

set.seed(1)
boxplot_point_df <- boxplot_df %>%
  mutate(
    cell_class_index = as.numeric(.data$cell_class),
    platform_offset = if_else(.data$platform == "3-end RNA-seq", -0.19, 0.19),
    x_point = .data$cell_class_index + .data$platform_offset + runif(n(), -0.09, 0.09),
    point_colour = scales::alpha(platform_palette[.data$platform], 0.32)
  )

boxplot_point_raster <- ragg::agg_capture(width = 5.2, height = 5.2, units = "in", res = 300, background = NA)
grid::grid.newpage()
grid::pushViewport(grid::viewport(xscale = boxplot_x_limits, yscale = boxplot_y_limits, clip = "on"))
grid::grid.points(
  x = grid::unit(boxplot_point_df$x_point, "native"),
  y = grid::unit(boxplot_point_df$score, "native"),
  pch = 16,
  size = grid::unit(0.65, "mm"),
  gp = grid::gpar(col = boxplot_point_df$point_colour)
)
grid::popViewport()
boxplot_point_image <- boxplot_point_raster()
grDevices::dev.off()

p_boxplot <- ggplot(boxplot_df, aes(x = .data$cell_class, y = .data$score, fill = .data$platform)) +
  annotation_custom(
    grid::rasterGrob(boxplot_point_image, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"), interpolate = TRUE),
    xmin = -Inf,
    xmax = Inf,
    ymin = -Inf,
    ymax = Inf
  ) +
  geom_boxplot(
    position = position_dodge(width = 0.76),
    width = 0.58,
    outlier.shape = NA,
    linewidth = 0.22,
    colour = "#2B2B2B",
    alpha = 0.18
  ) +
  scale_x_discrete(expand = expansion(add = 0.6)) +
  scale_y_continuous(limits = boxplot_y_limits, expand = expansion(mult = 0)) +
  scale_fill_manual(values = platform_palette) +
  labs(
    x = "Coarse cell class",
    y = metric_labels[[primary_score_column]],
    fill = NULL,
    title = sprintf("Cell perturbation %s by shared coarse cell class", str_to_lower(primary_metric_label)),
    subtitle = "Top 6 shared cell classes by retained cell count; zero-score cells were removed"
  ) +
  theme_nature(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 8.9),
    axis.text.y = element_text(size = 11.05),
    panel.grid.major.y = element_line(colour = "#E8E8E8", linewidth = 0.22),
    legend.position = "top"
  )
save_plot_png_pdf(file.path(plot_dir, "06_cell_level_score_boxplot_by_cellclass.png"), p_boxplot, width = 5.2, height = 5.98)

distance_plot_df <- distribution_distance_summary %>%
  mutate(
    cell_class = factor(.data$cell_class, levels = rev(cell_class_order)),
    metric_label = factor(.data$metric_label, levels = metric_labels[score_columns]),
    label = if_else(is.na(.data$ks_statistic), "NA", sprintf("%.2f", .data$ks_statistic))
  )

p_distance <- ggplot(distance_plot_df, aes(x = .data$metric_label, y = .data$cell_class, fill = .data$ks_statistic)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = .data$label), colour = "#2B2B2B", size = 2.6, fontface = "bold") +
  scale_fill_gradient(low = "#F7F7F7", high = "#D95F02", limits = c(0, 1), oob = scales::squish, na.value = "#E5E5E5") +
  labs(x = NULL, y = NULL, fill = "KS statistic", title = "Cell-level distribution differences between platforms") +
  theme_nature(base_size = 9) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank())
save_plot_png_pdf(file.path(plot_dir, "06_distribution_distance_heatmap.png"), p_distance, width = 7.8, height = max(4.4, 0.34 * length(unique(distance_plot_df$cell_class)) + 1.5))

primary_mean_cor <- profile_correlation_summary %>% filter(.data$metric == primary_score_column, .data$profile_stat == "mean_score")
primary_median_cor <- profile_correlation_summary %>% filter(.data$metric == primary_score_column, .data$profile_stat == "median_score")
retained_counts <- cell_scores %>% count(.data$platform, name = "n_retained_cells")
total_counts <- cell_scores_all %>% count(.data$platform, name = "n_total_cells")
shared_cell_classes <- intersect(
  unique(as.character(cell_scores$cell_class[cell_scores$platform == "3-end RNA-seq"])),
  unique(as.character(cell_scores$cell_class[cell_scores$platform == "Full-length RNA-seq"]))
)
largest_distribution_diffs <- distribution_distance_summary %>%
  filter(.data$metric == primary_score_column) %>%
  arrange(desc(.data$ks_statistic)) %>%
  slice_head(n = 5) %>%
  transmute(line = sprintf("  - %s: KS=%.3f; median delta full-minus-3'=%.3f", .data$cell_class, .data$ks_statistic, .data$median_delta_full_minus_three)) %>%
  pull(.data$line)
if (length(largest_distribution_diffs) == 0) {
  largest_distribution_diffs <- "  - No distribution comparisons available."
}

summary_lines <- c(
  "Cell perturbation score platform consistency summary",
  sprintf("Input directory: %s", input_dir),
  sprintf("3-end cell-score input: %s", three_end_score_path),
  sprintf("Full-length cell-score input: %s", full_length_score_path),
  sprintf("3-end marker input: %s", three_end_marker_path),
  sprintf("Full-length marker input: %s", full_length_marker_path),
  sprintf("Output directory: %s", output_dir),
  sprintf("Primary metric: %s", metric_labels[[primary_score_column]]),
  sprintf("Total cells: 3' RNA-seq n=%s; Full-length RNA-seq n=%s", format(total_counts$n_total_cells[match("3-end RNA-seq", total_counts$platform)], big.mark = ","), format(total_counts$n_total_cells[match("Full-length RNA-seq", total_counts$platform)], big.mark = ",")),
  sprintf("Retained nonzero-primary-score cells after marker-based cell-class filtering: 3' RNA-seq n=%s; Full-length RNA-seq n=%s", format(retained_counts$n_retained_cells[match("3-end RNA-seq", retained_counts$platform)], big.mark = ","), format(retained_counts$n_retained_cells[match("Full-length RNA-seq", retained_counts$platform)], big.mark = ",")),
  "Cells with primary metric equal to zero were removed before profile, distribution, and ranking analyses.",
  sprintf("Shared coarse cell classes: %d (%s)", length(shared_cell_classes), paste(sort(shared_cell_classes), collapse = ", ")),
  sprintf("Overlapping cell_id values across platforms: %d", length(shared_cell_ids)),
  "No direct cell-level pairing was performed; platform-specific clusters were first harmonized to marker-derived coarse cell classes.",
  sprintf("Primary metric mean-profile Pearson: %.3f; Spearman: %.3f", primary_mean_cor$pearson, primary_mean_cor$spearman),
  sprintf("Primary metric median-profile Pearson: %.3f; Spearman: %.3f", primary_median_cor$pearson, primary_median_cor$spearman),
  "Largest primary-metric distribution differences by KS statistic:",
  largest_distribution_diffs,
  sprintf("Tables: %s", table_dir),
  sprintf("Plots: %s", plot_dir)
)
writeLines(summary_lines, file.path(output_dir, "cell_perturbation_platform_consistency_summary.txt"), useBytes = TRUE)

message("Finished cell perturbation platform consistency analysis.")
message("Tables: ", table_dir)
message("Plots: ", plot_dir)
message("Summary: ", file.path(output_dir, "cell_perturbation_platform_consistency_summary.txt"))
