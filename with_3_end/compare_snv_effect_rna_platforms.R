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

input_dir <- if (length(args) >= 1) {
  args[[1]]
} else {
  "I:/极空间/文献插图/最终版/补充材料/snv_effect_by_cell_type"
}

output_dir <- if (length(args) >= 2) {
  args[[2]]
} else {
  file.path(input_dir, "rna_platform_consistency")
}

three_end_score_path <- env_or_default(
  "PRISMSNV_THREE_END_SNV_PATH",
  file.path(input_dir, "3-end_snv_perturbation_scores_by_celltype.csv")
)
full_length_score_path <- env_or_default(
  "PRISMSNV_FULL_LENGTH_SNV_PATH",
  file.path(input_dir, "full-length_snv_perturbation_scores_by_celltype.csv")
)
three_end_marker_path <- env_or_default(
  "PRISMSNV_THREE_END_MARKER_PATH",
  file.path(input_dir, "3-end_cell_cluster_marker_genes_top10.csv")
)
full_length_marker_path <- env_or_default(
  "PRISMSNV_FULL_LENGTH_MARKER_PATH",
  file.path(input_dir, "full-length_cell_cluster_marker_genes_top10.csv")
)

score_columns <- "perturb_euclidean_distance"
primary_score_column <- "perturb_euclidean_distance"
drop_cell_classes <- c("Unassigned", "Low information")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
table_dir <- file.path(output_dir, "tables")
plot_dir <- file.path(output_dir, "plots")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

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
  "NK cell" = c(
    "GNLY", "NKG7", "KLRD1", "KLRK1", "CTSW", "GZMB", "IL2RB", "STAT4", "SYTL3"
  ),
  "Mast cell" = c(
    "CPA3", "TPSB2", "TPSAB1", "MS4A2", "KIT", "HDC", "HPGDS", "HPGD",
    "SLC18A2", "SLC24A3", "IL18R1", "IL1RL1", "GATA2", "VWA5A"
  ),
  "Other immune" = c(
    "PTPRC", "CD52", "SRGN", "SAMSN1", "ARHGAP15", "IRF7", "IL3RA", "RHEX",
    "SEL1L3", "SLC15A4", "PAG1", "STK4", "RUBCN", "LPXN"
  ),
  "Melanocyte" = c(
    "MLANA", "PMEL", "DCT", "TYR", "TYRP1", "MITF", "TRPM1", "KIT", "FMN1",
    "LINC01320", "GCNT2", "INPP4B", "PLEKHA5", "SOX5", "CAPN3"
  ),
  "Adnexal epithelial" = c(
    "KRT7", "KRT8", "KRT18", "KRT19", "CLDN4", "ELF3", "AZGP1", "PIP", "DCD",
    "SCGB1B2P", "MUCL1", "SLC12A2", "TSPAN8", "ESRRG", "PLA2R1", "ZG16B", "NDRG2"
  ),
  "Schwann" = c(
    "SOX10", "ERBB3", "PMP22", "MPZ", "PLP1", "CDH19", "GPM6B", "NTM", "L1CAM",
    "NRXN1", "NRXN3", "CADM1", "GAS7", "FRMD5"
  ),
  "Sebocyte" = c(
    "FASN", "FADS2", "APOE", "APOC1", "ACSBG1", "ACSL1", "ALOX15B", "HMGCS1",
    "FDPS", "THRSP", "CYB5A", "MGST1"
  ),
  "Low information" = c(
    "RPS3", "RPS6", "RPS7", "RPS18", "RPS21", "RPS23", "RPS27", "RPL10A", "RPL12",
    "RPL17", "RPL23A", "RPL31", "RPL32", "RPL34", "RPL35A", "RPL36", "RPL36A",
    "RPL39", "MT-CO1", "MT-CO2", "MT-CO3", "MT-ND1", "MT-ND2", "HSPA1A", "HSPA1B",
    "UBC", "GAS5"
  )
)

cell_class_order <- setdiff(names(marker_panels), drop_cell_classes)
metric_labels <- c(
  "perturb_euclidean_distance" = "Perturbation Euclidean distance"
)
primary_metric_label <- unname(metric_labels[[primary_score_column]])

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

check_file_exists <- function(path) {
  if (!file.exists(path)) {
    stop("Required input file does not exist: ", path, call. = FALSE)
  }
}

walk(
  c(three_end_score_path, full_length_score_path, three_end_marker_path, full_length_marker_path),
  check_file_exists
)

obsolete_plot_files <- list.files(
  plot_dir,
  pattern = "^(0[4-9]|[1-9][0-9])_.*\\.(png|pdf)$",
  full.names = TRUE
)
if (length(obsolete_plot_files) > 0) {
  unlink(obsolete_plot_files, force = TRUE)
}

obsolete_table_files <- list.files(
  table_dir,
  pattern = paste0(
    "^(coarse_cellclass_gene_scores|exact_snv_platform_correlation_summary|",
    "gene_level_platform_correlation_summary|shared_gene_cellclass_pairs_(long|wide)|",
    "shared_snv_cellclass_pairs_(long|wide)|top[0-9]+_(gene_overlap_summary|",
    "rank_consistent_genes_by_cellclass|rank_consistent_snvs_by_cellclass|snv_overlap_summary))\\.csv$"
  ),
  full.names = TRUE
)
if (length(obsolete_table_files) > 0) {
  unlink(obsolete_table_files, force = TRUE)
}

obsolete_enrichment_dir <- file.path(output_dir, "enrichment")
if (dir.exists(obsolete_enrichment_dir)) {
  unlink(obsolete_enrichment_dir, recursive = TRUE, force = TRUE)
}

normalize_gene <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    str_replace("\\..*$", "") %>%
    str_to_upper()
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) {
    y
  } else {
    x
  }
}

safe_numeric_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  mean(x)
}

first_non_empty <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & str_trim(x) != "" & str_trim(x) != "."]
  if (length(x) == 0) {
    return(NA_character_)
  }
  x[[1]]
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

theme_nature <- function(base_size = 10) {
  axis_col <- "#2B2B2B"
  bg_col <- "white"
  strip_col <- "#F4F4F4"

  theme_bw(base_size = base_size) %+replace%
    theme(
      text = element_text(colour = axis_col, family = "sans"),
      plot.background = element_rect(fill = bg_col, colour = NA),
      panel.background = element_rect(fill = bg_col, colour = NA),
      panel.border = element_rect(fill = NA, colour = axis_col, linewidth = 0.4),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = axis_col, linewidth = 0.35),
      axis.ticks = element_line(colour = axis_col, linewidth = 0.35),
      axis.text = element_text(colour = axis_col, size = base_size * 0.82),
      axis.title = element_text(colour = axis_col, face = "bold", size = base_size * 0.95),
      plot.title = element_text(colour = axis_col, face = "bold", hjust = 0, size = base_size * 1.05),
      plot.subtitle = element_text(colour = axis_col, hjust = 0, size = base_size * 0.86),
      legend.background = element_rect(fill = bg_col, colour = NA),
      legend.key = element_rect(fill = bg_col, colour = NA),
      legend.key.size = unit(0.55, "lines"),
      legend.title = element_text(colour = axis_col, face = "bold", size = base_size * 0.82),
      legend.text = element_text(colour = axis_col, size = base_size * 0.74),
      strip.background = element_rect(fill = strip_col, colour = axis_col, linewidth = 0.3),
      strip.text = element_text(colour = axis_col, face = "bold", size = base_size * 0.78),
      plot.margin = margin(5, 6, 5, 5)
    )
}

save_plot_png_pdf <- function(filename_png, plot, width, height, dpi = 300, limitsize = FALSE) {
  ggsave(
    filename = filename_png,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white",
    limitsize = limitsize
  )
  ggsave(
    filename = sub("\\.png$", ".pdf", filename_png, ignore.case = TRUE),
    plot = plot,
    width = width,
    height = height,
    device = grDevices::pdf,
    bg = "white",
    limitsize = limitsize
  )
}

read_marker_table <- function(path, platform_label) {
  marker_df <- read_csv(
    path,
    col_types = cols(.default = col_character()),
    show_col_types = FALSE
  )
  names(marker_df) <- str_to_lower(names(marker_df))
  required_cols <- c("group", "names")
  missing_cols <- setdiff(required_cols, names(marker_df))
  if (length(missing_cols) > 0) {
    stop(
      "Marker file is missing required columns: ",
      path,
      " | ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
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
      tibble(
        cell_class = cell_class,
        marker_hits = length(matched),
        matched_marker_genes = paste(matched, collapse = ";")
      )
    }
  ) %>%
    mutate(panel_priority = match(.data$cell_class, names(marker_panels))) %>%
    arrange(desc(.data$marker_hits), .data$panel_priority)

  top <- score_df %>% slice(1)
  second_hits <- if (nrow(score_df) >= 2) score_df$marker_hits[[2]] else 0L

  if (is.na(top$marker_hits) || top$marker_hits <= 0) {
    return(tibble(
      cell_class = "Unassigned",
      marker_hits = 0L,
      second_marker_hits = second_hits,
      assignment_margin = 0L,
      matched_marker_genes = NA_character_
    ))
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
    summarise(
      n_marker_genes = n_distinct(.data$gene),
      marker_genes = paste(unique(.data$gene), collapse = ";"),
      .groups = "drop"
    ) %>%
    mutate(assignment = map(str_split(.data$marker_genes, ";"), assign_one_celltype)) %>%
    unnest("assignment")
}

read_score_table <- function(path, platform_label, assignment_df) {
  score_df <- read_csv(
    path,
    col_types = cols(.default = col_character()),
    show_col_types = FALSE
  )
  required_cols <- c("celltype", "SNV", score_columns)
  missing_cols <- setdiff(required_cols, names(score_df))
  if (length(missing_cols) > 0) {
    stop(
      "Score file is missing required columns: ",
      path,
      " | ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  score_df %>%
    mutate(
      platform = platform_label,
      celltype = as.character(.data$celltype),
      SNV = as.character(.data$SNV),
      across(all_of(score_columns), ~ suppressWarnings(as.numeric(.x)))
    ) %>%
    left_join(
      assignment_df %>%
        select(
          "platform",
          "celltype",
          "cell_class",
          "marker_hits",
          "matched_marker_genes"
        ),
      by = c("platform", "celltype")
    ) %>%
    mutate(cell_class = if_else(is.na(.data$cell_class), "Unassigned", .data$cell_class)) %>%
    filter(!.data$cell_class %in% drop_cell_classes) %>%
    group_by(.data$platform, .data$cell_class, .data$SNV) %>%
    summarise(
      perturb_euclidean_distance = safe_numeric_mean(.data$perturb_euclidean_distance),
      n_source_clusters = n_distinct(.data$celltype),
      source_celltypes = paste(sort(unique(.data$celltype)), collapse = ";"),
      gene_symbol = first_non_empty(.data$Gene),
      variant_func = first_non_empty(.data$Func),
      .groups = "drop"
    )
}


make_cellclass_summary <- function(score_df) {
  score_df %>%
    group_by(.data$platform, .data$cell_class) %>%
    summarise(
      n_unique_snvs = n_distinct(.data$SNV),
      n_snv_cellclass_rows = n(),
      source_clusters = paste(sort(unique(unlist(str_split(.data$source_celltypes, ";")))), collapse = ";"),
      .groups = "drop"
    ) %>%
    arrange(.data$platform, factor(.data$cell_class, levels = cell_class_order))
}

make_cellclass_profile <- function(score_df) {
  score_df %>%
    pivot_longer(
      cols = all_of(score_columns),
      names_to = "metric",
      values_to = "score"
    ) %>%
    filter(is.finite(.data$score)) %>%
    group_by(.data$platform, .data$cell_class, .data$metric) %>%
    summarise(
      n_snv_rows = n(),
      n_unique_snvs = n_distinct(.data$SNV),
      mean_score = mean(.data$score),
      median_score = median(.data$score),
      p90_score = as.numeric(quantile(.data$score, 0.90, names = FALSE)),
      p95_score = as.numeric(quantile(.data$score, 0.95, names = FALSE)),
      top50_mean_score = mean(head(sort(.data$score, decreasing = TRUE), min(50L, n()))),
      .groups = "drop"
    ) %>%
    mutate(
      cell_class = factor(.data$cell_class, levels = cell_class_order),
      metric_label = metric_labels[.data$metric]
    ) %>%
    arrange(.data$metric, .data$cell_class, .data$platform)
}

make_cellclass_profile_wide <- function(profile_df) {
  inner_join(
    profile_df %>% filter(.data$platform == "3-end RNA-seq"),
    profile_df %>% filter(.data$platform == "Full-length RNA-seq"),
    by = c("cell_class", "metric", "metric_label"),
    suffix = c("_three_end", "_full_length")
  )
}

summarise_cellclass_profile_correlations <- function(profile_wide_df) {
  profile_stats <- c("mean_score", "median_score", "p90_score", "p95_score", "top50_mean_score")

  map_dfr(
    score_columns,
    function(metric_name) {
      metric_df <- profile_wide_df %>% filter(.data$metric == metric_name)
      map_dfr(
        profile_stats,
        function(stat_name) {
          three_col <- paste0(stat_name, "_three_end")
          full_col <- paste0(stat_name, "_full_length")
          tibble(
            metric = metric_name,
            metric_label = metric_labels[[metric_name]],
            profile_stat = stat_name,
            profile_stat_label = recode(
              stat_name,
              "mean_score" = "Mean",
              "median_score" = "Median",
              "p90_score" = "P90",
              "p95_score" = "P95",
              "top50_mean_score" = "Top 50 mean"
            ),
            n_cell_classes = nrow(metric_df),
            pearson = safe_cor(metric_df[[three_col]], metric_df[[full_col]], "pearson", min_n = 3L),
            spearman = safe_cor(metric_df[[three_col]], metric_df[[full_col]], "spearman", min_n = 3L)
          )
        }
      )
    }
  )
}


three_end_markers <- read_marker_table(three_end_marker_path, "3-end RNA-seq")
full_length_markers <- read_marker_table(full_length_marker_path, "Full-length RNA-seq")
marker_df <- bind_rows(three_end_markers, full_length_markers)
celltype_assignment <- assign_celltypes(marker_df)

write_csv(celltype_assignment, file.path(table_dir, "celltype_marker_assignments.csv"))

three_end_scores <- read_score_table(three_end_score_path, "3-end RNA-seq", celltype_assignment)
full_length_scores <- read_score_table(full_length_score_path, "Full-length RNA-seq", celltype_assignment)
combined_scores <- bind_rows(three_end_scores, full_length_scores) %>%
  mutate(cell_class = factor(.data$cell_class, levels = cell_class_order))

write_csv(make_cellclass_summary(combined_scores), file.path(table_dir, "cellclass_score_coverage.csv"))
write_csv(combined_scores, file.path(table_dir, "coarse_cellclass_snv_scores.csv"))

cellclass_profile <- make_cellclass_profile(combined_scores)
cellclass_profile_wide <- make_cellclass_profile_wide(cellclass_profile)
cellclass_profile_correlation <- summarise_cellclass_profile_correlations(cellclass_profile_wide)

write_csv(cellclass_profile, file.path(table_dir, "cellclass_score_profile.csv"))
write_csv(cellclass_profile_wide, file.path(table_dir, "cellclass_score_profile_wide.csv"))
write_csv(cellclass_profile_correlation, file.path(table_dir, "cellclass_profile_correlation_summary.csv"))

cellclass_primary_mean <- cellclass_profile_correlation %>%
  filter(.data$metric == primary_score_column, .data$profile_stat == "mean_score")
cellclass_primary_median <- cellclass_profile_correlation %>%
  filter(.data$metric == primary_score_column, .data$profile_stat == "median_score")

notes <- c(
  "Cross-platform SNV perturbation score consistency summary",
  sprintf("Input directory: %s", input_dir),
  sprintf("3-end SNV input: %s", three_end_score_path),
  sprintf("Full-length SNV input: %s", full_length_score_path),
  sprintf("3-end marker input: %s", three_end_marker_path),
  sprintf("Full-length marker input: %s", full_length_marker_path),
  sprintf("Output directory: %s", output_dir),
  sprintf("Coarse cell classes retained: %s", paste(sort(unique(as.character(combined_scores$cell_class))), collapse = ", ")),
  sprintf(
    "Cell-class profile %s mean-score Pearson: %.3f; Spearman: %.3f; cell classes: %s",
    metric_labels[[primary_score_column]],
    cellclass_primary_mean$pearson,
    cellclass_primary_mean$spearman,
    cellclass_primary_mean$n_cell_classes
  ),
  sprintf(
    "Cell-class profile %s median-score Pearson: %.3f; Spearman: %.3f; cell classes: %s",
    metric_labels[[primary_score_column]],
    cellclass_primary_median$pearson,
    cellclass_primary_median$spearman,
    cellclass_primary_median$n_cell_classes
  ),
  "The analysis is restricted to Euclidean-distance coarse cell-class profiles.",
  "Only plots 01-03 and their required profile tables are generated."
)
writeLines(notes, file.path(output_dir, "consistency_summary.txt"), useBytes = TRUE)

assignment_plot_df <- celltype_assignment %>%
  mutate(
    celltype_numeric = suppressWarnings(as.integer(.data$celltype)),
    celltype_plot = factor(.data$celltype, levels = unique(.data$celltype[order(.data$platform, .data$celltype_numeric)])),
    cell_class = factor(.data$cell_class, levels = c(cell_class_order, drop_cell_classes))
  )

p_assignment <- ggplot(
  assignment_plot_df,
  aes(x = .data$celltype_plot, y = .data$cell_class, size = .data$marker_hits, fill = .data$platform)
) +
  geom_point(shape = 21, colour = "#2B2B2B", alpha = 0.85, stroke = 0.25) +
  facet_wrap(vars(.data$platform), scales = "free_x", nrow = 1) +
  scale_size_continuous(range = c(1.5, 5.5), breaks = c(1, 3, 5, 10), limits = c(0, NA)) +
  scale_fill_manual(values = c("3-end RNA-seq" = "#4E79A7", "Full-length RNA-seq" = "#E15759")) +
  labs(
    x = "Original cluster",
    y = "Coarse cell class",
    size = "Marker hits",
    fill = NULL,
    title = "Marker-based coarse cell-class assignment"
  ) +
  theme_nature(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.y = element_line(colour = "#E8E8E8", linewidth = 0.22),
    legend.position = "right"
  )

save_plot_png_pdf(
  file.path(plot_dir, "01_marker_based_cellclass_assignment.png"),
  p_assignment,
  width = 9.0,
  height = 4.8
)

cellclass_profile_plot_df <- cellclass_profile_wide %>%
  filter(.data$metric == primary_score_column) %>%
  transmute(
    cell_class = factor(.data$cell_class, levels = cell_class_order),
    cell_class_label = recode(
      as.character(.data$cell_class),
      "Adnexal epithelial" = "Adnexal epi.",
      "Perivascular smooth muscle" = "Perivascular",
      "Myeloid/dendritic" = "Myeloid/DC",
      "Plasma cell" = "Plasma",
      "Other immune" = "Other imm.",
      .default = as.character(.data$cell_class)
    ),
    score_three_end = .data$mean_score_three_end,
    score_full_length = .data$mean_score_full_length
  ) %>%
  mutate(
    rank_three_end = min_rank(desc(.data$score_three_end)),
    rank_full_length = min_rank(desc(.data$score_full_length))
  ) %>%
  mutate(
    label_x = case_when(
      as.character(.data$cell_class) == "Plasma cell" ~ .data$rank_three_end + 0.35,
      as.character(.data$cell_class) == "Adnexal epithelial" ~ .data$rank_three_end + 0.35,
      as.character(.data$cell_class) == "Keratinocyte" ~ .data$rank_three_end + 0.55,
      as.character(.data$cell_class) == "Mast cell" ~ .data$rank_three_end - 0.55,
      as.character(.data$cell_class) == "Fibroblast" ~ .data$rank_three_end - 0.55,
      as.character(.data$cell_class) == "Myeloid/dendritic" ~ .data$rank_three_end + 0.55,
      as.character(.data$cell_class) == "Perivascular smooth muscle" ~ .data$rank_three_end + 0.55,
      as.character(.data$cell_class) == "Endothelial" ~ .data$rank_three_end + 0.55,
      as.character(.data$cell_class) == "T cell" ~ .data$rank_three_end - 0.35,
      TRUE ~ .data$rank_three_end + 0.35
    ),
    label_y = case_when(
      as.character(.data$cell_class) == "Plasma cell" ~ .data$rank_full_length + 0.25,
      as.character(.data$cell_class) == "Adnexal epithelial" ~ .data$rank_full_length - 0.25,
      as.character(.data$cell_class) == "Keratinocyte" ~ .data$rank_full_length - 0.35,
      as.character(.data$cell_class) == "Mast cell" ~ .data$rank_full_length - 0.35,
      as.character(.data$cell_class) == "Fibroblast" ~ .data$rank_full_length + 0.35,
      as.character(.data$cell_class) == "Myeloid/dendritic" ~ .data$rank_full_length - 0.35,
      as.character(.data$cell_class) == "Perivascular smooth muscle" ~ .data$rank_full_length + 0.35,
      as.character(.data$cell_class) == "Endothelial" ~ .data$rank_full_length - 0.35,
      as.character(.data$cell_class) == "T cell" ~ .data$rank_full_length + 0.35,
      TRUE ~ .data$rank_full_length + 0.35
    )
  )
cellclass_rank_max <- max(
  cellclass_profile_plot_df$rank_three_end,
  cellclass_profile_plot_df$rank_full_length,
  na.rm = TRUE
)

p_cellclass_profile <- ggplot(
  cellclass_profile_plot_df,
  aes(x = .data$rank_three_end, y = .data$rank_full_length, colour = .data$cell_class)
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#666666", linewidth = 0.4) +
  geom_point(size = 9.5, alpha = 0.92) +
  geom_text(
    aes(x = .data$label_x, y = .data$label_y, label = .data$cell_class_label),
    size = 7.35,
    colour = "#2B2B2B",
    check_overlap = FALSE,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = nature_palette, drop = FALSE) +
  scale_x_reverse(
    breaks = seq_len(cellclass_rank_max),
    limits = c(cellclass_rank_max + 0.8, 0.2),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_reverse(
    breaks = seq_len(cellclass_rank_max),
    limits = c(cellclass_rank_max + 0.8, 0.2),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_equal(clip = "off") +
  labs(
    x = sprintf("3' RNA-seq rank (1 = highest mean %s)", str_to_lower(primary_metric_label)),
    y = sprintf("Full-length RNA-seq rank (1 = highest mean %s)", str_to_lower(primary_metric_label)),
    colour = "Cell class",
    title = "Coarse cell-class perturbation ranks are concordant",
    subtitle = sprintf(
      "Across retained cell classes: Spearman rank correlation = %.2f",
      cellclass_primary_mean$spearman
    )
  ) +
  theme_nature(base_size = 10) +
  theme(
    legend.position = "none",
    panel.grid.major = element_line(colour = "#E8E8E8", linewidth = 0.22),
    plot.margin = margin(5, 18, 5, 8)
  )

save_plot_png_pdf(
  file.path(plot_dir, "02_cellclass_profile_mean_score_consistency.png"),
  p_cellclass_profile,
  width = 7.0,
  height = 5.2
)

profile_snv_counts <- combined_scores %>%
  group_by(.data$platform) %>%
  summarise(n_unique_snvs = n_distinct(.data$SNV), .groups = "drop")
profile_snv_caption <- sprintf(
  "SNVs included after cell-class filtering: 3' RNA-seq n=%s; Full-length RNA-seq n=%s.",
  format(profile_snv_counts$n_unique_snvs[match("3-end RNA-seq", profile_snv_counts$platform)] %||% 0L, big.mark = ","),
  format(profile_snv_counts$n_unique_snvs[match("Full-length RNA-seq", profile_snv_counts$platform)] %||% 0L, big.mark = ",")
)

profile_cor_heatmap_df <- cellclass_profile_correlation %>%
  mutate(
    profile_stat_label = factor(
      .data$profile_stat_label,
      levels = c("Mean", "Median", "P90", "P95", "Top 50 mean")
    ),
    metric_label = factor(.data$metric_label, levels = metric_labels[score_columns]),
    label = if_else(is.na(.data$pearson), "NA", sprintf("%.2f", .data$pearson))
  )

p_profile_cor <- ggplot(
  profile_cor_heatmap_df,
  aes(x = .data$profile_stat_label, y = .data$metric_label, fill = .data$pearson)
) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = .data$label), colour = "#2B2B2B", size = 4.5, fontface = "bold") +
  scale_fill_gradient2(
    low = "#3B75AF",
    mid = "#F7F7F7",
    high = "#D95F02",
    midpoint = 0,
    limits = c(-1, 1),
    oob = scales::squish,
    na.value = "#E5E5E5"
  ) +
  labs(
    x = "Cell-class profile statistic",
    y = NULL,
    fill = "Pearson",
    title = "Cell-class profile correlations across sequencing strategies",
    caption = profile_snv_caption
  ) +
  theme_nature(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank())

save_plot_png_pdf(
  file.path(plot_dir, "03_cellclass_profile_correlation_heatmap.png"),
  p_profile_cor,
  width = 6.2,
  height = 3.6
)

message("Finished cross-platform consistency analysis.")
message("Tables: ", table_dir)
message("Plots: ", plot_dir)
message("Summary: ", file.path(output_dir, "consistency_summary.txt"))
