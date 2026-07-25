#!/usr/bin/env Rscript

# RStudio-friendly script: edit the constants below, then run `source()`.
#
# This script reads `snv_perturbation_scores_by_celltype.csv`, ranks SNVs by
# per-SNV max(score_euclidean) (tie-break: mean score, then distinct cell type
# count), keeps the top 50 SNVs, and draws the circular ring body first:
# - outer labels: SNV IDs
# - outer numeric ring: perturbation score lollipops
# - AlphaGenome boxplot rings: RNA-seq and ATAC-seq absolute quantile summaries
# - inner categorical rings: dominant Func, CLNSIG, VEP impact, and
#   MutationTaster prediction when present
# - separate top-left-style summary panel for later manual composition

# ==============================
# User Config (edit these)
# ==============================
# Run the corresponding upstream analysis in advance to generate this result.
INPUT_CSV <- "/path/to/input/snv_perturbation_scores_by_celltype.csv"
OUTPUT_DIR <- NA_character_ # Default: <input_dir>/top50_snv_circular
SCORE_COLUMN <- "score_cosine_sim"
TOP_N <- 100L
OUTPUT_STEM <- "top50_snv_circular"

# Run the AlphaGenome analysis in advance to generate these results.
ALPHAGENOME_BOXPLOT_STATS_CSV <- "/path/to/alphagenome/alphagenome_strict_skin_boxplot_stats.csv"
ALPHAGENOME_ATAC_METADATA_CSV <- "/path/to/alphagenome/alphagenome_atac_strict_skin_metadata.csv"
# Run snpEff in advance to generate this result.
SNPEFF_IMPACT_SUMMARY_TSV <- "/path/to/annotations/snpeff.impact_summary.tsv"
# Run RegulationSpotter in advance to generate these results.
REGULATIONSPOTTER_ALL_VARS_TSV <- "/path/to/regulationspotter/all_vars.tsv"
REGULATIONSPOTTER_SPOTTER_TSV <- "/path/to/regulationspotter/all_results_spotter.tsv"
REGULATIONSPOTTER_TASTER_TSV <- "/path/to/regulationspotter/all_results_taster.tsv"
# Run VEP in advance to generate this result.
VEP_ANNOTATION_CSV <- "/path/to/annotations/vep_output.csv"
# Run MutationTaster in advance to generate this result.
MUTATIONTASTER_ANNOTATION_TSV <- "/path/to/annotations/MutationTaster.tsv"

PLOT_WIDTH <- 11
PLOT_HEIGHT <- 11
SUMMARY_PANEL_WIDTH <- 10
SUMMARY_PANEL_HEIGHT <- 5.5
PLOT_DPI <- 600
PLOT_TITLE <- NA_character_ # Keep NA for a clean ring-only panel.
SHOW_LEGENDS <- FALSE # First draw the circular body; switch to TRUE if needed.
SCORE_RING_UPPER_QUANTILE <- 0.99
SUMMARY_PANEL_LEGEND_PANEL_WIDTH <- 5.0
SUMMARY_PANEL_LEGEND_RIGHT_MARGIN <- 8.2
SUMMARY_PANEL_LEGEND_X_OFFSETS <- c(0.85, 3.10)
SUMMARY_PANEL_LEGEND_TEXT_CEX <- 0.52
SUMMARY_PANEL_LEGEND_Y_INTERSP <- 0.82
SUMMARY_PANEL_LEGEND_BLOCK_GAP <- 0.02
SUMMARY_PANEL_NON_TOP_POINT_COLOR <- "#BDBDBD"
SUMMARY_PANEL_NUMERIC_BOX_WIDTH <- 0.4
SUMMARY_PANEL_NUMERIC_POINT_WIDTH <- 0.14
SUMMARY_PANEL_CATEGORY_BAR_WIDTH <- 0.60
SUMMARY_PANEL_NUMERIC_COLUMN_WIDTH_FACTOR <- 2.8
SUMMARY_PANEL_CATEGORY_LABEL_CEX <- 1.04
SHOW_OUTER_RING_VALUE_LABELS <- TRUE
OUTER_RING_VALUE_AXIS_SIDE <- "left"
OUTER_RING_VALUE_LABEL_CEX <- 0.42
OUTER_RING_VALUE_LABEL_COLOR <- "#4D4D4D"

START_DEGREE <- 90
SECTOR_GAP_DEGREE <- 1
FINAL_GAP_DEGREE <- 90 # Reserve the upper-left quadrant for composite panels.
LABEL_CEX <- 0.3
FIGURE_FONT_FAMILY <- "sans"

SCORE_RING_COLOR <- "#7B4AB7"
SCORE_RING_BASELINE_COLOR <- "#D9C7F1"
ALPHAGENOME_RNA_BOX_COLOR <- "#2B8CBE"
ALPHAGENOME_ATAC_BOX_COLOR <- "#238B45"
UNKNOWN_FUNC_COLOR <- "#BDBDBD"
UNKNOWN_CATEGORY_COLOR <- "#D9D9D9"
VEP_IMPACT_COLORS <- c(
  HIGH = "#B2182B",
  MODERATE = "#EF8A62",
  LOW = "#67A9CF",
  MODIFIER = "#D9D9D9",
  Unknown = UNKNOWN_CATEGORY_COLOR
)
MUTATIONTASTER_PREDICTION_COLORS <- c(
  "Deleterious (ClinVar)" = "#B2182B",
  "Deleterious (fs/PTC)" = "#D6604D",
  Deleterious = "#F4A582",
  "Benign (auto)" = "#92C5DE",
  Benign = "#2166AC",
  Unknown = UNKNOWN_CATEGORY_COLOR
)
RING_BORDER_COLOR <- "#222222"
RING_BACKGROUND_COLOR <- "#FFFFFF"
LABEL_TRACK_HEIGHT <- 0.25
NUMERIC_TRACK_HEIGHT <- 0.11
CATEGORY_TRACK_HEIGHT <- NUMERIC_TRACK_HEIGHT / 2
SCORE_RING_TRACK_HEIGHT <- NUMERIC_TRACK_HEIGHT
ALPHAGENOME_RNA_TRACK_HEIGHT <- NUMERIC_TRACK_HEIGHT
ALPHAGENOME_ATAC_TRACK_HEIGHT <- NUMERIC_TRACK_HEIGHT
SNPEFF_IMPACT_TRACK_HEIGHT <- CATEGORY_TRACK_HEIGHT
FUNC_TRACK_HEIGHT <- CATEGORY_TRACK_HEIGHT
VEP_IMPACT_TRACK_HEIGHT <- CATEGORY_TRACK_HEIGHT
MUTATIONTASTER_TRACK_HEIGHT <- CATEGORY_TRACK_HEIGHT
ALPHAGENOME_BOX_WIDTH_FRACTION <- 0.52

require_package <- function(pkg, install_hint = NULL) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    return(invisible(TRUE))
  }
  if (is.null(install_hint)) {
    install_hint <- pkg
  }
  stop(
    sprintf(
      "Package '%s' is required. Install it first, for example with install.packages('%s').",
      pkg,
      install_hint
    ),
    call. = FALSE
  )
}

normalize_path_string <- function(path) {
  if (is.na(path) || !nzchar(path)) {
    return(path)
  }
  path <- trimws(as.character(path))
  path <- gsub("\\\\", "/", path)
  path <- sub("^([A-Za-z]):$", "\\1:/", path)
  path <- sub("^([A-Za-z]):/{2,}", "\\1:/", path)
  path
}

ensure_dir <- function(path) {
  path <- normalize_path_string(path)
  if (is.na(path) || !nzchar(path)) {
    stop("Output directory path is empty.", call. = FALSE)
  }
  if (!dir.exists(path)) {
    ok <- dir.create(path, recursive = TRUE, showWarnings = FALSE)
    if (!ok && !dir.exists(path)) {
      stop(sprintf("Failed to create directory: %s", path), call. = FALSE)
    }
  }
  if (!dir.exists(path)) {
    stop(sprintf("Directory does not exist: %s", path), call. = FALSE)
  }
  path
}

check_input_csv <- function(input_csv) {
  input_csv <- normalize_path_string(input_csv)
  if (!file.exists(input_csv)) {
    stop(sprintf("INPUT_CSV not found: %s", input_csv), call. = FALSE)
  }
  info <- file.info(input_csv)
  if (is.na(info$isdir) || info$isdir) {
    stop(sprintf("INPUT_CSV must be a file, got directory: %s", input_csv), call. = FALSE)
  }
  readable <- suppressWarnings(file.access(input_csv, 4) == 0)
  if (!isTRUE(readable)) {
    stop(sprintf("INPUT_CSV is not readable: %s", input_csv), call. = FALSE)
  }
  input_csv
}

check_optional_csv <- function(input_csv, label) {
  input_csv <- normalize_path_string(input_csv)
  if (is.na(input_csv) || !nzchar(input_csv)) {
    return(NA_character_)
  }
  if (!file.exists(input_csv)) {
    stop(sprintf("%s not found: %s", label, input_csv), call. = FALSE)
  }
  info <- file.info(input_csv)
  if (is.na(info$isdir) || info$isdir) {
    stop(sprintf("%s must be a file, got directory: %s", label, input_csv), call. = FALSE)
  }
  readable <- suppressWarnings(file.access(input_csv, 4) == 0)
  if (!isTRUE(readable)) {
    stop(sprintf("%s is not readable: %s", label, input_csv), call. = FALSE)
  }
  input_csv
}

check_output_dir <- function(output_dir) {
  output_dir <- ensure_dir(output_dir)
  writable <- suppressWarnings(file.access(output_dir, 2) == 0)
  if (!isTRUE(writable)) {
    stop(sprintf("OUTPUT_DIR is not writable: %s", output_dir), call. = FALSE)
  }

  test_path <- file.path(output_dir, sprintf(".vespa_write_test_%s.tmp", as.integer(Sys.time())))
  ok <- FALSE
  err <- NULL
  tryCatch(
    {
      ok <- file.create(test_path)
      if (isTRUE(ok)) {
        unlink(test_path)
      }
    },
    error = function(e) {
      err <<- e
    }
  )
  if (!isTRUE(ok)) {
    msg <- if (!is.null(err)) conditionMessage(err) else "file.create returned FALSE"
    stop(sprintf("OUTPUT_DIR write test failed: %s (%s)", output_dir, msg), call. = FALSE)
  }

  output_dir
}

resolve_output_dir <- function(input_path, output_dir) {
  if (!is.na(output_dir) && nzchar(output_dir)) {
    return(normalize_path_string(output_dir))
  }
  file.path(dirname(normalize_path_string(input_path)), "top50_snv_circular")
}

validate_opts <- function(opts) {
  if (!is.character(opts$input) || length(opts$input) != 1 || !nzchar(opts$input)) {
    stop("INPUT_CSV must be a non-empty string.", call. = FALSE)
  }
  if (!is.character(opts$output_stem) || length(opts$output_stem) != 1 || !nzchar(opts$output_stem)) {
    stop("OUTPUT_STEM must be a non-empty string.", call. = FALSE)
  }
  if (!is.character(opts$score_column) || length(opts$score_column) != 1 || !nzchar(opts$score_column)) {
    stop("SCORE_COLUMN must be a non-empty string.", call. = FALSE)
  }
  if (!is.finite(opts$top_n) || opts$top_n <= 0) {
    stop("TOP_N must be a positive integer.", call. = FALSE)
  }
  if (!is.finite(opts$plot_width) || opts$plot_width <= 0) {
    stop("PLOT_WIDTH must be positive.", call. = FALSE)
  }
  if (!is.finite(opts$plot_height) || opts$plot_height <= 0) {
    stop("PLOT_HEIGHT must be positive.", call. = FALSE)
  }
  if (!is.finite(opts$plot_dpi) || opts$plot_dpi <= 0) {
    stop("PLOT_DPI must be positive.", call. = FALSE)
  }
  if (!is.finite(opts$score_ring_upper_quantile) ||
      opts$score_ring_upper_quantile <= 0 ||
      opts$score_ring_upper_quantile > 1) {
    stop("SCORE_RING_UPPER_QUANTILE must be in (0, 1].", call. = FALSE)
  }
  invisible(TRUE)
}

attach_required_packages <- function() {
  required_packages <- c("circlize", "dplyr", "readr")
  for (pkg in required_packages) {
    require_package(pkg)
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  invisible(TRUE)
}

sanitize_category <- function(values, default = "Unknown") {
  values <- trimws(as.character(values))
  values[is.na(values) | !nzchar(values) | values == "." | values == "NA"] <- default
  values
}

normalize_vep_impact <- function(values) {
  values <- toupper(sanitize_category(values, default = "Unknown"))
  values[values == "UNKNOWN"] <- "Unknown"
  values
}

sanitize_binary_flag <- function(values) {
  values <- suppressWarnings(as.numeric(trimws(as.character(values))))
  values[!is.finite(values)] <- 0
  as.integer(values > 0)
}

map_regulationspotter_prediction <- function(values) {
  values <- suppressWarnings(as.integer(as.character(values)))
  labels <- rep("Unknown", length(values))
  labels[values == 1] <- "Disease-causing"
  labels[values == 21] <- "Functional (strong)"
  labels[values == 22] <- "Functional"
  labels[values == 28] <- "Non-functional"
  labels[values == 29] <- "Non-functional (strong)"
  labels[values == 99] <- "Known polymorphism"
  labels
}

map_regulationtaster_prediction <- function(values) {
  values <- suppressWarnings(as.integer(as.character(values)))
  labels <- rep("Unknown", length(values))
  labels[values == 1] <- "Disease-causing (ClinVar 5)"
  labels[values == 2] <- "Disease-causing (ClinVar 4)"
  labels[values == 3] <- "Disease-causing"
  labels[values == 95] <- "Polymorphism"
  labels[values == 99] <- "Known polymorphism"
  labels
}

complement_dna_base <- function(values) {
  values <- toupper(trimws(as.character(values)))
  out <- values
  out[values == "A"] <- "T"
  out[values == "T"] <- "A"
  out[values == "C"] <- "G"
  out[values == "G"] <- "C"
  invalid_mask <- is.na(values) | !nzchar(values) | !(values %in% c("A", "C", "G", "T"))
  out[invalid_mask] <- NA_character_
  out
}

read_regulationspotter_variant_map <- function(all_vars_tsv) {
  all_vars_tsv <- check_optional_csv(all_vars_tsv, "REGULATIONSPOTTER_ALL_VARS_TSV")
  if (is.na(all_vars_tsv) || !nzchar(all_vars_tsv)) {
    return(NULL)
  }

  message(sprintf("[INFO] Reading RegulationSpotter variant map TSV: %s", all_vars_tsv))
  variant_df <- readr::read_tsv(
    file = all_vars_tsv,
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "minimal",
    col_types = readr::cols(.default = readr::col_character())
  )
  variant_df <- as.data.frame(variant_df, stringsAsFactors = FALSE)
  required_columns <- c("var_number", "chromosome", "position", "ref", "alt")
  missing_columns <- setdiff(required_columns, colnames(variant_df))
  if (length(missing_columns) > 0) {
    stop(
      sprintf("REGULATIONSPOTTER_ALL_VARS_TSV is missing required columns: %s", paste(missing_columns, collapse = ", ")),
      call. = FALSE
    )
  }

  variant_df <- data.frame(
    var_number = suppressWarnings(as.integer(variant_df$var_number)),
    SNV = build_snv_id_from_parts(variant_df$chromosome, variant_df$position, variant_df$ref, variant_df$alt),
    SNV_reverse_complement = build_snv_id_from_parts(
      variant_df$chromosome,
      variant_df$position,
      complement_dna_base(variant_df$ref),
      complement_dna_base(variant_df$alt)
    ),
    stringsAsFactors = FALSE
  )
  variant_df <- variant_df[is.finite(variant_df$var_number), , drop = FALSE]

  primary_df <- data.frame(
    var_number = variant_df$var_number,
    SNV = variant_df$SNV,
    snv_match_type = "direct",
    snv_match_priority = 1L,
    stringsAsFactors = FALSE
  )
  reverse_df <- data.frame(
    var_number = variant_df$var_number,
    SNV = variant_df$SNV_reverse_complement,
    snv_match_type = "reverse_complement",
    snv_match_priority = 2L,
    stringsAsFactors = FALSE
  )
  variant_df <- rbind(primary_df, reverse_df)
  variant_df <- variant_df[!is.na(variant_df$SNV) & nzchar(variant_df$SNV), , drop = FALSE]
  variant_df <- variant_df[!duplicated(variant_df), , drop = FALSE]
  variant_df
}

build_snv_id_from_parts <- function(chr, position, ref, alt) {
  chr <- trimws(as.character(chr))
  position <- trimws(as.character(position))
  ref <- toupper(trimws(as.character(ref)))
  alt <- toupper(trimws(as.character(alt)))
  chr <- sub("^chr", "", chr, ignore.case = TRUE)
  chr[chr %in% c("M", "MT", "MITO")] <- "M"
  invalid_mask <- is.na(chr) | !nzchar(chr) |
    is.na(position) | !nzchar(position) |
    is.na(ref) | !nzchar(ref) |
    is.na(alt) | !nzchar(alt) |
    chr %in% c(".", "NA") |
    position %in% c(".", "NA") |
    ref %in% c(".", "NA") |
    alt %in% c(".", "NA")
  snv <- paste0("chr", chr, ":", position, "_", ref, ">", alt)
  snv[invalid_mask] <- NA_character_
  snv
}

extract_snv_chromosome <- function(snv) {
  chrom <- sub(":.*$", "", trimws(as.character(snv)))
  chrom <- sub("^chr", "", chrom, ignore.case = TRUE)
  toupper(chrom)
}

is_allowed_plot_snv <- function(snv) {
  extract_snv_chromosome(snv) %in% c(as.character(1:12), "X", "Y")
}

rescale_to_unit <- function(values, upper_quantile = 1) {
  values <- as.numeric(values)
  finite_values <- values[is.finite(values)]
  scaled <- rep(NA_real_, length(values))
  if (length(finite_values) == 0) {
    return(scaled)
  }
  upper_quantile <- as.numeric(upper_quantile)[1]
  if (!is.finite(upper_quantile) || upper_quantile <= 0 || upper_quantile > 1) {
    stop("upper_quantile must be in (0, 1].", call. = FALSE)
  }
  min_value <- min(finite_values)
  max_value <- if (upper_quantile < 1) {
    stats::quantile(finite_values, probs = upper_quantile, na.rm = TRUE, names = FALSE)
  } else {
    max(finite_values)
  }
  if (!is.finite(max_value) || max_value <= min_value) {
    scaled[is.finite(values)] <- 1
    return(scaled)
  }
  clipped_values <- pmin(values[is.finite(values)], max_value)
  scaled[is.finite(values)] <- (clipped_values - min_value) / (max_value - min_value)
  scaled
}

prepare_input_table <- function(df, score_column) {
  required_columns <- c("celltype", "SNV", "Func", score_column)
  missing_columns <- setdiff(required_columns, colnames(df))
  if (length(missing_columns) > 0) {
    stop(
      sprintf("Input CSV is missing required columns: %s", paste(missing_columns, collapse = ", ")),
      call. = FALSE
    )
  }

  df$celltype <- trimws(as.character(df$celltype))
  df$SNV <- trimws(as.character(df$SNV))
  df$Func <- sanitize_category(df$Func, default = "Unknown")
  if ("ExonicFunc" %in% colnames(df)) {
    df$ExonicFunc <- sanitize_category(df$ExonicFunc, default = "Unknown")
  }
  if ("CLNSIG" %in% colnames(df)) {
    df$CLNSIG <- sanitize_category(df$CLNSIG, default = "Unknown")
  }
  df$plot_score <- suppressWarnings(as.numeric(df[[score_column]]))

  invalid_mask <- is.na(df$SNV) | !nzchar(df$SNV) |
    is.na(df$celltype) | !nzchar(df$celltype) |
    !is.finite(df$plot_score)

  if (any(invalid_mask)) {
    message(sprintf(
      "[WARN] Removing %d rows with missing/invalid SNV, celltype, or %s.",
      sum(invalid_mask),
      score_column
    ))
    df <- df[!invalid_mask, , drop = FALSE]
  }

  if (nrow(df) == 0) {
    stop("No valid rows remain after input cleaning.", call. = FALSE)
  }

  allowed_snv_mask <- is_allowed_plot_snv(df$SNV)
  if (any(!allowed_snv_mask)) {
    message(sprintf(
      "[INFO] Removing %d rows outside chr1-12, chrX, and chrY before SNV ranking.",
      sum(!allowed_snv_mask)
    ))
    df <- df[allowed_snv_mask, , drop = FALSE]
  }

  if (nrow(df) == 0) {
    stop("No valid rows remain after chromosome filtering; kept only chr1-12, chrX, and chrY.", call. = FALSE)
  }

  df
}

build_dominant_category_summary <- function(df, category_col, output_col) {
  category_summary <- df |>
    dplyr::group_by(SNV, category_value = .data[[category_col]]) |>
    dplyr::summarise(
      category_row_count = dplyr::n(),
      category_top_score = max(plot_score),
      .groups = "drop"
    ) |>
    dplyr::arrange(SNV, dplyr::desc(category_row_count), dplyr::desc(category_top_score), category_value) |>
    dplyr::group_by(SNV) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup()

  out <- data.frame(
    SNV = category_summary$SNV,
    category_value = sanitize_category(category_summary$category_value, default = "Unknown"),
    stringsAsFactors = FALSE
  )
  colnames(out)[2] <- output_col
  out
}

build_ranked_snv_summary <- function(df, score_upper_quantile = 1) {
  snv_summary <- df |>
    dplyr::group_by(SNV) |>
    dplyr::summarise(
      max_score = max(plot_score),
      mean_score = mean(plot_score),
      distinct_celltype_count = dplyr::n_distinct(celltype),
      total_rows = dplyr::n(),
      .groups = "drop"
    )

  ranked_summary <- snv_summary
  category_specs <- c(
    dominant_func = "Func"
  )
  for (output_col in names(category_specs)) {
    input_col <- unname(category_specs[[output_col]])
    if (input_col %in% colnames(df)) {
      ranked_summary <- ranked_summary |>
        dplyr::left_join(
          build_dominant_category_summary(df, input_col, output_col),
          by = "SNV"
        )
    }
  }

  ranked_summary <- ranked_summary |>
    dplyr::arrange(dplyr::desc(max_score), dplyr::desc(mean_score), dplyr::desc(distinct_celltype_count), SNV) |>
    dplyr::mutate(rank = dplyr::row_number())

  if (nrow(ranked_summary) == 0) {
    stop("No SNVs available after aggregation.", call. = FALSE)
  }

  ranked_summary |>
    dplyr::mutate(
      score_norm = rescale_to_unit(max_score, upper_quantile = score_upper_quantile)
    )
}

build_top_snv_summary <- function(df, top_n, score_upper_quantile = 1) {
  ranked_summary <- build_ranked_snv_summary(df, score_upper_quantile = score_upper_quantile)
  top_n <- min(as.integer(top_n), nrow(ranked_summary))
  ranked_summary |>
    dplyr::slice_head(n = top_n)
}

make_category_palette <- function(values, unknown_color = UNKNOWN_CATEGORY_COLOR, palette = "Set 3", fixed_colors = NULL) {
  values <- sanitize_category(values, default = "Unknown")
  category_levels <- sort(unique(values))
  if (length(category_levels) == 0) {
    return(setNames(character(0), character(0)))
  }
  if (!is.null(fixed_colors) && length(fixed_colors) > 0) {
    fixed_colors <- fixed_colors[!is.na(names(fixed_colors)) & nzchar(names(fixed_colors))]
    fixed_levels <- intersect(names(fixed_colors), category_levels)
    category_levels <- c(fixed_levels, setdiff(category_levels, fixed_levels))
  }
  category_colors <- setNames(grDevices::hcl.colors(length(category_levels), palette = palette), category_levels)
  if (!is.null(fixed_colors) && length(fixed_colors) > 0) {
    fixed_levels <- intersect(names(fixed_colors), names(category_colors))
    category_colors[fixed_levels] <- fixed_colors[fixed_levels]
  }
  if ("Unknown" %in% names(category_colors)) {
    category_colors["Unknown"] <- unknown_color
  }
  category_colors
}

build_category_tracks <- function(summary_df, opts) {
  track_specs <- list(
    list(
      column = "snpeff_impact",
      label = "snpEff impact",
      track_height = opts$snpeff_impact_track_height,
      unknown_color = UNKNOWN_CATEGORY_COLOR,
      palette = "Set 2",
      fixed_colors = VEP_IMPACT_COLORS
    ),
    list(
      column = "dominant_func",
      label = "Func",
      track_height = opts$func_track_height,
      unknown_color = UNKNOWN_FUNC_COLOR,
      palette = "Set 3"
    ),
    list(
      column = "vep_impact",
      label = "VEP impact",
      track_height = opts$vep_impact_track_height,
      unknown_color = UNKNOWN_CATEGORY_COLOR,
      palette = "Set 2",
      fixed_colors = VEP_IMPACT_COLORS
    ),
    list(
      column = "mutationtaster_prediction",
      label = "MutationTaster prediction",
      track_height = opts$mutationtaster_track_height,
      unknown_color = UNKNOWN_CATEGORY_COLOR,
      palette = "Set 2",
      fixed_colors = MUTATIONTASTER_PREDICTION_COLORS
    ),
    list(
      column = "dominant_clinsig",
      label = "ClinSig",
      track_height = opts$clinsig_track_height,
      unknown_color = UNKNOWN_CATEGORY_COLOR,
      palette = "Set 2"
    )
  )

  Filter(
    function(spec) {
      spec$column %in% colnames(summary_df) &&
        any(sanitize_category(summary_df[[spec$column]], default = "Unknown") != "Unknown")
    },
    track_specs
  )
}

append_snpeff_impact_summary <- function(summary_df, snpeff_impact_summary_tsv) {
  snpeff_impact_summary_tsv <- check_optional_csv(snpeff_impact_summary_tsv, "SNPEFF_IMPACT_SUMMARY_TSV")
  if (is.na(snpeff_impact_summary_tsv) || !nzchar(snpeff_impact_summary_tsv)) {
    return(summary_df)
  }

  message(sprintf("[INFO] Reading snpEff impact summary TSV: %s", snpeff_impact_summary_tsv))
  snpeff_df <- readr::read_tsv(
    file = snpeff_impact_summary_tsv,
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "minimal",
    col_types = readr::cols(.default = readr::col_character())
  )
  snpeff_df <- as.data.frame(snpeff_df, stringsAsFactors = FALSE)

  required_columns <- c("snv_id", "impact")
  missing_columns <- setdiff(required_columns, colnames(snpeff_df))
  if (length(missing_columns) > 0) {
    stop(
      sprintf("SNPEFF_IMPACT_SUMMARY_TSV is missing required columns: %s", paste(missing_columns, collapse = ", ")),
      call. = FALSE
    )
  }

  snpeff_df <- data.frame(
    SNV = trimws(as.character(snpeff_df$snv_id)),
    snpeff_impact = normalize_vep_impact(snpeff_df$impact),
    stringsAsFactors = FALSE
  )
  snpeff_df <- snpeff_df[
    !is.na(snpeff_df$SNV) & nzchar(snpeff_df$SNV) & snpeff_df$SNV %in% summary_df$SNV,
    ,
    drop = FALSE
  ]
  snpeff_df <- snpeff_df[!duplicated(snpeff_df$SNV), , drop = FALSE]

  if (nrow(snpeff_df) == 0) {
    summary_df$snpeff_impact <- "Unknown"
    message(sprintf("[WARN] snpEff impact summary matched 0/%d plotted SNVs.", nrow(summary_df)))
    return(summary_df)
  }

  merged_df <- dplyr::left_join(summary_df, snpeff_df, by = "SNV")
  matched_count <- sum(!is.na(merged_df$snpeff_impact))
  merged_df$snpeff_impact <- sanitize_category(merged_df$snpeff_impact, default = "Unknown")
  if (matched_count < nrow(merged_df)) {
    message(sprintf(
      "[WARN] snpEff impact summary matched %d/%d plotted SNVs.",
      matched_count,
      nrow(merged_df)
    ))
  }

  merged_df
}

append_regulationspotter_summary <- function(summary_df, all_vars_tsv, spotter_tsv) {
  variant_map <- read_regulationspotter_variant_map(all_vars_tsv)
  spotter_tsv <- check_optional_csv(spotter_tsv, "REGULATIONSPOTTER_SPOTTER_TSV")
  if (is.null(variant_map) || is.na(spotter_tsv) || !nzchar(spotter_tsv)) {
    return(summary_df)
  }

  message(sprintf("[INFO] Reading RegulationSpotter spotter TSV: %s", spotter_tsv))
  spotter_df <- readr::read_tsv(
    file = spotter_tsv,
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "minimal",
    col_types = readr::cols(.default = readr::col_character())
  )
  spotter_df <- as.data.frame(spotter_df, stringsAsFactors = FALSE)
  required_columns <- c(
    "var_number", "prediction", "promoter", "h3k4me3", "dhs", "tfbs",
    "enhancer", "polymerase", "interactions", "score", "cadd_scaled", "phylop", "phastcons"
  )
  missing_columns <- setdiff(required_columns, colnames(spotter_df))
  if (length(missing_columns) > 0) {
    stop(
      sprintf("REGULATIONSPOTTER_SPOTTER_TSV is missing required columns: %s", paste(missing_columns, collapse = ", ")),
      call. = FALSE
    )
  }

  spotter_df <- data.frame(
    var_number = suppressWarnings(as.integer(spotter_df$var_number)),
    regspotter_prediction = map_regulationspotter_prediction(spotter_df$prediction),
    regspotter_score = suppressWarnings(as.numeric(spotter_df$score)),
    regspotter_promoter = sanitize_category(spotter_df$promoter, default = "None"),
    regspotter_h3k4me3 = sanitize_category(spotter_df$h3k4me3, default = "None"),
    regspotter_dhs = sanitize_binary_flag(spotter_df$dhs),
    regspotter_tfbs = sanitize_binary_flag(spotter_df$tfbs),
    regspotter_enhancer = sanitize_binary_flag(spotter_df$enhancer),
    regspotter_polymerase = sanitize_binary_flag(spotter_df$polymerase),
    regspotter_interactions = sanitize_binary_flag(spotter_df$interactions),
    regspotter_cadd_scaled = suppressWarnings(as.numeric(spotter_df$cadd_scaled)),
    regspotter_phylop = suppressWarnings(as.numeric(spotter_df$phylop)),
    regspotter_phastcons = suppressWarnings(as.numeric(spotter_df$phastcons)),
    stringsAsFactors = FALSE
  )
  spotter_df <- spotter_df[is.finite(spotter_df$var_number), , drop = FALSE]
  spotter_df <- spotter_df[!duplicated(spotter_df$var_number), , drop = FALSE]
  spotter_df <- dplyr::left_join(spotter_df, variant_map, by = "var_number")
  spotter_df <- spotter_df[!is.na(spotter_df$SNV) & nzchar(spotter_df$SNV) & spotter_df$SNV %in% summary_df$SNV, , drop = FALSE]
  if (nrow(spotter_df) > 0) {
    direct_matches <- sum(spotter_df$snv_match_type == "direct", na.rm = TRUE)
    reverse_matches <- sum(spotter_df$snv_match_type == "reverse_complement", na.rm = TRUE)
    message(sprintf("[INFO] RegulationSpotter spotter join matches: direct=%d, reverse_complement=%d.", direct_matches, reverse_matches))
  }
  spotter_df <- spotter_df[order(spotter_df$var_number, spotter_df$snv_match_priority), , drop = FALSE]
  spotter_df <- spotter_df[!duplicated(spotter_df$var_number), , drop = FALSE]
  spotter_df <- spotter_df[!duplicated(spotter_df$SNV), , drop = FALSE]

  if (nrow(spotter_df) == 0) {
    message(sprintf("[WARN] RegulationSpotter spotter TSV matched 0/%d plotted SNVs.", nrow(summary_df)))
    return(summary_df)
  }

  spotter_df$regspotter_score_norm <- rescale_to_unit(spotter_df$regspotter_score)
  join_df <- spotter_df[, setdiff(colnames(spotter_df), "var_number"), drop = FALSE]
  merged_df <- dplyr::left_join(summary_df, join_df, by = "SNV")
  matched_count <- sum(!is.na(merged_df$regspotter_prediction))
  if (matched_count < nrow(merged_df)) {
    message(sprintf(
      "[WARN] RegulationSpotter spotter TSV matched %d/%d plotted SNVs.",
      matched_count,
      nrow(merged_df)
    ))
  }

  merged_df
}

append_regulationtaster_summary <- function(summary_df, all_vars_tsv, taster_tsv) {
  variant_map <- read_regulationspotter_variant_map(all_vars_tsv)
  taster_tsv <- check_optional_csv(taster_tsv, "REGULATIONSPOTTER_TASTER_TSV")
  if (is.null(variant_map) || is.na(taster_tsv) || !nzchar(taster_tsv)) {
    return(summary_df)
  }

  message(sprintf("[INFO] Reading RegulationSpotter taster TSV: %s", taster_tsv))
  taster_df <- readr::read_tsv(
    file = taster_tsv,
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "minimal",
    col_types = readr::cols(.default = readr::col_character())
  )
  taster_df <- as.data.frame(taster_df, stringsAsFactors = FALSE)
  required_columns <- c(
    "var_number", "prediction", "confidence", "distance_splice_site", "splicing",
    "frameshift", "prematurestop", "splice_site", "start_atg", "poly_a", "kozak"
  )
  missing_columns <- setdiff(required_columns, colnames(taster_df))
  if (length(missing_columns) > 0) {
    stop(
      sprintf("REGULATIONSPOTTER_TASTER_TSV is missing required columns: %s", paste(missing_columns, collapse = ", ")),
      call. = FALSE
    )
  }

  prediction_rank <- c(
    "1" = 5,
    "2" = 4,
    "3" = 3,
    "95" = 2,
    "99" = 1,
    "0" = 0
  )
  taster_df <- data.frame(
    var_number = suppressWarnings(as.integer(taster_df$var_number)),
    prediction_code = trimws(as.character(taster_df$prediction)),
    confidence = suppressWarnings(as.numeric(taster_df$confidence)),
    distance_splice_site = suppressWarnings(as.numeric(taster_df$distance_splice_site)),
    any_splicing = sanitize_binary_flag(taster_df$splicing),
    any_frameshift = sanitize_binary_flag(taster_df$frameshift),
    any_prematurestop = sanitize_binary_flag(taster_df$prematurestop),
    any_splice_site = sanitize_binary_flag(taster_df$splice_site),
    any_start_atg = sanitize_binary_flag(taster_df$start_atg),
    any_poly_a = sanitize_binary_flag(taster_df$poly_a),
    any_kozak = sanitize_binary_flag(taster_df$kozak),
    stringsAsFactors = FALSE
  )
  taster_df <- taster_df[is.finite(taster_df$var_number), , drop = FALSE]
  taster_df$prediction_rank <- unname(prediction_rank[taster_df$prediction_code])
  taster_df$prediction_rank[!is.finite(taster_df$prediction_rank)] <- 0

  summarised_df <- taster_df |>
    dplyr::group_by(var_number) |>
    dplyr::summarise(
      regulationtaster_prediction_code = prediction_code[which.max(prediction_rank)],
      regulationtaster_prediction = map_regulationtaster_prediction(prediction_code[which.max(prediction_rank)]),
      regulationtaster_prediction_mixed = dplyr::n_distinct(prediction_code[!is.na(prediction_code) & nzchar(prediction_code)]) > 1,
      regulationtaster_max_confidence = if (any(is.finite(confidence))) max(confidence, na.rm = TRUE) else NA_real_,
      regulationtaster_min_distance_splice_site = if (any(is.finite(distance_splice_site))) min(abs(distance_splice_site[is.finite(distance_splice_site)]), na.rm = TRUE) else NA_real_,
      regulationtaster_any_splicing = as.integer(any(any_splicing > 0, na.rm = TRUE)),
      regulationtaster_any_frameshift = as.integer(any(any_frameshift > 0, na.rm = TRUE)),
      regulationtaster_any_prematurestop = as.integer(any(any_prematurestop > 0, na.rm = TRUE)),
      regulationtaster_any_splice_site = as.integer(any(any_splice_site > 0, na.rm = TRUE)),
      regulationtaster_any_start_atg = as.integer(any(any_start_atg > 0, na.rm = TRUE)),
      regulationtaster_any_poly_a = as.integer(any(any_poly_a > 0, na.rm = TRUE)),
      regulationtaster_any_kozak = as.integer(any(any_kozak > 0, na.rm = TRUE)),
      .groups = "drop"
    )
  summarised_df <- as.data.frame(summarised_df, stringsAsFactors = FALSE)
  summarised_df <- dplyr::left_join(summarised_df, variant_map, by = "var_number")
  summarised_df <- summarised_df[!is.na(summarised_df$SNV) & nzchar(summarised_df$SNV) & summarised_df$SNV %in% summary_df$SNV, , drop = FALSE]
  if (nrow(summarised_df) > 0) {
    direct_matches <- sum(summarised_df$snv_match_type == "direct", na.rm = TRUE)
    reverse_matches <- sum(summarised_df$snv_match_type == "reverse_complement", na.rm = TRUE)
    message(sprintf("[INFO] RegulationSpotter taster join matches: direct=%d, reverse_complement=%d.", direct_matches, reverse_matches))
  }
  summarised_df <- summarised_df[order(summarised_df$var_number, summarised_df$snv_match_priority), , drop = FALSE]
  summarised_df <- summarised_df[!duplicated(summarised_df$var_number), , drop = FALSE]
  summarised_df <- summarised_df[!duplicated(summarised_df$SNV), , drop = FALSE]

  if (nrow(summarised_df) == 0) {
    message(sprintf("[WARN] RegulationSpotter taster TSV matched 0/%d plotted SNVs.", nrow(summary_df)))
    return(summary_df)
  }

  join_df <- summarised_df[, setdiff(colnames(summarised_df), "var_number"), drop = FALSE]
  merged_df <- dplyr::left_join(summary_df, join_df, by = "SNV")
  matched_count <- sum(!is.na(merged_df$regulationtaster_prediction))
  if (matched_count < nrow(merged_df)) {
    message(sprintf(
      "[WARN] RegulationSpotter taster TSV matched %d/%d plotted SNVs.",
      matched_count,
      nrow(merged_df)
    ))
  }

  merged_df
}

append_vep_impact_summary <- function(summary_df, vep_annotation_csv) {
  vep_annotation_csv <- check_optional_csv(vep_annotation_csv, "VEP_ANNOTATION_CSV")
  if (is.na(vep_annotation_csv) || !nzchar(vep_annotation_csv)) {
    return(summary_df)
  }

  message(sprintf("[INFO] Reading VEP annotation CSV: %s", vep_annotation_csv))
  vep_df <- readr::read_csv(
    file = vep_annotation_csv,
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "minimal",
    col_types = readr::cols(.default = readr::col_character())
  )
  vep_df <- as.data.frame(vep_df, stringsAsFactors = FALSE)

  required_columns <- c("variant", "impact")
  missing_columns <- setdiff(required_columns, colnames(vep_df))
  if (length(missing_columns) > 0) {
    stop(
      sprintf("VEP_ANNOTATION_CSV is missing required columns: %s", paste(missing_columns, collapse = ", ")),
      call. = FALSE
    )
  }

  impact_rank <- c(HIGH = 4, MODERATE = 3, LOW = 2, MODIFIER = 1, Unknown = 0)
  vep_df <- data.frame(
    SNV = trimws(as.character(vep_df$variant)),
    vep_impact = normalize_vep_impact(vep_df$impact),
    stringsAsFactors = FALSE
  )
  vep_df <- vep_df[!is.na(vep_df$SNV) & nzchar(vep_df$SNV) & vep_df$SNV %in% summary_df$SNV, , drop = FALSE]
  vep_df$impact_rank <- unname(impact_rank[vep_df$vep_impact])
  vep_df$impact_rank[!is.finite(vep_df$impact_rank)] <- 0

  if (nrow(vep_df) == 0) {
    summary_df$vep_impact <- "Unknown"
    message(sprintf("[WARN] VEP annotation matched 0/%d plotted SNVs.", nrow(summary_df)))
    return(summary_df)
  }

  impact_summary <- vep_df |>
    dplyr::arrange(SNV, dplyr::desc(impact_rank), vep_impact) |>
    dplyr::group_by(SNV) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup() |>
    dplyr::select(SNV, vep_impact)

  merged_df <- dplyr::left_join(summary_df, impact_summary, by = "SNV")
  matched_count <- sum(!is.na(merged_df$vep_impact))
  merged_df$vep_impact <- sanitize_category(merged_df$vep_impact, default = "Unknown")
  if (matched_count < nrow(merged_df)) {
    message(sprintf(
      "[WARN] VEP annotation matched %d/%d plotted SNVs.",
      matched_count,
      nrow(merged_df)
    ))
  }

  merged_df
}

append_mutationtaster_prediction_summary <- function(summary_df, mutationtaster_annotation_tsv) {
  mutationtaster_annotation_tsv <- check_optional_csv(mutationtaster_annotation_tsv, "MUTATIONTASTER_ANNOTATION_TSV")
  if (is.na(mutationtaster_annotation_tsv) || !nzchar(mutationtaster_annotation_tsv)) {
    return(summary_df)
  }

  message(sprintf("[INFO] Reading MutationTaster annotation TSV: %s", mutationtaster_annotation_tsv))
  required_columns <- c("Chr", "Position", "Prediction", "Ref", "Alt")
  mt_df <- readr::read_tsv(
    file = mutationtaster_annotation_tsv,
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "minimal",
    col_types = readr::cols(.default = readr::col_character()),
    col_select = dplyr::all_of(required_columns)
  )
  mt_df <- as.data.frame(mt_df, stringsAsFactors = FALSE)
  missing_columns <- setdiff(required_columns, colnames(mt_df))
  if (length(missing_columns) > 0) {
    stop(
      sprintf("MUTATIONTASTER_ANNOTATION_TSV is missing required columns: %s", paste(missing_columns, collapse = ", ")),
      call. = FALSE
    )
  }

  prediction_rank <- c(
    "Deleterious (ClinVar)" = 5,
    "Deleterious (fs/PTC)" = 4,
    Deleterious = 3,
    "Benign (auto)" = 2,
    Benign = 1,
    Unknown = 0
  )
  mt_df <- data.frame(
    SNV = build_snv_id_from_parts(mt_df$Chr, mt_df$Position, mt_df$Ref, mt_df$Alt),
    mutationtaster_prediction = sanitize_category(mt_df$Prediction, default = "Unknown"),
    stringsAsFactors = FALSE
  )
  mt_df <- mt_df[
    !is.na(mt_df$SNV) & nzchar(mt_df$SNV) & mt_df$SNV %in% summary_df$SNV,
    ,
    drop = FALSE
  ]
  mt_df$prediction_rank <- unname(prediction_rank[mt_df$mutationtaster_prediction])
  mt_df$prediction_rank[!is.finite(mt_df$prediction_rank)] <- 0

  if (nrow(mt_df) == 0) {
    summary_df$mutationtaster_prediction <- "Unknown"
    message(sprintf("[WARN] MutationTaster annotation matched 0/%d plotted SNVs.", nrow(summary_df)))
    return(summary_df)
  }

  prediction_summary <- mt_df |>
    dplyr::arrange(SNV, dplyr::desc(prediction_rank), mutationtaster_prediction) |>
    dplyr::group_by(SNV) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup() |>
    dplyr::select(SNV, mutationtaster_prediction)

  merged_df <- dplyr::left_join(summary_df, prediction_summary, by = "SNV")
  matched_count <- sum(!is.na(merged_df$mutationtaster_prediction))
  merged_df$mutationtaster_prediction <- sanitize_category(merged_df$mutationtaster_prediction, default = "Unknown")
  if (matched_count < nrow(merged_df)) {
    message(sprintf(
      "[WARN] MutationTaster annotation matched %d/%d plotted SNVs.",
      matched_count,
      nrow(merged_df)
    ))
  }

  merged_df
}

read_alphagenome_atac_metadata <- function(metadata_csv) {
  metadata_csv <- check_optional_csv(metadata_csv, "ALPHAGENOME_ATAC_METADATA_CSV")
  if (is.na(metadata_csv) || !nzchar(metadata_csv)) {
    return(NULL)
  }

  metadata_df <- readr::read_csv(
    file = metadata_csv,
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "minimal"
  )
  metadata_df <- as.data.frame(metadata_df, stringsAsFactors = FALSE)
  required_columns <- c("type", "name", "rows")
  missing_columns <- setdiff(required_columns, colnames(metadata_df))
  if (length(missing_columns) > 0) {
    stop(
      sprintf("ALPHAGENOME_ATAC_METADATA_CSV is missing required columns: %s", paste(missing_columns, collapse = ", ")),
      call. = FALSE
    )
  }

  track_name <- metadata_df$name[match("track_name", metadata_df$type)]
  biosample_name <- metadata_df$name[match("biosample_name", metadata_df$type)]
  message(sprintf(
    "[INFO] Loaded AlphaGenome ATAC metadata: track=%s; biosample=%s.",
    ifelse(is.na(track_name), "unknown", track_name),
    ifelse(is.na(biosample_name), "unknown", biosample_name)
  ))

  metadata_df
}

boxplot_stat_columns <- function(prefix) {
  paste0(prefix, c("_min", "_q1", "_median", "_q3", "_max"))
}

append_alphagenome_boxplot_summary <- function(summary_df, boxplot_stats_csv) {
  boxplot_stats_csv <- check_optional_csv(boxplot_stats_csv, "ALPHAGENOME_BOXPLOT_STATS_CSV")
  if (is.na(boxplot_stats_csv) || !nzchar(boxplot_stats_csv)) {
    return(summary_df)
  }

  message(sprintf("[INFO] Reading AlphaGenome boxplot stats: %s", boxplot_stats_csv))
  stats_df <- readr::read_csv(
    file = boxplot_stats_csv,
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "minimal"
  )
  stats_df <- as.data.frame(stats_df, stringsAsFactors = FALSE)
  if (!"input_snv" %in% colnames(stats_df)) {
    stop("ALPHAGENOME_BOXPLOT_STATS_CSV is missing required column: input_snv", call. = FALSE)
  }

  track_specs <- list(
    list(prefix = "rna_abs_quantile", output_prefix = "ag_rna"),
    list(prefix = "atac_abs_quantile", output_prefix = "ag_atac")
  )

  merged_df <- summary_df
  for (spec in track_specs) {
    required_columns <- boxplot_stat_columns(spec$prefix)
    missing_columns <- setdiff(required_columns, colnames(stats_df))
    if (length(missing_columns) > 0) {
      stop(
        sprintf(
          "ALPHAGENOME_BOXPLOT_STATS_CSV is missing required columns for %s: %s",
          spec$prefix,
          paste(missing_columns, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    join_df <- stats_df[, c("input_snv", required_columns), drop = FALSE]
    colnames(join_df) <- c("SNV", paste0(spec$output_prefix, c("_min", "_q1", "_median", "_q3", "_max")))
    join_df$SNV <- trimws(as.character(join_df$SNV))
    join_df <- join_df[!is.na(join_df$SNV) & nzchar(join_df$SNV), , drop = FALSE]
    join_df <- join_df[!duplicated(join_df$SNV), , drop = FALSE]
    for (value_col in setdiff(colnames(join_df), "SNV")) {
      join_df[[value_col]] <- suppressWarnings(as.numeric(join_df[[value_col]]))
    }

    merged_df <- dplyr::left_join(merged_df, join_df, by = "SNV")
  }

  matched_count <- sum(!is.na(merged_df$ag_rna_median) | !is.na(merged_df$ag_atac_median))
  if (matched_count < nrow(merged_df)) {
    message(sprintf(
      "[WARN] AlphaGenome boxplot stats matched %d/%d summary SNVs.",
      matched_count,
      nrow(merged_df)
    ))
  }

  merged_df
}

build_boxplot_tracks <- function(summary_df, opts) {
  track_specs <- list(
    list(
      prefix = "ag_rna",
      label = "AlphaGenome RNA abs quantile",
      track_height = opts$alphagenome_rna_track_height,
      line_color = ALPHAGENOME_RNA_BOX_COLOR,
      fill_color = grDevices::adjustcolor(ALPHAGENOME_RNA_BOX_COLOR, alpha.f = 0.32)
    ),
    list(
      prefix = "ag_atac",
      label = "AlphaGenome ATAC abs quantile",
      track_height = opts$alphagenome_atac_track_height,
      line_color = ALPHAGENOME_ATAC_BOX_COLOR,
      fill_color = grDevices::adjustcolor(ALPHAGENOME_ATAC_BOX_COLOR, alpha.f = 0.32)
    )
  )

  Filter(
    function(spec) {
      required_columns <- paste0(spec$prefix, c("_min", "_q1", "_median", "_q3", "_max"))
      all(required_columns %in% colnames(summary_df)) &&
        any(is.finite(suppressWarnings(as.numeric(summary_df[[paste0(spec$prefix, "_median")]]))))
    },
    track_specs
  )
}

draw_boxplot_track <- function(summary_df, spec, opts) {
  stat_cols <- paste0(spec$prefix, c("_min", "_q1", "_median", "_q3", "_max"))
  circlize::circos.trackPlotRegion(
    factors = summary_df$SNV,
    ylim = c(0, 1),
    bg.col = opts$ring_background_color,
    bg.border = opts$ring_border_color,
    bg.lwd = 0.45,
    track.height = spec$track_height,
    panel.fun = function(x, y) {
      snv <- circlize::CELL_META$sector.index
      row_idx <- match(snv, summary_df$SNV)
      values <- suppressWarnings(as.numeric(unlist(summary_df[row_idx, stat_cols], use.names = FALSE)))
      names(values) <- c("min", "q1", "median", "q3", "max")
      if (any(!is.finite(values))) {
        return(invisible(NULL))
      }
      values <- pmin(pmax(values, 0), 1)
      x_center <- circlize::CELL_META$xcenter
      x_limits <- circlize::CELL_META$xlim
      half_width <- diff(x_limits) * opts$alphagenome_box_width_fraction / 2
      x_left <- x_center - half_width
      x_right <- x_center + half_width

      circlize::circos.segments(
        x0 = x_center,
        y0 = values["min"],
        x1 = x_center,
        y1 = values["max"],
        col = spec$line_color,
        lwd = 0.75
      )
      circlize::circos.segments(
        x0 = x_left,
        y0 = values["min"],
        x1 = x_right,
        y1 = values["min"],
        col = spec$line_color,
        lwd = 0.65
      )
      circlize::circos.segments(
        x0 = x_left,
        y0 = values["max"],
        x1 = x_right,
        y1 = values["max"],
        col = spec$line_color,
        lwd = 0.65
      )
      if (values["q3"] > values["q1"]) {
        circlize::circos.rect(
          xleft = x_left,
          ybottom = values["q1"],
          xright = x_right,
          ytop = values["q3"],
          col = spec$fill_color,
          border = spec$line_color,
          lwd = 0.7
        )
      }
      circlize::circos.segments(
        x0 = x_left,
        y0 = values["median"],
        x1 = x_right,
        y1 = values["median"],
        col = spec$line_color,
        lwd = 0.95
      )
      if (values["q3"] <= values["q1"]) {
        circlize::circos.points(
          x = x_center,
          y = values["median"],
          pch = 16,
          cex = 0.28,
          col = spec$line_color
        )
      }
    }
  )
}

draw_lollipop_track <- function(summary_df, value_col, point_color, baseline_color, opts) {
  circlize::circos.trackPlotRegion(
    factors = summary_df$SNV,
    ylim = c(0, 1),
    bg.col = opts$ring_background_color,
    bg.border = opts$ring_border_color,
    bg.lwd = 0.45,
    track.height = opts$score_ring_track_height,
    panel.fun = function(x, y) {
      snv <- circlize::CELL_META$sector.index
      row_idx <- match(snv, summary_df$SNV)
      x_center <- circlize::CELL_META$xcenter
      y_value <- summary_df[[value_col]][row_idx]
      if (!is.finite(y_value)) {
        y_value <- 0
      }
      x_limits <- circlize::CELL_META$xlim

      circlize::circos.lines(
        x = x_limits,
        y = c(0, 0),
        col = baseline_color,
        lwd = 0.85,
        straight = TRUE
      )
      circlize::circos.segments(
        x0 = x_center,
        y0 = 0,
        x1 = x_center,
        y1 = y_value,
        col = point_color,
        lwd = 1.15
      )
      circlize::circos.points(
        x = x_center,
        y = y_value,
        pch = 16,
        cex = 0.42,
        col = point_color
      )
    }
  )
}

draw_category_band_track <- function(summary_df, value_col, category_colors, track_height, opts) {
  circlize::circos.trackPlotRegion(
    factors = summary_df$SNV,
    ylim = c(0, 1),
    bg.col = opts$ring_background_color,
    bg.border = opts$ring_border_color,
    bg.lwd = 0.45,
    track.height = track_height,
    panel.fun = function(x, y) {
      snv <- circlize::CELL_META$sector.index
      row_idx <- match(snv, summary_df$SNV)
      category_value <- sanitize_category(summary_df[[value_col]][row_idx], default = "Unknown")
      x_limits <- circlize::CELL_META$xlim
      fill_color <- category_colors[[category_value]]
      if (is.null(fill_color) || is.na(fill_color)) {
        fill_color <- UNKNOWN_CATEGORY_COLOR
      }

      circlize::circos.rect(
        xleft = x_limits[1],
        ybottom = 0,
        xright = x_limits[2],
        ytop = 1,
        col = fill_color,
        border = opts$ring_border_color,
        lwd = 0.32
      )
    }
  )
}

draw_label_track <- function(summary_df, opts) {
  circlize::circos.trackPlotRegion(
    factors = summary_df$SNV,
    ylim = c(0, 1),
    bg.border = NA,
    track.height = opts$label_track_height,
    panel.fun = function(x, y) {
      circlize::circos.text(
        x = circlize::CELL_META$xcenter,
        y = 0.08,
        labels = circlize::CELL_META$sector.index,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5),
        cex = opts$label_cex,
        col = "#222222"
      )
    }
  )
}

draw_circular_plot <- function(summary_df, opts, numeric_summary_df = summary_df, png_path = NULL, pdf_path = NULL) {
  category_tracks <- build_category_tracks(summary_df, opts = opts)
  boxplot_tracks <- build_boxplot_tracks(summary_df, opts = opts)
  category_palettes <- setNames(
    lapply(
      category_tracks,
      function(spec) {
        make_category_palette(
          values = summary_df[[spec$column]],
          unknown_color = spec$unknown_color,
          palette = spec$palette,
          fixed_colors = spec$fixed_colors
        )
      }
    ),
    vapply(category_tracks, function(spec) spec$column, character(1))
  )

  plot_once <- function() {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    circlize::circos.clear()
    on.exit(circlize::circos.clear(), add = TRUE)

    graphics::par(
      mar = if (isTRUE(opts$show_legends)) c(1.0, 1.0, 2.2, 8.0) else c(0.8, 0.8, 0.8, 0.8),
      family = opts$font_family,
      xpd = NA
    )

    gap_after <- if (nrow(summary_df) > 1) {
      c(rep(opts$sector_gap_degree, nrow(summary_df) - 1L), opts$final_gap_degree)
    } else {
      opts$final_gap_degree
    }

    circlize::circos.par(
      start.degree = opts$start_degree,
      gap.after = gap_after,
      cell.padding = c(0, 0, 0, 0),
      track.margin = c(0.004, 0.004),
      points.overflow.warning = FALSE,
      canvas.xlim = if (isTRUE(opts$show_legends)) c(-1.45, 1.85) else c(-1.22, 1.22),
      canvas.ylim = if (isTRUE(opts$show_legends)) c(-1.55, 1.55) else c(-1.22, 1.22)
    )

    circlize::circos.initialize(
      factors = summary_df$SNV,
      xlim = cbind(rep(0, nrow(summary_df)), rep(1, nrow(summary_df)))
    )

    # Outermost track: SNV labels.
    draw_label_track(summary_df = summary_df, opts = opts)

    # Numeric ring: per-SNV max perturbation score.
    draw_lollipop_track(
      summary_df = summary_df,
      value_col = "score_norm",
      point_color = opts$score_ring_color,
      baseline_color = opts$score_ring_baseline_color,
      opts = opts
    )

    # AlphaGenome rings: per-SNV min/Q1/median/Q3/max summaries.
    for (spec in boxplot_tracks) {
      draw_boxplot_track(summary_df = summary_df, spec = spec, opts = opts)
    }

    draw_outer_ring_value_labels(
      summary_df = summary_df,
      numeric_summary_df = numeric_summary_df,
      boxplot_tracks = boxplot_tracks,
      opts = opts
    )

    # Inner annotation rings: dominant categorical annotations per SNV.
    for (spec in category_tracks) {
        draw_category_band_track(
          summary_df = summary_df,
          value_col = spec$column,
          category_colors = category_palettes[[spec$column]],
          track_height = spec$track_height,
          opts = opts
        )
      }

    if (!is.na(opts$plot_title) && nzchar(opts$plot_title)) {
      graphics::title(main = opts$plot_title, cex.main = 1.02, line = 0.4)
    }

    if (isTRUE(opts$show_legends)) {
      func_colors <- category_palettes[["dominant_func"]]
      if (!is.null(func_colors) && length(func_colors) > 0) {
        graphics::legend(
          x = "right",
          inset = c(-0.28, 0.02),
          title = "Func",
          legend = names(func_colors),
          fill = unname(func_colors),
          border = NA,
          bty = "n",
          cex = 0.72,
          pt.cex = 1.2,
          y.intersp = 0.9,
          ncol = if (length(func_colors) > 12) 2 else 1,
          xpd = NA
        )
      }

      graphics::legend(
        x = "left",
        inset = c(-0.12, 0.02),
        legend = c(
          sprintf("max %s", opts$score_column),
          vapply(boxplot_tracks, function(spec) spec$label, character(1)),
          vapply(category_tracks, function(spec) spec$label, character(1))
        ),
        col = c(
          opts$score_ring_color,
          vapply(boxplot_tracks, function(spec) spec$line_color, character(1)),
          rep(NA_character_, length(category_tracks))
        ),
        lty = c(1, rep(1, length(boxplot_tracks)), rep(NA, length(category_tracks))),
        lwd = c(2, rep(1.2, length(boxplot_tracks)), rep(NA, length(category_tracks))),
        pch = c(16, rep(NA, length(boxplot_tracks)), rep(15, length(category_tracks))),
        pt.cex = c(1, rep(NA, length(boxplot_tracks)), rep(1.2, length(category_tracks))),
        bty = "n",
        cex = 0.70,
        text.col = "#222222",
        xpd = NA,
        y.intersp = 1.0
      )
    }
  }

  save_plot <- function(file_path, device) {
    device(file_path)
    on.exit(grDevices::dev.off(), add = TRUE)
    plot_once()
    invisible(TRUE)
  }

  if (!is.null(png_path)) {
    save_plot(
      png_path,
      function(file_path) {
        grDevices::png(
          filename = file_path,
          width = opts$plot_width,
          height = opts$plot_height,
          units = "in",
          res = opts$plot_dpi,
          bg = "white"
        )
      }
    )
  }

  if (!is.null(pdf_path)) {
    save_plot(
      pdf_path,
      function(file_path) {
        grDevices::pdf(
          file = file_path,
          width = opts$plot_width,
          height = opts$plot_height,
          useDingbats = FALSE,
          bg = "white"
        )
      }
    )
  }
}

contrast_text_color <- function(fill_color) {
  rgb <- grDevices::col2rgb(fill_color) / 255
  luminance <- 0.299 * rgb[1, ] + 0.587 * rgb[2, ] + 0.114 * rgb[3, ]
  ifelse(luminance < 0.55, "white", "#222222")
}

format_ring_axis_label <- function(values) {
  vapply(
    values,
    function(value) {
      if (!is.finite(value)) {
        return("NA")
      }
      formatC(signif(value, 3), format = "fg", flag = "#")
    },
    character(1)
  )
}

build_outer_ring_value_label_specs <- function(numeric_summary_df, boxplot_tracks, opts) {
  specs <- list()
  normalized_positions <- c(0, 0.5, 1)
  normalized_labels <- c("0", "0.5", "1")

  score_values <- suppressWarnings(as.numeric(numeric_summary_df$max_score))
  score_values <- score_values[is.finite(score_values)]
  if (length(score_values) > 0) {
    specs[[length(specs) + 1L]] <- list(
      track_index = 2L,
      positions = normalized_positions,
      labels = normalized_labels
    )
  }

  if (length(boxplot_tracks) > 0) {
    for (idx in seq_along(boxplot_tracks)) {
      spec <- boxplot_tracks[[idx]]
      stat_cols <- paste0(spec$prefix, c("_min", "_q1", "_median", "_q3", "_max"))
      if (!all(stat_cols %in% colnames(numeric_summary_df))) {
        next
      }
      ring_values <- suppressWarnings(as.numeric(unlist(numeric_summary_df[, stat_cols, drop = FALSE], use.names = FALSE)))
      ring_values <- ring_values[is.finite(ring_values)]
      if (length(ring_values) == 0) {
        next
      }
      specs[[length(specs) + 1L]] <- list(
        track_index = 2L + idx,
        positions = normalized_positions,
        labels = normalized_labels
      )
    }
  }

  specs
}

draw_outer_ring_value_labels <- function(summary_df, numeric_summary_df, boxplot_tracks, opts) {
  if (!isTRUE(opts$show_outer_ring_value_labels) || nrow(summary_df) == 0) {
    return(invisible(NULL))
  }

  label_specs <- build_outer_ring_value_label_specs(numeric_summary_df, boxplot_tracks, opts)
  if (length(label_specs) == 0) {
    return(invisible(NULL))
  }

  anchor_sector <- summary_df$SNV[[1]]
  for (spec in label_specs) {
    circlize::circos.yaxis(
      side = opts$outer_ring_value_axis_side,
      at = spec$positions,
      labels = spec$labels,
      sector.index = anchor_sector,
      track.index = spec$track_index,
      labels.cex = opts$outer_ring_value_label_cex,
      col = opts$outer_ring_value_label_color,
      labels.col = opts$outer_ring_value_label_color
    )
  }

  invisible(TRUE)
}

draw_summary_numeric_column <- function(values, x_center, color, title, opts, column_width = 1, highlight_mask = NULL, background_point_color = "#BDBDBD") {
  values <- suppressWarnings(as.numeric(values))
  if (is.null(highlight_mask)) {
    highlight_mask <- rep(TRUE, length(values))
  }
  highlight_mask <- as.logical(highlight_mask)
  if (length(highlight_mask) != length(values)) {
    highlight_mask <- rep(FALSE, length(values))
  }
  valid_mask <- is.finite(values)
  values <- values[valid_mask]
  highlight_mask <- highlight_mask[valid_mask]
  values <- pmin(pmax(values, 0), 1)
  if (length(values) == 0) {
    graphics::text(x_center, 0.5, "NA", cex = 0.72, col = "#666666")
    graphics::text(x_center, 1.08, title, cex = 0.66, font = 2, xpd = NA)
    return(invisible(NULL))
  }

  box_width <- opts$summary_panel_numeric_box_width * column_width
  point_width <- opts$summary_panel_numeric_point_width * column_width
  q_values <- stats::quantile(values, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE, names = FALSE)
  names(q_values) <- c("min", "q1", "median", "q3", "max")
  point_x <- x_center + stats::runif(length(values), -point_width, point_width)

  if (any(!highlight_mask)) {
    graphics::points(
      x = point_x[!highlight_mask],
      y = values[!highlight_mask],
      pch = 16,
      cex = 0.50,
      col = grDevices::adjustcolor(background_point_color, alpha.f = 0.75)
    )
  }
  if (any(highlight_mask)) {
    graphics::points(
      x = point_x[highlight_mask],
      y = values[highlight_mask],
      pch = 16,
      cex = 0.58,
      col = grDevices::adjustcolor(color, alpha.f = 0.80)
    )
  }
  graphics::segments(x_center, q_values["min"], x_center, q_values["max"], col = "#333333", lwd = 0.7)
  graphics::segments(x_center - box_width / 3, q_values["min"], x_center + box_width / 3, q_values["min"], col = "#333333", lwd = 0.7)
  graphics::segments(x_center - box_width / 3, q_values["max"], x_center + box_width / 3, q_values["max"], col = "#333333", lwd = 0.7)
  graphics::rect(
    x_center - box_width / 2,
    q_values["q1"],
    x_center + box_width / 2,
    q_values["q3"],
    col = grDevices::adjustcolor("white", alpha.f = 0.82),
    border = "#333333",
    lwd = 0.8
  )
  graphics::segments(
    x_center - box_width / 2,
    q_values["median"],
    x_center + box_width / 2,
    q_values["median"],
    col = "#333333",
    lwd = 1.2
  )
  label_values <- q_values[c("q1", "median", "q3")]
  graphics::text(
    x = x_center,
    y = label_values,
    labels = sprintf("%.2f", label_values),
    cex = 0.48,
    col = "#777777",
    font = 2
  )
  graphics::text(x_center, 1.08, title, cex = 0.58, font = 2, xpd = NA)
  invisible(NULL)
}

draw_summary_category_column <- function(values, x_center, category_colors, title, opts, drop_unknown = FALSE) {
  values <- sanitize_category(values, default = "Unknown")
  if (isTRUE(drop_unknown)) {
    values <- values[values != "Unknown"]
  }
  if (length(values) == 0) {
    graphics::text(x_center, 0.5, "NA", cex = 0.72, col = "#666666")
    graphics::text(x_center, 1.08, title, cex = 0.58, font = 2, xpd = NA)
    return(invisible(NULL))
  }
  category_levels <- c(intersect(names(category_colors), unique(values)), setdiff(sort(unique(values)), names(category_colors)))
  if (length(category_levels) == 0) {
    graphics::text(x_center, 0.5, "NA", cex = 0.72, col = "#666666")
    graphics::text(x_center, 1.08, title, cex = 0.58, font = 2, xpd = NA)
    return(invisible(NULL))
  }
  counts <- table(factor(values, levels = category_levels))
  counts <- counts[counts > 0]
  fractions <- as.numeric(counts) / sum(counts)
  names(fractions) <- names(counts)
  fractions <- sort(fractions, decreasing = TRUE)

  bar_width <- opts$summary_panel_category_bar_width
  bottom <- 0
  for (category_value in names(fractions)) {
    top <- bottom + fractions[[category_value]]
    fill_color <- category_colors[[category_value]]
    if (is.null(fill_color) || is.na(fill_color)) {
      fill_color <- UNKNOWN_CATEGORY_COLOR
    }
    graphics::rect(
      x_center - bar_width / 2,
      bottom,
      x_center + bar_width / 2,
      top,
      col = fill_color,
      border = "#333333",
      lwd = 0.55
    )
    if (fractions[[category_value]] >= 0.08) {
      graphics::text(
        x = x_center,
        y = (bottom + top) / 2,
        labels = sprintf("%.1f", 100 * fractions[[category_value]]),
        cex = opts$summary_panel_category_label_cex,
        col = contrast_text_color(fill_color),
        font = 2,
        srt = 90
      )
    }
    bottom <- top
  }
  graphics::text(x_center, 1.08, title, cex = 0.58, font = 2, xpd = NA)
  invisible(NULL)
}

summary_panel_category_title <- function(spec) {
  if (identical(spec$column, "snpeff_impact")) {
    return("snpEff\nimpact\n(%)")
  }
  if (identical(spec$column, "dominant_func")) {
    return("Func\n(%)")
  }
  if (identical(spec$column, "vep_impact")) {
    return("VEP\nimpact\n(%)")
  }
  if (identical(spec$column, "mutationtaster_prediction")) {
    return("MT\nprediction\n(%)")
  }
  if (identical(spec$column, "dominant_clinsig")) {
    return("ClinSig\n(%)")
  }
  paste0(spec$label, "\n(%)")
}

legend_label_from_multiline <- function(label) {
  label <- gsub("[\r\n]+", " ", as.character(label))
  trimws(gsub("\\s+", " ", label))
}

draw_summary_panel_legends <- function(numeric_specs, category_tracks, category_palettes, total_columns, opts) {
  legend_x_positions <- total_columns + opts$summary_panel_legend_x_offsets
  legend_y_top <- 1.12
  legend_y <- legend_y_top
  legend_col_idx <- 1L

  draw_legend_block <- function(...) {
    legend_args <- list(...)
    common_args <- list(
      x = legend_x_positions[legend_col_idx],
      y = legend_y,
      xjust = 0,
      yjust = 1,
      bty = "n",
      cex = opts$summary_panel_legend_text_cex,
      y.intersp = opts$summary_panel_legend_y_intersp,
      xpd = NA,
      text.col = "#222222"
    )
    legend_info <- do.call(graphics::legend, c(common_args, legend_args, list(plot = FALSE)))
    if ((legend_y - legend_info$rect$h) < 0.05 && legend_col_idx < length(legend_x_positions)) {
      legend_col_idx <<- legend_col_idx + 1L
      legend_y <<- legend_y_top
      common_args$x <- legend_x_positions[legend_col_idx]
      common_args$y <- legend_y
      legend_info <- do.call(graphics::legend, c(common_args, legend_args, list(plot = FALSE)))
    }
    do.call(graphics::legend, c(common_args, legend_args))
    legend_y <<- legend_y - legend_info$rect$h - opts$summary_panel_legend_block_gap
    invisible(legend_info)
  }

  if (length(numeric_specs) > 0) {
    draw_legend_block(
      title = "Numeric metrics",
      legend = vapply(numeric_specs, function(spec) legend_label_from_multiline(spec$title), character(1)),
      col = vapply(numeric_specs, function(spec) spec$color, character(1)),
      pch = 16,
      pt.cex = 0.85
    )
  }

  for (spec in category_tracks) {
    category_colors <- category_palettes[[spec$column]]
    if (identical(spec$column, "mutationtaster_prediction")) {
      category_colors <- category_colors[names(category_colors) != "Unknown"]
    }
    if (is.null(category_colors) || length(category_colors) == 0) {
      next
    }
    draw_legend_block(
      title = spec$label,
      legend = names(category_colors),
      fill = unname(category_colors),
      border = "#333333",
      ncol = if (length(category_colors) > 10) 2 else 1
    )
  }

  invisible(TRUE)
}

draw_summary_panel <- function(summary_df, opts, numeric_summary_df = summary_df, png_path = NULL, pdf_path = NULL) {
  category_tracks <- build_category_tracks(summary_df, opts = opts)
  category_palettes <- setNames(
    lapply(
      category_tracks,
      function(spec) {
        make_category_palette(
          values = summary_df[[spec$column]],
          unknown_color = spec$unknown_color,
          palette = spec$palette,
          fixed_colors = spec$fixed_colors
        )
      }
    ),
    vapply(category_tracks, function(spec) spec$column, character(1))
  )

  numeric_specs <- list(
    list(column = "score_norm", title = "VESPA\nscore", color = opts$score_ring_color),
    list(column = "ag_rna_median", title = "AlphaGenome\nRNA", color = ALPHAGENOME_RNA_BOX_COLOR),
    list(column = "ag_atac_median", title = "AlphaGenome\nATAC", color = ALPHAGENOME_ATAC_BOX_COLOR)
  )
  numeric_specs <- Filter(
    function(spec) {
      spec$column %in% colnames(numeric_summary_df) &&
        any(is.finite(suppressWarnings(as.numeric(numeric_summary_df[[spec$column]]))))
    },
    numeric_specs
  )

  numeric_column_width <- opts$summary_panel_numeric_column_width_factor
  category_column_width <- 1
  total_columns <- length(numeric_specs) * numeric_column_width + length(category_tracks) * category_column_width
  if (total_columns == 0) {
    return(invisible(FALSE))
  }

  plot_once <- function() {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(
      mar = c(2.0, 3.0, 2.5, opts$summary_panel_legend_right_margin),
      family = opts$font_family,
      xpd = NA,
      mgp = c(1.6, 0.35, 0)
    )
    graphics::plot.new()
    graphics::plot.window(
      xlim = c(0.45, total_columns + 0.55 + opts$summary_panel_legend_panel_width),
      ylim = c(0, 1.14),
      xaxs = "i",
      yaxs = "i"
    )
    graphics::abline(h = c(0, 0.25, 0.5, 0.75, 1.0), col = "#E5E5E5", lwd = 0.55)
    graphics::axis(2, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"), las = 1, cex.axis = 0.62, tck = -0.015)
    graphics::mtext("Scaled value / fraction", side = 2, line = 1.55, cex = 0.66)

    set.seed(1)
    numeric_highlight_mask <- numeric_summary_df$SNV %in% summary_df$SNV
    x_pos <- 0.45
    for (spec in numeric_specs) {
      x_center <- x_pos + numeric_column_width / 2
      draw_summary_numeric_column(
        values = numeric_summary_df[[spec$column]],
        x_center = x_center,
        color = spec$color,
        title = spec$title,
        opts = opts,
        column_width = numeric_column_width,
        highlight_mask = numeric_highlight_mask,
        background_point_color = opts$summary_panel_non_top_point_color
      )
      x_pos <- x_pos + numeric_column_width
    }

    for (spec in category_tracks) {
      x_center <- x_pos + category_column_width / 2
      draw_summary_category_column(
        values = summary_df[[spec$column]],
        x_center = x_center,
        category_colors = category_palettes[[spec$column]],
        title = summary_panel_category_title(spec),
        opts = opts,
        drop_unknown = identical(spec$column, "mutationtaster_prediction")
      )
      x_pos <- x_pos + category_column_width
    }

    graphics::rect(
      xleft = 0.45,
      ybottom = 0,
      xright = total_columns + 0.55,
      ytop = 1.14,
      border = "#333333",
      lwd = 0.8
    )
    draw_summary_panel_legends(
      numeric_specs = numeric_specs,
      category_tracks = category_tracks,
      category_palettes = category_palettes,
      total_columns = total_columns,
      opts = opts
    )
  }

  save_plot <- function(file_path, device) {
    device(file_path)
    on.exit(grDevices::dev.off(), add = TRUE)
    plot_once()
    invisible(TRUE)
  }

  if (!is.null(png_path)) {
    save_plot(
      png_path,
      function(file_path) {
        grDevices::png(
          filename = file_path,
          width = opts$summary_panel_width,
          height = opts$summary_panel_height,
          units = "in",
          res = opts$plot_dpi,
          bg = "white"
        )
      }
    )
  }

  if (!is.null(pdf_path)) {
    save_plot(
      pdf_path,
      function(file_path) {
        grDevices::pdf(
          file = file_path,
          width = opts$summary_panel_width,
          height = opts$summary_panel_height,
          useDingbats = FALSE,
          bg = "white"
        )
      }
    )
  }

  invisible(TRUE)
}

run_top50_snv_circular_plot <- function(opts) {
  validate_opts(opts)

  opts$input <- check_input_csv(opts$input)
  opts$output_dir <- resolve_output_dir(opts$input, opts$output_dir)
  opts$output_dir <- check_output_dir(opts$output_dir)

  attach_required_packages()

  message(sprintf("[INFO] Reading input CSV: %s", opts$input))
  input_df <- readr::read_csv(
    file = opts$input,
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "minimal"
  )
  input_df <- as.data.frame(input_df, stringsAsFactors = FALSE)
  input_df <- prepare_input_table(input_df, opts$score_column)

  message(sprintf("[INFO] Loaded %d valid rows across %d unique SNVs.", nrow(input_df), length(unique(input_df$SNV))))

  ranked_summary <- build_ranked_snv_summary(
    input_df,
    score_upper_quantile = opts$score_ring_upper_quantile
  )
  if (nrow(ranked_summary) < opts$top_n) {
    message(sprintf("[INFO] Only %d unique SNVs available; plotting all of them.", nrow(ranked_summary)))
  }
  invisible(read_alphagenome_atac_metadata(opts$alphagenome_atac_metadata_csv))
  full_numeric_summary <- append_alphagenome_boxplot_summary(
    summary_df = ranked_summary,
    boxplot_stats_csv = opts$alphagenome_boxplot_stats_csv
  )
  top_summary <- build_top_snv_summary(
    input_df,
    top_n = opts$top_n,
    score_upper_quantile = opts$score_ring_upper_quantile
  )
  top_summary <- full_numeric_summary[match(top_summary$SNV, full_numeric_summary$SNV), , drop = FALSE]
  top_summary <- append_regulationspotter_summary(
    summary_df = top_summary,
    all_vars_tsv = opts$regulationspotter_all_vars_tsv,
    spotter_tsv = opts$regulationspotter_spotter_tsv
  )
  top_summary <- append_regulationtaster_summary(
    summary_df = top_summary,
    all_vars_tsv = opts$regulationspotter_all_vars_tsv,
    taster_tsv = opts$regulationspotter_taster_tsv
  )
  top_summary <- append_snpeff_impact_summary(
    summary_df = top_summary,
    snpeff_impact_summary_tsv = opts$snpeff_impact_summary_tsv
  )
  top_summary <- append_vep_impact_summary(
    summary_df = top_summary,
    vep_annotation_csv = opts$vep_annotation_csv
  )
  top_summary <- append_mutationtaster_prediction_summary(
    summary_df = top_summary,
    mutationtaster_annotation_tsv = opts$mutationtaster_annotation_tsv
  )

  summary_csv <- file.path(opts$output_dir, paste0(opts$output_stem, "_summary.csv"))
  png_path <- file.path(opts$output_dir, paste0(opts$output_stem, ".png"))
  pdf_path <- file.path(opts$output_dir, paste0(opts$output_stem, ".pdf"))
  summary_panel_png_path <- file.path(opts$output_dir, paste0(opts$output_stem, "_summary_panel.png"))
  summary_panel_pdf_path <- file.path(opts$output_dir, paste0(opts$output_stem, "_summary_panel.pdf"))

  readr::write_csv(top_summary, summary_csv)
  message(sprintf("[INFO] Saved top-SNV summary to %s", summary_csv))

  draw_circular_plot(
    summary_df = top_summary,
    opts = opts,
    numeric_summary_df = full_numeric_summary,
    png_path = png_path,
    pdf_path = pdf_path
  )
  draw_summary_panel(
    summary_df = top_summary,
    opts = opts,
    numeric_summary_df = full_numeric_summary,
    png_path = summary_panel_png_path,
    pdf_path = summary_panel_pdf_path
  )

  message(sprintf("[INFO] Saved PNG plot to %s", png_path))
  message(sprintf("[INFO] Saved PDF plot to %s", pdf_path))
  message(sprintf("[INFO] Saved summary panel PNG plot to %s", summary_panel_png_path))
  message(sprintf("[INFO] Saved summary panel PDF plot to %s", summary_panel_pdf_path))
  message(sprintf("[INFO] Completed circular plot generation under %s", opts$output_dir))
}

opts <- list(
  input = INPUT_CSV,
  output_dir = OUTPUT_DIR,
  score_column = SCORE_COLUMN,
  top_n = as.integer(TOP_N),
  output_stem = OUTPUT_STEM,
  alphagenome_boxplot_stats_csv = ALPHAGENOME_BOXPLOT_STATS_CSV,
  alphagenome_atac_metadata_csv = ALPHAGENOME_ATAC_METADATA_CSV,
  snpeff_impact_summary_tsv = SNPEFF_IMPACT_SUMMARY_TSV,
  regulationspotter_all_vars_tsv = REGULATIONSPOTTER_ALL_VARS_TSV,
  regulationspotter_spotter_tsv = REGULATIONSPOTTER_SPOTTER_TSV,
  regulationspotter_taster_tsv = REGULATIONSPOTTER_TASTER_TSV,
  vep_annotation_csv = VEP_ANNOTATION_CSV,
  mutationtaster_annotation_tsv = MUTATIONTASTER_ANNOTATION_TSV,
  plot_width = PLOT_WIDTH,
  plot_height = PLOT_HEIGHT,
  summary_panel_width = SUMMARY_PANEL_WIDTH,
  summary_panel_height = SUMMARY_PANEL_HEIGHT,
  summary_panel_legend_panel_width = SUMMARY_PANEL_LEGEND_PANEL_WIDTH,
  summary_panel_legend_right_margin = SUMMARY_PANEL_LEGEND_RIGHT_MARGIN,
  summary_panel_legend_x_offsets = SUMMARY_PANEL_LEGEND_X_OFFSETS,
  summary_panel_legend_text_cex = SUMMARY_PANEL_LEGEND_TEXT_CEX,
  summary_panel_legend_y_intersp = SUMMARY_PANEL_LEGEND_Y_INTERSP,
  summary_panel_legend_block_gap = SUMMARY_PANEL_LEGEND_BLOCK_GAP,
  summary_panel_non_top_point_color = SUMMARY_PANEL_NON_TOP_POINT_COLOR,
  summary_panel_numeric_box_width = SUMMARY_PANEL_NUMERIC_BOX_WIDTH,
  summary_panel_numeric_point_width = SUMMARY_PANEL_NUMERIC_POINT_WIDTH,
  summary_panel_category_bar_width = SUMMARY_PANEL_CATEGORY_BAR_WIDTH,
  summary_panel_numeric_column_width_factor = SUMMARY_PANEL_NUMERIC_COLUMN_WIDTH_FACTOR,
  summary_panel_category_label_cex = SUMMARY_PANEL_CATEGORY_LABEL_CEX,
  show_outer_ring_value_labels = SHOW_OUTER_RING_VALUE_LABELS,
  outer_ring_value_axis_side = OUTER_RING_VALUE_AXIS_SIDE,
  outer_ring_value_label_cex = OUTER_RING_VALUE_LABEL_CEX,
  outer_ring_value_label_color = OUTER_RING_VALUE_LABEL_COLOR,
  score_ring_upper_quantile = SCORE_RING_UPPER_QUANTILE,
  plot_dpi = PLOT_DPI,
  plot_title = PLOT_TITLE,
  show_legends = SHOW_LEGENDS,
  start_degree = START_DEGREE,
  sector_gap_degree = SECTOR_GAP_DEGREE,
  final_gap_degree = FINAL_GAP_DEGREE,
  label_cex = LABEL_CEX,
  font_family = FIGURE_FONT_FAMILY,
  score_ring_color = SCORE_RING_COLOR,
  score_ring_baseline_color = SCORE_RING_BASELINE_COLOR,
  ring_border_color = RING_BORDER_COLOR,
  ring_background_color = RING_BACKGROUND_COLOR,
  label_track_height = LABEL_TRACK_HEIGHT,
  score_ring_track_height = SCORE_RING_TRACK_HEIGHT,
  alphagenome_rna_track_height = ALPHAGENOME_RNA_TRACK_HEIGHT,
  alphagenome_atac_track_height = ALPHAGENOME_ATAC_TRACK_HEIGHT,
  snpeff_impact_track_height = SNPEFF_IMPACT_TRACK_HEIGHT,
  func_track_height = FUNC_TRACK_HEIGHT,
  vep_impact_track_height = VEP_IMPACT_TRACK_HEIGHT,
  mutationtaster_track_height = MUTATIONTASTER_TRACK_HEIGHT,
  alphagenome_box_width_fraction = ALPHAGENOME_BOX_WIDTH_FRACTION
)

run_top50_snv_circular_plot(opts)
