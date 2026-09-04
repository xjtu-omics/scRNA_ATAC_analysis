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

script_dir <- tryCatch(
  dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]))),
  error = function(e) getwd()
)

three_end_path <- if (length(args) >= 1) {
  args[[1]]
} else {
  Sys.getenv(
    "PRISMSNV_THREE_END_PATIENT_SNV_PATH",
    unset = file.path(script_dir, "input", "3-end_private_snv_scores_by_patient_celltype.csv")
  )
}

full_length_path <- if (length(args) >= 2) {
  args[[2]]
} else {
  Sys.getenv(
    "PRISMSNV_FULL_LENGTH_PATIENT_SNV_PATH",
    unset = file.path(script_dir, "input", "full-length_private_snv_scores_by_patient_celltype.csv")
  )
}

output_dir <- if (length(args) >= 3) {
  args[[3]]
} else {
  file.path(script_dir, "private_snv_gene_celltype_recurrence")
}

n_bootstrap <- if (length(args) >= 4) as.integer(args[[4]]) else 2000L
n_null <- if (length(args) >= 5) as.integer(args[[5]]) else 2000L
seed <- if (length(args) >= 6) as.integer(args[[6]]) else 20260901L
mode <- if (length(args) >= 7) tolower(args[[7]]) else "auto"
min_patients_per_cohort <- if (length(args) >= 8) as.integer(args[[8]]) else 1L
match_k <- if (length(args) >= 9) as.integer(args[[9]]) else 10L

if (!mode %in% c("auto", "strict", "audit")) {
  stop("Mode must be one of: auto, strict, audit.", call. = FALSE)
}
if (!is.finite(n_bootstrap) || n_bootstrap < 1) {
  stop("n_bootstrap must be a positive integer.", call. = FALSE)
}
if (!is.finite(n_null) || n_null < 1) {
  stop("n_null must be a positive integer.", call. = FALSE)
}
if (!is.finite(min_patients_per_cohort) || min_patients_per_cohort < 1) {
  stop("min_patients_per_cohort must be a positive integer.", call. = FALSE)
}
if (!is.finite(match_k) || match_k < 1) {
  stop("match_k must be a positive integer.", call. = FALSE)
}

set.seed(seed)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
table_dir <- file.path(output_dir, "tables")
plot_dir <- file.path(output_dir, "plots")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

field_aliases <- list(
  patient_id = c("patient_id", "patient", "donor_id", "donor", "individual_id", "subject_id"),
  cell_type = c("cell_type", "celltype", "cell_class", "cell_cluster"),
  snv = c("snv", "SNV", "snv_id", "variant", "variant_id"),
  gene = c("gene", "Gene", "gene_symbol", "symbol"),
  score = c(
    "perturb_euclidean_distance",
    "score_euclidean",
    "perturbation_score_euclidean",
    "euclidean_distance",
    "score"
  ),
  gene_expression = c(
    "gene_expression",
    "mean_gene_expression",
    "expression_mean",
    "mean_expression",
    "avg_expression"
  ),
  gene_length = c("gene_length", "gene_length_bp", "length_bp"),
  callable_bases = c("callable_bases", "callable_bp", "n_callable_bases"),
  covered_cell_count = c(
    "covered_cell_count",
    "n_covered_cells",
    "coverage_cell_count",
    "n_cells_covered"
  ),
  carrier_cell_count = c(
    "carrier_cell_count",
    "n_carrier_cells",
    "alt_cell_count",
    "n_alt_cells"
  ),
  sample_id = c("sample_id", "sample", "library_id", "library")
)

basic_fields <- c("cell_type", "snv", "gene", "score")
strict_fields <- c(
  "patient_id",
  basic_fields,
  "gene_expression",
  "gene_length",
  "callable_bases",
  "covered_cell_count",
  "carrier_cell_count"
)
matching_covariates <- c(
  "gene_expression",
  "gene_length",
  "callable_bases",
  "n_snvs",
  "covered_cell_count",
  "carrier_cell_count"
)

check_file_exists <- function(path) {
  if (!file.exists(path)) {
    stop("Required input file does not exist: ", path, call. = FALSE)
  }
}

check_file_exists(three_end_path)
check_file_exists(full_length_path)

first_alias <- function(columns, aliases) {
  hit <- aliases[aliases %in% columns]
  if (length(hit) == 0) NA_character_ else hit[[1]]
}

make_field_map <- function(columns) {
  tibble(
    canonical_field = names(field_aliases),
    source_field = map_chr(field_aliases, ~ first_alias(columns, .x)),
    available = !is.na(.data$source_field)
  )
}

normalize_gene <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    str_replace("\\..*$", "") %>%
    str_to_upper()
}

split_gene_annotation <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  str_split(x, "[;,|]")
}

normalize_snv <- function(x) {
  original <- as.character(x)
  cleaned <- original %>%
    str_trim() %>%
    str_to_upper() %>%
    str_replace_all("\\s+", "") %>%
    str_replace("^CHR", "")

  parsed <- str_match(cleaned, "^([^:]+):([0-9]+)[_:]([ACGTN]+)>([ACGTN]+)$")
  canonical <- ifelse(
    !is.na(parsed[, 1]),
    paste0("CHR", parsed[, 2], ":", parsed[, 3], ":", parsed[, 4], ">", parsed[, 5]),
    cleaned
  )
  canonical[canonical == ""] <- NA_character_
  canonical
}

safe_spearman <- function(x, y, min_n = 3L) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < min_n || length(unique(x[keep])) < 2 || length(unique(y[keep])) < 2) {
    return(NA_real_)
  }
  suppressWarnings(cor(x[keep], y[keep], method = "spearman"))
}

safe_spearman_p <- function(x, y, min_n = 3L) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < min_n || length(unique(x[keep])) < 2 || length(unique(y[keep])) < 2) {
    return(NA_real_)
  }
  suppressWarnings(cor.test(x[keep], y[keep], method = "spearman", exact = FALSE)$p.value)
}

safe_quantile <- function(x, probability) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  as.numeric(quantile(x, probability, names = FALSE, type = 7))
}

percentile_rank_high <- function(x) {
  valid <- is.finite(x)
  out <- rep(NA_real_, length(x))
  if (!any(valid)) return(out)
  n <- sum(valid)
  if (n == 1) {
    out[valid] <- 1
  } else {
    out[valid] <- (rank(x[valid], ties.method = "average") - 1) / (n - 1)
  }
  out
}

standardize_input <- function(path, platform_label) {
  raw <- read_csv(path, col_types = cols(.default = col_character()), show_col_types = FALSE)
  mapping <- make_field_map(names(raw))
  names(mapping)[names(mapping) == "source_field"] <- paste0("source_field_", platform_label)
  names(mapping)[names(mapping) == "available"] <- paste0("available_", platform_label)

  source_lookup <- setNames(
    mapping[[paste0("source_field_", platform_label)]],
    mapping$canonical_field
  )

  output <- tibble(.source_row = seq_len(nrow(raw)))
  for (field in names(field_aliases)) {
    source_name <- source_lookup[[field]]
    if (!is.na(source_name)) {
      output[[field]] <- raw[[source_name]]
    } else {
      output[[field]] <- NA_character_
    }
  }

  output <- output %>%
    mutate(
      platform = platform_label,
      patient_id = str_trim(as.character(.data$patient_id)),
      sample_id = str_trim(as.character(.data$sample_id)),
      cell_type = str_trim(as.character(.data$cell_type)),
      snv_original = str_trim(as.character(.data$snv)),
      snv = normalize_snv(.data$snv_original),
      gene_list = split_gene_annotation(.data$gene),
      score = suppressWarnings(as.numeric(.data$score)),
      gene_expression = suppressWarnings(as.numeric(.data$gene_expression)),
      gene_length = suppressWarnings(as.numeric(.data$gene_length)),
      callable_bases = suppressWarnings(as.numeric(.data$callable_bases)),
      covered_cell_count = suppressWarnings(as.numeric(.data$covered_cell_count)),
      carrier_cell_count = suppressWarnings(as.numeric(.data$carrier_cell_count))
    ) %>%
    select(-"gene") %>%
    unnest_longer("gene_list", values_to = "gene") %>%
    mutate(gene = normalize_gene(.data$gene)) %>%
    filter(
      !is.na(.data$snv),
      !is.na(.data$gene),
      .data$gene != "",
      !is.na(.data$cell_type),
      .data$cell_type != "",
      is.finite(.data$score)
    )

  list(data = output, mapping = mapping, raw_columns = names(raw), raw_rows = nrow(raw))
}

three <- standardize_input(three_end_path, "three_end")
full <- standardize_input(full_length_path, "full_length")

field_audit <- full_join(three$mapping, full$mapping, by = "canonical_field") %>%
  mutate(
    required_for = case_when(
      .data$canonical_field %in% basic_fields ~ "descriptive_and_strict",
      .data$canonical_field %in% strict_fields ~ "strict",
      TRUE ~ "optional"
    ),
    accepted_aliases = map_chr(.data$canonical_field, ~ paste(field_aliases[[.x]], collapse = ";"))
  ) %>%
  select(
    "canonical_field",
    "required_for",
    "accepted_aliases",
    "source_field_three_end",
    "available_three_end",
    "source_field_full_length",
    "available_full_length"
  )

write_csv(field_audit, file.path(table_dir, "input_schema_audit.csv"))

combined_basic <- bind_rows(three$data, full$data)

cohort_basic_summary <- combined_basic %>%
  group_by(.data$platform) %>%
  summarise(
    input_rows_after_gene_expansion = n(),
    n_unique_snvs = n_distinct(.data$snv),
    n_unique_genes = n_distinct(.data$gene),
    n_gene_celltype_pairs = n_distinct(paste(.data$gene, .data$cell_type, sep = "||")),
    .groups = "drop"
  )

snvs_three <- sort(unique(three$data$snv))
snvs_full <- sort(unique(full$data$snv))
genes_three <- sort(unique(three$data$gene))
genes_full <- sort(unique(full$data$gene))
shared_snvs <- intersect(snvs_three, snvs_full)
shared_genes <- intersect(genes_three, genes_full)

snv_overlap_summary <- tibble(
  n_three_end = length(snvs_three),
  n_full_length = length(snvs_full),
  n_shared = length(shared_snvs),
  union_size = length(union(snvs_three, snvs_full)),
  jaccard = ifelse(length(union(snvs_three, snvs_full)) > 0, length(shared_snvs) / length(union(snvs_three, snvs_full)), NA_real_),
  overlap_fraction_three_end = ifelse(length(snvs_three) > 0, length(shared_snvs) / length(snvs_three), NA_real_),
  overlap_fraction_full_length = ifelse(length(snvs_full) > 0, length(shared_snvs) / length(snvs_full), NA_real_)
)

gene_overlap_summary <- tibble(
  n_three_end = length(genes_three),
  n_full_length = length(genes_full),
  n_shared = length(shared_genes),
  union_size = length(union(genes_three, genes_full)),
  jaccard = ifelse(length(union(genes_three, genes_full)) > 0, length(shared_genes) / length(union(genes_three, genes_full)), NA_real_),
  overlap_fraction_three_end = ifelse(length(genes_three) > 0, length(shared_genes) / length(genes_three), NA_real_),
  overlap_fraction_full_length = ifelse(length(genes_full) > 0, length(shared_genes) / length(genes_full), NA_real_)
)

write_csv(cohort_basic_summary, file.path(table_dir, "descriptive_cohort_summary.csv"))
write_csv(snv_overlap_summary, file.path(table_dir, "descriptive_exact_snv_overlap_summary.csv"))
write_csv(tibble(snv = shared_snvs), file.path(table_dir, "descriptive_shared_exact_snvs.csv"))
write_csv(gene_overlap_summary, file.path(table_dir, "descriptive_gene_overlap_summary.csv"))
write_csv(tibble(gene = shared_genes), file.path(table_dir, "descriptive_shared_genes.csv"))

strict_missing <- field_audit %>%
  filter(
    .data$canonical_field %in% strict_fields,
    !coalesce(.data$available_three_end, FALSE) | !coalesce(.data$available_full_length, FALSE)
  )

patients_three <- sort(unique(three$data$patient_id[!is.na(three$data$patient_id) & three$data$patient_id != ""]))
patients_full <- sort(unique(full$data$patient_id[!is.na(full$data$patient_id) & full$data$patient_id != ""]))
patient_overlap <- intersect(patients_three, patients_full)

patient_overlap_summary <- tibble(
  n_three_end = length(patients_three),
  n_full_length = length(patients_full),
  n_shared = length(patient_overlap),
  patient_sets_disjoint = length(patients_three) > 0 && length(patients_full) > 0 && length(patient_overlap) == 0,
  three_end_patients = paste(patients_three, collapse = ";"),
  full_length_patients = paste(patients_full, collapse = ";"),
  shared_patients = paste(patient_overlap, collapse = ";")
)
write_csv(patient_overlap_summary, file.path(table_dir, "patient_overlap_summary.csv"))

write_status <- function(status, detail_lines) {
  lines <- c(
    paste0("status=", status),
    paste0("three_end_input=", normalizePath(three_end_path, winslash = "/", mustWork = FALSE)),
    paste0("full_length_input=", normalizePath(full_length_path, winslash = "/", mustWork = FALSE)),
    paste0("output_dir=", normalizePath(output_dir, winslash = "/", mustWork = FALSE)),
    paste0("seed=", seed),
    paste0("n_bootstrap=", n_bootstrap),
    paste0("n_null=", n_null),
    paste0("min_patients_per_cohort=", min_patients_per_cohort),
    detail_lines
  )
  writeLines(lines, file.path(output_dir, "analysis_status.txt"), useBytes = TRUE)
}

if (mode == "audit" || nrow(strict_missing) > 0) {
  missing_text <- if (nrow(strict_missing) > 0) {
    paste(strict_missing$canonical_field, collapse = ",")
  } else {
    "none"
  }

  status <- if (nrow(strict_missing) > 0) {
    "BLOCKED_MISSING_PATIENT_LEVEL_FIELDS"
  } else {
    "AUDIT_ONLY"
  }

  write_status(
    status,
    c(
      paste0("missing_strict_fields=", missing_text),
      paste0("descriptive_exact_snv_overlap=", length(shared_snvs)),
      paste0("descriptive_exact_snv_jaccard=", format(snv_overlap_summary$jaccard, scientific = TRUE, digits = 8)),
      paste0("descriptive_shared_genes=", length(shared_genes)),
      "patient_bootstrap_run=FALSE",
      "matched_null_run=FALSE",
      "note=Patient bootstrap and matched-null inference require patient-level SNV rows and all six matching covariates."
    )
  )

  message("Descriptive audit outputs were written to: ", output_dir)
  if (nrow(strict_missing) > 0) {
    message("Strict inference was not run. Missing fields: ", missing_text)
  }

  if (mode == "strict" && nrow(strict_missing) > 0) {
    quit(save = "no", status = 2L)
  }
  quit(save = "no", status = 0L)
}

if (length(patients_three) == 0 || length(patients_full) == 0) {
  stop("Patient identifiers are present as columns but contain no usable values.", call. = FALSE)
}
if (length(patient_overlap) > 0) {
  stop(
    "Patient cohorts are not disjoint. Shared patient IDs: ",
    paste(patient_overlap, collapse = ", "),
    call. = FALSE
  )
}

strict_data <- combined_basic %>%
  filter(
    !is.na(.data$patient_id),
    .data$patient_id != "",
    is.finite(.data$gene_expression),
    is.finite(.data$gene_length),
    is.finite(.data$callable_bases),
    is.finite(.data$covered_cell_count),
    is.finite(.data$carrier_cell_count),
    .data$gene_expression >= 0,
    .data$gene_length > 0,
    .data$callable_bases >= 0,
    .data$covered_cell_count >= 0,
    .data$carrier_cell_count >= 0
  )

if (nrow(strict_data) == 0) {
  stop("No rows remain after validating patient IDs and matching covariates.", call. = FALSE)
}

private_snv_audit <- strict_data %>%
  distinct(.data$platform, .data$patient_id, .data$snv) %>%
  group_by(.data$platform, .data$snv) %>%
  summarise(
    n_patients = n_distinct(.data$patient_id),
    patient_ids = paste(sort(unique(.data$patient_id)), collapse = ";"),
    is_private_within_cohort = .data$n_patients == 1,
    .groups = "drop"
  )

write_csv(private_snv_audit, file.path(table_dir, "private_snv_audit.csv"))

private_data <- strict_data %>%
  inner_join(
    private_snv_audit %>% filter(.data$is_private_within_cohort) %>% select("platform", "snv"),
    by = c("platform", "snv")
  )

private_cohort_summary <- private_data %>%
  group_by(.data$platform) %>%
  summarise(
    n_patients = n_distinct(.data$patient_id),
    n_private_snvs = n_distinct(.data$snv),
    n_private_snv_genes = n_distinct(.data$gene),
    .groups = "drop"
  )
write_csv(private_cohort_summary, file.path(table_dir, "private_cohort_summary.csv"))

private_snvs_three <- private_data %>% filter(.data$platform == "three_end") %>% pull(.data$snv) %>% unique()
private_snvs_full <- private_data %>% filter(.data$platform == "full_length") %>% pull(.data$snv) %>% unique()
private_shared_snvs <- sort(intersect(private_snvs_three, private_snvs_full))
private_snv_union <- union(private_snvs_three, private_snvs_full)

private_exact_overlap <- tibble(
  n_three_end_private_snvs = length(private_snvs_three),
  n_full_length_private_snvs = length(private_snvs_full),
  n_shared_exact_private_snvs = length(private_shared_snvs),
  union_size = length(private_snv_union),
  jaccard = ifelse(length(private_snv_union) > 0, length(private_shared_snvs) / length(private_snv_union), NA_real_)
)
write_csv(private_exact_overlap, file.path(table_dir, "private_exact_snv_overlap_summary.csv"))
write_csv(tibble(snv = private_shared_snvs), file.path(table_dir, "private_shared_exact_snvs.csv"))

private_genes_three <- private_data %>% filter(.data$platform == "three_end") %>% pull(.data$gene) %>% unique()
private_genes_full <- private_data %>% filter(.data$platform == "full_length") %>% pull(.data$gene) %>% unique()
common_private_genes <- sort(intersect(private_genes_three, private_genes_full))

common_gene_audit <- private_data %>%
  filter(.data$gene %in% common_private_genes) %>%
  group_by(.data$platform, .data$gene) %>%
  summarise(
    n_patients = n_distinct(.data$patient_id),
    n_private_snvs = n_distinct(.data$snv),
    patient_ids = paste(sort(unique(.data$patient_id)), collapse = ";"),
    private_snvs = paste(sort(unique(.data$snv)), collapse = ";"),
    .groups = "drop"
  )
write_csv(common_gene_audit, file.path(table_dir, "common_private_gene_audit.csv"))

patient_snv_scores <- private_data %>%
  group_by(.data$platform, .data$patient_id, .data$cell_type, .data$gene, .data$snv) %>%
  summarise(
    snv_score = median(.data$score, na.rm = TRUE),
    gene_expression = median(.data$gene_expression, na.rm = TRUE),
    gene_length = median(.data$gene_length, na.rm = TRUE),
    callable_bases = median(.data$callable_bases, na.rm = TRUE),
    covered_cell_count = median(.data$covered_cell_count, na.rm = TRUE),
    carrier_cell_count = median(.data$carrier_cell_count, na.rm = TRUE),
    n_source_rows = n(),
    .groups = "drop"
  )

patient_gene_celltype <- patient_snv_scores %>%
  group_by(.data$platform, .data$patient_id, .data$cell_type, .data$gene) %>%
  summarise(
    gene_celltype_score = median(.data$snv_score, na.rm = TRUE),
    n_snvs = n_distinct(.data$snv),
    gene_expression = median(.data$gene_expression, na.rm = TRUE),
    gene_length = median(.data$gene_length, na.rm = TRUE),
    callable_bases = median(.data$callable_bases, na.rm = TRUE),
    covered_cell_count = sum(.data$covered_cell_count, na.rm = TRUE),
    carrier_cell_count = sum(.data$carrier_cell_count, na.rm = TRUE),
    private_snvs = paste(sort(unique(.data$snv)), collapse = ";"),
    .groups = "drop"
  ) %>%
  group_by(.data$platform, .data$patient_id, .data$cell_type) %>%
  mutate(gene_celltype_rank = percentile_rank_high(.data$gene_celltype_score)) %>%
  ungroup()

write_csv(patient_snv_scores, file.path(table_dir, "patient_snv_celltype_scores.csv"))
write_csv(patient_gene_celltype, file.path(table_dir, "patient_gene_celltype_scores.csv"))

cohort_gene_celltype <- patient_gene_celltype %>%
  group_by(.data$platform, .data$cell_type, .data$gene) %>%
  summarise(
    n_patients = n_distinct(.data$patient_id),
    cohort_rank_score = median(.data$gene_celltype_rank, na.rm = TRUE),
    cohort_raw_score = median(.data$gene_celltype_score, na.rm = TRUE),
    gene_expression = median(.data$gene_expression, na.rm = TRUE),
    gene_length = median(.data$gene_length, na.rm = TRUE),
    callable_bases = median(.data$callable_bases, na.rm = TRUE),
    n_snvs = sum(.data$n_snvs, na.rm = TRUE),
    covered_cell_count = sum(.data$covered_cell_count, na.rm = TRUE),
    carrier_cell_count = sum(.data$carrier_cell_count, na.rm = TRUE),
    patient_ids = paste(sort(unique(.data$patient_id)), collapse = ";"),
    .groups = "drop"
  )

write_csv(cohort_gene_celltype, file.path(table_dir, "cohort_gene_celltype_scores.csv"))

cohort_three <- cohort_gene_celltype %>%
  filter(
    .data$platform == "three_end",
    .data$gene %in% common_private_genes,
    .data$n_patients >= min_patients_per_cohort
  ) %>%
  select(-"platform")
cohort_full <- cohort_gene_celltype %>%
  filter(
    .data$platform == "full_length",
    .data$gene %in% common_private_genes,
    .data$n_patients >= min_patients_per_cohort
  ) %>%
  select(-"platform")

cohort_full_null_candidates <- cohort_gene_celltype %>%
  filter(
    .data$platform == "full_length",
    .data$n_patients >= min_patients_per_cohort
  ) %>%
  select(-"platform")

shared_pairs <- inner_join(
  cohort_three,
  cohort_full,
  by = c("cell_type", "gene"),
  suffix = c("_three_end", "_full_length")
) %>%
  filter(
    .data$gene %in% common_private_genes,
    is.finite(.data$cohort_rank_score_three_end),
    is.finite(.data$cohort_rank_score_full_length)
  )

if (nrow(shared_pairs) < 3) {
  stop(
    "Fewer than three shared gene-cell-type pairs pass the patient-level filters. ",
    "Lower min_patients_per_cohort only for an explicitly descriptive sensitivity analysis.",
    call. = FALSE
  )
}

write_csv(shared_pairs, file.path(table_dir, "shared_private_gene_celltype_pairs.csv"))

observed_correlations <- bind_rows(
  tibble(
    scope = "all_pairs",
    cell_type = "All",
    n_pairs = nrow(shared_pairs),
    spearman_rho = safe_spearman(shared_pairs$cohort_rank_score_three_end, shared_pairs$cohort_rank_score_full_length),
    p_value = safe_spearman_p(shared_pairs$cohort_rank_score_three_end, shared_pairs$cohort_rank_score_full_length)
  ),
  shared_pairs %>%
    group_by(.data$cell_type) %>%
    summarise(
      scope = "within_cell_type",
      n_pairs = n(),
      spearman_rho = safe_spearman(.data$cohort_rank_score_three_end, .data$cohort_rank_score_full_length),
      p_value = safe_spearman_p(.data$cohort_rank_score_three_end, .data$cohort_rank_score_full_length),
      .groups = "drop"
    ) %>%
    select("scope", "cell_type", "n_pairs", "spearman_rho", "p_value")
)

aggregate_bootstrap_cohort <- function(patient_df, sampled_patients) {
  weights <- tibble(patient_id = sampled_patients) %>%
    count(.data$patient_id, name = "bootstrap_weight")

  patient_df %>%
    inner_join(weights, by = "patient_id") %>%
    uncount(.data$bootstrap_weight) %>%
    group_by(.data$cell_type, .data$gene) %>%
    summarise(cohort_rank_score = median(.data$gene_celltype_rank, na.rm = TRUE), .groups = "drop")
}

bootstrap_one <- function(iteration) {
  sample_three <- sample(patients_three, length(patients_three), replace = TRUE)
  sample_full <- sample(patients_full, length(patients_full), replace = TRUE)

  boot_three <- aggregate_bootstrap_cohort(
    patient_gene_celltype %>% filter(.data$platform == "three_end"),
    sample_three
  ) %>%
    filter(.data$gene %in% common_private_genes)
  boot_full <- aggregate_bootstrap_cohort(
    patient_gene_celltype %>% filter(.data$platform == "full_length"),
    sample_full
  ) %>%
    filter(.data$gene %in% common_private_genes)

  boot_pairs <- shared_pairs %>%
    select("cell_type", "gene") %>%
    inner_join(boot_three, by = c("cell_type", "gene")) %>%
    rename(score_three_end = "cohort_rank_score") %>%
    inner_join(boot_full, by = c("cell_type", "gene")) %>%
    rename(score_full_length = "cohort_rank_score")

  bind_rows(
    tibble(
      iteration = iteration,
      scope = "all_pairs",
      cell_type = "All",
      n_pairs = nrow(boot_pairs),
      spearman_rho = safe_spearman(boot_pairs$score_three_end, boot_pairs$score_full_length)
    ),
    boot_pairs %>%
      group_by(.data$cell_type) %>%
      summarise(
        iteration = iteration,
        scope = "within_cell_type",
        n_pairs = n(),
        spearman_rho = safe_spearman(.data$score_three_end, .data$score_full_length),
        .groups = "drop"
      ) %>%
      select("iteration", "scope", "cell_type", "n_pairs", "spearman_rho")
  )
}

bootstrap_correlations <- map_dfr(seq_len(n_bootstrap), bootstrap_one)
write_csv(bootstrap_correlations, file.path(table_dir, "patient_bootstrap_spearman.csv"))

bootstrap_ci <- bootstrap_correlations %>%
  group_by(.data$scope, .data$cell_type) %>%
  summarise(
    n_bootstrap_requested = n_bootstrap,
    n_bootstrap_valid = sum(is.finite(.data$spearman_rho)),
    valid_fraction = .data$n_bootstrap_valid / .data$n_bootstrap_requested,
    ci_low = safe_quantile(.data$spearman_rho, 0.025),
    ci_median = safe_quantile(.data$spearman_rho, 0.5),
    ci_high = safe_quantile(.data$spearman_rho, 0.975),
    .groups = "drop"
  )

correlation_summary <- observed_correlations %>%
  left_join(bootstrap_ci, by = c("scope", "cell_type"))
write_csv(correlation_summary, file.path(table_dir, "patient_bootstrap_correlation_summary.csv"))

prepare_matching_space <- function(target_df, candidate_df) {
  target <- target_df %>%
    transmute(
      cohort = "three_end",
      cell_type = .data$cell_type,
      gene = .data$gene,
      score = .data$cohort_rank_score_three_end,
      gene_expression = .data$gene_expression_three_end,
      gene_length = .data$gene_length_three_end,
      callable_bases = .data$callable_bases_three_end,
      n_snvs = .data$n_snvs_three_end,
      covered_cell_count = .data$covered_cell_count_three_end,
      carrier_cell_count = .data$carrier_cell_count_three_end
    )

  candidates <- candidate_df %>%
    transmute(
      cohort = "full_length",
      cell_type = .data$cell_type,
      gene = .data$gene,
      score = .data$cohort_rank_score,
      across(all_of(matching_covariates), identity)
    )

  combined <- bind_rows(target, candidates)
  for (covariate in matching_covariates) {
    transformed <- log1p(pmax(combined[[covariate]], 0))
    center <- mean(transformed, na.rm = TRUE)
    spread <- sd(transformed, na.rm = TRUE)
    if (!is.finite(spread) || spread == 0) spread <- 1
    combined[[paste0(covariate, "_z")]] <- (transformed - center) / spread
  }

  list(
    target = combined %>% filter(.data$cohort == "three_end"),
    candidates = combined %>% filter(.data$cohort == "full_length")
  )
}

candidate_full <- cohort_full_null_candidates %>%
  filter(if_all(all_of(c("cohort_rank_score", matching_covariates)), is.finite))
matching_space <- prepare_matching_space(shared_pairs, candidate_full)
z_columns <- paste0(matching_covariates, "_z")

match_one_null <- function(iteration) {
  target <- matching_space$target
  candidates <- matching_space$candidates
  target_order <- sample(seq_len(nrow(target)))
  used_candidate_rows <- integer(0)
  matched_rows <- vector("list", nrow(target))

  for (step in seq_along(target_order)) {
    target_index <- target_order[[step]]
    target_row <- target[target_index, , drop = FALSE]
    candidate_indices <- which(
      candidates$cell_type == target_row$cell_type &
        candidates$gene != target_row$gene &
        !seq_len(nrow(candidates)) %in% used_candidate_rows
    )

    reused <- FALSE
    if (length(candidate_indices) == 0) {
      candidate_indices <- which(
        candidates$cell_type == target_row$cell_type & candidates$gene != target_row$gene
      )
      reused <- TRUE
    }
    if (length(candidate_indices) == 0) next

    target_vector <- as.numeric(target_row[1, z_columns, drop = TRUE])
    candidate_matrix <- as.matrix(candidates[candidate_indices, z_columns, drop = FALSE])
    distances <- sqrt(rowSums((sweep(candidate_matrix, 2, target_vector, "-"))^2))
    nearest_order <- order(distances)
    nearest_order <- head(nearest_order, min(match_k, length(nearest_order)))
    nearest_indices <- candidate_indices[nearest_order]
    nearest_distances <- distances[nearest_order]
    temperature <- median(nearest_distances[is.finite(nearest_distances)])
    if (!is.finite(temperature) || temperature <= 0) temperature <- 1
    probabilities <- exp(-nearest_distances / temperature)
    probabilities <- probabilities / sum(probabilities)
    chosen_position <- sample.int(length(nearest_indices), 1, prob = probabilities)
    chosen <- nearest_indices[[chosen_position]]
    if (!reused) used_candidate_rows <- c(used_candidate_rows, chosen)

    balance <- as.numeric(candidates[chosen, z_columns, drop = TRUE]) - target_vector
    matched_rows[[target_index]] <- tibble(
      target_gene = target_row$gene,
      matched_gene = candidates$gene[[chosen]],
      cell_type = target_row$cell_type,
      score_three_end = target_row$score,
      score_full_length = candidates$score[[chosen]],
      match_distance = distances[match(chosen, candidate_indices)],
      reused_candidate = reused,
      mean_abs_standardized_difference = mean(abs(balance), na.rm = TRUE)
    )
  }

  matched <- bind_rows(matched_rows)
  tibble(
    iteration = iteration,
    n_pairs = nrow(matched),
    spearman_rho = safe_spearman(matched$score_three_end, matched$score_full_length),
    median_match_distance = median(matched$match_distance, na.rm = TRUE),
    mean_abs_standardized_difference = mean(matched$mean_abs_standardized_difference, na.rm = TRUE),
    reused_candidate_fraction = mean(matched$reused_candidate, na.rm = TRUE)
  )
}

matched_null <- map_dfr(seq_len(n_null), match_one_null)
write_csv(matched_null, file.path(table_dir, "matched_covariate_null_spearman.csv"))

observed_rho <- observed_correlations %>%
  filter(.data$scope == "all_pairs") %>%
  pull(.data$spearman_rho)
empirical_p <- if (is.finite(observed_rho)) {
  (1 + sum(matched_null$spearman_rho >= observed_rho, na.rm = TRUE)) /
    (1 + sum(is.finite(matched_null$spearman_rho)))
} else {
  NA_real_
}

matched_null_summary <- tibble(
  observed_spearman_rho = observed_rho,
  n_null_requested = n_null,
  n_null_valid = sum(is.finite(matched_null$spearman_rho)),
  null_median = safe_quantile(matched_null$spearman_rho, 0.5),
  null_q025 = safe_quantile(matched_null$spearman_rho, 0.025),
  null_q975 = safe_quantile(matched_null$spearman_rho, 0.975),
  empirical_p_upper = empirical_p,
  median_match_distance = median(matched_null$median_match_distance, na.rm = TRUE),
  median_mean_abs_standardized_difference = median(matched_null$mean_abs_standardized_difference, na.rm = TRUE),
  median_reused_candidate_fraction = median(matched_null$reused_candidate_fraction, na.rm = TRUE)
)
write_csv(matched_null_summary, file.path(table_dir, "matched_covariate_null_summary.csv"))

celltype_palette <- setNames(
  rep(c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#6A3D9A", "#666666"), length.out = length(unique(shared_pairs$cell_type))),
  sort(unique(shared_pairs$cell_type))
)

theme_cell <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 7, colour = "black"),
      legend.text = element_text(size = 6.5),
      legend.title = element_blank(),
      plot.title = element_text(size = 9, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 7.5, colour = "#333333", hjust = 0),
      strip.background = element_rect(fill = "#F1F1F1", colour = NA),
      strip.text = element_text(size = 7, face = "bold"),
      panel.spacing = grid::unit(5, "pt")
    )
}

save_publication_plot <- function(plot, filename_stem, width, height) {
  ggsave(file.path(plot_dir, paste0(filename_stem, ".pdf")), plot, width = width, height = height, device = grDevices::cairo_pdf, bg = "white")
  ggsave(file.path(plot_dir, paste0(filename_stem, ".svg")), plot, width = width, height = height, device = grDevices::svg, bg = "white")
  ggsave(file.path(plot_dir, paste0(filename_stem, ".png")), plot, width = width, height = height, dpi = 600, bg = "white")
  ggsave(
    file.path(plot_dir, paste0(filename_stem, ".tiff")),
    plot,
    width = width,
    height = height,
    dpi = 600,
    compression = "lzw",
    bg = "white"
  )
}

independence_plot_df <- bind_rows(
  private_cohort_summary %>%
    transmute(
      category = "Patients",
      platform = .data$platform,
      value = .data$n_patients
    ),
  tibble(
    category = "Private SNVs",
    platform = c("three_end", "full_length"),
    value = c(length(private_snvs_three), length(private_snvs_full))
  ),
  tibble(
    category = "Private-SNV genes",
    platform = c("three_end", "full_length"),
    value = c(length(private_genes_three), length(private_genes_full))
  )
) %>%
  mutate(
    platform = recode(.data$platform, three_end = "3' RNA-seq", full_length = "Full-length RNA-seq"),
    category = factor(.data$category, levels = c("Patients", "Private SNVs", "Private-SNV genes"))
  )

p_independence <- ggplot(independence_plot_df, aes(x = .data$category, y = .data$value, fill = .data$platform)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64) +
  geom_text(aes(label = comma(.data$value)), position = position_dodge(width = 0.72), vjust = -0.25, size = 2.2) +
  scale_fill_manual(values = c("3' RNA-seq" = "#0072B2", "Full-length RNA-seq" = "#D55E00")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    x = NULL,
    y = "Count",
    title = "Independent cohorts and private-variant overlap",
    subtitle = sprintf(
      "Shared patients: %d; shared exact private SNVs: %d; SNV Jaccard = %.3g",
      length(patient_overlap),
      length(private_shared_snvs),
      private_exact_overlap$jaccard
    ),
    fill = NULL
  ) +
  theme_cell() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "top")

p_scatter <- ggplot(
  shared_pairs,
  aes(
    x = .data$cohort_rank_score_three_end,
    y = .data$cohort_rank_score_full_length,
    colour = .data$cell_type
  )
) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.45, linetype = "dashed", colour = "#888888") +
  geom_point(size = 1.7, alpha = 0.78) +
  scale_colour_manual(values = celltype_palette) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    x = "3' RNA-seq cohort rank score",
    y = "Full-length RNA-seq cohort rank score",
    title = "Independent private SNVs converge at gene-by-cell-type level",
    subtitle = sprintf(
      "Each point is one gene-by-cell-type pair; Spearman rho = %.2f (patient bootstrap 95%% CI %.2f to %.2f); n = %d pairs",
      observed_rho,
      correlation_summary$ci_low[correlation_summary$scope == "all_pairs"],
      correlation_summary$ci_high[correlation_summary$scope == "all_pairs"],
      nrow(shared_pairs)
    ),
    colour = NULL
  ) +
  theme_cell() +
  theme(legend.position = "right")

forest_df <- correlation_summary %>%
  filter(.data$scope == "within_cell_type", is.finite(.data$spearman_rho)) %>%
  arrange(.data$spearman_rho) %>%
  mutate(cell_type = factor(.data$cell_type, levels = .data$cell_type))

p_forest <- ggplot(forest_df, aes(y = .data$cell_type, x = .data$spearman_rho)) +
  geom_vline(xintercept = 0, linewidth = 0.45, linetype = "dashed", colour = "#888888") +
  geom_errorbarh(aes(xmin = .data$ci_low, xmax = .data$ci_high), height = 0, linewidth = 0.65, colour = "#4D4D4D") +
  geom_point(aes(colour = .data$cell_type), size = 2.2) +
  scale_colour_manual(values = celltype_palette, guide = "none") +
  coord_cartesian(xlim = c(-1, 1)) +
  labs(
    x = "Spearman rho (patient bootstrap 95% CI)",
    y = NULL,
    title = "Convergence within matched cell types",
    subtitle = "Patients, rather than cells or SNVs, are resampled"
  ) +
  theme_cell()

p_null <- ggplot(matched_null %>% filter(is.finite(.data$spearman_rho)), aes(x = .data$spearman_rho)) +
  geom_histogram(aes(y = after_stat(density)), bins = 45, fill = "#B8C6D1", colour = "white", linewidth = 0.2) +
  geom_density(linewidth = 0.75, colour = "#4D4D4D") +
  geom_vline(xintercept = observed_rho, colour = "#D55E00", linewidth = 0.9) +
  labs(
    x = "Spearman rho under matched-gene null",
    y = "Density",
    title = "Observed convergence exceeds a six-covariate matched null",
    subtitle = sprintf(
      "Matched for expression, gene length, callable bases, SNV count, coverage and carrier cells; empirical P = %.3g",
      empirical_p
    )
  ) +
  theme_cell()

save_publication_plot(p_independence, "01_cohort_independence_and_overlap", 7.01, 3.7)
save_publication_plot(p_scatter, "02_gene_celltype_rank_convergence", 7.01, 5.5)
save_publication_plot(p_forest, "03_patient_bootstrap_celltype_correlations", 5.6, max(3.8, 0.28 * nrow(forest_df) + 1.7))
save_publication_plot(p_null, "04_six_covariate_matched_null", 5.6, 4.0)

write_status(
  "COMPLETE",
  c(
    paste0("patients_three_end=", paste(patients_three, collapse = ";")),
    paste0("patients_full_length=", paste(patients_full, collapse = ";")),
    paste0("shared_patients=", length(patient_overlap)),
    paste0("private_snvs_three_end=", length(private_snvs_three)),
    paste0("private_snvs_full_length=", length(private_snvs_full)),
    paste0("shared_exact_private_snvs=", length(private_shared_snvs)),
    paste0("private_snv_jaccard=", format(private_exact_overlap$jaccard, scientific = TRUE, digits = 8)),
    paste0("common_private_genes=", length(common_private_genes)),
    paste0("shared_gene_celltype_pairs=", nrow(shared_pairs)),
    paste0("observed_spearman_rho=", format(observed_rho, digits = 8)),
    paste0("bootstrap_ci_low=", format(correlation_summary$ci_low[correlation_summary$scope == "all_pairs"], digits = 8)),
    paste0("bootstrap_ci_high=", format(correlation_summary$ci_high[correlation_summary$scope == "all_pairs"], digits = 8)),
    paste0("matched_null_empirical_p=", format(empirical_p, scientific = TRUE, digits = 8)),
    "patient_bootstrap_run=TRUE",
    "matched_null_run=TRUE"
  )
)

message("Analysis completed. Outputs: ", output_dir)
