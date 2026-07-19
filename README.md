# scRNA and scATAC Benchmark Scripts

Benchmark scripts for comparing single-cell RNA-seq and single-cell ATAC-seq technologies. The Python scripts summarize quality metrics, capture efficiency, biological signal, and composite scoring across sequencing platforms, while the R scripts assess RNA-seq platform consistency (3-end vs. full-length) using SNV-effect and cell-perturbation evidence. All scripts export summary tables and publication-style figures.

## 📁 Repository Layout

```
scRNA_ATAC_benchmark/
├── with_other_platform/          # Python cross-platform benchmark scripts
│   ├── benchmark_scRNA_metrics_compute.py
│   └── benchmark_scATAC_metrics_compute.py
└── with_3_end/                   # R 3-end vs. full-length consistency scripts
    ├── compare_snv_effect_rna_platforms.R
    ├── compare_cell_perturbation_rna_platforms.R
    └── combine_rna_platform_consistency.R
```

## 🐍 Python Scripts (`with_other_platform/`)

### `benchmark_scRNA_metrics_compute.py`

Notebook-derived scRNA-seq benchmark that compares many technologies (e.g. SPLiT-seq, 10x Multiome, SmartSeq2, SeekGene, Microwell-seq, Drop-seq, 10x v2/v3, 10x snRNA) with **SeekGene** highlighted as the reference.

- **Inputs**: per-technology matrices declared in the `DATASETS` list, supporting `h5ad` and 10x `h5` formats. Optional GTF is used to map Ensembl IDs to gene symbols for robust marker/MT detection. Counts are auto-detected across common layers (`counts`, `raw`, `umi`, …).
- **Metrics computed**:
  - Gene/UMI complexity and library saturation across UMI targets (`UMI_TARGETS = [500, 1000, 2000, 5000, 10000]`).
  - Marker-based cell-type labeling using a brain marker panel (`Neuron_EX/IN`, `Oligo`, `Astro`, `Micro`, `OPC`, `Endo`, `Peri`).
  - Biological separability via silhouette score (marker-first, Leiden fallback).
  - Cell-type composition, lncRNA capture (median lncRNA genes per cell), and reproducibility.
  - A normalized, weighted **composite score** ranking technologies overall.
- **Statistics/plots**: Kruskal–Wallis tests and significance bars against SeekGene; one figure per metric.
- **Config knobs**: `MAX_CELLS_PER_DS = 8000`, random seed `1234`, `OUTDIR = "benchmark_out"`.
- **Outputs**: summary tables and publication-oriented **PDF** figures under `benchmark_out/` (all `savefig` calls are forced to vector PDF, font type 42).

### `benchmark_scATAC_metrics_compute.py`

All-in-one scATAC-seq QC and benchmark built on **SnapATAC2**, robust to new/old API signatures and HDF5 file-lock issues.

- **Inputs** (top-of-file config): `FRAG` (fragments `.tsv.gz`), optional `PEAKS` (BED) and `TENXH5` (for barcode whitelist), and `GENOME` (`mm10 | mm39 | hg38 | hg19`).
- **Workflow**: fragment import with barcode whitelist → TSS enrichment (TSSe) → FRiP after peak sanitization → peak-matrix `n_peaks` (auto fallback to tile-matrix `n_tiles`) → nucleosome/fragment-length signal → doublet scoring (Scrublet on tile matrix).
- **Config knobs**: `DO_NUCLEOSOME = True`, `DOUBLETS_THRESHOLD = 0.5`, `SCRUBLET_JOBS = 4`, `OUTDIR = "qc_out"`.
- **Outputs**: fragment statistics, QC panels, filtering summaries, and final scoring reports written under `qc_out*` directories.

## 📊 R Scripts (`with_3_end/`)

These scripts compare two RNA-seq strategies — **3-end** and **full-length** — projected onto shared, marker-derived coarse cell classes. Each script takes an optional input directory and output directory as command-line arguments (or environment overrides) and otherwise uses built-in defaults.

### `compare_snv_effect_rna_platforms.R`

Compares **SNV-effect (perturbation) profiles** across the two platforms at the coarse cell-class level.

- **Inputs**: `3-end_snv_perturbation_scores_by_celltype.csv`, `full-length_snv_perturbation_scores_by_celltype.csv`, and matching `*_cell_cluster_marker_genes_top10.csv` marker tables.
- **Does**: marker-based coarse cell-class assignment, per-cell-class SNV score profiles, and cross-platform (Spearman) profile correlations.
- **Outputs** (`rna_platform_consistency/`): tables including `celltype_marker_assignments.csv`, `coarse_cellclass_snv_scores.csv`, `cellclass_score_profile[_wide].csv`, `cellclass_profile_correlation_summary.csv`, plus assignment, rank-concordance, and correlation figures.

### `compare_cell_perturbation_rna_platforms.R`

Compares **cell-level perturbation scores** (primary metric: `perturbation_score_euclidean`) across the two platforms.

- **Inputs**: `3-end_cell_perturbation_scores.csv`, `full-length_cell_perturbation_scores.csv`, and the top-10 marker gene tables.
- **Does**: marker-based cell-class assignment, perturbation score profiles, coverage/cluster-count summaries, rank agreement, and distribution-distance comparisons.
- **Outputs** (`cell_perturbation_rna_platform_consistency/`): tables including `cell_perturbation_scores_with_cellclass.csv`, `cell_perturbation_score_profile[_wide].csv`, `cell_perturbation_profile_correlation_summary.csv`, `cell_perturbation_distribution_distance_summary.csv`, `top{N}_cells_by_metric.csv`, plus assignment, rank, and correlation-heatmap figures. Requires the `ragg` package.

### `combine_rna_platform_consistency.R`

Integrates the two orthogonal analyses above into a joint consistency assessment.

- **Inputs**: the profile/correlation/distribution CSVs produced by the two comparison scripts.
- **Does**: builds joint cell-class consistency profiles, joint correlation and consistency-score summaries, rank-concordance across evidence layers, and a supplementary-evidence summary.
- **Outputs** (`combined_rna_platform_consistency/`): tables including `joint_cellclass_consistency_profile.csv`, `joint_profile_correlation_summary.csv`, `joint_consistency_score_summary.csv`, `joint_cellclass_rank_concordance.csv`, `joint_supplementary_evidence_summary.csv`, plus overview scatter, rank, and rank-heatmap figures.

## 🧰 Requirements

**Python** (scientific environment): `numpy`, `pandas`, `matplotlib`, `scipy`, `scikit-learn`, `anndata`, `scanpy`, `h5py`, and `snapatac2` (scATAC workflow).

**R**: `ggplot2`, `dplyr`, `tidyr`, `readr`, `stringr`, `purrr`, `scales`, and `ragg` (used by `compare_cell_perturbation_rna_platforms.R`). Each R script verifies required packages up front and stops with an install message if any are missing.

## ⚙️ Before Running

The Python scripts contain hard-coded working directories and dataset paths, for example:

```python
os.chdir("/groups/adv2105_gp/yichen/Yi/multi/scRNA-out")
```

Review and update these paths first, and check dataset-specific configuration such as `DATASETS`, `FRAG`, `PEAKS`, `TENXH5`, `GENOME`, `OUTDIR`, and filtering thresholds. Note the scripts are notebook-derived and execute top-level code on import.

The R scripts default their input/output directories to paths under `I:/极空间/...`; pass explicit arguments (or set the corresponding environment variables) to point at your own data.

## ▶️ Running the Workflows

Python:

```bash
python with_other_platform/benchmark_scRNA_metrics_compute.py
python with_other_platform/benchmark_scATAC_metrics_compute.py
```

R (arguments optional; defaults used when omitted):

```bash
Rscript with_3_end/compare_snv_effect_rna_platforms.R [input_dir] [output_dir]
Rscript with_3_end/compare_cell_perturbation_rna_platforms.R [input_dir] [output_dir]
Rscript with_3_end/combine_rna_platform_consistency.R [base_dir] [output_dir]
```

Run the two `compare_*` R scripts before `combine_rna_platform_consistency.R`, since the combine step consumes their output tables.

## 📈 Outputs

Typical outputs are summary tables (CSV), QC reports, and publication-oriented figures (vector PDF from Python; PNG + PDF from R). The Python scRNA workflow writes to `benchmark_out/`, the scATAC workflow to `qc_out*` directories, and each R script writes its tables and figures into its own consistency output directory.
