# scRNA and scATAC Benchmark Scripts

Benchmark scripts for comparing single-cell RNA-seq and single-cell ATAC-seq technologies. The Python scripts summarize quality metrics, capture efficiency, biological signal, and scoring results, while the R scripts compare RNA-seq platform consistency using SNV-effect and cell-perturbation evidence, then export tables and publication-style figures.

## 📁 Repository Contents

### Python scripts (`with_other_platform/`)

- `benchmark_scRNA_metrics_compute.py`  
  Runs the scRNA-seq benchmark workflow and writes summary metrics and figures. Evaluates transcriptome capture and downstream biological signal across technologies, including gene/UMI complexity, marker detection, cell-type composition, lncRNA capture, reproducibility, and composite scoring. Outputs are written under `benchmark_out/`.

- `benchmark_scATAC_metrics_compute.py`  
  Runs the scATAC-seq QC and benchmark workflow with SnapATAC2. Evaluates chromatin accessibility data quality and benchmark performance, including fragment statistics, TSSe, FRiP, nucleosome signal, doublet scoring, filtering summaries, and final scoring reports. Outputs are written under `qc_out*` directories.

### R scripts (`with_3_end/`)

- `compare_snv_effect_rna_platforms.R`  
  Compares SNV-effect profiles across RNA-seq platforms at the coarse cell-class level. Performs marker-based cell-class assignment, builds per-cell-class SNV score profiles, and computes cross-platform profile correlations. Outputs to `rna_platform_consistency/`.

- `compare_cell_perturbation_rna_platforms.R`  
  Compares cell-level perturbation scores across RNA-seq platforms. Performs marker-based cell-class assignment, builds perturbation score profiles, and evaluates rank agreement and distribution distances between platforms. Outputs to `cell_perturbation_rna_platform_consistency/`.

- `combine_rna_platform_consistency.R`  
  Combines the two orthogonal analyses (SNV-effect and cell-perturbation) into joint consistency profiles, rank-concordance summaries, and supplementary evidence tables and figures. Outputs to `combined_rna_platform_consistency/`.

## 🧰 Requirements

**Python** (scientific environment): `numpy`, `pandas`, `matplotlib`, `scipy`, `scikit-learn`, `anndata`, `scanpy`, `h5py`, and `snapatac2` for the scATAC workflow.

**R**: `ggplot2`, `dplyr`, `tidyr`, `readr`, `stringr`, `purrr`, `scales`, and `ragg` (for `compare_cell_perturbation_rna_platforms.R`). Each R script checks for missing packages and stops with an install message if any are absent.

## ⚙️ Before Running

The Python scripts contain hard-coded working directories and dataset paths, for example:

```python
os.chdir("/groups/adv2105_gp/yichen/Yi/multi/scRNA-out")
```

Review and update these paths before execution, and check dataset-specific configuration such as `DATASETS`, `FRAG`, `PEAKS`, `TENXH5`, `GENOME`, `OUTDIR`, and filtering thresholds.

The R scripts accept an input directory and an output directory as command-line arguments and fall back to built-in defaults otherwise. Confirm these paths before running.

## ▶️ Running the Workflows

Python:

```bash
python with_other_platform/benchmark_scRNA_metrics_compute.py
python with_other_platform/benchmark_scATAC_metrics_compute.py
```

R (arguments are optional; defaults are used when omitted):

```bash
Rscript with_3_end/compare_snv_effect_rna_platforms.R [input_dir] [output_dir]
Rscript with_3_end/compare_cell_perturbation_rna_platforms.R [input_dir] [output_dir]
Rscript with_3_end/combine_rna_platform_consistency.R [base_dir] [output_dir]
```

## 📊 Outputs

Typical outputs include summary tables (CSV), QC reports, benchmark plots, and publication-oriented PNG/PDF figures. The Python scRNA workflow writes to `benchmark_out/`, the scATAC workflow uses `qc_out*` directories, and the R scripts write per-analysis tables and figures into their respective output directories.
