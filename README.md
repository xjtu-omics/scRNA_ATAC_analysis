# scRNA and scATAC Benchmark Scripts

This repository contains standalone Python scripts for benchmarking single-cell RNA-seq and single-cell ATAC-seq datasets. The analyses compare multiple single-cell technologies, summarize quality metrics, and generate benchmark tables and publication-style figures.

## What This Repository Does

The scRNA workflow evaluates transcriptome capture and downstream biological signal across technologies, including gene/UMI complexity, marker detection, cell-type composition, lncRNA capture, reproducibility, and composite scoring.

The scATAC workflow evaluates chromatin accessibility data quality and benchmark performance, including fragment statistics, TSSe, FRiP, nucleosome signal, doublet scoring, filtering summaries, and final scoring reports.

## Repository Contents

- `benchmark_scRNA_metrics_compute.py`  
  Runs the scRNA-seq benchmark workflow and writes summary metrics and figures.

- `benchmark_scATAC_metrics_compute.py`  
  Runs the scATAC-seq QC and benchmark workflow with SnapATAC2.

## Requirements

The scripts expect a scientific Python environment with commonly used single-cell analysis packages:

- `numpy`
- `pandas`
- `matplotlib`
- `scipy`
- `scikit-learn`
- `anndata`
- `scanpy`
- `h5py`
- `snapatac2` for the scATAC workflow

No package metadata or environment file is currently included.

## Before Running

Both scripts contain hard-coded working directories and dataset paths, for example:

```python
os.chdir("/groups/adv2105_gp/yichen/Yi/multi/scRNA-out")
```

Review and update these paths before execution. Also check dataset-specific configuration values such as `DATASETS`, `FRAG`, `PEAKS`, `TENXH5`, `GENOME`, `OUTDIR`, and filtering thresholds.

## Running the Workflows

After confirming paths and dependencies, run:

```bash
python benchmark_scRNA_metrics_compute.py
python benchmark_scATAC_metrics_compute.py
```

These scripts can be long-running and may require large input matrices or fragment files. The scATAC script includes multiple notebook-derived workflow blocks, so inspect the active configuration before running the full file.

## Outputs

Typical outputs include summary tables, QC reports, benchmark plots, and publication-oriented PDF figures. The scRNA workflow writes to `benchmark_out/`. The scATAC workflow uses `qc_out*` directories and, in later blocks, absolute output paths under `/groups/...`.
