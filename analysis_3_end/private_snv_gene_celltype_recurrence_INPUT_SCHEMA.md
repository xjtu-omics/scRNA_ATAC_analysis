# Patient-level input schema

`analyze_private_snv_gene_celltype_recurrence.R` requires one CSV per cohort. Each row must represent one patient-SNV-cell-type observation before cohort-level aggregation.

Required columns:

| Column | Meaning |
|---|---|
| `patient_id` | Globally unique patient identifier. Identifiers must not be reused between cohorts. |
| `cell_type` | Harmonized cell type or cell class. |
| `snv` | Exact SNV identifier, preferably `chr:position_REF>ALT`. |
| `gene` | Gene symbol hit by the SNV. Multiple symbols may be separated by `;`, `,`, or `|`. |
| `score` | Per-patient, per-SNV, per-cell-type PrismSNV Euclidean displacement score. |
| `gene_expression` | Mean expression of the gene in the corresponding patient and cell type. |
| `gene_length` | Gene length in base pairs. Use the same annotation release for both cohorts. |
| `callable_bases` | Number of callable bases for the gene in the corresponding patient and cell type. |
| `covered_cell_count` | Number of cells with REF or ALT coverage at the SNV. |
| `carrier_cell_count` | Number of ALT-supporting cells at the SNV. |

Optional column:

| Column | Meaning |
|---|---|
| `sample_id` | Sample or library identifier used for provenance. |

The script derives `n_snvs` after patient-level aggregation. It verifies that each retained SNV occurs in exactly one patient within each cohort, verifies that patient sets are disjoint, and reports exact SNV overlap between cohorts.

Example command:

```powershell
& "C:\path\to\Rscript.exe" `
  "analysis_3_end\analyze_private_snv_gene_celltype_recurrence.R" `
  "path\to\3-end_private_snv_scores_by_patient_celltype.csv" `
  "path\to\full-length_private_snv_scores_by_patient_celltype.csv" `
  "analysis_3_end\private_snv_gene_celltype_recurrence" `
  2000 2000 20260901 strict 1 10
```

The two input tables currently used by the original Figure 4 workflow do not contain `patient_id`, expression, callable-base, coverage, or carrier-cell fields. They are suitable only for the descriptive overlap audit and cannot support patient bootstrap or matched-null inference.
