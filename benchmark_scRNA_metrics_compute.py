#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib, matplotlib.pyplot as plt
matplotlib.rcParams['savefig.format'] = 'pdf'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

__orig_savefig = plt.savefig
def __savefig_pdf(fname, *args, **kwargs):
    import os
    s = str(fname)
    if not s.lower().endswith('.pdf'):
        s = os.path.splitext(s)[0] + '.pdf'
    return __orig_savefig(s, *args, **kwargs)
plt.savefig = __savefig_pdf
print('[config] Matplotlib configured: all savefig() will write PDF files.')
try:
    import scanpy as sc
    sc.settings.file_format_figs = 'pdf'
    print("[config] scanpy: file_format_figs set to 'pdf'")
except Exception:
    pass


# In[2]:


import os
os.chdir("/groups/adv2105_gp/yichen/Yi/multi/scRNA-out")


# 
# # Single-cell Technology Benchmark — **Full + lncRNA & Composite Score**
# 
# This notebook compares multiple single-cell technologies and highlights **SeekGene**, with new additions:
# - **lncRNA capture ability**: number of lncRNA genes detected per cell (median);
# - **Composite score**: normalize multiple metrics and compute a weighted sum to see which technology is best overall.
#
# > Plotting: matplotlib only; one plot per cell; no specific colors set.
# 

# In[3]:


# ==== Setup & imports ====
# !pip install --user scanpy anndata scipy numpy pandas h5py scikit-learn

import os, math, warnings, random, itertools
from dataclasses import dataclass
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread
import anndata as ad
import scanpy as sc
from sklearn.metrics import silhouette_score
import matplotlib
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore", category=UserWarning)
np.random.seed(1234); random.seed(1234)

OUTDIR = "benchmark_out"
os.makedirs(OUTDIR, exist_ok=True)

print("scanpy:", sc.__version__)
print("anndata:", ad.__version__)
print("numpy:", np.__version__)
print("pandas:", pd.__version__)
print("matplotlib:", matplotlib.__version__)


# In[4]:


# === Gene symbol mapping from GTF (for marker/MT detection robustness) ===
def map_varnames_to_symbols_from_gtf(adata, gtf_path, inplace=True, uppercase=True):
    """
    If var_names are Ensembl IDs (or mixed), use the GTF to map them to gene symbols (gene_name).
    - inplace=True: modify adata.var_names directly; otherwise return a copy.
    - uppercase=True: convert to uppercase to match the marker list (usually uppercase).
    """
    import gzip, re
    mp = {}
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    try:
        with opener(gtf_path, "rt") as fh:
            for line in fh:
                if not line or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9 or parts[2] != "gene":
                    continue
                attrs = parts[8]
                gid_m = re.search(r'gene_id "([^"]+)"', attrs)
                gnm_m = re.search(r'gene_name "([^"]+)"', attrs)
                if gid_m and gnm_m:
                    mp[gid_m.group(1)] = gnm_m.group(1)
    except FileNotFoundError:
        print(f"[warn] GTF not found: {gtf_path}; skip mapping")
        return adata

    A = adata if inplace else adata.copy()
    idx = A.var_names.astype(str)
    is_ens = idx.str.startswith(("ENSG", "ENSMUSG", "ENS"))
    # Only attempt mapping for Ensembl-like entries; keep other names as their original values
    mapped = idx.where(~is_ens, idx.map(mp))
    mapped = mapped.fillna(idx)  # keep original name for those that cannot be mapped
    if uppercase:
        mapped = mapped.str.upper()
    A.var["gene_symbols"] = mapped
    A.var_names = mapped
    A.var_names_make_unique()
    return A


# In[5]:


# ==== Dataset registry (paths as provided) ====
@dataclass
class DS:
    name: str
    tech: str
    path: str
    loader: str   # "h5ad" | "10x_h5"
    prefer_counts: bool = True

DATASETS = [
    DS("SPLiT-seq",      "SPLiT-seq", "/groups/adv2105_gp/yichen/Yi/multi/split_seq/merged_splitseq.h5ad", "h5ad",  False),
    DS("10x_multiome", "10x Multiome", "/groups/adv2105_gp/yichen/Yi/multi/10x_multiome/M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_filtered_feature_bc_matrix.h5", "10x_h5", True),
    DS("SSv2_nonMy",     "SmartSeq2", "/groups/adv2105_gp/yichen/Yi/multi/smartseqv2_GSE132042_brain_non_myeloid/smartseqv2_GSE132042_brain_non_myeloid.h5ad", "h5ad", False),
    DS("SSv2_my",        "SmartSeq2", "/groups/adv2105_gp/yichen/Yi/multi/smartseqv2_GSE132042_brain_myeloid/smartseqv2_GSE132042_brain_myeloid.h5ad", "h5ad", False),
    DS("SeekGene",       "SeekGene",  "/groups/adv2105_gp/yichen/Yi/multi/seekgene/filtered_feature_bc_matrix.h5", "10x_h5", True),
    DS("Microwell-seq",  "Microwell", "/groups/adv2105_gp/yichen/Yi/multi/microwell-seq/GSE153562_RAW_extracted/brain/GSM_brain_DGE.merged_by_intersection.h5ad", "h5ad", False),
    DS("Drop-seq",       "Drop-seq",  "/groups/adv2105_gp/yichen/Yi/multi/drop_seq/GSE116470.Dropseq.brain.all_regions.raw.h5ad", "h5ad", True),
    DS("10x v3",         "10x v3",    "/groups/adv2105_gp/yichen/Yi/multi/10xv3/10xv3.h5ad", "h5ad", True),
    DS("10x v2",         "10x v2",    "/groups/adv2105_gp/yichen/Yi/multi/10xv2_scRNA/10xv2_scRNA.h5ad", "h5ad", True),
    DS("10x snRNA",      "10x snRNA", "/groups/adv2105_gp/yichen/Yi/multi/10x_snRNA/filtered_gene_bc_matrices/mm10/mm10.h5ad", "h5ad", True),
]

# Control the sample size (for the full set, set the upper limit very large)
MAX_CELLS_PER_DS = 8000
UMI_TARGETS = [500, 1000, 2000, 5000, 10000]
SEEKGENE_TECH = "SeekGene"


# In[6]:


# ==== Counts-aware loading & helpers ====
USE_RAW_IF_PRESENT = True
LAYER_CANDIDATES = ["counts", "raw", "umi", "UMI", "X_counts"]

def looks_like_counts(X, int_frac_thresh=0.95, q95_thresh=3):
    if sp.issparse(X):
        data = X.data
        if data.size == 0:
            return False
    else:
        data = np.asarray(X).ravel()
        if data.size == 0:
            return False
    frac_int = (np.abs(data - np.round(data)) < 1e-6).sum() / data.size
    q95 = np.percentile(data, 95)
    return (frac_int >= int_frac_thresh) and (q95 >= q95_thresh)

def load_any(ds):
    if ds.loader == "h5ad":
        a = sc.read_h5ad(ds.path)
        counts_src = "X"
        if not looks_like_counts(a.X):
            if USE_RAW_IF_PRESENT and getattr(a, "raw", None) is not None:
                try:
                    raw = a.raw.to_adata()
                    if looks_like_counts(raw.X):
                        a = raw
                        counts_src = "raw.X"
                except Exception:
                    pass
        if counts_src == "X" and not looks_like_counts(a.X):
            for key in LAYER_CANDIDATES:
                if key in getattr(a, "layers", {}):
                    Xcand = a.layers[key]
                    Xcand = Xcand.tocsr() if sp.issparse(Xcand) else sp.csr_matrix(Xcand)
                    if looks_like_counts(Xcand):
                        a.X = Xcand
                        counts_src = f"layers['{key}']"
                        break
    elif ds.loader == "10x_h5":
        a = sc.read_10x_h5(ds.path)
        # If this is 10x Multiome, keep only Gene Expression features
        if "feature_types" in a.var.columns:
            ft = a.var["feature_types"].astype(str)
            if ft.isin(["Gene Expression", "GeneExpression", "gex"]).any():
                a = a[:, ft.isin(["Gene Expression", "GeneExpression", "gex"])].copy()
                counts_src = "10x_h5:GeneExpression"
            else:
                counts_src = "10x_h5"
        else:
            counts_src = "10x_h5"
    else:
        raise ValueError(ds.loader)

    a.X = a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)
    a.var_names = a.var_names.astype(str)
    a.var_names_make_unique()
    a.obs["tech"] = ds.tech
    a.obs["dataset"] = ds.name
    # Automatic gene name mapping (when Ensembl IDs are detected)
    try:
        if 'GTF_PATH' in globals():
            idx = a.var_names.astype(str)
            if idx.str.startswith(('ENSG','ENSMUSG','ENS')).any():
                a = map_varnames_to_symbols_from_gtf(a, GTF_PATH, inplace=True, uppercase=True)
                a.uns['gene_symbol_mapped'] = True
    except Exception as _e:
        print('[warn] gene symbol mapping skipped:', _e)

    a.uns["counts_source"] = counts_src
    return a

def is_counts_matrix(adata):
    return looks_like_counts(adata.X)

def nnz_per_row(X):
    return np.diff(X.indptr) if sp.issparse(X) else (X>0).sum(axis=1).A1

def libsize_per_row(X):
    return (np.asarray(X.sum(axis=1)).ravel() if sp.issparse(X) else X.sum(axis=1))

def norm_log1p_cpm(adata):
    a = adata.copy()
    if is_counts_matrix(a):
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
    return a


# In[7]:


# ==== Marker sets & robust label ====
MARKERS = {
    "Neuron_EX": ["SLC17A7","VGLUT1","SLC17A6","TBR1","SATB2","RBFOX3","MAP2","SNAP25","TUBB3"],
    "Neuron_IN": ["GAD1","GAD2","DLX1","DLX2","RELN","SLC6A1","PVALB","SST","VIP"],
    "Oligo":     ["MBP","MOG","PLP1","MOBP","MAL","MAG","CLDN11"],
    "Astro":     ["ALDH1L1","SLC1A3","AQP4","GLUL","GJA1","GFAP"],
    "Micro":     ["CX3CR1","P2RY12","AIF1","TYROBP","CSF1R","SPI1"],
    "OPC":       ["PDGFRA","CSPG4","SOX10","OLIG1","OLIG2"],
    "Endo":      ["KDR","PECAM1","CLDN5","FLT1","RGS5"],
    "Peri":      ["PDGFRB","RGS5","KCNJ8","ABCC9"],
}

def marker_label(adata_norm, markers, min_score=0.2):
    genes_upper = np.asarray([g.upper() for g in adata_norm.var_names])
    M = {}
    for ct, lst in markers.items():
        idx_list = []
        for g in lst:
            hit = np.where(genes_upper == g.upper())[0]
            if hit.size > 0:
                idx_list.extend(hit.tolist())
        if len(idx_list) == 0:
            M[ct] = np.zeros(adata_norm.n_obs)
        else:
            Xsub = adata_norm[:, idx_list].X
            score = (np.asarray(Xsub.mean(axis=1)).ravel() if sp.issparse(Xsub) else Xsub.mean(axis=1))
            M[ct] = score
    Mdf = pd.DataFrame(M, index=adata_norm.obs_names)
    max_ct = Mdf.idxmax(axis=1)
    max_val = Mdf.max(axis=1)
    label = max_ct.where(max_val >= min_score, other="Unknown")
    return label, Mdf

def compute_silhouette(adata_norm, labels, n_pcs=50, min_cells=100, min_classes=2):
    return silhouette_marker_or_leiden(adata_norm, labels, n_pcs=n_pcs, min_cells=min_cells, min_classes=min_classes)

    idx = labels[labels!="Unknown"].index
    if len(idx) < 100 or len(set(labels.loc[idx])) < 2:
        return np.nan
    sc.pp.pca(adata_norm, n_comps=min(n_pcs, adata_norm.n_vars-1))
    Xp = adata_norm.obsm["X_pca"][adata_norm.obs_names.get_indexer(idx), :]
    try:
        return silhouette_score(Xp, labels.loc[idx].values, metric="euclidean")
    except Exception:
        return np.nan


# In[8]:


# === Silhouette with marker-first, Leiden-fallback ===
from sklearn.metrics import silhouette_score

def silhouette_marker_or_leiden(
    adata_norm,
    labels_marker,
    n_pcs=50,
    min_cells=100,       # minimum number of valid cells
    min_classes=2        # minimum number of classes
):
    # First try marker labels
    idx = labels_marker[labels_marker != "Unknown"].index
    if len(idx) >= min_cells and labels_marker.loc[idx].nunique() >= min_classes:
        sc.pp.pca(adata_norm, n_comps=min(n_pcs, max(2, adata_norm.n_vars-1)))
        Xp = adata_norm.obsm["X_pca"][adata_norm.obs_names.get_indexer(idx), :]
        try:
            return float(silhouette_score(Xp, labels_marker.loc[idx].values, metric="euclidean"))
        except Exception:
            pass

    # Fallback: Leiden unsupervised clustering
    sc.pp.pca(adata_norm, n_comps=min(n_pcs, max(2, adata_norm.n_vars-1)))
    sc.pp.neighbors(adata_norm)
    sc.tl.leiden(adata_norm, key_added="leiden_fallback", resolution=1.0)
    lab2 = adata_norm.obs["leiden_fallback"]
    if lab2.nunique() < min_classes or adata_norm.n_obs < min_cells:
        return np.nan
    Xp = adata_norm.obsm["X_pca"]
    try:
        return float(silhouette_score(Xp, lab2.values, metric="euclidean"))
    except Exception:
        return np.nan


# In[9]:


# ==== Main evaluation loop ====
from collections import defaultdict

rows_metrics = []
comp_rows = []
rng = np.random.RandomState(2024)

for ds in DATASETS:
    if not os.path.exists(ds.path):
        print(f"[WARN] Missing: {ds.path} -> skip")
        continue
    print(f"\n[LOAD] {ds.name} | {ds.tech} | {ds.path}")
    a = load_any(ds)

    if a.n_obs > MAX_CELLS_PER_DS:
        idx = np.random.choice(a.n_obs, size=MAX_CELLS_PER_DS, replace=False)
        a = a[idx, :].copy()
        a.obs_names_make_unique()

    n_genes_cell = nnz_per_row(a.X)
    zeros_cell = 1.0 - (n_genes_cell / a.n_vars)
    lib = libsize_per_row(a.X)
    mt_mask = np.array([str(g).upper().startswith("MT-") for g in a.var_names])
    pct_mt = (np.asarray(a.X[:, mt_mask].sum(axis=1)).ravel() / np.maximum(lib, 1)) * 100.0 if mt_mask.any() else np.zeros(a.n_obs)

    nnz_gene = (np.asarray(a.X.astype(bool).sum(axis=0)).ravel() if sp.issparse(a.X) else (a.X>0).sum(axis=0).A1)
    zeros_gene = 1.0 - (nnz_gene / a.n_obs)

    median_genes = float(np.median(n_genes_cell))
    median_umi = float(np.median(lib))
    mean_zeros_cell = float(np.mean(zeros_cell))
    mean_zeros_gene = float(np.mean(zeros_gene))
    mt_median = float(np.median(pct_mt))

    a_norm = norm_log1p_cpm(a)
    reps = []
    for _ in range(5):
        perm = np.random.permutation(a_norm.n_obs)
        half = a_norm.n_obs // 2
        idxA, idxB = perm[:half], perm[half:]
        gmeanA = np.asarray(a_norm[idxA, :].X.mean(axis=0)).ravel()
        gmeanB = np.asarray(a_norm[idxB, :].X.mean(axis=0)).ravel()
        ok = (gmeanA > 0) | (gmeanB > 0)
        reps.append(np.corrcoef(gmeanA[ok], gmeanB[ok])[0,1] if ok.sum()>=50 else np.nan)
    rep_mean = float(np.nanmean(reps)); rep_sd = float(np.nanstd(reps))

    labels, _ = marker_label(a_norm, MARKERS, min_score=0.2)
    sil = compute_silhouette(a_norm, labels)

    do_complexity = is_counts_matrix(a)
    if do_complexity:
        for t in UMI_TARGETS:
            rows_, cols_, data_ = [], [], []
            sum_per_cell = np.asarray(a.X.sum(axis=1)).ravel()
            for i in range(a.n_obs):
                start, end = a.X.indptr[i], a.X.indptr[i+1]
                if start==end: 
                    continue
                idx = a.X.indices[start:end]
                val = a.X.data[start:end].astype(int)
                L = sum_per_cell[i]
                p = 1.0 if L <= 0 else min(1.0, t / L)
                thin = np.random.binomial(val, p)
                keep = thin > 0
                if keep.any():
                    rows_.extend([i]*keep.sum()); cols_.extend(idx[keep]); data_.extend(thin[keep])
            Xthin = sp.csr_matrix((np.array(data_), (np.array(rows_), np.array(cols_))), shape=a.X.shape)
            gdet = nnz_per_row(Xthin)
            comp_rows.append(dict(dataset=ds.name, tech=ds.tech, umi=t, median_genes=float(np.median(gdet))))
    else:
        print("  [INFO] Skip complexity curve (no counts available)")

    rows_metrics.append(dict(
        dataset=ds.name, tech=ds.tech, n_cells=a.n_obs, n_genes=a.n_vars,
        median_genes_per_cell=median_genes, median_umi_per_cell=median_umi,
        mean_zeros_cell=mean_zeros_cell, mean_zeros_gene=mean_zeros_gene,
        mt_pct_median=mt_median, reproducibility_r_mean=rep_mean, reproducibility_r_sd=rep_sd,
        silhouette_marker=(np.nan if sil is None else sil),
        counts_available=do_complexity,
        counts_source=a.uns.get("counts_source", "X")
    ))

df = pd.DataFrame(rows_metrics)
complexity_long = pd.DataFrame(comp_rows)
df.to_csv(os.path.join(OUTDIR, "summary_metrics.csv"), index=False)
if not complexity_long.empty:
    complexity_long.to_csv(os.path.join(OUTDIR, "complexity_long.csv"), index=False)

display(df)


# In[10]:


# ==== Stats helpers ====
from scipy.stats import kruskal, mannwhitneyu

def p_adjust_bh(pvals):
    p = np.asarray(pvals, dtype=float)
    n = p.size
    order = np.argsort(p)
    p_sorted = p[order]
    ranks = np.arange(1, n+1)
    adj_sorted = p_sorted * n / ranks
    adj_sorted = np.minimum.accumulate(adj_sorted[::-1])[::-1]
    adj = np.empty_like(adj_sorted)
    adj[order] = adj_sorted
    return np.clip(adj, 0, 1)

def stats_summary(df_in, metric):
    d = df_in[['tech', metric]].dropna()
    grouped = {t: g[metric].values for t, g in d.groupby('tech')}
    techs = list(grouped.keys())
    groups = [grouped[t] for t in techs]
    H, p = (kruskal(*groups) if len(groups) >= 2 else (np.nan, np.nan))
    rows = []
    for a, b in itertools.combinations(techs, 2):
        x, y = grouped[a], grouped[b]
        if len(x) >= 1 and len(y) >= 1:
            U, p_raw = mannwhitneyu(x, y, alternative='two-sided')
            rows.append(dict(A=a, B=b, p_raw=p_raw,
                             A_median=np.median(x), B_median=np.median(y)))
    res = pd.DataFrame(rows)
    if not res.empty:
        res['p_fdr'] = p_adjust_bh(res['p_raw'].values)
    return H, p, res

def star(p):
    if np.isnan(p): return 'ns'
    if p < 1e-4: return '****'
    if p < 1e-3: return '***'
    if p < 1e-2: return '**'
    if p < 0.05: return '*'
    return 'ns'


# In[11]:


# ==== Plot 1: Sensitivity (median genes/cell) with sig vs SeekGene ====
def bar_with_sig_vs_seekgene(df_in, metric, title, out_prefix):
    d = df_in[['tech', metric]].dropna()
    stat_med = d.groupby('tech')[metric].median().sort_values(ascending=False)
    techs = stat_med.index.tolist()
    vals = stat_med.values.tolist()
    fig, ax = plt.subplots(figsize=(9,4))
    ax.bar(range(len(techs)), vals)
    ax.set_xticks(range(len(techs)))
    ax.set_xticklabels(techs, rotation=45, ha='right')
    ax.set_ylabel(metric)
    ax.set_title(title)
    if 'SeekGene' in techs:
        H, p_global, pair = stats_summary(df_in, metric)
        if pair is not None and not pair.empty:
            sg = 'SeekGene'
            y_max = max(vals) if len(vals)>0 else 1.0
            y = y_max * 1.05
            step = (y_max * 0.04) if y_max>0 else 0.1
            rows = pair[(pair['A']==sg) | (pair['B']==sg)]
            for _, r in rows.iterrows():
                other = r['B'] if r['A']==sg else r['A']
                if other not in techs:
                    continue
                i_sg = techs.index(sg); i_oth = techs.index(other)
                x1, x2 = min(i_sg, i_oth), max(i_sg, i_oth)
                ax.plot([x1, x1, x2, x2], [y, y+step*0.2, y+step*0.2, y], linewidth=1)
                ax.text((x1+x2)/2, y+step*0.25, star(r['p_fdr']), ha='center', va='bottom')
                y += step
            if not np.isnan(p_global):
                ax.text(0.01, 0.95, f"Kruskal p={p_global:.2g}", transform=ax.transAxes, va='top')
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, f"{out_prefix}.png"), dpi=150)
    plt.show()

bar_with_sig_vs_seekgene(df, "median_genes_per_cell", "Sensitivity (median genes/cell)", "median_genes_with_sig")


# In[12]:


# ==== Plot 2: Reproducibility ====
fig, ax = plt.subplots(figsize=(9,4))
techs = df["tech"].unique().tolist()
data_box = [df[df.tech==t]["reproducibility_r_mean"].dropna().values for t in techs]
ax.boxplot(data_box, labels=techs, showfliers=False)
ax.set_xticklabels(techs, rotation=45, ha='right')
ax.set_ylabel("Split-half gene-mean Pearson r")
ax.set_title("Reproducibility (higher is better)")
fig.tight_layout()
fig.savefig(os.path.join(OUTDIR, "reproducibility_box.png"), dpi=150)
plt.show()

H,p,pairs = stats_summary(df, "reproducibility_r_mean")
print("Kruskal (reproducibility):", H, p)
display(pairs.head(10))


# In[13]:


# ==== Plot 3: Silhouette ====
bar_with_sig_vs_seekgene(df, "silhouette_marker", "Biological separability (silhouette)", "silhouette_with_sig")
H,p,pairs = stats_summary(df, "silhouette_marker")
print("Kruskal (silhouette):", H, p)
display(pairs.head(10))


# In[14]:


# ==== Plot 4: Zeros per cell ====
fig, ax = plt.subplots(figsize=(9,4))
techs = df["tech"].unique().tolist()
data_box = [df[df.tech==t]["mean_zeros_cell"].dropna().values for t in techs]
ax.violinplot(data_box, showmeans=True, showmedians=True)
ax.set_xticks(range(1, len(techs)+1))
ax.set_xticklabels(techs, rotation=45, ha='right')
ax.set_ylabel("Zeros per cell (fraction)")
ax.set_title("Sparsity (lower is better)")
fig.tight_layout()
fig.savefig(os.path.join(OUTDIR, "zeros_per_cell_violin.png"), dpi=150)
plt.show()

H,p,pairs = stats_summary(df, "mean_zeros_cell")
print("Kruskal (zeros/cell):", H, p)
display(pairs.head(10))


# In[15]:


# ==== Plot 5: Complexity curves (tech-level median ± SEM) ====
if 'complexity_long' in globals() and not complexity_long.empty:
    g = (complexity_long
         .groupby(["tech","umi"])["median_genes"]
         .agg(["median","mean","std","count"])
         .reset_index())
    g["sem"] = g["std"] / np.sqrt(g["count"].clip(lower=1))

    plt.figure(figsize=(7,5))
    for tech in sorted(g["tech"].unique()):
        sub = g[g.tech==tech].sort_values("umi")
        plt.errorbar(sub["umi"].values, sub["median"].values, yerr=sub["sem"].values,
                     marker='o', label=tech)
    plt.xlabel("Target UMIs per cell (downsampled)")
    plt.ylabel("Median genes detected")
    plt.title("Complexity curves (median ± SEM)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, "complexity_curves_tech_agg.png"), dpi=150)
    plt.show()
else:
    print("No complexity_long available; counts datasets may be missing.")


# In[16]:


# ==== Plot 6: Cell-type composition ====
def compute_composition(DATASETS, MARKERS, max_cells=6000):
    counts = []
    for ds in DATASETS:
        if not os.path.exists(ds.path):
            continue
        a = load_any(ds)
        if a.n_obs > max_cells:
            idx = np.random.choice(a.n_obs, size=max_cells, replace=False)
            a = a[idx,:].copy()
        labels, _ = marker_label(norm_log1p_cpm(a), MARKERS, min_score=0.2)
        tab = labels.value_counts().to_dict()
        tab['tech'] = ds.tech
        counts.append(tab)
    comp = pd.DataFrame(counts).fillna(0).groupby('tech').sum()
    comp = comp.div(comp.sum(axis=1), axis=0)
    return comp

comp = compute_composition(DATASETS, MARKERS)
comp.to_csv(os.path.join(OUTDIR, "composition_by_tech.csv"))

fig, ax = plt.subplots(figsize=(10,4))
bottom = np.zeros(comp.shape[0])
x = np.arange(comp.shape[0])
for ct in comp.columns:
    ax.bar(x, comp[ct].values, bottom=bottom, label=ct)
    bottom += comp[ct].values
ax.set_xticks(x)
ax.set_xticklabels(comp.index, rotation=45, ha='right')
ax.set_ylabel("Proportion")
ax.set_title("Cell-type composition (marker-based)")
ax.legend(ncol=4, fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(OUTDIR, "celltype_composition.png"), dpi=150)
plt.show()


# In[17]:


# ==== Plot 7: Tech-level gene-mean correlation heatmap ====
def tech_gene_means(DATASETS):
    gm = {}
    for ds in DATASETS:
        if not os.path.exists(ds.path): 
            continue
        a = load_any(ds)
        a_norm = norm_log1p_cpm(a)
        gmean = np.asarray(a_norm.X.mean(axis=0)).ravel()
        gm.setdefault(ds.tech, []).append(pd.Series(gmean, index=a_norm.var_names))
    if not gm:
        return pd.DataFrame()
    inter = None
    for lst in gm.values():
        for s in lst:
            inter = set(s.index) if inter is None else (inter & set(s.index))
    if not inter:
        return pd.DataFrame()
    genes = sorted(list(inter))
    tech_mean = {}
    for tech, lst in gm.items():
        aligned = [s.reindex(genes) for s in lst]
        tech_mean[tech] = pd.concat(aligned, axis=1).mean(axis=1)
    return pd.DataFrame(tech_mean)

G = tech_gene_means(DATASETS)
if not G.empty:
    C = G.corr(method='pearson')
    fig, ax = plt.subplots(figsize=(6,5))
    im = ax.imshow(C.values, aspect='auto')
    ax.set_xticks(range(C.shape[1])); ax.set_xticklabels(C.columns, rotation=45, ha='right')
    ax.set_yticks(range(C.shape[0])); ax.set_yticklabels(C.index)
    ax.set_title("Tech-level gene-mean correlation")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "tech_correlation_heatmap.png"), dpi=150)
    plt.show()
else:
    print("Correlation heatmap skipped (no overlapping genes).")


# In[18]:


# ==== Plot 8: Distribution of per-dataset median genes per cell ====
fig, ax = plt.subplots(figsize=(9,4))
order = df.groupby("tech")["median_genes_per_cell"].median().sort_values(ascending=False).index.tolist()
for tech in order:
    vals = df[df.tech==tech]["median_genes_per_cell"].dropna().values
    if vals.size == 0:
        continue
    ax.hist(vals, bins=10, histtype='step', label=tech)
ax.set_xlabel("Median genes per cell (per-dataset)")
ax.set_ylabel("Number of datasets")
ax.set_title("Distribution of per-dataset median genes per cell")
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(OUTDIR, "genes_per_cell_hist.png"), dpi=150)
plt.show()


# In[19]:


# ==== Plot 9: Genes vs UMIs scatter (per dataset) ====
def genes_umis_scatter(ds, max_cells=5000):
    if not os.path.exists(ds.path): 
        return
    a = load_any(ds)
    if a.n_obs > max_cells:
        idx = np.random.choice(a.n_obs, size=max_cells, replace=False)
        a = a[idx,:].copy()
    g = nnz_per_row(a.X)
    u = libsize_per_row(a.X)
    fig, ax = plt.subplots(figsize=(5,4))
    ax.scatter(u, g, s=3, alpha=0.5)
    ax.set_xlabel("UMIs per cell")
    ax.set_ylabel("Genes detected per cell")
    ax.set_title(f"{ds.tech}: genes vs UMIs (sampled)")
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, f"scatter_{ds.name}.png"), dpi=150)
    plt.show()

for ds in DATASETS:
    genes_umis_scatter(ds, max_cells=5000)


# In[20]:


# ==== Plot 10: mt% and per-gene zeros ====
fig, axes = plt.subplots(1,2, figsize=(12,4))

order_mt = df.groupby("tech")["mt_pct_median"].median().sort_values().index.tolist()
vals_mt = [df[df.tech==t]["mt_pct_median"].median() for t in order_mt]
axes[0].bar(range(len(order_mt)), vals_mt)
axes[0].set_xticks(range(len(order_mt)))
axes[0].set_xticklabels(order_mt, rotation=45, ha='right')
axes[0].set_ylabel("Median mt% (per-dataset)")
axes[0].set_title("Mitochondrial content")

order_zg = df.groupby("tech")["mean_zeros_gene"].median().sort_values().index.tolist()
vals_zg = [df[df.tech==t]["mean_zeros_gene"].median() for t in order_zg]
axes[1].bar(range(len(order_zg)), vals_zg)
axes[1].set_xticks(range(len(order_zg)))
axes[1].set_xticklabels(order_zg, rotation=45, ha='right')
axes[1].set_ylabel("Zeros per gene (fraction)")
axes[1].set_title("Sparsity per gene")

plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "mt_and_zerosgene.png"), dpi=150)
plt.show()


# In[21]:


# ==== Plot 11: Marker-panel detection heatmap ====
def marker_panel_detection(DATASETS, MARKERS, max_cells=6000):
    pertech = {}
    for ds in DATASETS:
        if not os.path.exists(ds.path):
            continue
        a = load_any(ds)
        if a.n_obs > max_cells:
            idx = np.random.choice(a.n_obs, size=max_cells, replace=False)
            a = a[idx,:].copy()
        X = a.X
        genes_upper = np.asarray([g.upper() for g in a.var_names])
        panel_rates = {}
        for ct, lst in MARKERS.items():
            gene_rates = []
            for g in lst:
                hit = np.where(genes_upper == g.upper())[0]
                if hit.size == 0: 
                    continue
                Xsub = X[:, hit]
                detected = (np.asarray((Xsub.sum(axis=1) > 0)).ravel() if sp.issparse(Xsub) else (Xsub.sum(axis=1) > 0))
                gene_rates.append(float(detected.mean()))
            panel_rates[ct] = (float(np.mean(gene_rates)) if gene_rates else np.nan)
        pertech.setdefault(ds.tech, []).append(panel_rates)
    rows = {}
    for tech, lst in pertech.items():
        df_ = pd.DataFrame(lst)
        rows[tech] = df_.mean(axis=0)
    return pd.DataFrame(rows).T

panel = marker_panel_detection(DATASETS, MARKERS)
panel.to_csv(os.path.join(OUTDIR, "marker_panel_detection.csv"))
if not panel.empty:
    fig, ax = plt.subplots(figsize=(8,4))
    im = ax.imshow(panel.values, aspect='auto')
    ax.set_xticks(range(panel.shape[1])); ax.set_xticklabels(panel.columns, rotation=45, ha='right')
    ax.set_yticks(range(panel.shape[0])); ax.set_yticklabels(panel.index)
    ax.set_title("Marker-panel detection fraction")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "marker_panel_detection_heatmap.png"), dpi=150)
    plt.show()
else:
    print("Marker-panel detection heatmap skipped (no data).")


# In[22]:


# ==== Plot 12: HVG counts ====
rows = []
for ds in DATASETS:
    if not os.path.exists(ds.path): 
        continue
    a = load_any(ds)
    a_norm = norm_log1p_cpm(a)
    try:
        sc.pp.highly_variable_genes(a_norm, n_top_genes=2000, flavor="seurat_v3", subset=False)
        hvg_count = int(a_norm.var['highly_variable'].sum())
    except Exception:
        hvg_count = np.nan
    rows.append(dict(dataset=ds.name, tech=ds.tech, hvg_count=hvg_count))
df_hvg = pd.DataFrame(rows)
df_hvg.to_csv(os.path.join(OUTDIR, "hvg_counts.csv"), index=False)

fig, ax = plt.subplots(figsize=(9,4))
order = df_hvg.groupby("tech")["hvg_count"].median().sort_values(ascending=False).index.tolist()
vals = [df_hvg[df_hvg.tech==t]["hvg_count"].median() for t in order]
ax.bar(range(len(order)), vals)
ax.set_xticks(range(len(order))); ax.set_xticklabels(order, rotation=45, ha='right')
ax.set_ylabel("HVG count (median per tech)")
ax.set_title("Highly variable genes (n_top_genes=2000)")
fig.tight_layout()
fig.savefig(os.path.join(OUTDIR, "hvg_counts_by_tech.png"), dpi=150)
plt.show()


# In[23]:


# ==== Plot 13: Average cells per sample ====
tech_sample = (
    df.groupby("tech")
      .agg(n_datasets=("dataset","nunique"),
           total_cells=("n_cells","sum"),
           avg_cells_per_dataset=("n_cells","mean"),
           median_cells_per_dataset=("n_cells","median"),
           sd_cells_per_dataset=("n_cells","std"))
      .reset_index()
      .sort_values("avg_cells_per_dataset", ascending=False)
)
display(tech_sample)

plt.figure(figsize=(8,4))
xt = tech_sample["tech"].tolist()
yv = tech_sample["avg_cells_per_dataset"].tolist()
plt.bar(range(len(xt)), yv)
plt.xticks(range(len(xt)), xt, rotation=45, ha='right')
plt.ylabel("Average cells per sample")
plt.title("Average cells per sample by technology")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "avg_cells_per_sample_by_tech.png"), dpi=150)
plt.show()


# In[24]:


# ==== Plot 14: Sample-level aggregation bars + significance vs SeekGene ====
def bar_sample_mean_with_sig(df_in, metric, title, out_prefix):
    stat_mean = df_in.groupby("tech")[metric].mean().dropna().sort_values(ascending=False)
    techs = stat_mean.index.tolist()
    vals = stat_mean.values.tolist()
    plt.figure(figsize=(9,4))
    plt.bar(range(len(techs)), vals)
    plt.xticks(range(len(techs)), techs, rotation=45, ha='right')
    plt.ylabel(metric)
    plt.title(title)
    if "SeekGene" in techs:
        H, p_global, pair = stats_summary(df_in, metric)
        rows = pair[(pair['A']=="SeekGene") | (pair['B']=="SeekGene")] if pair is not None and not pair.empty else pd.DataFrame()
        if not rows.empty:
            y_max = max(vals) if len(vals)>0 else 1.0
            y = y_max * 1.05
            step = (y_max * 0.04) if y_max>0 else 0.1
            for _, r in rows.iterrows():
                other = r['B'] if r['A']=="SeekGene" else r['A']
                if other not in techs: 
                    continue
                i_sg = techs.index("SeekGene")
                i_oth = techs.index(other)
                x1, x2 = min(i_sg, i_oth), max(i_sg, i_oth)
                plt.plot([x1, x1, x2, x2], [y, y+step*0.2, y+step*0.2, y], linewidth=1)
                plt.text((x1+x2)/2, y+step*0.25, star(r["p_fdr"]), ha='center', va='bottom')
                y += step
            if not np.isnan(p_global):
                plt.text(0.01, 0.95, f"Kruskal p={p_global:.2g}", transform=plt.gca().transAxes, va='top')
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, f"{out_prefix}.png"), dpi=150)
    plt.show()

bar_sample_mean_with_sig(df, "median_genes_per_cell", "Sensitivity (sample-mean of per-dataset medians)", "samplemean_genes_per_cell")
bar_sample_mean_with_sig(df, "reproducibility_r_mean", "Reproducibility (sample-mean)", "samplemean_reproducibility")
bar_sample_mean_with_sig(df, "silhouette_marker", "Biological separability (sample-mean)", "samplemean_silhouette")


# 
# ## New Figure A — lncRNA capture count
# Prefer the biotype column included in the matrix; if absent, set `GTF_PATH` to read GTF annotation and identify lncRNA.

# In[25]:


# ---- A1. Optional GTF annotation ----
GTF_PATH = None  # e.g.: "/path/to/gencode.vM31.annotation.gtf" (mouse)

LNC_KEYWORDS = [
    "lncrna", "lincRNA".lower(), "antisense", "processed_transcript",
    "sense_intronic", "sense_overlapping", "macro_lncRNA".lower(),
    "3prime_overlapping_ncRNA".lower(), "bidirectional_promoter_lncRNA".lower(),
    "non_coding", "noncoding", "long_noncoding", "lnc"
]

def load_lnc_from_gtf(gtf_path):
    import gzip, re
    if gtf_path is None or not os.path.exists(gtf_path):
        return set(), set()
    opener = gzip.open if gtf_path.endswith(".gz") else open
    gene_ids_lnc, gene_names_lnc = set(), set()
    with opener(gtf_path, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attrs = parts[8]
            def get_attr(key):
                m = re.search(rf'{key} "(.*?)"', attrs)
                return m.group(1) if m else None
            gid = get_attr("gene_id")
            gname = get_attr("gene_name")
            gtype = get_attr("gene_biotype") or get_attr("gene_type") or ""
            gt_low = (gtype or "").lower()
            is_lnc = any(k in gt_low for k in LNC_KEYWORDS)
            if is_lnc:
                if gid: gene_ids_lnc.add(gid)
                if gname: gene_names_lnc.add(gname)
    return gene_ids_lnc, gene_names_lnc

GTF_LNC_IDS, GTF_LNC_NAMES = load_lnc_from_gtf(GTF_PATH)
print("GTF lncRNA ids:", len(GTF_LNC_IDS), "names:", len(GTF_LNC_NAMES))


# In[26]:


# ---- A2. Compute lnc capture per dataset ----
def var_lnc_mask(adata):
    cols = [c for c in adata.var.columns if any(k in c.lower() for k in ["biotype","gene_type","type","biotypes"])]
    mask = None
    for col in cols:
        vals = adata.var[col].astype(str).str.lower()
        m = vals.apply(lambda x: any(k in x for k in LNC_KEYWORDS))
        if m.any():
            mask = m.values
            break
    if mask is None or not np.any(mask):
        mask = np.zeros(adata.n_vars, dtype=bool)
        if "gene_id" in adata.var.columns and len(GTF_LNC_IDS)>0:
            ids = adata.var["gene_id"].astype(str).values
            mask |= np.array([g in GTF_LNC_IDS for g in ids])
        if not np.any(mask) and len(GTF_LNC_NAMES)>0:
            names = adata.var_names.astype(str).values
            mask |= np.array([g in GTF_LNC_NAMES for g in names])
    return mask

lnc_rows = []
for ds in DATASETS:
    if not os.path.exists(ds.path): 
        continue
    a = load_any(ds)
    m = var_lnc_mask(a)
    n_lnc_in_matrix = int(m.sum())
    if n_lnc_in_matrix == 0:
        print(f"[INFO] {ds.name}: no lnc labels; skip.")
        lnc_rows.append(dict(dataset=ds.name, tech=ds.tech, lnc_genes_in_matrix=0,
                             lnc_median_genes_per_cell=np.nan, lnc_median_frac=np.nan))
        continue
    Xlnc = a[:, m].X
    gdet = (np.diff(Xlnc.indptr) if sp.issparse(Xlnc) else (Xlnc>0).sum(axis=1).A1)
    lnc_rows.append(dict(
        dataset=ds.name, tech=ds.tech, lnc_genes_in_matrix=n_lnc_in_matrix,
        lnc_median_genes_per_cell=float(np.median(gdet)),
        lnc_median_frac=float(np.median(gdet / np.maximum(1, nnz_per_row(a.X))))
    ))

df_lnc = pd.DataFrame(lnc_rows)
df_lnc.to_csv(os.path.join(OUTDIR, "lnc_capture_by_dataset.csv"), index=False)
display(df_lnc)


# In[27]:


# ---- A3. Plot lnc capture by technology ----
agg = df_lnc.groupby("tech")["lnc_median_genes_per_cell"].median().dropna().sort_values(ascending=False)
plt.figure(figsize=(8,4))
plt.bar(range(len(agg.index)), agg.values)
plt.xticks(range(len(agg.index)), agg.index, rotation=45, ha='right')
plt.ylabel("Median lncRNA genes per cell (per-technology median)")
plt.title("lncRNA capture")
if "SeekGene" in agg.index:
    i = list(agg.index).index("SeekGene")
    plt.text(i, agg.values[i]*1.03 if not np.isnan(agg.values[i]) else 0, "SeekGene", ha='center', va='bottom', fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "lnc_capture_by_tech.png"), dpi=150)
plt.show()


# 
# ## New Figure B — Composite Score
#
# Normalize key metrics (Min-Max), adjust by direction, and take an equally weighted average to output the total score for each technology:
# - Higher is better: `median_genes_per_cell`, `reproducibility_r_mean`, `silhouette_marker`, `lnc_median_genes_per_cell`, `complexity_auc`
# - Lower is better: `mean_zeros_cell`, `mean_zeros_gene`, `mt_pct_median`
# You can modify `WEIGHTS` to change the weights.
# 

# In[28]:


# ---- B1. Complexity AUC by technology ----
def complexity_auc_by_tech(comp_long):
    if comp_long is None or comp_long.empty:
        return pd.Series(dtype=float)
    aucs = {}
    for tech, sub in comp_long.groupby("tech"):
        sub = sub.sort_values("umi")
        if sub.empty:
            continue
        x = sub["umi"].values.astype(float)
        y = sub["median_genes"].values.astype(float)
        aucs[tech] = np.trapz(y, x)
    return pd.Series(aucs, name="complexity_auc")

auc_s = complexity_auc_by_tech(complexity_long if 'complexity_long' in globals() else None)
display(auc_s.sort_values(ascending=False))


# In[29]:


# ---- B2. Assemble per-technology table ----
pertech = (df.groupby("tech")[[
    "median_genes_per_cell","reproducibility_r_mean","silhouette_marker",
    "mean_zeros_cell","mean_zeros_gene","mt_pct_median"
]].median())

if 'df_lnc' in globals() and not df_lnc.empty:
    lnc_agg = df_lnc.groupby("tech")["lnc_median_genes_per_cell"].median()
    pertech = pertech.join(lnc_agg, how="left")
else:
    pertech["lnc_median_genes_per_cell"] = np.nan

if 'auc_s' in globals() and not auc_s.empty:
    pertech = pertech.join(auc_s, how="left")
else:
    pertech["complexity_auc"] = np.nan

display(pertech)
pertech.to_csv(os.path.join(OUTDIR, "per_technology_raw_metrics.csv"))


# In[30]:


# ---- B3. Normalize, weight, and score ----
WEIGHTS = {
    "median_genes_per_cell": 1.0,
    "reproducibility_r_mean": 1.0,
    "silhouette_marker": 1.0,
    "lnc_median_genes_per_cell": 1.0,
    "complexity_auc": 1.0,
    "mean_zeros_cell": 1.0,
    "mean_zeros_gene": 1.0,
    "mt_pct_median": 1.0
}

higher_better = {
    "median_genes_per_cell": True,
    "reproducibility_r_mean": True,
    "silhouette_marker": True,
    "lnc_median_genes_per_cell": True,
    "complexity_auc": True,
    "mean_zeros_cell": False,
    "mean_zeros_gene": False,
    "mt_pct_median": False
}

def minmax_norm(s):
    s = s.copy()
    if s.notna().sum() <= 1:
        return pd.Series(np.zeros_like(s), index=s.index)
    vmin, vmax = s.min(skipna=True), s.max(skipna=True)
    if np.isclose(vmin, vmax):
        return pd.Series(np.ones_like(s), index=s.index)
    return (s - vmin) / (vmax - vmin)

norm_cols = {}
for col, hb in higher_better.items():
    if col not in pertech.columns:
        continue
    s = pertech[col]
    n = minmax_norm(s)
    if not hb:
        n = 1 - n
    norm_cols[col] = n

norm_df = pd.DataFrame(norm_cols)
scores = {}
for tech, row in norm_df.iterrows():
    w_on = {k:w for k,w in WEIGHTS.items() if (k in row.index and not pd.isna(row[k]))}
    if not w_on:
        scores[tech] = np.nan
        continue
    s = sum(row[k]*w for k,w in w_on.items()) / sum(w_on.values())
    scores[tech] = float(s)

score_s = pd.Series(scores, name="composite_score").sort_values(ascending=False)
display(score_s)
score_s.to_csv(os.path.join(OUTDIR, "composite_scores.csv"))


# In[31]:


# ---- B4. Plot composite scores ----
plt.figure(figsize=(9,4))
plt.bar(range(len(score_s.index)), score_s.values)
plt.xticks(range(len(score_s.index)), score_s.index, rotation=45, ha='right')
plt.ylabel("Composite score (0–1)")
plt.title("Overall technology score (equal weights)")
if "SeekGene" in score_s.index:
    i = list(score_s.index).index("SeekGene")
    plt.text(i, score_s.values[i]*1.03 if not np.isnan(score_s.values[i]) else 0, "SeekGene", ha='center', va='bottom', fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "composite_score_by_tech.png"), dpi=150)
plt.show()

if "SeekGene" in score_s.index:
    rank = list(score_s.index).index("SeekGene") + 1
    print(f"SeekGene rank: {rank} / {len(score_s)}  |  score = {score_s['SeekGene']:.3f}")
else:
    print("SeekGene not found in the ranking.")


# In[32]:


# ==== Quick summary ====
def tech_stat(col, higher_is_better=True):
    s = df.groupby("tech")[col].median().dropna()
    if s.empty:
        return None, None
    best = s.idxmax() if higher_is_better else s.idxmin()
    return s, best

sens, best_sens = tech_stat("median_genes_per_cell", True)
rep , best_rep  = tech_stat("reproducibility_r_mean", True)
silh, best_silh = tech_stat("silhouette_marker", True)

print("========= QUICK SUMMARY =========")
if sens is not None:
    print(f"[Sensitivity] best: {best_sens}  (median genes/cell = {sens[best_sens]:.0f})")
if rep is not None:
    print(f"[Reproducibility] best: {best_rep}  (split-half r = {rep[best_rep]:.3f})")
if silh is not None:
    print(f"[Separability] best: {best_silh}  (silhouette = {silh[best_silh]:.3f})")
print("=================================")
display(df[df.tech=="SeekGene"].sort_values("dataset"))

