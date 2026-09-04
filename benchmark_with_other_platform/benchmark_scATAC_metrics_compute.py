#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
os.chdir("/path/to/scATAC/brain")


# In[ ]:


# ===================== scATAC-seq QC (SnapATAC2) — robust all-in-one with panels & doublet =====================
# Features: HDF5 lock fix; import_fragments+whitelist; auto-clean peaks then FRiP (compatible with new/old signatures);
# optional peak matrix->n_peaks (on failure auto fallback to tile matrix->n_tiles); fraglen fallback; doublet run on tile;
# PyDataFrameElem-compatible export and filter; all figures written into a single PDF.

# ---------------- Environment fix (NFS/HDF5) ----------------
import os, sys, re, gzip, time, gc, warnings, logging
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
os.environ.setdefault("HDF5_DISABLE_VERSION_CHECK", "2")  # optional

# ---------------- Parameters (modify as needed) ----------------
FRAG   = "M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_atac_fragments.tsv.gz"
PEAKS  = "M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_atac_peaks.bed"   # if absent: None
TENXH5 = "M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_filtered_feature_bc_matrix.h5"  # if absent: None
GENOME = "mm10"   # choose: mm10 | mm39 | hg38 | hg19

OUTDIR = "qc_out"
PREFIX = "snap_nb"

# Filter thresholds
MIN_COUNTS, MAX_COUNTS = 1000, 200000
MIN_TSSE,   MIN_FRIP   = 6.0, 0.20
MAX_MITO,   MAX_DUP    = 0.20, 0.80

# Optional features
DO_NUCLEOSOME       = True       # mono-/nucleosome-free fragment ratio
NUC_MAX_LINES       = None       # debug rate limit; None=full
WRITE_PDF           = True       # uniform PDF export
REBUILD_PEAK_MATRIX = True       # build peak matrix (on failure auto fallback to tile)
DOUBLETS_THRESHOLD  = 0.5        # doublet_score threshold
SCRUBLET_JOBS       = 4          # scrublet parallel jobs

# Raw backend priority: try hdf5 first, fall back to zarr
BACKEND_PREF = "hdf5"

# ---------------- Dependencies ----------------
import numpy as np, pandas as pd, h5py
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
warnings.filterwarnings("ignore", message=".*_import_from_c.*")  # silence pyo3-polars old-interface warning
try:
    import snapatac2 as snap
except Exception as e:
    print("ERROR: snapatac2 is required. Suggestion: conda create -n snap2 -y python=3.10 snapatac2 h5py pandas matplotlib", file=sys.stderr)
    raise

# ---------------- Utility functions ----------------
def ensure_dir(p): os.makedirs(p, exist_ok=True); return p

def map_genome(label: str):
    label=(label or "").strip().lower()
    if label in ["mm10","mouse","grcm38"]: return snap.genome.mm10
    if label in ["mm39","grcm39"]:         return snap.genome.mm39
    if label in ["hg38","grch38","human"]: return snap.genome.hg38
    if label in ["hg19","grch37"]:         return snap.genome.hg19
    raise ValueError(f"Unsupported genome: {label}")

def read_tenx_barcodes_from_h5(h5_path: str):
    if not h5_path or (not os.path.exists(h5_path)): return None
    try:
        with h5py.File(h5_path, "r") as f:
            if "matrix" in f and "barcodes" in f["matrix"]:
                ds=f["matrix/barcodes"]; 
                return set(x.decode("utf-8") if isinstance(x,bytes) else str(x) for x in ds[()])
            for key in ["barcodes","obs_names","barcodes_names"]:
                if key in f:
                    ds=f[key]; 
                    return set(x.decode("utf-8") if isinstance(x,bytes) else str(x) for x in ds[()])
    except Exception as e:
        print(f"[WARN] Failed to read 10x H5 barcodes: {e}")
    return None

def sanitize_bed(infile: str, outfile: str):
    """Ensure a standard 3-column BED; remove #/track/browser and invalid lines."""
    if infile is None or (not os.path.exists(infile)):
        raise FileNotFoundError(f"Peaks file not found: {infile}")
    op_in = gzip.open if infile.endswith(".gz") else open
    bad = kept = 0
    with op_in(infile, "rt") as fin, open(outfile, "w") as fout:
        for ln in fin:
            s=ln.strip()
            if (not s) or s.startswith("#") or s.startswith("track") or s.startswith("browser"):
                bad+=1; continue
            parts = re.split(r"\s+", s)
            if len(parts) < 3:
                bad+=1; continue
            try:
                a = int(parts[1]); b = int(parts[2])
            except Exception:
                bad+=1; continue
            if b <= a:
                bad+=1; continue
            print(f"{parts[0]}\t{a}\t{b}", file=fout); kept += 1
    print(f"[Clean BED] kept={kept:,} dropped={bad:,}  ->  {outfile}")

def compute_nucleosome_signal_by_streaming(frag_path, keep_barcodes=None,
                                           nfr_max=147, mono_min=147, mono_max=294, max_lines=None):
    nfr, mono = {}, {}
    i=0
    with gzip.open(frag_path, "rt") as fh:
        for line in fh:
            if max_lines is not None and i>=max_lines: break
            i+=1
            if not line or line[0]=="#": continue
            parts=line.rstrip("\n").split("\t")
            if len(parts)<4: continue
            try: L = int(parts[2]) - int(parts[1])
            except Exception: continue
            bc = parts[3]
            if keep_barcodes is not None and bc not in keep_barcodes: continue
            if L < nfr_max:               nfr[bc]  = nfr.get(bc,0)  + 1
            elif mono_min <= L < mono_max: mono[bc] = mono.get(bc,0) + 1
    ks = (keep_barcodes if keep_barcodes is not None else set(list(nfr.keys())+list(mono.keys())))
    return {bc: (mono.get(bc,0)/(nfr.get(bc,0) if nfr.get(bc,0)>0 else 1.0)) for bc in ks}

def barcode_knee_plot(counts, title):
    vals = np.sort(np.asarray(counts))[::-1]; ranks = np.arange(1,len(vals)+1)
    plt.figure(); plt.plot(ranks, vals); plt.xscale("log"); plt.yscale("log")
    plt.xlabel("barcode rank (log10)"); plt.ylabel("n_fragment (log10)"); plt.title(title); plt.tight_layout()

def hist_plot(x, title, xlabel):
    v=np.asarray(x, dtype=float)
    plt.figure(); plt.hist(v[~np.isnan(v)], bins=60)
    plt.title(title); plt.xlabel(xlabel); plt.ylabel("count"); plt.tight_layout()

def scatter_plot(x,y,title,xlabel,ylabel,logx=False,logy=False):
    plt.figure(); plt.scatter(x,y,s=3,alpha=0.6)
    if logx: plt.xscale("log")
    if logy: plt.yscale("log")
    plt.title(title); plt.xlabel(xlabel); plt.ylabel(ylabel); plt.tight_layout()

def close_all_handles(globs: dict):
    for k,v in list(globs.items()):
        if hasattr(v, "close"):
            try: v.close()
            except Exception: pass
    gc.collect()

def unique_path(base_dir, base_name, backend):
    ts = time.strftime("%Y%m%d-%H%M%S"); pid = os.getpid()
    return os.path.join(base_dir, f"{base_name}.raw.{ts}.{pid}.{('h5ad' if backend=='hdf5' else 'zarr')}")

def _safe_row_nnz(X):
    try:
        import scipy.sparse as sp
        if sp.issparse(X): return np.asarray(X.getnnz(axis=1)).ravel()
        Xn = np.asarray(X)
        if Xn.ndim == 2: return np.count_nonzero(Xn, axis=1)
    except Exception:
        pass
    return None

def frag_len_histogram(frag_path, max_len=800, step=1):
    """Streaming fallback: returns x(1..max_len), y(percent), overflow(count of >max_len)"""
    bins = np.zeros(max_len + 1, dtype=np.int64)  # index 0 collects overflow
    with gzip.open(frag_path, "rt") as fh:
        for ln in fh:
            if not ln or ln[0] == "#": continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 3: continue
            try: L = int(p[2]) - int(p[1])
            except Exception: continue
            if 0 < L <= max_len: bins[L] += 1
            else: bins[0] += 1
    idx = np.arange(1, max_len+1, step)
    val = bins[1:max_len+1:step].astype(float)
    val = val/val.sum() if val.sum()>0 else val
    return idx, val, int(bins[0])

# ---------------- Main flow: import + basic QC ----------------
ensure_dir(OUTDIR)
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
genome_obj = map_genome(GENOME)

barcodes_keep = read_tenx_barcodes_from_h5(TENXH5)
if barcodes_keep is not None:
    print(f"[INFO] Read {len(barcodes_keep):,} barcodes from 10x H5; will be used as whitelist.")

close_all_handles(globals())
backend   = BACKEND_PREF
raw_path  = unique_path(OUTDIR, PREFIX, backend)

print("[INFO] Importing fragments and computing n_fragment / frac_dup / frac_mito ...")
try:
    adata = snap.pp.import_fragments(
        FRAG, chrom_sizes=genome_obj, file=raw_path,
        sorted_by_barcode=False, whitelist=list(barcodes_keep) if barcodes_keep is not None else None,
        backend=("hdf5" if backend=="hdf5" else "zarr")
    )
except RuntimeError as e:
    print(f"[WARN] HDF5 import failed ({e}); retrying with Zarr backend.")
    backend  = "zarr"
    raw_path = unique_path(OUTDIR, PREFIX, backend)
    adata    = snap.pp.import_fragments(
        FRAG, chrom_sizes=genome_obj, file=raw_path,
        sorted_by_barcode=False, whitelist=list(barcodes_keep) if barcodes_keep is not None else None,
        backend="zarr"
    )

print("[INFO] Computing TSSe ...")
snap.metrics.tsse(adata, genome_obj)
print("[INFO] TSSe done")

# ---------------- FRiP (compute after cleaning peaks) ----------------
frip_added = False
PEAKS_CLEAN = None
if PEAKS and os.path.exists(PEAKS):
    PEAKS_CLEAN = os.path.join(OUTDIR, f"{PREFIX}.peaks.clean.bed")
    sanitize_bed(PEAKS, PEAKS_CLEAN)
    print("[INFO] Computing FRiP ...")
    try:
        snap.metrics.frip(adata, {"frip": PEAKS_CLEAN})     # new signature
        frip_added = True
    except TypeError:
        snap.metrics.frip(adata, PEAKS_CLEAN)               # old signature
        if "frip" not in adata.obs:
            for c in ["FRiP","Frip","FRIP"]:
                if c in adata.obs: adata.obs["frip"]=adata.obs[c]; break
        frip_added = True
    except Exception as e:
        print(f"[WARN] FRiP computation failed: {e}")
print("[INFO] FRiP done" if frip_added else "[INFO] FRiP not computed (no peaks or failed)")

# ---------------- peak matrix -> n_peaks (fall back to tile on failure) ----------------
peakA = None
if REBUILD_PEAK_MATRIX and PEAKS and os.path.exists(PEAKS_CLEAN if PEAKS_CLEAN else PEAKS):
    try:
        peak_matrix_path = os.path.join(OUTDIR, f"{PREFIX}.peak_matrix.{('h5ad' if backend=='hdf5' else 'zarr')}")
        peakA = snap.pp.make_peak_matrix(adata, file=peak_matrix_path, peak_file=(PEAKS_CLEAN or PEAKS))
        nnz = _safe_row_nnz(peakA.X)
        if nnz is not None:
            adata.obs["n_peaks"] = nnz
            print("[INFO] Wrote adata.obs['n_peaks']")
    except Exception as e:
        print(f"[WARN] Failed to build peak matrix: {e}")

tileA = None
if "n_peaks" not in adata.obs:
    try:
        print("[INFO] Building tile(5kb) matrix as fallback ...")
        tile_path = os.path.join(OUTDIR, f"{PREFIX}.tile_5kb.zarr")
        tileA = snap.pp.make_tile_matrix(adata, file=tile_path, bin_size=5000)
        X = tileA.X
        import scipy.sparse as sp
        n_tiles = (X > 0).sum(axis=1).A1 if sp.issparse(X) else np.count_nonzero(np.asarray(X), axis=1)
        adata.obs["n_tiles"] = np.asarray(n_tiles, dtype=float)
        print("[INFO] Wrote adata.obs['n_tiles'] (tile-based)")
    except Exception as e:
        print(f"[WARN] Tile matrix failed: {e}")

# ---------------- nucleosome signal (optional) ----------------
if DO_NUCLEOSOME:
    print("[INFO] Computing nucleosome signal (mono/NFR); may be slow ...")
    bc_set = set(adata.obs_names)
    nuc_signal = compute_nucleosome_signal_by_streaming(FRAG, keep_barcodes=bc_set, max_lines=NUC_MAX_LINES)
    adata.obs["nucleosome_signal"] = pd.Series(nuc_signal).reindex(adata.obs_names).astype(float)
    print("[INFO] nucleosome_signal added to adata.obs")

# ---------------- Export metrics (PyDataFrameElem-compatible) ----------------
ensure_dir(OUTDIR)
metrics_csv = os.path.join(OUTDIR, f"{PREFIX}.qc_metrics.csv")
wanted = ["n_fragment","frac_dup","frac_mito","tsse","frip","n_peaks","n_tiles","nucleosome_signal"]
avail  = [k for k in wanted if k in adata.obs]
df = pd.DataFrame(index=adata.obs_names)
for k in avail: df[k] = np.asarray(adata.obs[k])
df.to_csv(metrics_csv); print("[INFO] Metrics CSV:", metrics_csv)

# ---------------- Filter & subset export ----------------
n_fragment = np.asarray(adata.obs["n_fragment"], dtype=float)
tsse       = np.asarray(adata.obs["tsse"],       dtype=float)
frac_mito  = np.asarray(adata.obs["frac_mito"],  dtype=float)
frac_dup   = np.asarray(adata.obs["frac_dup"],   dtype=float)
mask = (
    (n_fragment >= MIN_COUNTS) &
    (n_fragment <= MAX_COUNTS) &
    (tsse       >= MIN_TSSE)   &
    (frac_mito  <= MAX_MITO)   &
    (frac_dup   <= MAX_DUP)
)
if ("frip" in adata.obs) and (PEAKS and os.path.exists(PEAKS)):
    frip = np.asarray(adata.obs["frip"], dtype=float)
    mask = mask & (frip >= MIN_FRIP)
adata.obs["qc_pass"] = mask

idx = np.where(mask)[0]
obs_names_list = list(adata.obs_names)
pass_barcodes  = [obs_names_list[i] for i in idx]

filtered_h5ad = os.path.join(OUTDIR, f"{PREFIX}.filtered.h5ad")
try:
    _ = adata.subset(obs_indices=idx, inplace=False, out=filtered_h5ad, backend="hdf5")
    print("[INFO] Filtered h5ad:", filtered_h5ad, "(subset/out)")
except Exception as e:
    print(f"[WARN] subset(out=...) failed; re-importing using whitelist: {e}")
    tmp = os.path.join(OUTDIR, f"{PREFIX}.filtered.tmp.h5ad")
    adata_f = snap.pp.import_fragments(
        FRAG, chrom_sizes=genome_obj, file=tmp, sorted_by_barcode=False, whitelist=pass_barcodes
    )
    snap.metrics.tsse(adata_f, genome_obj)
    if PEAKS and os.path.exists(PEAKS_CLEAN or PEAKS):
        try:    snap.metrics.frip(adata_f, {"frip": (PEAKS_CLEAN or PEAKS)})
        except: 
            try: snap.metrics.frip(adata_f, (PEAKS_CLEAN or PEAKS))
            except: pass
    os.replace(tmp, filtered_h5ad)
    if hasattr(adata_f, "close"):
        try: adata_f.close()
        except Exception: pass
    print("[INFO] Filtered h5ad:", filtered_h5ad, "(re-import whitelist)")

filtered_barcodes_txt = os.path.join(OUTDIR, f"{PREFIX}.filtered_barcodes.txt")
pd.Series(pass_barcodes).to_csv(filtered_barcodes_txt, index=False, header=False)
print("[INFO] Passing barcodes:", filtered_barcodes_txt)

# ---------------- Scrublet doublet (run on tile or peak) ----------------
dbl_key = None
srcA_for_dbl = None
if tileA is not None:
    srcA_for_dbl = tileA
elif peakA is not None:
    srcA_for_dbl = peakA

if srcA_for_dbl is None:
    try:
        # Last-resort fallback: temporarily build a tile matrix just for doublet
        print("[INFO] Temporarily building tile matrix for doublet ...")
        tmp_tile = os.path.join(OUTDIR, f"{PREFIX}.tile_5kb.doublet.zarr")
        srcA_for_dbl = snap.pp.make_tile_matrix(adata, file=tmp_tile, bin_size=5000)
    except Exception as e:
        print(f"[WARN] Failed to build doublet matrix: {e}")

if srcA_for_dbl is not None:
    try:
        snap.pp.scrublet(srcA_for_dbl, n_jobs=SCRUBLET_JOBS)  # writes srcA_for_dbl.obs['doublet_score']
        s = pd.Series(np.asarray(srcA_for_dbl.obs["doublet_score"]), index=list(srcA_for_dbl.obs_names), name="doublet_score")
        adata.obs["doublet_score"] = s.reindex(adata.obs_names).astype(float).values
        adata.obs["doublet"] = np.asarray(adata.obs["doublet_score"], dtype=float) > float(DOUBLETS_THRESHOLD)
        dbl_key = "doublet_score"
        print("[INFO] doublet_score written; threshold >", float(DOUBLETS_THRESHOLD), " marked as doublet")
    except Exception as e:
        print(f"[WARN] scrublet failed: {e}")
else:
    print("[INFO] No doublet input matrix built (skipping scrublet)")

# ---------------- Unified PDF report (all pages) ----------------
pdf_path = os.path.join(OUTDIR, f"{PREFIX}.qc_report.pdf") if WRITE_PDF else None
pdf = PdfPages(pdf_path) if WRITE_PDF else None
def savefig_or_inline():
    if pdf is not None: pdf.savefig(); plt.close()
    else: plt.show()

# 1) Knee
barcode_knee_plot(adata.obs["n_fragment"].to_numpy(), "Barcode knee (n_fragment)"); savefig_or_inline()

# 2) Fragment length distribution (official or fallback)
ok_fsd = True
try:
    fsd = snap.metrics.frag_size_distr(adata)
    ok_fsd = fsd is not None
except Exception:
    ok_fsd = False
if ok_fsd:
    x = np.arange(len(fsd)); y = np.array(fsd, dtype=float); x = x[1:]; y = y[1:]
    plt.figure(); plt.plot(x, y); plt.yscale("log")
    plt.xlabel("fragment length (bp)"); plt.ylabel("count (log)")
    plt.title("Fragment length distribution"); plt.tight_layout(); savefig_or_inline()
    # Percentage version (closer to common reports)
    y_pct = y / y.sum() if y.sum()>0 else y
    plt.figure(); plt.plot(x, y_pct)
    plt.xlabel("Size of Fragments (bp)"); plt.ylabel("Fragments (%)")
    plt.title("Fragment Size Distribution (percent)"); plt.tight_layout(); savefig_or_inline()
else:
    print("[INFO] Official frag_size_distr unavailable; using fallback streaming version")
    x, y, overflow = frag_len_histogram(FRAG, max_len=800, step=1)
    plt.figure(); plt.plot(x, y)
    plt.xlabel("Fragment size (bp)"); plt.ylabel("Fraction")
    plt.title(f"Fragment length distribution (fallback); overflow>{800}bp: {overflow:,}")
    plt.tight_layout(); savefig_or_inline()

# 3) TSSe vs n_fragment (scatter)
scatter_plot(np.asarray(adata.obs["n_fragment"]), np.asarray(adata.obs["tsse"]),
             "TSSe vs n_fragment", "n_fragment", "TSSe", logx=True); savefig_or_inline()

# 4) TSSe vs log10(nFrags) density + threshold
nf = np.asarray(adata.obs["n_fragment"], dtype=float)
ts = np.asarray(adata.obs["tsse"], dtype=float)
log_nf = np.log10(np.clip(nf, 1, None))
plt.figure(figsize=(6,6))
plt.hist2d(log_nf, ts, bins=200, cmap="viridis"); cbar=plt.colorbar(); cbar.set_label("density")
plt.axvline(np.log10(max(MIN_COUNTS,1)), ls="--", c="k"); plt.axhline(MIN_TSSE, ls="--", c="k")
med_nf = float(np.median(nf)); med_ts = float(np.median(ts))
n_pass = int(np.sum(np.asarray(adata.obs["qc_pass"], dtype=bool)))
txt = (f"{PREFIX}\n"
       f"nCells Pass Filter = {n_pass}\n"
       f"Median Frags = {med_nf:.1f}\n"
       f"Median TSS Enrichment = {med_ts:.4f}")
plt.text(log_nf.min()+0.05, ts.max()*0.95, txt, va="top")
plt.xlabel("Log10 (Unique Fragments)"); plt.ylabel("TSS Enrichment")
plt.title("TSSe vs Log10(Unique Fragments) — density"); plt.tight_layout(); savefig_or_inline()

# 5) Various histograms
if "frip" in adata.obs: hist_plot(np.asarray(adata.obs["frip"]), "Histogram: FRiP", "FRiP"); savefig_or_inline()
hist_plot(np.asarray(adata.obs["n_fragment"]), "Histogram: n_fragment", "n_fragment"); savefig_or_inline()
hist_plot(np.asarray(adata.obs["tsse"]),      "Histogram: TSSe",       "TSSe");       savefig_or_inline()
if "frac_mito" in adata.obs: hist_plot(np.asarray(adata.obs["frac_mito"]), "Histogram: frac_mito", "frac_mito"); savefig_or_inline()
if "frac_dup"  in adata.obs: hist_plot(np.asarray(adata.obs["frac_dup"]),  "Histogram: frac_dup",  "frac_dup");  savefig_or_inline()
if "n_peaks"   in adata.obs: hist_plot(np.asarray(adata.obs["n_peaks"]),   "Histogram: n_peaks per cell", "n_peaks"); savefig_or_inline()
if "n_tiles"   in adata.obs: hist_plot(np.asarray(adata.obs["n_tiles"]),   "Histogram: n_tiles per cell", "n_tiles"); savefig_or_inline()

# 6) Association scatter
if "frip" in adata.obs:
    scatter_plot(np.asarray(adata.obs["n_fragment"]), np.asarray(adata.obs["frip"]),
                 "FRiP vs n_fragment", "n_fragment","FRiP", logx=True); savefig_or_inline()
    scatter_plot(np.asarray(adata.obs["frip"]), np.asarray(adata.obs["tsse"]),
                 "TSSe vs FRiP", "FRiP","TSSe"); savefig_or_inline()

if "frac_mito" in adata.obs:
    scatter_plot(np.asarray(adata.obs["n_fragment"]), np.asarray(adata.obs["frac_mito"]),
                 "Mito fraction vs n_fragment","n_fragment","frac_mito", logx=True); savefig_or_inline()
if "frac_dup" in adata.obs:
    scatter_plot(np.asarray(adata.obs["n_fragment"]), np.asarray(adata.obs["frac_dup"]),
                 "Duplicate fraction vs n_fragment","n_fragment","frac_dup", logx=True); savefig_or_inline()

# 7) Pass/fail (TSSe vs n_fragment)
try:
    x = np.asarray(adata.obs["n_fragment"], dtype=float)
    y = np.asarray(adata.obs["tsse"],       dtype=float)
    keep = np.asarray(adata.obs["qc_pass"], dtype=bool)
    plt.figure()
    plt.scatter(x[~keep], y[~keep], s=3, alpha=0.4, label="fail")
    plt.scatter(x[keep],  y[keep],  s=3, alpha=0.6, label="pass")
    plt.xscale("log"); plt.xlabel("n_fragment"); plt.ylabel("TSSe")
    plt.title("QC pass/fail (TSSe vs n_fragment)"); plt.legend(); plt.tight_layout(); savefig_or_inline()
except Exception as e:
    print("[WARN] pass/fail visualization failed:", e)

# 8) Threshold-line histograms
for col, thr, label in [("tsse", MIN_TSSE, "TSSe >= MIN_TSSE"),
                        ("frac_mito", MAX_MITO, "frac_mito <= MAX_MITO"),
                        ("frac_dup",  MAX_DUP,  "frac_dup <= MAX_DUP")]:
    if col in adata.obs:
        v = np.asarray(adata.obs[col], dtype=float)
        plt.figure(); plt.hist(v[~np.isnan(v)], bins=60); plt.axvline(thr)
        plt.title(f"Threshold: {label}"); plt.xlabel(col); plt.ylabel("count"); plt.tight_layout(); savefig_or_inline()

if ("frip" in adata.obs) and (PEAKS and os.path.exists(PEAKS)):
    v = np.asarray(adata.obs["frip"], dtype=float)
    plt.figure(); plt.hist(v[~np.isnan(v)], bins=60); plt.axvline(MIN_FRIP)
    plt.title("Threshold: FRiP >= MIN_FRIP"); plt.xlabel("frip"); plt.ylabel("count"); plt.tight_layout(); savefig_or_inline()

# 9) Doublet-related plots (if successful)
if dbl_key and dbl_key in adata.obs:
    v = np.asarray(adata.obs[dbl_key], dtype=float); thr = float(DOUBLETS_THRESHOLD)
    plt.figure(); plt.hist(v[~np.isnan(v)], bins=60); plt.axvline(thr, ls="--", c="r")
    plt.xlabel("doublet_score"); plt.ylabel("count"); plt.title("Histogram: doublet_score")
    plt.tight_layout(); savefig_or_inline()

    plt.figure(); plt.scatter(np.asarray(adata.obs["n_fragment"]), v, s=3, alpha=0.5)
    plt.xscale("log"); plt.xlabel("n_fragment"); plt.ylabel("doublet_score")
    plt.title("doublet_score vs n_fragment"); plt.tight_layout(); savefig_or_inline()

    plt.figure(); plt.scatter(np.asarray(adata.obs["tsse"]), v, s=3, alpha=0.5)
    plt.xlabel("TSSe"); plt.ylabel("doublet_score"); plt.title("doublet_score vs TSSe")
    plt.tight_layout(); savefig_or_inline()

    plt.figure()
    x = np.asarray(adata.obs["n_fragment"], dtype=float)
    y = np.asarray(adata.obs["tsse"], dtype=float)
    kk = np.asarray(adata.obs["doublet"], dtype=bool)
    plt.scatter(x[~kk], y[~kk], s=3, alpha=0.3, label="singlet")
    plt.scatter(x[ kk], y[ kk], s=5, alpha=0.6, label="doublet")
    plt.xscale("log"); plt.xlabel("n_fragment"); plt.ylabel("TSSe")
    plt.title("TSSe vs n_fragment (doublet overlay)"); plt.legend()
    plt.tight_layout(); savefig_or_inline()

if WRITE_PDF and pdf is not None:
    pdf.close(); print("[INFO] Unified QC PDF:", pdf_path)

# Close any open handles
close_all_handles(globals())
# ===============================================================================================================


# In[3]:


# ===================== scATAC-seq QC (SnapATAC2) — robust all-in-one with panels & doublet =====================
# Features: HDF5 lock fix; import_fragments+whitelist; auto-clean peaks then FRiP (compatible with new/old signatures);
# optional peak matrix->n_peaks (on failure auto fallback to tile matrix->n_tiles); fraglen fallback; doublet run on tile;
# PyDataFrameElem-compatible export and filter; all figures written into a single PDF.

# ---------------- Environment fix (NFS/HDF5) ----------------
import os, sys, re, gzip, time, gc, warnings, logging
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
os.environ.setdefault("HDF5_DISABLE_VERSION_CHECK", "2")  # optional

# ---------------- Parameters (modify as needed) ----------------
FRAG   = "atac_v1_E18_brain_flash_5k_fragments.tsv.gz"
PEAKS  = "atac_v1_E18_brain_flash_5k_peaks.bed"   # if absent: None
TENXH5 = "atac_v1_E18_brain_flash_5k_filtered_peak_bc_matrix.h5"  # if absent: None
GENOME = "mm10"   # choose: mm10 | mm39 | hg38 | hg19

OUTDIR = "qc_out_atac"
PREFIX = "E18_brain_scATAC"

# Filter thresholds
MIN_COUNTS, MAX_COUNTS = 1000, 200000
MIN_TSSE,   MIN_FRIP   = 6.0, 0.20
MAX_MITO,   MAX_DUP    = 0.20, 0.80

# Optional features
DO_NUCLEOSOME       = True       # mono-/nucleosome-free fragment ratio
NUC_MAX_LINES       = None       # debug rate limit; None=full
WRITE_PDF           = True       # uniform PDF export
REBUILD_PEAK_MATRIX = True       # build peak matrix (on failure auto fallback to tile)
DOUBLETS_THRESHOLD  = 0.5        # doublet_score threshold
SCRUBLET_JOBS       = 4          # scrublet parallel jobs

# Raw backend priority: try hdf5 first, fall back to zarr
BACKEND_PREF = "hdf5"

# ---------------- Dependencies ----------------
import numpy as np, pandas as pd, h5py
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
warnings.filterwarnings("ignore", message=".*_import_from_c.*")  # silence pyo3-polars old-interface warning
try:
    import snapatac2 as snap
except Exception as e:
    print("ERROR: snapatac2 is required. Suggestion: conda create -n snap2 -y python=3.10 snapatac2 h5py pandas matplotlib", file=sys.stderr)
    raise

# ---------------- Utility functions ----------------
def ensure_dir(p): os.makedirs(p, exist_ok=True); return p

def map_genome(label: str):
    label=(label or "").strip().lower()
    if label in ["mm10","mouse","grcm38"]: return snap.genome.mm10
    if label in ["mm39","grcm39"]:         return snap.genome.mm39
    if label in ["hg38","grch38","human"]: return snap.genome.hg38
    if label in ["hg19","grch37"]:         return snap.genome.hg19
    raise ValueError(f"Unsupported genome: {label}")

def read_tenx_barcodes_from_h5(h5_path: str):
    if not h5_path or (not os.path.exists(h5_path)): return None
    try:
        with h5py.File(h5_path, "r") as f:
            if "matrix" in f and "barcodes" in f["matrix"]:
                ds=f["matrix/barcodes"]; 
                return set(x.decode("utf-8") if isinstance(x,bytes) else str(x) for x in ds[()])
            for key in ["barcodes","obs_names","barcodes_names"]:
                if key in f:
                    ds=f[key]; 
                    return set(x.decode("utf-8") if isinstance(x,bytes) else str(x) for x in ds[()])
    except Exception as e:
        print(f"[WARN] Failed to read 10x H5 barcodes: {e}")
    return None

def sanitize_bed(infile: str, outfile: str):
    """Ensure a standard 3-column BED; remove #/track/browser and invalid lines."""
    if infile is None or (not os.path.exists(infile)):
        raise FileNotFoundError(f"Peaks file not found: {infile}")
    op_in = gzip.open if infile.endswith(".gz") else open
    bad = kept = 0
    with op_in(infile, "rt") as fin, open(outfile, "w") as fout:
        for ln in fin:
            s=ln.strip()
            if (not s) or s.startswith("#") or s.startswith("track") or s.startswith("browser"):
                bad+=1; continue
            parts = re.split(r"\s+", s)
            if len(parts) < 3:
                bad+=1; continue
            try:
                a = int(parts[1]); b = int(parts[2])
            except Exception:
                bad+=1; continue
            if b <= a:
                bad+=1; continue
            print(f"{parts[0]}\t{a}\t{b}", file=fout); kept += 1
    print(f"[Clean BED] kept={kept:,} dropped={bad:,}  ->  {outfile}")

def compute_nucleosome_signal_by_streaming(frag_path, keep_barcodes=None,
                                           nfr_max=147, mono_min=147, mono_max=294, max_lines=None):
    nfr, mono = {}, {}
    i=0
    with gzip.open(frag_path, "rt") as fh:
        for line in fh:
            if max_lines is not None and i>=max_lines: break
            i+=1
            if not line or line[0]=="#": continue
            parts=line.rstrip("\n").split("\t")
            if len(parts)<4: continue
            try: L = int(parts[2]) - int(parts[1])
            except Exception: continue
            bc = parts[3]
            if keep_barcodes is not None and bc not in keep_barcodes: continue
            if L < nfr_max:               nfr[bc]  = nfr.get(bc,0)  + 1
            elif mono_min <= L < mono_max: mono[bc] = mono.get(bc,0) + 1
    ks = (keep_barcodes if keep_barcodes is not None else set(list(nfr.keys())+list(mono.keys())))
    return {bc: (mono.get(bc,0)/(nfr.get(bc,0) if nfr.get(bc,0)>0 else 1.0)) for bc in ks}

def barcode_knee_plot(counts, title):
    vals = np.sort(np.asarray(counts))[::-1]; ranks = np.arange(1,len(vals)+1)
    plt.figure(); plt.plot(ranks, vals); plt.xscale("log"); plt.yscale("log")
    plt.xlabel("barcode rank (log10)"); plt.ylabel("n_fragment (log10)"); plt.title(title); plt.tight_layout()

def hist_plot(x, title, xlabel):
    v=np.asarray(x, dtype=float)
    plt.figure(); plt.hist(v[~np.isnan(v)], bins=60)
    plt.title(title); plt.xlabel(xlabel); plt.ylabel("count"); plt.tight_layout()

def scatter_plot(x,y,title,xlabel,ylabel,logx=False,logy=False):
    plt.figure(); plt.scatter(x,y,s=3,alpha=0.6)
    if logx: plt.xscale("log")
    if logy: plt.yscale("log")
    plt.title(title); plt.xlabel(xlabel); plt.ylabel(ylabel); plt.tight_layout()

def close_all_handles(globs: dict):
    for k,v in list(globs.items()):
        if hasattr(v, "close"):
            try: v.close()
            except Exception: pass
    gc.collect()

def unique_path(base_dir, base_name, backend):
    ts = time.strftime("%Y%m%d-%H%M%S"); pid = os.getpid()
    return os.path.join(base_dir, f"{base_name}.raw.{ts}.{pid}.{('h5ad' if backend=='hdf5' else 'zarr')}")

def _safe_row_nnz(X):
    try:
        import scipy.sparse as sp
        if sp.issparse(X): return np.asarray(X.getnnz(axis=1)).ravel()
        Xn = np.asarray(X)
        if Xn.ndim == 2: return np.count_nonzero(Xn, axis=1)
    except Exception:
        pass
    return None

def frag_len_histogram(frag_path, max_len=800, step=1):
    """Streaming fallback: returns x(1..max_len), y(percent), overflow(count of >max_len)"""
    bins = np.zeros(max_len + 1, dtype=np.int64)  # index 0 collects overflow
    with gzip.open(frag_path, "rt") as fh:
        for ln in fh:
            if not ln or ln[0] == "#": continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 3: continue
            try: L = int(p[2]) - int(p[1])
            except Exception: continue
            if 0 < L <= max_len: bins[L] += 1
            else: bins[0] += 1
    idx = np.arange(1, max_len+1, step)
    val = bins[1:max_len+1:step].astype(float)
    val = val/val.sum() if val.sum()>0 else val
    return idx, val, int(bins[0])

# ---------------- Main flow: import + basic QC ----------------
ensure_dir(OUTDIR)
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
genome_obj = map_genome(GENOME)

barcodes_keep = read_tenx_barcodes_from_h5(TENXH5)
if barcodes_keep is not None:
    print(f"[INFO] Read {len(barcodes_keep):,} barcodes from 10x H5; will be used as whitelist.")

close_all_handles(globals())
backend   = BACKEND_PREF
raw_path  = unique_path(OUTDIR, PREFIX, backend)

print("[INFO] Importing fragments and computing n_fragment / frac_dup / frac_mito ...")
try:
    adata = snap.pp.import_fragments(
        FRAG, chrom_sizes=genome_obj, file=raw_path,
        sorted_by_barcode=False, whitelist=list(barcodes_keep) if barcodes_keep is not None else None,
        backend=("hdf5" if backend=="hdf5" else "zarr")
    )
except RuntimeError as e:
    print(f"[WARN] HDF5 import failed ({e}); retrying with Zarr backend.")
    backend  = "zarr"
    raw_path = unique_path(OUTDIR, PREFIX, backend)
    adata    = snap.pp.import_fragments(
        FRAG, chrom_sizes=genome_obj, file=raw_path,
        sorted_by_barcode=False, whitelist=list(barcodes_keep) if barcodes_keep is not None else None,
        backend="zarr"
    )

print("[INFO] Computing TSSe ...")
snap.metrics.tsse(adata, genome_obj)
print("[INFO] TSSe done")

# ---------------- FRiP (compute after cleaning peaks) ----------------
frip_added = False
PEAKS_CLEAN = None
if PEAKS and os.path.exists(PEAKS):
    PEAKS_CLEAN = os.path.join(OUTDIR, f"{PREFIX}.peaks.clean.bed")
    sanitize_bed(PEAKS, PEAKS_CLEAN)
    print("[INFO] Computing FRiP ...")
    try:
        snap.metrics.frip(adata, {"frip": PEAKS_CLEAN})     # new signature
        frip_added = True
    except TypeError:
        snap.metrics.frip(adata, PEAKS_CLEAN)               # old signature
        if "frip" not in adata.obs:
            for c in ["FRiP","Frip","FRIP"]:
                if c in adata.obs: adata.obs["frip"]=adata.obs[c]; break
        frip_added = True
    except Exception as e:
        print(f"[WARN] FRiP computation failed: {e}")
print("[INFO] FRiP done" if frip_added else "[INFO] FRiP not computed (no peaks or failed)")

# ---------------- peak matrix -> n_peaks (fall back to tile on failure) ----------------
peakA = None
if REBUILD_PEAK_MATRIX and PEAKS and os.path.exists(PEAKS_CLEAN if PEAKS_CLEAN else PEAKS):
    try:
        peak_matrix_path = os.path.join(OUTDIR, f"{PREFIX}.peak_matrix.{('h5ad' if backend=='hdf5' else 'zarr')}")
        peakA = snap.pp.make_peak_matrix(adata, file=peak_matrix_path, peak_file=(PEAKS_CLEAN or PEAKS))
        nnz = _safe_row_nnz(peakA.X)
        if nnz is not None:
            adata.obs["n_peaks"] = nnz
            print("[INFO] Wrote adata.obs['n_peaks']")
    except Exception as e:
        print(f"[WARN] Failed to build peak matrix: {e}")

tileA = None
if "n_peaks" not in adata.obs:
    try:
        print("[INFO] Building tile(5kb) matrix as fallback ...")
        tile_path = os.path.join(OUTDIR, f"{PREFIX}.tile_5kb.zarr")
        tileA = snap.pp.make_tile_matrix(adata, file=tile_path, bin_size=5000)
        X = tileA.X
        import scipy.sparse as sp
        n_tiles = (X > 0).sum(axis=1).A1 if sp.issparse(X) else np.count_nonzero(np.asarray(X), axis=1)
        adata.obs["n_tiles"] = np.asarray(n_tiles, dtype=float)
        print("[INFO] Wrote adata.obs['n_tiles'] (tile-based)")
    except Exception as e:
        print(f"[WARN] Tile matrix failed: {e}")

# ---------------- nucleosome signal (optional) ----------------
if DO_NUCLEOSOME:
    print("[INFO] Computing nucleosome signal (mono/NFR); may be slow ...")
    bc_set = set(adata.obs_names)
    nuc_signal = compute_nucleosome_signal_by_streaming(FRAG, keep_barcodes=bc_set, max_lines=NUC_MAX_LINES)
    adata.obs["nucleosome_signal"] = pd.Series(nuc_signal).reindex(adata.obs_names).astype(float)
    print("[INFO] nucleosome_signal added to adata.obs")

# ---------------- Export metrics (PyDataFrameElem-compatible) ----------------
ensure_dir(OUTDIR)
metrics_csv = os.path.join(OUTDIR, f"{PREFIX}.qc_metrics.csv")
wanted = ["n_fragment","frac_dup","frac_mito","tsse","frip","n_peaks","n_tiles","nucleosome_signal"]
avail  = [k for k in wanted if k in adata.obs]
df = pd.DataFrame(index=adata.obs_names)
for k in avail: df[k] = np.asarray(adata.obs[k])
df.to_csv(metrics_csv); print("[INFO] Metrics CSV:", metrics_csv)

# ---------------- Filter & subset export ----------------
n_fragment = np.asarray(adata.obs["n_fragment"], dtype=float)
tsse       = np.asarray(adata.obs["tsse"],       dtype=float)
frac_mito  = np.asarray(adata.obs["frac_mito"],  dtype=float)
frac_dup   = np.asarray(adata.obs["frac_dup"],   dtype=float)
mask = (
    (n_fragment >= MIN_COUNTS) &
    (n_fragment <= MAX_COUNTS) &
    (tsse       >= MIN_TSSE)   &
    (frac_mito  <= MAX_MITO)   &
    (frac_dup   <= MAX_DUP)
)
if ("frip" in adata.obs) and (PEAKS and os.path.exists(PEAKS)):
    frip = np.asarray(adata.obs["frip"], dtype=float)
    mask = mask & (frip >= MIN_FRIP)
adata.obs["qc_pass"] = mask

idx = np.where(mask)[0]
obs_names_list = list(adata.obs_names)
pass_barcodes  = [obs_names_list[i] for i in idx]

filtered_h5ad = os.path.join(OUTDIR, f"{PREFIX}.filtered.h5ad")
try:
    _ = adata.subset(obs_indices=idx, inplace=False, out=filtered_h5ad, backend="hdf5")
    print("[INFO] Filtered h5ad:", filtered_h5ad, "(subset/out)")
except Exception as e:
    print(f"[WARN] subset(out=...) failed; re-importing using whitelist: {e}")
    tmp = os.path.join(OUTDIR, f"{PREFIX}.filtered.tmp.h5ad")
    adata_f = snap.pp.import_fragments(
        FRAG, chrom_sizes=genome_obj, file=tmp, sorted_by_barcode=False, whitelist=pass_barcodes
    )
    snap.metrics.tsse(adata_f, genome_obj)
    if PEAKS and os.path.exists(PEAKS_CLEAN or PEAKS):
        try:    snap.metrics.frip(adata_f, {"frip": (PEAKS_CLEAN or PEAKS)})
        except: 
            try: snap.metrics.frip(adata_f, (PEAKS_CLEAN or PEAKS))
            except: pass
    os.replace(tmp, filtered_h5ad)
    if hasattr(adata_f, "close"):
        try: adata_f.close()
        except Exception: pass
    print("[INFO] Filtered h5ad:", filtered_h5ad, "(re-import whitelist)")

filtered_barcodes_txt = os.path.join(OUTDIR, f"{PREFIX}.filtered_barcodes.txt")
pd.Series(pass_barcodes).to_csv(filtered_barcodes_txt, index=False, header=False)
print("[INFO] Passing barcodes:", filtered_barcodes_txt)

# ---------------- Scrublet doublet (run on tile or peak) ----------------
dbl_key = None
srcA_for_dbl = None
if tileA is not None:
    srcA_for_dbl = tileA
elif peakA is not None:
    srcA_for_dbl = peakA

if srcA_for_dbl is None:
    try:
        # Last-resort fallback: temporarily build a tile matrix just for doublet
        print("[INFO] Temporarily building tile matrix for doublet ...")
        tmp_tile = os.path.join(OUTDIR, f"{PREFIX}.tile_5kb.doublet.zarr")
        srcA_for_dbl = snap.pp.make_tile_matrix(adata, file=tmp_tile, bin_size=5000)
    except Exception as e:
        print(f"[WARN] Failed to build doublet matrix: {e}")

if srcA_for_dbl is not None:
    try:
        snap.pp.scrublet(srcA_for_dbl, n_jobs=SCRUBLET_JOBS)  # writes srcA_for_dbl.obs['doublet_score']
        s = pd.Series(np.asarray(srcA_for_dbl.obs["doublet_score"]), index=list(srcA_for_dbl.obs_names), name="doublet_score")
        adata.obs["doublet_score"] = s.reindex(adata.obs_names).astype(float).values
        adata.obs["doublet"] = np.asarray(adata.obs["doublet_score"], dtype=float) > float(DOUBLETS_THRESHOLD)
        dbl_key = "doublet_score"
        print("[INFO] doublet_score written; threshold >", float(DOUBLETS_THRESHOLD), " marked as doublet")
    except Exception as e:
        print(f"[WARN] scrublet failed: {e}")
else:
    print("[INFO] No doublet input matrix built (skipping scrublet)")

# ---------------- Unified PDF report (all pages) ----------------
pdf_path = os.path.join(OUTDIR, f"{PREFIX}.qc_report.pdf") if WRITE_PDF else None
pdf = PdfPages(pdf_path) if WRITE_PDF else None
def savefig_or_inline():
    if pdf is not None: pdf.savefig(); plt.close()
    else: plt.show()

# 1) Knee
barcode_knee_plot(adata.obs["n_fragment"].to_numpy(), "Barcode knee (n_fragment)"); savefig_or_inline()

# 2) Fragment length distribution (official or fallback)
ok_fsd = True
try:
    fsd = snap.metrics.frag_size_distr(adata)
    ok_fsd = fsd is not None
except Exception:
    ok_fsd = False
if ok_fsd:
    x = np.arange(len(fsd)); y = np.array(fsd, dtype=float); x = x[1:]; y = y[1:]
    plt.figure(); plt.plot(x, y); plt.yscale("log")
    plt.xlabel("fragment length (bp)"); plt.ylabel("count (log)")
    plt.title("Fragment length distribution"); plt.tight_layout(); savefig_or_inline()
    # Percentage version (closer to common reports)
    y_pct = y / y.sum() if y.sum()>0 else y
    plt.figure(); plt.plot(x, y_pct)
    plt.xlabel("Size of Fragments (bp)"); plt.ylabel("Fragments (%)")
    plt.title("Fragment Size Distribution (percent)"); plt.tight_layout(); savefig_or_inline()
else:
    print("[INFO] Official frag_size_distr unavailable; using fallback streaming version")
    x, y, overflow = frag_len_histogram(FRAG, max_len=800, step=1)
    plt.figure(); plt.plot(x, y)
    plt.xlabel("Fragment size (bp)"); plt.ylabel("Fraction")
    plt.title(f"Fragment length distribution (fallback); overflow>{800}bp: {overflow:,}")
    plt.tight_layout(); savefig_or_inline()

# 3) TSSe vs n_fragment (scatter)
scatter_plot(np.asarray(adata.obs["n_fragment"]), np.asarray(adata.obs["tsse"]),
             "TSSe vs n_fragment", "n_fragment", "TSSe", logx=True); savefig_or_inline()

# 4) TSSe vs log10(nFrags) density + threshold
nf = np.asarray(adata.obs["n_fragment"], dtype=float)
ts = np.asarray(adata.obs["tsse"], dtype=float)
log_nf = np.log10(np.clip(nf, 1, None))
plt.figure(figsize=(6,6))
plt.hist2d(log_nf, ts, bins=200, cmap="viridis"); cbar=plt.colorbar(); cbar.set_label("density")
plt.axvline(np.log10(max(MIN_COUNTS,1)), ls="--", c="k"); plt.axhline(MIN_TSSE, ls="--", c="k")
med_nf = float(np.median(nf)); med_ts = float(np.median(ts))
n_pass = int(np.sum(np.asarray(adata.obs["qc_pass"], dtype=bool)))
txt = (f"{PREFIX}\n"
       f"nCells Pass Filter = {n_pass}\n"
       f"Median Frags = {med_nf:.1f}\n"
       f"Median TSS Enrichment = {med_ts:.4f}")
plt.text(log_nf.min()+0.05, ts.max()*0.95, txt, va="top")
plt.xlabel("Log10 (Unique Fragments)"); plt.ylabel("TSS Enrichment")
plt.title("TSSe vs Log10(Unique Fragments) — density"); plt.tight_layout(); savefig_or_inline()

# 5) Various histograms
if "frip" in adata.obs: hist_plot(np.asarray(adata.obs["frip"]), "Histogram: FRiP", "FRiP"); savefig_or_inline()
hist_plot(np.asarray(adata.obs["n_fragment"]), "Histogram: n_fragment", "n_fragment"); savefig_or_inline()
hist_plot(np.asarray(adata.obs["tsse"]),      "Histogram: TSSe",       "TSSe");       savefig_or_inline()
if "frac_mito" in adata.obs: hist_plot(np.asarray(adata.obs["frac_mito"]), "Histogram: frac_mito", "frac_mito"); savefig_or_inline()
if "frac_dup"  in adata.obs: hist_plot(np.asarray(adata.obs["frac_dup"]),  "Histogram: frac_dup",  "frac_dup");  savefig_or_inline()
if "n_peaks"   in adata.obs: hist_plot(np.asarray(adata.obs["n_peaks"]),   "Histogram: n_peaks per cell", "n_peaks"); savefig_or_inline()
if "n_tiles"   in adata.obs: hist_plot(np.asarray(adata.obs["n_tiles"]),   "Histogram: n_tiles per cell", "n_tiles"); savefig_or_inline()

# 6) Association scatter
if "frip" in adata.obs:
    scatter_plot(np.asarray(adata.obs["n_fragment"]), np.asarray(adata.obs["frip"]),
                 "FRiP vs n_fragment", "n_fragment","FRiP", logx=True); savefig_or_inline()
    scatter_plot(np.asarray(adata.obs["frip"]), np.asarray(adata.obs["tsse"]),
                 "TSSe vs FRiP", "FRiP","TSSe"); savefig_or_inline()

if "frac_mito" in adata.obs:
    scatter_plot(np.asarray(adata.obs["n_fragment"]), np.asarray(adata.obs["frac_mito"]),
                 "Mito fraction vs n_fragment","n_fragment","frac_mito", logx=True); savefig_or_inline()
if "frac_dup" in adata.obs:
    scatter_plot(np.asarray(adata.obs["n_fragment"]), np.asarray(adata.obs["frac_dup"]),
                 "Duplicate fraction vs n_fragment","n_fragment","frac_dup", logx=True); savefig_or_inline()

# 7) Pass/fail (TSSe vs n_fragment)
try:
    x = np.asarray(adata.obs["n_fragment"], dtype=float)
    y = np.asarray(adata.obs["tsse"],       dtype=float)
    keep = np.asarray(adata.obs["qc_pass"], dtype=bool)
    plt.figure()
    plt.scatter(x[~keep], y[~keep], s=3, alpha=0.4, label="fail")
    plt.scatter(x[keep],  y[keep],  s=3, alpha=0.6, label="pass")
    plt.xscale("log"); plt.xlabel("n_fragment"); plt.ylabel("TSSe")
    plt.title("QC pass/fail (TSSe vs n_fragment)"); plt.legend(); plt.tight_layout(); savefig_or_inline()
except Exception as e:
    print("[WARN] pass/fail visualization failed:", e)

# 8) Threshold-line histograms
for col, thr, label in [("tsse", MIN_TSSE, "TSSe >= MIN_TSSE"),
                        ("frac_mito", MAX_MITO, "frac_mito <= MAX_MITO"),
                        ("frac_dup",  MAX_DUP,  "frac_dup <= MAX_DUP")]:
    if col in adata.obs:
        v = np.asarray(adata.obs[col], dtype=float)
        plt.figure(); plt.hist(v[~np.isnan(v)], bins=60); plt.axvline(thr)
        plt.title(f"Threshold: {label}"); plt.xlabel(col); plt.ylabel("count"); plt.tight_layout(); savefig_or_inline()

if ("frip" in adata.obs) and (PEAKS and os.path.exists(PEAKS)):
    v = np.asarray(adata.obs["frip"], dtype=float)
    plt.figure(); plt.hist(v[~np.isnan(v)], bins=60); plt.axvline(MIN_FRIP)
    plt.title("Threshold: FRiP >= MIN_FRIP"); plt.xlabel("frip"); plt.ylabel("count"); plt.tight_layout(); savefig_or_inline()

# 9) Doublet-related plots (if successful)
if dbl_key and dbl_key in adata.obs:
    v = np.asarray(adata.obs[dbl_key], dtype=float); thr = float(DOUBLETS_THRESHOLD)
    plt.figure(); plt.hist(v[~np.isnan(v)], bins=60); plt.axvline(thr, ls="--", c="r")
    plt.xlabel("doublet_score"); plt.ylabel("count"); plt.title("Histogram: doublet_score")
    plt.tight_layout(); savefig_or_inline()

    plt.figure(); plt.scatter(np.asarray(adata.obs["n_fragment"]), v, s=3, alpha=0.5)
    plt.xscale("log"); plt.xlabel("n_fragment"); plt.ylabel("doublet_score")
    plt.title("doublet_score vs n_fragment"); plt.tight_layout(); savefig_or_inline()

    plt.figure(); plt.scatter(np.asarray(adata.obs["tsse"]), v, s=3, alpha=0.5)
    plt.xlabel("TSSe"); plt.ylabel("doublet_score"); plt.title("doublet_score vs TSSe")
    plt.tight_layout(); savefig_or_inline()

    plt.figure()
    x = np.asarray(adata.obs["n_fragment"], dtype=float)
    y = np.asarray(adata.obs["tsse"], dtype=float)
    kk = np.asarray(adata.obs["doublet"], dtype=bool)
    plt.scatter(x[~kk], y[~kk], s=3, alpha=0.3, label="singlet")
    plt.scatter(x[ kk], y[ kk], s=5, alpha=0.6, label="doublet")
    plt.xscale("log"); plt.xlabel("n_fragment"); plt.ylabel("TSSe")
    plt.title("TSSe vs n_fragment (doublet overlay)"); plt.legend()
    plt.tight_layout(); savefig_or_inline()

if WRITE_PDF and pdf is not None:
    pdf.close(); print("[INFO] Unified QC PDF:", pdf_path)

# Close any open handles
close_all_handles(globals())
# ===============================================================================================================


# In[4]:


import os
os.chdir("/path/to/scATAC/brain/txci-atac")


# In[6]:


# ===================== scATAC-seq QC (SnapATAC2) — robust all-in-one with panels & doublet =====================
# For GSM7852211_mm10.* files (fragments + counts sparse matrix rows/cols)
# Features: HDF5 lock fix; import_fragments+whitelist; auto-build peaks from rows and clean before FRiP;
# optional peak matrix->n_peaks (on failure auto fallback to tile matrix->n_tiles); fraglen fallback; doublet run on tile;
# PyDataFrameElem-compatible export and filter; all figures written into a single PDF.

# ---------------- Environment fix (NFS/HDF5) ----------------
import os, sys, re, gzip, time, gc, warnings, logging
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
os.environ.setdefault("HDF5_DISABLE_VERSION_CHECK", "2")  # optional

# ---------------- Parameters (modify as needed) ----------------
FRAG        = "GSM7852211_mm10.merged.fragments.tsv.gz"
COUNTS_ROWS = "GSM7852211_mm10.counts.sparseMatrix.rows.txt.gz"  # may describe peaks (chr:start-end, etc.)
COUNTS_COLS = "GSM7852211_mm10.counts.sparseMatrix.cols.txt.gz"  # barcodes
PEAKS       = None  # if the rows above are parseable, generated automatically here; otherwise specify a BED manually
TENXH5      = None  # this dataset is not 10x H5; build whitelist from COUNTS_COLS
GENOME      = "mm10"   # choose: mm10 | mm39 | hg38 | hg19

OUTDIR = "qc_out_txci"
PREFIX = "GSM7852211_scATAC"

# Filter thresholds
MIN_COUNTS, MAX_COUNTS = 1000, 200000
MIN_TSSE,   MIN_FRIP   = 6.0, 0.20
MAX_MITO,   MAX_DUP    = 0.20, 0.80

# Optional features
DO_NUCLEOSOME       = True
NUC_MAX_LINES       = None
WRITE_PDF           = True
REBUILD_PEAK_MATRIX = True
DOUBLETS_THRESHOLD  = 0.5
SCRUBLET_JOBS       = 4

BACKEND_PREF = "hdf5"  # Raw backend priority: try hdf5 first, fall back to zarr

# ---------------- Dependencies ----------------
import numpy as np, pandas as pd, h5py
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
warnings.filterwarnings("ignore", message=".*_import_from_c.*")
try:
    import snapatac2 as snap
except Exception as e:
    print("ERROR: snapatac2 is required. Suggestion: conda create -n snap2 -y python=3.10 snapatac2 h5py pandas matplotlib", file=sys.stderr)
    raise

# ---------------- Utility functions ----------------
def ensure_dir(p): os.makedirs(p, exist_ok=True); return p

def map_genome(label: str):
    label=(label or "").strip().lower()
    if label in ["mm10","mouse","grcm38"]: return snap.genome.mm10
    if label in ["mm39","grcm39"]:         return snap.genome.mm39
    if label in ["hg38","grch38","human"]: return snap.genome.hg38
    if label in ["hg19","grch37"]:         return snap.genome.hg19
    raise ValueError(f"Unsupported genome: {label}")

def sanitize_bed(infile: str, outfile: str):
    """Ensure a standard 3-column BED; remove #/track/browser and invalid lines."""
    if infile is None or (not os.path.exists(infile)):
        raise FileNotFoundError(f"Peaks file not found: {infile}")
    op_in = gzip.open if infile.endswith(".gz") else open
    bad = kept = 0
    with op_in(infile, "rt") as fin, open(outfile, "w") as fout:
        for ln in fin:
            s=ln.strip()
            if (not s) or s.startswith("#") or s.startswith("track") or s.startswith("browser"):
                bad+=1; continue
            parts = re.split(r"\s+", s)
            if len(parts) < 3:
                bad+=1; continue
            try:
                a = int(parts[1]); b = int(parts[2])
            except Exception:
                bad+=1; continue
            if b <= a:
                bad+=1; continue
            print(f"{parts[0]}\t{a}\t{b}", file=fout); kept += 1
    print(f"[Clean BED] kept={kept:,} dropped={bad:,}  ->  {outfile}")

def try_build_peaks_from_rows(rows_path: str, out_bed: str):
    """Build BED from rows text: supports common formats such as 'chr:start-end', 'chr\tstart\tend', 'chr_start_end'"""
    if rows_path is None or (not os.path.exists(rows_path)):
        return None
    cnt = kept = 0
    with gzip.open(rows_path, "rt") if rows_path.endswith(".gz") else open(rows_path, "rt") as f, \
         open(out_bed, "w") as bed:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"): continue
            cnt += 1
            chrom=None; start=None; end=None
            # 1) chr:start-end
            m = re.match(r"^(chr\S+):(\d+)-(\d+)$", s)
            if m:
                chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))
            # 2) chr \t start \t end
            if chrom is None:
                p = re.split(r"[\t ]+", s)
                if len(p) >= 3 and p[1].isdigit() and p[2].isdigit() and p[0].startswith("chr"):
                    chrom, start, end = p[0], int(p[1]), int(p[2])
            # 3) chr_start_end
            if chrom is None:
                m2 = re.match(r"^(chr\S+)[_\-:](\d+)[_\-:](\d+)$", s)
                if m2:
                    chrom, start, end = m2.group(1), int(m2.group(2)), int(m2.group(3))
            if chrom is not None and end > start:
                bed.write(f"{chrom}\t{start}\t{end}\n"); kept += 1
    if kept > 0:
        print(f"[PEAKS] built from rows: {kept:,}/{cnt:,} -> {out_bed}")
        return out_bed
    print("[PEAKS] rows didn't look like peaks; skip.")
    return None

def read_barcodes_from_cols(cols_path: str):
    """Read barcodes from cols text (one per line)"""
    if cols_path is None or (not os.path.exists(cols_path)): return None
    bcodes=[]
    with gzip.open(cols_path, "rt") if cols_path.endswith(".gz") else open(cols_path, "rt") as f:
        for ln in f:
            s=ln.strip()
            if s: bcodes.append(s)
    print(f"[WHITELIST] loaded barcodes: {len(bcodes):,} from {cols_path}")
    return set(bcodes)

def compute_nucleosome_signal_by_streaming(frag_path, keep_barcodes=None,
                                           nfr_max=147, mono_min=147, mono_max=294, max_lines=None):
    nfr, mono = {}, {}
    i=0
    with gzip.open(frag_path, "rt") as fh:
        for line in fh:
            if max_lines is not None and i>=max_lines: break
            i+=1
            if not line or line[0]=="#": continue
            parts=line.rstrip("\n").split("\t")
            if len(parts)<4: continue
            try: L = int(parts[2]) - int(parts[1])
            except Exception: continue
            bc = parts[3]
            if keep_barcodes is not None and bc not in keep_barcodes: continue
            if L < nfr_max:               nfr[bc]  = nfr.get(bc,0)  + 1
            elif mono_min <= L < mono_max: mono[bc] = mono.get(bc,0) + 1
    ks = (keep_barcodes if keep_barcodes is not None else set(list(nfr.keys())+list(mono.keys())))
    return {bc: (mono.get(bc,0)/(nfr.get(bc,0) if nfr.get(bc,0)>0 else 1.0)) for bc in ks}

def barcode_knee_plot(counts, title):
    vals = np.sort(np.asarray(counts))[::-1]; ranks = np.arange(1,len(vals)+1)
    plt.figure(); plt.plot(ranks, vals); plt.xscale("log"); plt.yscale("log")
    plt.xlabel("barcode rank (log10)"); plt.ylabel("n_fragment (log10)"); plt.title(title); plt.tight_layout()

def hist_plot(x, title, xlabel):
    v=np.asarray(x, dtype=float)
    plt.figure(); plt.hist(v[~np.isnan(v)], bins=60)
    plt.title(title); plt.xlabel(xlabel); plt.ylabel("count"); plt.tight_layout()

def scatter_plot(x,y,title,xlabel,ylabel,logx=False,logy=False):
    plt.figure(); plt.scatter(x,y,s=3,alpha=0.6)
    if logx: plt.xscale("log")
    if logy: plt.yscale("log")
    plt.title(title); plt.xlabel(xlabel); plt.ylabel(ylabel); plt.tight_layout()

def close_all_handles(globs: dict):
    for k,v in list(globs.items()):
        if hasattr(v, "close"):
            try: v.close()
            except Exception: pass
    gc.collect()

def unique_path(base_dir, base_name, backend):
    ts = time.strftime("%Y%m%d-%H%M%S"); pid = os.getpid()
    return os.path.join(base_dir, f"{base_name}.raw.{ts}.{pid}.{('h5ad' if backend=='hdf5' else 'zarr')}")

def _safe_row_nnz(X):
    try:
        import scipy.sparse as sp
        if sp.issparse(X): return np.asarray(X.getnnz(axis=1)).ravel()
        Xn = np.asarray(X)
        if Xn.ndim == 2: return np.count_nonzero(Xn, axis=1)
    except Exception:
        pass
    return None

def frag_len_histogram(frag_path, max_len=800, step=1):
    """Streaming fallback: returns x(1..max_len), y(percent), overflow(count of >max_len)"""
    bins = np.zeros(max_len + 1, dtype=np.int64)  # index 0 collects overflow
    with gzip.open(frag_path, "rt") as fh:
        for ln in fh:
            if not ln or ln[0] == "#": continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 3: continue
            try: L = int(p[2]) - int(p[1])
            except Exception: continue
            if 0 < L <= max_len: bins[L] += 1
            else: bins[0] += 1
    idx = np.arange(1, max_len+1, step)
    val = bins[1:max_len+1:step].astype(float)
    val = val/val.sum() if val.sum()>0 else val
    return idx, val, int(bins[0])

# ---------------- Main flow: import + basic QC ----------------
ensure_dir(OUTDIR)
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
genome_obj = map_genome(GENOME)

# 1) barcodes (whitelist) from cols
barcodes_keep = read_barcodes_from_cols(COUNTS_COLS)
if barcodes_keep is not None:
    print(f"[INFO] Built whitelist from counts/cols: {len(barcodes_keep):,} entries")

# 2) peaks from rows (build BED if possible)
built_peaks = None
if PEAKS is None and COUNTS_ROWS and os.path.exists(COUNTS_ROWS):
    built_peaks = try_build_peaks_from_rows(COUNTS_ROWS, os.path.join(OUTDIR, f"{PREFIX}.peaks.from_rows.bed"))
    if built_peaks: 
        # Clean again to ensure 3 columns and valid coordinates
        clean_bed = os.path.join(OUTDIR, f"{PREFIX}.peaks.from_rows.clean.bed")
        sanitize_bed(built_peaks, clean_bed)
        PEAKS = clean_bed

# Import fragments (hdf5->zarr fallback; avoid file locking)
close_all_handles(globals())
backend   = BACKEND_PREF
raw_path  = unique_path(OUTDIR, PREFIX, backend)

print("[INFO] Importing fragments and computing n_fragment / frac_dup / frac_mito ...")
try:
    adata = snap.pp.import_fragments(
        FRAG, chrom_sizes=genome_obj, file=raw_path,
        sorted_by_barcode=False, whitelist=list(barcodes_keep) if barcodes_keep is not None else None,
        backend=("hdf5" if backend=="hdf5" else "zarr")
    )
except RuntimeError as e:
    print(f"[WARN] HDF5 import failed ({e}); retrying with Zarr backend.")
    backend  = "zarr"
    raw_path = unique_path(OUTDIR, PREFIX, backend)
    adata    = snap.pp.import_fragments(
        FRAG, chrom_sizes=genome_obj, file=raw_path,
        sorted_by_barcode=False, whitelist=list(barcodes_keep) if barcodes_keep is not None else None,
        backend="zarr"
    )

print("[INFO] Computing TSSe ...")
snap.metrics.tsse(adata, genome_obj)
print("[INFO] TSSe done")

# ---------------- FRiP (if PEAKS available) ----------------
frip_added = False
PEAKS_CLEAN = None
if PEAKS and os.path.exists(PEAKS):
    # Clean once more defensively
    PEAKS_CLEAN = os.path.join(OUTDIR, f"{PREFIX}.peaks.clean.bed")
    sanitize_bed(PEAKS, PEAKS_CLEAN)
    print("[INFO] Computing FRiP ...")
    try:
        snap.metrics.frip(adata, {"frip": PEAKS_CLEAN})     # new signature
        frip_added = True
    except TypeError:
        snap.metrics.frip(adata, PEAKS_CLEAN)               # old signature
        if "frip" not in adata.obs:
            for c in ["FRiP","Frip","FRIP"]:
                if c in adata.obs: adata.obs["frip"]=adata.obs[c]; break
        frip_added = True
    except Exception as e:
        print(f"[WARN] FRiP computation failed: {e}")
print("[INFO] FRiP done" if frip_added else "[INFO] FRiP not computed (no peaks or failed)")

# ---------------- peak matrix -> n_peaks (fall back to tile on failure) ----------------
peakA = None
if REBUILD_PEAK_MATRIX and PEAKS_CLEAN and os.path.exists(PEAKS_CLEAN):
    try:
        peak_matrix_path = os.path.join(OUTDIR, f"{PREFIX}.peak_matrix.{('h5ad' if backend=='hdf5' else 'zarr')}")
        peakA = snap.pp.make_peak_matrix(adata, file=peak_matrix_path, peak_file=PEAKS_CLEAN)
        nnz = _safe_row_nnz(peakA.X)
        if nnz is not None:
            adata.obs["n_peaks"] = nnz
            print("[INFO] Wrote adata.obs['n_peaks']")
    except Exception as e:
        print(f"[WARN] Failed to build peak matrix: {e}")

tileA = None
if "n_peaks" not in adata.obs:
    try:
        print("[INFO] Building tile(5kb) matrix as fallback ...")
        tile_path = os.path.join(OUTDIR, f"{PREFIX}.tile_5kb.zarr")
        tileA = snap.pp.make_tile_matrix(adata, file=tile_path, bin_size=5000)
        X = tileA.X
        import scipy.sparse as sp
        n_tiles = (X > 0).sum(axis=1).A1 if sp.issparse(X) else np.count_nonzero(np.asarray(X), axis=1)
        adata.obs["n_tiles"] = np.asarray(n_tiles, dtype=float)
        print("[INFO] Wrote adata.obs['n_tiles'] (tile-based)")
    except Exception as e:
        print(f"[WARN] Tile matrix failed: {e}")

# ---------------- nucleosome signal (optional) ----------------
if DO_NUCLEOSOME:
    print("[INFO] Computing nucleosome signal (mono/NFR); may be slow ...")
    bc_set = set(adata.obs_names)
    nuc_signal = compute_nucleosome_signal_by_streaming(FRAG, keep_barcodes=bc_set, max_lines=NUC_MAX_LINES)
    adata.obs["nucleosome_signal"] = pd.Series(nuc_signal).reindex(adata.obs_names).astype(float)
    print("[INFO] nucleosome_signal added to adata.obs")

# ---------------- Export metrics (PyDataFrameElem-compatible) ----------------
ensure_dir(OUTDIR)
metrics_csv = os.path.join(OUTDIR, f"{PREFIX}.qc_metrics.csv")
wanted = ["n_fragment","frac_dup","frac_mito","tsse","frip","n_peaks","n_tiles","nucleosome_signal"]
avail  = [k for k in wanted if k in adata.obs]
df = pd.DataFrame(index=adata.obs_names)
for k in avail: df[k] = np.asarray(adata.obs[k])
df.to_csv(metrics_csv); print("[INFO] Metrics CSV:", metrics_csv)

# ---------------- Filter & subset export ----------------
n_fragment = np.asarray(adata.obs["n_fragment"], dtype=float)
tsse       = np.asarray(adata.obs["tsse"],       dtype=float)
frac_mito  = np.asarray(adata.obs["frac_mito"],  dtype=float)
frac_dup   = np.asarray(adata.obs["frac_dup"],   dtype=float)
mask = (
    (n_fragment >= MIN_COUNTS) &
    (n_fragment <= MAX_COUNTS) &
    (tsse       >= MIN_TSSE)   &
    (frac_mito  <= MAX_MITO)   &
    (frac_dup   <= MAX_DUP)
)
if ("frip" in adata.obs) and (PEAKS_CLEAN and os.path.exists(PEAKS_CLEAN)):
    frip = np.asarray(adata.obs["frip"], dtype=float)
    mask = mask & (frip >= MIN_FRIP)
adata.obs["qc_pass"] = mask

idx = np.where(mask)[0]
obs_names_list = list(adata.obs_names)
pass_barcodes  = [obs_names_list[i] for i in idx]

filtered_h5ad = os.path.join(OUTDIR, f"{PREFIX}.filtered.h5ad")
try:
    _ = adata.subset(obs_indices=idx, inplace=False, out=filtered_h5ad, backend="hdf5")
    print("[INFO] Filtered h5ad:", filtered_h5ad, "(subset/out)")
except Exception as e:
    print(f"[WARN] subset(out=...) failed; re-importing using whitelist: {e}")
    tmp = os.path.join(OUTDIR, f"{PREFIX}.filtered.tmp.h5ad")
    adata_f = snap.pp.import_fragments(
        FRAG, chrom_sizes=genome_obj, file=tmp, sorted_by_barcode=False, whitelist=pass_barcodes
    )
    snap.metrics.tsse(adata_f, genome_obj)
    if PEAKS_CLEAN and os.path.exists(PEAKS_CLEAN):
        try:    snap.metrics.frip(adata_f, {"frip": PEAKS_CLEAN})
        except: 
            try: snap.metrics.frip(adata_f, PEAKS_CLEAN)
            except: pass
    os.replace(tmp, filtered_h5ad)
    if hasattr(adata_f, "close"):
        try: adata_f.close()
        except Exception: pass
    print("[INFO] Filtered h5ad:", filtered_h5ad, "(re-import whitelist)")

filtered_barcodes_txt = os.path.join(OUTDIR, f"{PREFIX}.filtered_barcodes.txt")
pd.Series(pass_barcodes).to_csv(filtered_barcodes_txt, index=False, header=False)
print("[INFO] Passing barcodes:", filtered_barcodes_txt)

# ---------------- Scrublet doublet (run on tile or peak) ----------------
dbl_key = None
srcA_for_dbl = None
if tileA is not None:
    srcA_for_dbl = tileA
elif peakA is not None:
    srcA_for_dbl = peakA

if srcA_for_dbl is None:
    try:
        print("[INFO] Temporarily building tile matrix for doublet ...")
        tmp_tile = os.path.join(OUTDIR, f"{PREFIX}.tile_5kb.doublet.zarr")
        srcA_for_dbl = snap.pp.make_tile_matrix(adata, file=tmp_tile, bin_size=5000)
    except Exception as e:
        print(f"[WARN] Failed to build doublet matrix: {e}")

if srcA_for_dbl is not None:
    try:
        snap.pp.scrublet(srcA_for_dbl, n_jobs=SCRUBLET_JOBS)  # writes srcA_for_dbl.obs['doublet_score']
        s = pd.Series(np.asarray(srcA_for_dbl.obs["doublet_score"]), index=list(srcA_for_dbl.obs_names), name="doublet_score")
        adata.obs["doublet_score"] = s.reindex(adata.obs_names).astype(float).values
        adata.obs["doublet"] = np.asarray(adata.obs["doublet_score"], dtype=float) > float(DOUBLETS_THRESHOLD)
        dbl_key = "doublet_score"
        print("[INFO] doublet_score written; threshold >", float(DOUBLETS_THRESHOLD), " marked as doublet")
    except Exception as e:
        print(f"[WARN] scrublet failed: {e}")
else:
    print("[INFO] No doublet input matrix built (skipping scrublet)")

# ---------------- Unified PDF report (all pages) ----------------
pdf_path = os.path.join(OUTDIR, f"{PREFIX}.qc_report.pdf") if WRITE_PDF else None
pdf = PdfPages(pdf_path) if WRITE_PDF else None
def savefig_or_inline():
    if pdf is not None: pdf.savefig(); plt.close()
    else: plt.show()

# 1) Knee
barcode_knee_plot(adata.obs["n_fragment"].to_numpy(), "Barcode knee (n_fragment)"); savefig_or_inline()

# 2) Fragment length distribution (official or fallback)
ok_fsd = True
try:
    fsd = snap.metrics.frag_size_distr(adata)
    ok_fsd = fsd is not None
except Exception:
    ok_fsd = False
if ok_fsd:
    x = np.arange(len(fsd)); y = np.array(fsd, dtype=float); x = x[1:]; y = y[1:]
    plt.figure(); plt.plot(x, y); plt.yscale("log")
    plt.xlabel("fragment length (bp)"); plt.ylabel("count (log)")
    plt.title("Fragment length distribution"); plt.tight_layout(); savefig_or_inline()
    # Percentage version
    y_pct = y / y.sum() if y.sum()>0 else y
    plt.figure(); plt.plot(x, y_pct)
    plt.xlabel("Size of Fragments (bp)"); plt.ylabel("Fragments (%)")
    plt.title("Fragment Size Distribution (percent)"); plt.tight_layout(); savefig_or_inline()
else:
    print("[INFO] Official frag_size_distr unavailable; using fallback streaming version")
    x, y, overflow = frag_len_histogram(FRAG, max_len=800, step=1)
    plt.figure(); plt.plot(x, y)
    plt.xlabel("Fragment size (bp)"); plt.ylabel("Fraction")
    plt.title(f"Fragment length distribution (fallback); overflow>{800}bp: {overflow:,}")
    plt.tight_layout(); savefig_or_inline()

# 3) TSSe vs n_fragment (scatter)
scatter_plot(np.asarray(adata.obs["n_fragment"]), np.asarray(adata.obs["tsse"]),
             "TSSe vs n_fragment", "n_fragment", "TSSe", logx=True); savefig_or_inline()

# 4) TSSe vs log10(nFrags) density + threshold
nf = np.asarray(adata.obs["n_fragment"], dtype=float)
ts = np.asarray(adata.obs["tsse"], dtype=float)
log_nf = np.log10(np.clip(nf, 1, None))
plt.figure(figsize=(6,6))
plt.hist2d(log_nf, ts, bins=200, cmap="viridis"); cbar=plt.colorbar(); cbar.set_label("density")
plt.axvline(np.log10(max(MIN_COUNTS,1)), ls="--", c="k"); plt.axhline(MIN_TSSE, ls="--", c="k")
med_nf = float(np.median(nf)); med_ts = float(np.median(ts))
n_pass = int(np.sum(np.asarray(adata.obs["qc_pass"], dtype=bool)))
txt = (f"{PREFIX}\n"
       f"nCells Pass Filter = {n_pass}\n"
       f"Median Frags = {med_nf:.1f}\n"
       f"Median TSS Enrichment = {med_ts:.4f}")
plt.text(log_nf.min()+0.05, ts.max()*0.95, txt, va="top")
plt.xlabel("Log10 (Unique Fragments)"); plt.ylabel("TSS Enrichment")
plt.title("TSSe vs Log10(Unique Fragments) — density"); plt.tight_layout(); savefig_or_inline()

# 5) Various histograms
if "frip" in adata.obs: hist_plot(np.asarray(adata.obs["frip"]), "Histogram: FRiP", "FRiP"); savefig_or_inline()
hist_plot(np.asarray(adata.obs["n_fragment"]), "Histogram: n_fragment", "n_fragment"); savefig_or_inline()
hist_plot(np.asarray(adata.obs["tsse"]),      "Histogram: TSSe",       "TSSe");       savefig_or_inline()
if "frac_mito" in adata.obs: hist_plot(np.asarray(adata.obs["frac_mito"]), "Histogram: frac_mito", "frac_mito"); savefig_or_inline()
if "frac_dup"  in adata.obs: hist_plot(np.asarray(adata.obs["frac_dup"]),  "Histogram: frac_dup",  "frac_dup");  savefig_or_inline()
if "n_peaks"   in adata.obs: hist_plot(np.asarray(adata.obs["n_peaks"]),   "Histogram: n_peaks per cell", "n_peaks"); savefig_or_inline()
if "n_tiles"   in adata.obs: hist_plot(np.asarray(adata.obs["n_tiles"]),   "Histogram: n_tiles per cell", "n_tiles"); savefig_or_inline()

# 6) Association scatter
if "frip" in adata.obs:
    scatter_plot(np.asarray(adata.obs["n_fragment"]), np.asarray(adata.obs["frip"]),
                 "FRiP vs n_fragment", "n_fragment","FRiP", logx=True); savefig_or_inline()
    scatter_plot(np.asarray(adata.obs["frip"]), np.asarray(adata.obs["tsse"]),
                 "TSSe vs FRiP", "FRiP","TSSe"); savefig_or_inline()

if "frac_mito" in adata.obs:
    scatter_plot(np.asarray(adata.obs["n_fragment"]), np.asarray(adata.obs["frac_mito"]),
                 "Mito fraction vs n_fragment","n_fragment","frac_mito", logx=True); savefig_or_inline()
if "frac_dup" in adata.obs:
    scatter_plot(np.asarray(adata.obs["n_fragment"]), np.asarray(adata.obs["frac_dup"]),
                 "Duplicate fraction vs n_fragment","n_fragment","frac_dup", logx=True); savefig_or_inline()

# 7) Pass/fail (TSSe vs n_fragment)
try:
    x = np.asarray(adata.obs["n_fragment"], dtype=float)
    y = np.asarray(adata.obs["tsse"],       dtype=float)
    keep = np.asarray(adata.obs["qc_pass"], dtype=bool)
    plt.figure()
    plt.scatter(x[~keep], y[~keep], s=3, alpha=0.4, label="fail")
    plt.scatter(x[keep],  y[keep],  s=3, alpha=0.6, label="pass")
    plt.xscale("log"); plt.xlabel("n_fragment"); plt.ylabel("TSSe")
    plt.title("QC pass/fail (TSSe vs n_fragment)"); plt.legend(); plt.tight_layout(); savefig_or_inline()
except Exception as e:
    print("[WARN] pass/fail visualization failed:", e)

# 8) Threshold-line histograms
for col, thr, label in [("tsse", MIN_TSSE, "TSSe >= MIN_TSSE"),
                        ("frac_mito", MAX_MITO, "frac_mito <= MAX_MITO"),
                        ("frac_dup",  MAX_DUP,  "frac_dup <= MAX_DUP")]:
    if col in adata.obs:
        v = np.asarray(adata.obs[col], dtype=float)
        plt.figure(); plt.hist(v[~np.isnan(v)], bins=60); plt.axvline(thr)
        plt.title(f"Threshold: {label}"); plt.xlabel(col); plt.ylabel("count"); plt.tight_layout(); savefig_or_inline()

if ("frip" in adata.obs) and PEAKS_CLEAN and os.path.exists(PEAKS_CLEAN):
    v = np.asarray(adata.obs["frip"], dtype=float)
    plt.figure(); plt.hist(v[~np.isnan(v)], bins=60); plt.axvline(MIN_FRIP)
    plt.title("Threshold: FRiP >= MIN_FRIP"); plt.xlabel("frip"); plt.ylabel("count"); plt.tight_layout(); savefig_or_inline()

# 9) Doublet-related plots (if successful)
if dbl_key and dbl_key in adata.obs:
    v = np.asarray(adata.obs[dbl_key], dtype=float); thr = float(DOUBLETS_THRESHOLD)
    plt.figure(); plt.hist(v[~np.isnan(v)], bins=60); plt.axvline(thr, ls="--", c="r")
    plt.xlabel("doublet_score"); plt.ylabel("count"); plt.title("Histogram: doublet_score")
    plt.tight_layout(); savefig_or_inline()

    plt.figure(); plt.scatter(np.asarray(adata.obs["n_fragment"]), v, s=3, alpha=0.5)
    plt.xscale("log"); plt.xlabel("n_fragment"); plt.ylabel("doublet_score")
    plt.title("doublet_score vs n_fragment"); plt.tight_layout(); savefig_or_inline()

    plt.figure(); plt.scatter(np.asarray(adata.obs["tsse"]), v, s=3, alpha=0.5)
    plt.xlabel("TSSe"); plt.ylabel("doublet_score"); plt.title("doublet_score vs TSSe")
    plt.tight_layout(); savefig_or_inline()

    plt.figure()
    x = np.asarray(adata.obs["n_fragment"], dtype=float)
    y = np.asarray(adata.obs["tsse"], dtype=float)
    kk = np.asarray(adata.obs["doublet"], dtype=bool)
    plt.scatter(x[~kk], y[~kk], s=3, alpha=0.3, label="singlet")
    plt.scatter(x[ kk], y[ kk], s=5, alpha=0.6, label="doublet")
    plt.xscale("log"); plt.xlabel("n_fragment"); plt.ylabel("TSSe")
    plt.title("TSSe vs n_fragment (doublet overlay)"); plt.legend()
    plt.tight_layout(); savefig_or_inline()

if WRITE_PDF and pdf is not None:
    pdf.close(); print("[INFO] Unified QC PDF:", pdf_path)

# Close any open handles
for obj in [locals().get('adata', None), locals().get('peakA', None), locals().get('tileA', None)]:
    if hasattr(obj, "close"):
        try: obj.close()
        except Exception: pass
# ===============================================================================================================


# In[8]:


# ===================== scATAC-seq QC (SnapATAC2) — robust all-in-one (txci-atac) =====================
# Temp/Cache/Outputs => /path/to/tmp
# Inputs in CWD:
#   - GSM7852211_mm10.merged.fragments.tsv.gz  (+ .tbi or .tbi.gz)
#   - GSM7852211_mm10.counts.sparseMatrix.rows.txt.gz  (peak-like rows)
#   - GSM7852211_mm10.counts.sparseMatrix.cols.txt.gz  (barcodes)

# ---------------- Env & params ----------------
import os, sys, re, gzip, time, gc, warnings, logging
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
os.environ.setdefault("HDF5_DISABLE_VERSION_CHECK", "2")  # optional

# <<< CHANGE HERE IF NEEDED >>>
FRAG        = "GSM7852211_mm10.merged.fragments.tsv.gz"
COUNTS_ROWS = "GSM7852211_mm10.counts.sparseMatrix.rows.txt.gz"
COUNTS_COLS = "GSM7852211_mm10.counts.sparseMatrix.cols.txt.gz"
GENOME      = "mm10"

# Force temp/cache/outputs to big scratch
BIG_TMP   = "/path/to/tmp"                 # temp + caches
OUTDIR    = "/path/to/txci_atac_qc"    # outputs
PREFIX    = "GSM7852211_scATAC"

# QC thresholds
MIN_COUNTS, MAX_COUNTS = 1000, 200000
MIN_TSSE,   MIN_FRIP   = 6.0,   0.20
MAX_MITO,   MAX_DUP    = 0.20,  0.80

# Options / tuning
DO_NUCLEOSOME       = True
NUC_MAX_LINES       = None
WRITE_PDF           = True
REBUILD_PEAK_MATRIX = False      # keep False to save disk; tile fallback enables doublet/complexity
DOUBLETS_THRESHOLD  = 0.5
SCRUBLET_JOBS       = 4

# Prefer Zarr; modest chunks; early min-frags filter
BACKEND_PREF       = "zarr"
CHUNK_SIZE         = 1_000_000
N_JOBS             = 4
MIN_FRAGS_IMPORT   = max(500, MIN_COUNTS)

# ---------------- Deps ----------------
import numpy as np, pandas as pd, h5py
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
warnings.filterwarnings("ignore", message=".*_import_from_c.*")
try:
    import snapatac2 as snap
except Exception as e:
    print("ERROR: snapatac2 not installed. Try: conda create -n snap2 -y python=3.10 snapatac2 h5py pandas matplotlib", file=sys.stderr)
    raise

# ---------------- Helpers ----------------
from pathlib import Path
def ensure_dir(p): Path(p).mkdir(parents=True, exist_ok=True); return p

def map_genome(label: str):
    ll=(label or "").strip().lower()
    if ll in ["mm10","mouse","grcm38"]: return snap.genome.mm10
    if ll in ["mm39","grcm39"]:         return snap.genome.mm39
    if ll in ["hg38","grch38","human"]: return snap.genome.hg38
    if ll in ["hg19","grch37"]:         return snap.genome.hg19
    raise ValueError(f"Unsupported genome: {label}")

def set_all_tmp(tmpdir: str):
    os.environ["TMPDIR"] = tmpdir
    os.environ["TMP"]    = tmpdir
    os.environ["TEMP"]   = tmpdir
    os.environ.setdefault("MPLCONFIGDIR", os.path.join(tmpdir, "mpl"))

def ensure_tbi_for_frag(frag_path: str):
    """Ensure .tbi exists; if only .tbi.gz exists, decompress to .tbi."""
    tbi = frag_path + ".tbi"
    tbi_gz = tbi + ".gz"
    if os.path.exists(tbi):
        return tbi
    if os.path.exists(tbi_gz):
        print(f"[INFO] Decompressing {tbi_gz} -> {tbi}")
        with gzip.open(tbi_gz, "rb") as fi, open(tbi, "wb") as fo:
            fo.write(fi.read())
        return tbi
    print(f"[WARN] No .tbi found for {frag_path} (SnapATAC2 may index to TMP).")
    return None

def sanitize_bed(infile: str, outfile: str):
    """Ensure 3-col BED, strip headers/bad lines."""
    if infile is None or (not os.path.exists(infile)): raise FileNotFoundError(infile)
    op_in = gzip.open if infile.endswith(".gz") else open
    bad=kept=0
    with op_in(infile, "rt") as fin, open(outfile, "w") as fo:
        for ln in fin:
            s=ln.strip()
            if (not s) or s.startswith("#") or s.startswith("track") or s.startswith("browser"):
                bad+=1; continue
            p=re.split(r"\s+", s)
            if len(p)<3: bad+=1; continue
            try: a=int(p[1]); b=int(p[2])
            except: bad+=1; continue
            if b<=a: bad+=1; continue
            fo.write(f"{p[0]}\t{a}\t{b}\n"); kept+=1
    print(f"[Clean BED] kept={kept:,} dropped={bad:,} -> {outfile}")
    return outfile

def try_build_peaks_from_rows(rows_path: str, out_bed: str):
    """Build peaks BED from rows: supports 'chr:start-end', 'chr start end', 'chr_start_end'."""
    if rows_path is None or (not os.path.exists(rows_path)): return None
    cnt=kept=0
    op = gzip.open if rows_path.endswith(".gz") else open
    with op(rows_path, "rt") as f, open(out_bed, "w") as bed:
        for line in f:
            s=line.strip()
            if not s or s.startswith("#"): continue
            cnt+=1
            chrom=start=end=None
            m=re.match(r"^(chr\S+):(\d+)-(\d+)$", s)
            if m: chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))
            if chrom is None:
                p=re.split(r"[\t ]+", s)
                if len(p)>=3 and p[1].isdigit() and p[2].isdigit() and p[0].startswith("chr"):
                    chrom, start, end = p[0], int(p[1]), int(p[2])
            if chrom is None:
                m2=re.match(r"^(chr\S+)[_\-:](\d+)[_\-:](\d+)$", s)
                if m2: chrom, start, end = m2.group(1), int(m2.group(2)), int(m2.group(3))
            if chrom is not None and end>start:
                bed.write(f"{chrom}\t{start}\t{end}\n"); kept+=1
    if kept>0:
        print(f"[PEAKS] built from rows: {kept:,}/{cnt:,} -> {out_bed}")
        return out_bed
    print("[PEAKS] rows not peak-like; skipping.")
    return None

def read_barcodes_from_cols(cols_path: str):
    if cols_path is None or (not os.path.exists(cols_path)): return None
    bcodes=[]
    op=gzip.open if cols_path.endswith(".gz") else open
    with op(cols_path, "rt") as f:
        for ln in f:
            s=ln.strip()
            if s: bcodes.append(s)
    print(f"[WHITELIST] loaded barcodes: {len(bcodes):,} from {cols_path}")
    return set(bcodes)

def compute_nucleosome_signal_by_streaming(frag_path, keep_barcodes=None,
                                           nfr_max=147, mono_min=147, mono_max=294, max_lines=None):
    nfr, mono = {}, {}
    i=0
    with gzip.open(frag_path, "rt") as fh:
        for line in fh:
            if max_lines is not None and i>=max_lines: break
            i+=1
            if not line or line[0]=="#": continue
            p=line.rstrip("\n").split("\t")
            if len(p)<4: continue
            try: L = int(p[2]) - int(p[1])
            except: continue
            bc=p[3]
            if keep_barcodes is not None and bc not in keep_barcodes: continue
            if L < nfr_max:               nfr[bc]  = nfr.get(bc,0)  + 1
            elif mono_min <= L < mono_max: mono[bc] = mono.get(bc,0) + 1
    ks = (keep_barcodes if keep_barcodes is not None else set(list(nfr.keys())+list(mono.keys())))
    return {bc: (mono.get(bc,0)/(nfr.get(bc,0) if nfr.get(bc,0)>0 else 1.0)) for bc in ks}

def barcode_knee_plot(counts, title):
    vals=np.sort(np.asarray(counts))[::-1]; ranks=np.arange(1,len(vals)+1)
    plt.figure(); plt.plot(ranks, vals); plt.xscale("log"); plt.yscale("log")
    plt.xlabel("barcode rank (log10)"); plt.ylabel("n_fragment (log10)")
    plt.title(title); plt.tight_layout()

def hist_plot(x, title, xlabel):
    v=np.asarray(x, dtype=float)
    plt.figure(); plt.hist(v[~np.isnan(v)], bins=60)
    plt.title(title); plt.xlabel(xlabel); plt.ylabel("count"); plt.tight_layout()

def scatter_plot(x,y,title,xlabel,ylabel,logx=False,logy=False):
    plt.figure(); plt.scatter(x,y,s=3,alpha=0.6)
    if logx: plt.xscale("log")
    if logy: plt.yscale("log")
    plt.title(title); plt.xlabel(xlabel); plt.ylabel(ylabel); plt.tight_layout()

def close_all_handles(globs: dict):
    for k,v in list(globs.items()):
        if hasattr(v,"close"):
            try: v.close()
            except Exception: pass
    gc.collect()

def unique_path(base_dir, base_name, backend):
    ts=time.strftime("%Y%m%d-%H%M%S"); pid=os.getpid()
    return os.path.join(base_dir, f"{base_name}.raw.{ts}.{pid}.{('h5ad' if backend=='hdf5' else 'zarr')}")

def _safe_row_nnz(X):
    try:
        import scipy.sparse as sp
        if sp.issparse(X): return np.asarray(X.getnnz(axis=1)).ravel()
        Xn=np.asarray(X)
        if Xn.ndim==2: return np.count_nonzero(Xn, axis=1)
    except Exception:
        pass
    return None

def frag_len_histogram(frag_path, max_len=800, step=1):
    bins=np.zeros(max_len+1, dtype=np.int64)  # 0 collects overflow
    with gzip.open(frag_path,"rt") as fh:
        for ln in fh:
            if not ln or ln[0]=="#": continue
            p=ln.rstrip("\n").split("\t")
            if len(p)<3: continue
            try: L=int(p[2])-int(p[1])
            except: continue
            if 0 < L <= max_len: bins[L]+=1
            else: bins[0]+=1
    idx=np.arange(1,max_len+1,step)
    val=bins[1:max_len+1:step].astype(float)
    val=val/val.sum() if val.sum()>0 else val
    return idx, val, int(bins[0])

# ---------------- Setup scratch & genome ----------------
ensure_dir(BIG_TMP); ensure_dir(OUTDIR); set_all_tmp(BIG_TMP)
genome_obj = map_genome(GENOME)

# Ensure .tbi present
ensure_tbi_for_frag(FRAG)

# Load whitelist from COLS
barcodes_keep = read_barcodes_from_cols(COUNTS_COLS)

# Build peaks from ROWS -> clean BED
PEAKS = None
if COUNTS_ROWS and os.path.exists(COUNTS_ROWS):
    built = try_build_peaks_from_rows(COUNTS_ROWS, os.path.join(OUTDIR, f"{PREFIX}.peaks.from_rows.bed"))
    if built:
        PEAKS = sanitize_bed(built, os.path.join(OUTDIR, f"{PREFIX}.peaks.from_rows.clean.bed"))

# ---------------- Import fragments (Zarr on scratch) ----------------
close_all_handles(globals())
backend   = BACKEND_PREF
raw_path  = unique_path(OUTDIR, PREFIX, backend)

print("[INFO] Importing fragments and computing n_fragment / frac_dup / frac_mito ...")
try:
    adata = snap.pp.import_fragments(
        FRAG,
        chrom_sizes=genome_obj,
        file=raw_path,
        sorted_by_barcode=False,
        whitelist=list(barcodes_keep) if barcodes_keep is not None else None,
        backend=("hdf5" if backend=="hdf5" else "zarr"),
        tempdir=BIG_TMP,
        chunk_size=CHUNK_SIZE,
        n_jobs=N_JOBS,
        min_num_fragments=MIN_FRAGS_IMPORT,
    )
except Exception as e:
    print(f"[WARN] import_fragments failed ({e}). Switching backend/file and retry …")
    backend = "hdf5" if backend=="zarr" else "zarr"
    raw_path = unique_path(OUTDIR, PREFIX, backend)
    adata = snap.pp.import_fragments(
        FRAG,
        chrom_sizes=genome_obj,
        file=raw_path,
        sorted_by_barcode=False,
        whitelist=list(barcodes_keep) if barcodes_keep is not None else None,
        backend=("hdf5" if backend=="hdf5" else "zarr"),
        tempdir=BIG_TMP,
        chunk_size=CHUNK_SIZE,
        n_jobs=N_JOBS,
        min_num_fragments=MIN_FRAGS_IMPORT,
    )
print("[INFO] Raw store:", raw_path)

# ---------------- TSSe with cache on scratch (fix NotFound) ----------------
from pathlib import Path as _Path
SNAP_CACHE = os.path.join(BIG_TMP, "snapatac2_cache")
_Path(SNAP_CACHE).mkdir(parents=True, exist_ok=True)
os.environ["XDG_CACHE_HOME"]      = SNAP_CACHE
os.environ["SNAPATAC2_CACHE_DIR"] = SNAP_CACHE

print("[INFO] Computing TSSe ...")
snap.metrics.tsse(adata, genome_obj)   # or: snap.metrics.tsse(adata, "/path/to/local.gtf/gff")
print("[INFO] TSSe done")

# ---------------- FRiP (if peaks available) ----------------
frip_added=False; PEAKS_CLEAN=None
if PEAKS and os.path.exists(PEAKS):
    PEAKS_CLEAN = os.path.join(OUTDIR, f"{PREFIX}.peaks.clean.bed")
    sanitize_bed(PEAKS, PEAKS_CLEAN)
    print("[INFO] Computing FRiP ...")
    try:
        snap.metrics.frip(adata, {"frip": PEAKS_CLEAN})
        frip_added=True
    except TypeError:
        snap.metrics.frip(adata, PEAKS_CLEAN)
        if "frip" not in adata.obs:
            for c in ["FRiP","Frip","FRIP"]:
                if c in adata.obs: adata.obs["frip"]=adata.obs[c]; break
        frip_added=True
    except Exception as e:
        print(f"[WARN] FRiP computation failed: {e}")
print("[INFO] FRiP done" if frip_added else "[INFO] FRiP not computed (no peaks or failed)")

# ---------------- peak matrix (optional; disabled to save disk) ----------------
peakA=None
if REBUILD_PEAK_MATRIX and PEAKS_CLEAN and os.path.exists(PEAKS_CLEAN):
    try:
        peak_matrix_path=os.path.join(OUTDIR, f"{PREFIX}.peak_matrix.{('h5ad' if backend=='hdf5' else 'zarr')}")
        peakA=snap.pp.make_peak_matrix(adata, file=peak_matrix_path, peak_file=PEAKS_CLEAN)
        nnz=_safe_row_nnz(peakA.X)
        if nnz is not None:
            adata.obs["n_peaks"]=nnz
            print("[INFO] n_peaks added")
    except Exception as e:
        print(f"[WARN] peak matrix failed: {e}")

# ---------------- tile matrix fallback (provides n_tiles & doublet input) ----------------
tileA=None
if "n_peaks" not in adata.obs:
    try:
        print("[INFO] Building 5kb tile matrix (fallback) …")
        tile_path=os.path.join(OUTDIR, f"{PREFIX}.tile_5kb.zarr")
        tileA=snap.pp.make_tile_matrix(adata, file=tile_path, bin_size=5000)
        X=tileA.X
        import scipy.sparse as sp
        n_tiles=(X>0).sum(axis=1).A1 if sp.issparse(X) else np.count_nonzero(np.asarray(X), axis=1)
        adata.obs["n_tiles"]=np.asarray(n_tiles, dtype=float)
        print("[INFO] n_tiles added")
    except Exception as e:
        print(f"[WARN] tile matrix failed: {e}")

# ---------------- nucleosome signal ----------------
if DO_NUCLEOSOME:
    print("[INFO] Computing nucleosome signal (mono/NFR) …")
    bc_set=set(adata.obs_names)
    nuc_signal=compute_nucleosome_signal_by_streaming(FRAG, keep_barcodes=bc_set, max_lines=NUC_MAX_LINES)
    adata.obs["nucleosome_signal"]=pd.Series(nuc_signal).reindex(adata.obs_names).astype(float)
    print("[INFO] nucleosome_signal added")

# ---------------- metrics CSV ----------------
ensure_dir(OUTDIR)
metrics_csv=os.path.join(OUTDIR, f"{PREFIX}.qc_metrics.csv")
wanted=["n_fragment","frac_dup","frac_mito","tsse","frip","n_peaks","n_tiles","nucleosome_signal"]
avail=[k for k in wanted if k in adata.obs]
df=pd.DataFrame(index=adata.obs_names)
for k in avail: df[k]=np.asarray(adata.obs[k])
df.to_csv(metrics_csv); print("[INFO] metrics CSV:", metrics_csv)

# ---------------- filtering & subset export ----------------
n_fragment=np.asarray(adata.obs["n_fragment"], dtype=float)
tsse      =np.asarray(adata.obs["tsse"],       dtype=float)
frac_mito =np.asarray(adata.obs["frac_mito"],  dtype=float)
frac_dup  =np.asarray(adata.obs["frac_dup"],   dtype=float)
mask=( (n_fragment>=MIN_COUNTS) & (n_fragment<=MAX_COUNTS) &
       (tsse>=MIN_TSSE) & (frac_mito<=MAX_MITO) & (frac_dup<=MAX_DUP) )
if ("frip" in adata.obs) and (PEAKS_CLEAN and os.path.exists(PEAKS_CLEAN)):
    frip=np.asarray(adata.obs["frip"], dtype=float)
    mask = mask & (frip>=MIN_FRIP)
adata.obs["qc_pass"]=mask

idx=np.where(mask)[0]
obs_list=list(adata.obs_names)
pass_barcodes=[obs_list[i] for i in idx]

filtered_h5ad=os.path.join(OUTDIR, f"{PREFIX}.filtered.h5ad")
try:
    _=adata.subset(obs_indices=idx, inplace=False, out=filtered_h5ad, backend="hdf5")
    print("[INFO] filtered:", filtered_h5ad, "(subset/out)")
except Exception as e:
    print(f"[WARN] subset(out=...) failed ({e}); re-import whitelist fallback …")
    tmp=os.path.join(OUTDIR, f"{PREFIX}.filtered.tmp.h5ad")
    adata_f=snap.pp.import_fragments(
        FRAM=FRAG, chrom_sizes=genome_obj, file=tmp,
        sorted_by_barcode=False, whitelist=pass_barcodes,
        backend="hdf5", tempdir=BIG_TMP, chunk_size=CHUNK_SIZE, n_jobs=N_JOBS,
        min_num_fragments=MIN_FRAGS_IMPORT,
    )
    snap.metrics.tsse(adata_f, genome_obj)
    if PEAKS_CLEAN and os.path.exists(PEAKS_CLEAN):
        try:    snap.metrics.frip(adata_f, {"frip": PEAKS_CLEAN})
        except: 
            try: snap.metrics.frip(adata_f, PEAKS_CLEAN)
            except: pass
    os.replace(tmp, filtered_h5ad)
    if hasattr(adata_f, "close"):
        try: adata_f.close()
        except Exception: pass
    print("[INFO] filtered:", filtered_h5ad, "(re-import whitelist)")

pd.Series(pass_barcodes).to_csv(os.path.join(OUTDIR, f"{PREFIX}.filtered_barcodes.txt"), index=False, header=False)

# ---------------- scrublet doublet (tile/peak) ----------------
dbl_key=None
srcA=tileA if tileA is not None else (peakA if peakA is not None else None)
if srcA is None:
    try:
        print("[INFO] Building temporary 5kb tile matrix for doublet …")
        tmp_tile=os.path.join(OUTDIR, f"{PREFIX}.tile_5kb.doublet.zarr")
        srcA=snap.pp.make_tile_matrix(adata, file=tmp_tile, bin_size=5000)
    except Exception as e:
        print(f("[WARN] temp tile for doublet failed: {e}"))

if srcA is not None:
    try:
        snap.pp.scrublet(srcA, n_jobs=SCRUBLET_JOBS)
        s=pd.Series(np.asarray(srcA.obs["doublet_score"]), index=list(srcA.obs_names), name="doublet_score")
        adata.obs["doublet_score"]=s.reindex(adata.obs_names).astype(float).values
        adata.obs["doublet"]=np.asarray(adata.obs["doublet_score"], dtype=float)>float(DOUBLETS_THRESHOLD)
        dbl_key="doublet_score"
        print("[INFO] doublet_score written; threshold >", float(DOUBLETS_THRESHOLD))
    except Exception as e:
        print(f"[WARN] scrublet failed: {e}")

# ---------------- Unified PDF report ----------------
pdf_path=os.path.join(OUTDIR, f"{PREFIX}.qc_report.pdf") if WRITE_PDF else None
pdf=PdfPages(pdf_path) if WRITE_PDF else None
def savefig_or_inline():
    if pdf is not None: pdf.savefig(); plt.close()
    else: plt.show()

# Knee
barcode_knee_plot(adata.obs["n_fragment"].to_numpy(), "Barcode knee (n_fragment)"); savefig_or_inline()

# Fragment length (official or fallback)
ok_fsd=True
try:
    fsd=snap.metrics.frag_size_distr(adata); ok_fsd=(fsd is not None)
except Exception: ok_fsd=False
if ok_fsd:
    x=np.arange(len(fsd)); y=np.array(fsd, dtype=float); x=x[1:]; y=y[1:]
    plt.figure(); plt.plot(x,y); plt.yscale("log")
    plt.xlabel("fragment length (bp)"); plt.ylabel("count (log)")
    plt.title("Fragment length distribution"); plt.tight_layout(); savefig_or_inline()
    yp=y/(y.sum() if y.sum()>0 else 1.0)
    plt.figure(); plt.plot(x,yp)
    plt.xlabel("Size of Fragments (bp)"); plt.ylabel("Fragments (%)")
    plt.title("Fragment Size Distribution (percent)"); plt.tight_layout(); savefig_or_inline()
else:
    x,y,overflow=frag_len_histogram(FRAG, max_len=800, step=1)
    plt.figure(); plt.plot(x,y)
    plt.xlabel("Fragment size (bp)"); plt.ylabel("Fraction")
    plt.title(f"Fragment length distribution (fallback); overflow>{800}bp: {overflow:,}")
    plt.tight_layout(); savefig_or_inline()

# TSSe vs n_fragment
scatter_plot(np.asarray(adata.obs["n_fragment"]), np.asarray(adata.obs["tsse"]),
             "TSSe vs n_fragment", "n_fragment", "TSSe", logx=True); savefig_or_inline()

# Density panel
nf=np.asarray(adata.obs["n_fragment"], dtype=float); ts=np.asarray(adata.obs["tsse"], dtype=float)
log_nf=np.log10(np.clip(nf,1,None))
plt.figure(figsize=(6,6))
plt.hist2d(log_nf, ts, bins=200, cmap="viridis"); cb=plt.colorbar(); cb.set_label("density")
plt.axvline(np.log10(max(MIN_COUNTS,1)), ls="--", c="k"); plt.axhline(MIN_TSSE, ls="--", c="k")
med_nf=float(np.median(nf)); med_ts=float(np.median(ts)); n_pass=int(np.sum(np.asarray(adata.obs["qc_pass"], dtype=bool)))
txt=(f"{PREFIX}\nPass = {n_pass}\nMedian Frags = {med_nf:.1f}\nMedian TSSe = {med_ts:.4f}")
plt.text(log_nf.min()+0.05, ts.max()*0.95, txt, va="top")
plt.xlabel("Log10 (Unique Fragments)"); plt.ylabel("TSS Enrichment")
plt.title("TSSe vs Log10(Unique Fragments) — density"); plt.tight_layout(); savefig_or_inline()

# Histograms
if "frip" in adata.obs: hist_plot(np.asarray(adata.obs["frip"]), "Histogram: FRiP", "FRiP"); savefig_or_inline()
hist_plot(np.asarray(adata.obs["n_fragment"]), "Histogram: n_fragment", "n_fragment"); savefig_or_inline()
hist_plot(np.asarray(adata.obs["tsse"]),      "Histogram: TSSe",       "TSSe");       savefig_or_inline()
if "frac_mito" in adata.obs: hist_plot(np.asarray(adata.obs["frac_mito"]), "Histogram: frac_mito", "frac_mito"); savefig_or_inline()
if "frac_dup"  in adata.obs: hist_plot(np.asarray(adata.obs["frac_dup"]),  "Histogram: frac_dup",  "frac_dup");  savefig_or_inline()
if "n_peaks"   in adata.obs: hist_plot(np.asarray(adata.obs["n_peaks"]),   "Histogram: n_peaks per cell", "n_peaks"); savefig_or_inline()
if "n_tiles"   in adata.obs: hist_plot(np.asarray(adata.obs["n_tiles"]),   "Histogram: n_tiles per cell", "n_tiles"); savefig_or_inline()

# Correlations
if "frip" in adata.obs:
    scatter_plot(np.asarray(adata.obs["n_fragment"]), np.asarray(adata.obs["frip"]),
                 "FRiP vs n_fragment", "n_fragment", "FRiP", logx=True); savefig_or_inline()
    scatter_plot(np.asarray(adata.obs["frip"]), np.asarray(adata.obs["tsse"]),
                 "TSSe vs FRiP", "FRiP", "TSSe"); savefig_or_inline()
if "frac_mito" in adata.obs:
    scatter_plot(np.asarray(adata.obs["n_fragment"]), np.asarray(adata.obs["frac_mito"]),
                 "Mito fraction vs n_fragment", "n_fragment", "frac_mito", logx=True); savefig_or_inline()
if "frac_dup" in adata.obs:
    scatter_plot(np.asarray(adata.obs["n_fragment"]), np.asarray(adata.obs["frac_dup"]),
                 "Duplicate fraction vs n_fragment", "n_fragment", "frac_dup", logx=True); savefig_or_inline()

# Pass/fail overlay
try:
    x=np.asarray(adata.obs["n_fragment"], dtype=float); y=np.asarray(adata.obs["tsse"], dtype=float); keep=np.asarray(adata.obs["qc_pass"], dtype=bool)
    plt.figure()
    plt.scatter(x[~keep], y[~keep], s=3, alpha=0.35, label="fail")
    plt.scatter(x[keep],  y[keep],  s=3, alpha=0.65, label="pass")
    plt.xscale("log"); plt.xlabel("n_fragment"); plt.ylabel("TSSe"); plt.title("QC pass/fail (TSSe vs n_fragment)")
    plt.legend(); plt.tight_layout(); savefig_or_inline()
except Exception as e:
    print("[WARN] pass/fail plot failed:", e)

# Threshold overlays
for col, thr, label in [("tsse", MIN_TSSE, "TSSe ≥ MIN_TSSE"),
                        ("frac_mito", MAX_MITO, "frac_mito ≤ MAX_MITO"),
                        ("frac_dup",  MAX_DUP,  "frac_dup ≤ MAX_DUP")]:
    if col in adata.obs:
        v=np.asarray(adata.obs[col], dtype=float)
        plt.figure(); plt.hist(v[~np.isnan(v)], bins=60); plt.axvline(thr)
        plt.title(f"Threshold: {label}"); plt.xlabel(col); plt.ylabel("count"); plt.tight_layout(); savefig_or_inline()

if ("frip" in adata.obs) and (PEAKS is not None):
    v=np.asarray(adata.obs["frip"], dtype=float)
    plt.figure(); plt.hist(v[~np.isnan(v)], bins=60); plt.axvline(MIN_FRIP)
    plt.title("Threshold: FRiP ≥ MIN_FRIP"); plt.xlabel("frip"); plt.ylabel("count"); plt.tight_layout(); savefig_or_inline()

# Doublet figs
if "doublet_score" in adata.obs:
    v=np.asarray(adata.obs["doublet_score"], dtype=float); thr=float(DOUBLETS_THRESHOLD)
    plt.figure(); plt.hist(v[~np.isnan(v)], bins=60); plt.axvline(thr, ls="--", c="r")
    plt.xlabel("doublet_score"); plt.ylabel("count"); plt.title("Histogram: doublet_score")
    plt.tight_layout(); savefig_or_inline()

    plt.figure(); plt.scatter(np.asarray(adata.obs["n_fragment"]), v, s=3, alpha=0.5)
    plt.xscale("log"); plt.xlabel("n_fragment"); plt.ylabel("doublet_score")
    plt.title("doublet_score vs n_fragment"); plt.tight_layout(); savefig_or_inline()

    plt.figure(); plt.scatter(np.asarray(adata.obs["tsse"]), v, s=3, alpha=0.5)
    plt.xlabel("TSSe"); plt.ylabel("doublet_score"); plt.title("doublet_score vs TSSe")
    plt.tight_layout(); savefig_or_inline()

    plt.figure()
    kk=np.asarray(adata.obs.get("doublet", np.zeros_like(v, dtype=bool)), dtype=bool)
    plt.scatter(np.asarray(adata.obs["n_fragment"])[~kk], np.asarray(adata.obs["tsse"])[~kk], s=3, alpha=0.3, label="singlet")
    plt.scatter(np.asarray(adata.obs["n_fragment"])[ kk], np.asarray(adata.obs["tsse"])[ kk], s=5, alpha=0.6, label="doublet")
    plt.xscale("log"); plt.xlabel("n_fragment"); plt.ylabel("TSSe")
    plt.title("TSSe vs n_fragment (doublet overlay)"); plt.legend()
    plt.tight_layout(); savefig_or_inline()

if WRITE_PDF and pdf is not None:
    pdf.close(); print("[INFO] Unified QC PDF:", pdf_path)

# Close handles
for obj in [locals().get('adata', None), locals().get('peakA', None), locals().get('tileA', None)]:
    if hasattr(obj, "close"):
        try: obj.close()
        except Exception: pass
# =====================================================================================================


# In[9]:


import os
os.chdir("/path/to/scATAC/brain/txci-atac")


# In[10]:


# ===================== scATAC-seq QC (SnapATAC2) — robust all-in-one (txci-atac) =====================
# Temp/Cache/Outputs => /path/to/tmp
# Inputs in CWD:
#   - GSM7852211_mm10.merged.fragments.tsv.gz  (+ .tbi or .tbi.gz)
#   - GSM7852211_mm10.counts.sparseMatrix.rows.txt.gz  (peak-like rows)
#   - GSM7852211_mm10.counts.sparseMatrix.cols.txt.gz  (barcodes)

# ---------------- Env & params ----------------
import os, sys, re, gzip, time, gc, warnings, logging
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
os.environ.setdefault("HDF5_DISABLE_VERSION_CHECK", "2")  # optional

# <<< CHANGE HERE IF NEEDED >>>
FRAG        = "GSM7852211_mm10.merged.fragments.tsv.gz"
COUNTS_ROWS = "GSM7852211_mm10.counts.sparseMatrix.rows.txt.gz"
COUNTS_COLS = "GSM7852211_mm10.counts.sparseMatrix.cols.txt.gz"
GENOME      = "mm10"

# Force temp/cache/outputs to big scratch
BIG_TMP   = "/path/to/tmp"                 # temp + caches
OUTDIR    = "/path/to/txci_atac_qc"    # outputs
PREFIX    = "GSM7852211_scATAC"

# QC thresholds
MIN_COUNTS, MAX_COUNTS = 1000, 200000
MIN_TSSE,   MIN_FRIP   = 6.0,   0.20
MAX_MITO,   MAX_DUP    = 0.20,  0.80

# Options / tuning
DO_NUCLEOSOME       = True
NUC_MAX_LINES       = None
WRITE_PDF           = True
REBUILD_PEAK_MATRIX = False      # keep False to save disk; tile fallback enables doublet/complexity
DOUBLETS_THRESHOLD  = 0.5
SCRUBLET_JOBS       = 4

# Prefer Zarr; modest chunks; early min-frags filter
BACKEND_PREF       = "zarr"
CHUNK_SIZE         = 1_000_000
N_JOBS             = 4
MIN_FRAGS_IMPORT   = max(500, MIN_COUNTS)

# ---------------- Deps ----------------
import numpy as np, pandas as pd, h5py
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
warnings.filterwarnings("ignore", message=".*_import_from_c.*")
try:
    import snapatac2 as snap
except Exception as e:
    print("ERROR: snapatac2 not installed. Try: conda create -n snap2 -y python=3.10 snapatac2 h5py pandas matplotlib", file=sys.stderr)
    raise

# ---------------- Helpers ----------------
from pathlib import Path
def ensure_dir(p): Path(p).mkdir(parents=True, exist_ok=True); return p

def map_genome(label: str):
    ll=(label or "").strip().lower()
    if ll in ["mm10","mouse","grcm38"]: return snap.genome.mm10
    if ll in ["mm39","grcm39"]:         return snap.genome.mm39
    if ll in ["hg38","grch38","human"]: return snap.genome.hg38
    if ll in ["hg19","grch37"]:         return snap.genome.hg19
    raise ValueError(f"Unsupported genome: {label}")

def set_all_tmp(tmpdir: str):
    os.environ["TMPDIR"] = tmpdir
    os.environ["TMP"]    = tmpdir
    os.environ["TEMP"]   = tmpdir
    os.environ.setdefault("MPLCONFIGDIR", os.path.join(tmpdir, "mpl"))

def ensure_tbi_for_frag(frag_path: str):
    """Ensure .tbi exists; if only .tbi.gz exists, decompress to .tbi."""
    tbi = frag_path + ".tbi"
    tbi_gz = tbi + ".gz"
    if os.path.exists(tbi):
        return tbi
    if os.path.exists(tbi_gz):
        print(f"[INFO] Decompressing {tbi_gz} -> {tbi}")
        with gzip.open(tbi_gz, "rb") as fi, open(tbi, "wb") as fo:
            fo.write(fi.read())
        return tbi
    print(f"[WARN] No .tbi found for {frag_path} (SnapATAC2 may index to TMP).")
    return None

def sanitize_bed(infile: str, outfile: str):
    """Ensure 3-col BED, strip headers/bad lines."""
    if infile is None or (not os.path.exists(infile)): raise FileNotFoundError(infile)
    op_in = gzip.open if infile.endswith(".gz") else open
    bad=kept=0
    with op_in(infile, "rt") as fin, open(outfile, "w") as fo:
        for ln in fin:
            s=ln.strip()
            if (not s) or s.startswith("#") or s.startswith("track") or s.startswith("browser"):
                bad+=1; continue
            p=re.split(r"\s+", s)
            if len(p)<3: bad+=1; continue
            try: a=int(p[1]); b=int(p[2])
            except: bad+=1; continue
            if b<=a: bad+=1; continue
            fo.write(f"{p[0]}\t{a}\t{b}\n"); kept+=1
    print(f"[Clean BED] kept={kept:,} dropped={bad:,} -> {outfile}")
    return outfile

def try_build_peaks_from_rows(rows_path: str, out_bed: str):
    """Build peaks BED from rows: supports 'chr:start-end', 'chr start end', 'chr_start_end'."""
    if rows_path is None or (not os.path.exists(rows_path)): return None
    cnt=kept=0
    op = gzip.open if rows_path.endswith(".gz") else open
    with op(rows_path, "rt") as f, open(out_bed, "w") as bed:
        for line in f:
            s=line.strip()
            if not s or s.startswith("#"): continue
            cnt+=1
            chrom=start=end=None
            m=re.match(r"^(chr\S+):(\d+)-(\d+)$", s)
            if m: chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))
            if chrom is None:
                p=re.split(r"[\t ]+", s)
                if len(p)>=3 and p[1].isdigit() and p[2].isdigit() and p[0].startswith("chr"):
                    chrom, start, end = p[0], int(p[1]), int(p[2])
            if chrom is None:
                m2=re.match(r"^(chr\S+)[_\-:](\d+)[_\-:](\d+)$", s)
                if m2: chrom, start, end = m2.group(1), int(m2.group(2)), int(m2.group(3))
            if chrom is not None and end>start:
                bed.write(f"{chrom}\t{start}\t{end}\n"); kept+=1
    if kept>0:
        print(f"[PEAKS] built from rows: {kept:,}/{cnt:,} -> {out_bed}")
        return out_bed
    print("[PEAKS] rows not peak-like; skipping.")
    return None

def read_barcodes_from_cols(cols_path: str):
    if cols_path is None or (not os.path.exists(cols_path)): return None
    bcodes=[]
    op=gzip.open if cols_path.endswith(".gz") else open
    with op(cols_path, "rt") as f:
        for ln in f:
            s=ln.strip()
            if s: bcodes.append(s)
    print(f"[WHITELIST] loaded barcodes: {len(bcodes):,} from {cols_path}")
    return set(bcodes)

def compute_nucleosome_signal_by_streaming(frag_path, keep_barcodes=None,
                                           nfr_max=147, mono_min=147, mono_max=294, max_lines=None):
    nfr, mono = {}, {}
    i=0
    with gzip.open(frag_path, "rt") as fh:
        for line in fh:
            if max_lines is not None and i>=max_lines: break
            i+=1
            if not line or line[0]=="#": continue
            p=line.rstrip("\n").split("\t")
            if len(p)<4: continue
            try: L = int(p[2]) - int(p[1])
            except: continue
            bc=p[3]
            if keep_barcodes is not None and bc not in keep_barcodes: continue
            if L < nfr_max:               nfr[bc]  = nfr.get(bc,0)  + 1
            elif mono_min <= L < mono_max: mono[bc] = mono.get(bc,0) + 1
    ks = (keep_barcodes if keep_barcodes is not None else set(list(nfr.keys())+list(mono.keys())))
    return {bc: (mono.get(bc,0)/(nfr.get(bc,0) if nfr.get(bc,0)>0 else 1.0)) for bc in ks}

def barcode_knee_plot(counts, title):
    vals=np.sort(np.asarray(counts))[::-1]; ranks=np.arange(1,len(vals)+1)
    plt.figure(); plt.plot(ranks, vals); plt.xscale("log"); plt.yscale("log")
    plt.xlabel("barcode rank (log10)"); plt.ylabel("n_fragment (log10)")
    plt.title(title); plt.tight_layout()

def hist_plot(x, title, xlabel):
    v=np.asarray(x, dtype=float)
    plt.figure(); plt.hist(v[~np.isnan(v)], bins=60)
    plt.title(title); plt.xlabel(xlabel); plt.ylabel("count"); plt.tight_layout()

def scatter_plot(x,y,title,xlabel,ylabel,logx=False,logy=False):
    plt.figure(); plt.scatter(x,y,s=3,alpha=0.6)
    if logx: plt.xscale("log")
    if logy: plt.yscale("log")
    plt.title(title); plt.xlabel(xlabel); plt.ylabel(ylabel); plt.tight_layout()

def close_all_handles(globs: dict):
    for k,v in list(globs.items()):
        if hasattr(v,"close"):
            try: v.close()
            except Exception: pass
    gc.collect()

def unique_path(base_dir, base_name, backend):
    ts=time.strftime("%Y%m%d-%H%M%S"); pid=os.getpid()
    return os.path.join(base_dir, f"{base_name}.raw.{ts}.{pid}.{('h5ad' if backend=='hdf5' else 'zarr')}")

def _safe_row_nnz(X):
    try:
        import scipy.sparse as sp
        if sp.issparse(X): return np.asarray(X.getnnz(axis=1)).ravel()
        Xn=np.asarray(X)
        if Xn.ndim==2: return np.count_nonzero(Xn, axis=1)
    except Exception:
        pass
    return None

def frag_len_histogram(frag_path, max_len=800, step=1):
    bins=np.zeros(max_len+1, dtype=np.int64)  # 0 collects overflow
    with gzip.open(frag_path,"rt") as fh:
        for ln in fh:
            if not ln or ln[0]=="#": continue
            p=ln.rstrip("\n").split("\t")
            if len(p)<3: continue
            try: L=int(p[2])-int(p[1])
            except: continue
            if 0 < L <= max_len: bins[L]+=1
            else: bins[0]+=1
    idx=np.arange(1,max_len+1,step)
    val=bins[1:max_len+1:step].astype(float)
    val=val/val.sum() if val.sum()>0 else val
    return idx, val, int(bins[0])

# ---------------- Setup scratch & genome ----------------
ensure_dir(BIG_TMP); ensure_dir(OUTDIR); set_all_tmp(BIG_TMP)
genome_obj = map_genome(GENOME)

# Ensure .tbi present
ensure_tbi_for_frag(FRAG)

# Load whitelist from COLS
barcodes_keep = read_barcodes_from_cols(COUNTS_COLS)

# Build peaks from ROWS -> clean BED
PEAKS = None
if COUNTS_ROWS and os.path.exists(COUNTS_ROWS):
    built = try_build_peaks_from_rows(COUNTS_ROWS, os.path.join(OUTDIR, f"{PREFIX}.peaks.from_rows.bed"))
    if built:
        PEAKS = sanitize_bed(built, os.path.join(OUTDIR, f"{PREFIX}.peaks.from_rows.clean.bed"))

# ---------------- Import fragments (Zarr on scratch) ----------------
close_all_handles(globals())
backend   = BACKEND_PREF
raw_path  = unique_path(OUTDIR, PREFIX, backend)

print("[INFO] Importing fragments and computing n_fragment / frac_dup / frac_mito ...")
try:
    adata = snap.pp.import_fragments(
        FRAG,
        chrom_sizes=genome_obj,
        file=raw_path,
        sorted_by_barcode=False,
        whitelist=list(barcodes_keep) if barcodes_keep is not None else None,
        backend=("hdf5" if backend=="hdf5" else "zarr"),
        tempdir=BIG_TMP,
        chunk_size=CHUNK_SIZE,
        n_jobs=N_JOBS,
        min_num_fragments=MIN_FRAGS_IMPORT,
    )
except Exception as e:
    print(f"[WARN] import_fragments failed ({e}). Switching backend/file and retry …")
    backend = "hdf5" if backend=="zarr" else "zarr"
    raw_path = unique_path(OUTDIR, PREFIX, backend)
    adata = snap.pp.import_fragments(
        FRAG,
        chrom_sizes=genome_obj,
        file=raw_path,
        sorted_by_barcode=False,
        whitelist=list(barcodes_keep) if barcodes_keep is not None else None,
        backend=("hdf5" if backend=="hdf5" else "zarr"),
        tempdir=BIG_TMP,
        chunk_size=CHUNK_SIZE,
        n_jobs=N_JOBS,
        min_num_fragments=MIN_FRAGS_IMPORT,
    )
print("[INFO] Raw store:", raw_path)

# ---------------- TSSe with cache on scratch (fix NotFound) ----------------
from pathlib import Path as _Path
SNAP_CACHE = os.path.join(BIG_TMP, "snapatac2_cache")
_Path(SNAP_CACHE).mkdir(parents=True, exist_ok=True)
os.environ["XDG_CACHE_HOME"]      = SNAP_CACHE
os.environ["SNAPATAC2_CACHE_DIR"] = SNAP_CACHE

print("[INFO] Computing TSSe ...")
snap.metrics.tsse(adata, genome_obj)   # or: snap.metrics.tsse(adata, "/path/to/local.gtf/gff")
print("[INFO] TSSe done")

# ---------------- FRiP (if peaks available) ----------------
frip_added=False; PEAKS_CLEAN=None
if PEAKS and os.path.exists(PEAKS):
    PEAKS_CLEAN = os.path.join(OUTDIR, f"{PREFIX}.peaks.clean.bed")
    sanitize_bed(PEAKS, PEAKS_CLEAN)
    print("[INFO] Computing FRiP ...")
    try:
        snap.metrics.frip(adata, {"frip": PEAKS_CLEAN})
        frip_added=True
    except TypeError:
        snap.metrics.frip(adata, PEAKS_CLEAN)
        if "frip" not in adata.obs:
            for c in ["FRiP","Frip","FRIP"]:
                if c in adata.obs: adata.obs["frip"]=adata.obs[c]; break
        frip_added=True
    except Exception as e:
        print(f"[WARN] FRiP computation failed: {e}")
print("[INFO] FRiP done" if frip_added else "[INFO] FRiP not computed (no peaks or failed)")

# ---------------- peak matrix (optional; disabled to save disk) ----------------
peakA=None
if REBUILD_PEAK_MATRIX and PEAKS_CLEAN and os.path.exists(PEAKS_CLEAN):
    try:
        peak_matrix_path=os.path.join(OUTDIR, f"{PREFIX}.peak_matrix.{('h5ad' if backend=='hdf5' else 'zarr')}")
        peakA=snap.pp.make_peak_matrix(adata, file=peak_matrix_path, peak_file=PEAKS_CLEAN)
        nnz=_safe_row_nnz(peakA.X)
        if nnz is not None:
            adata.obs["n_peaks"]=nnz
            print("[INFO] n_peaks added")
    except Exception as e:
        print(f"[WARN] peak matrix failed: {e}")

# ---------------- tile matrix fallback (provides n_tiles & doublet input) ----------------
tileA=None
if "n_peaks" not in adata.obs:
    try:
        print("[INFO] Building 5kb tile matrix (fallback) …")
        tile_path=os.path.join(OUTDIR, f"{PREFIX}.tile_5kb.zarr")
        tileA=snap.pp.make_tile_matrix(adata, file=tile_path, bin_size=5000)
        X=tileA.X
        import scipy.sparse as sp
        n_tiles=(X>0).sum(axis=1).A1 if sp.issparse(X) else np.count_nonzero(np.asarray(X), axis=1)
        adata.obs["n_tiles"]=np.asarray(n_tiles, dtype=float)
        print("[INFO] n_tiles added")
    except Exception as e:
        print(f"[WARN] tile matrix failed: {e}")

# ---------------- nucleosome signal ----------------
if DO_NUCLEOSOME:
    print("[INFO] Computing nucleosome signal (mono/NFR) …")
    bc_set=set(adata.obs_names)
    nuc_signal=compute_nucleosome_signal_by_streaming(FRAG, keep_barcodes=bc_set, max_lines=NUC_MAX_LINES)
    adata.obs["nucleosome_signal"]=pd.Series(nuc_signal).reindex(adata.obs_names).astype(float)
    print("[INFO] nucleosome_signal added")

# ---------------- metrics CSV ----------------
ensure_dir(OUTDIR)
metrics_csv=os.path.join(OUTDIR, f"{PREFIX}.qc_metrics.csv")
wanted=["n_fragment","frac_dup","frac_mito","tsse","frip","n_peaks","n_tiles","nucleosome_signal"]
avail=[k for k in wanted if k in adata.obs]
df=pd.DataFrame(index=adata.obs_names)
for k in avail: df[k]=np.asarray(adata.obs[k])
df.to_csv(metrics_csv); print("[INFO] metrics CSV:", metrics_csv)

# ---------------- filtering & subset export ----------------
n_fragment=np.asarray(adata.obs["n_fragment"], dtype=float)
tsse      =np.asarray(adata.obs["tsse"],       dtype=float)
frac_mito =np.asarray(adata.obs["frac_mito"],  dtype=float)
frac_dup  =np.asarray(adata.obs["frac_dup"],   dtype=float)
mask=( (n_fragment>=MIN_COUNTS) & (n_fragment<=MAX_COUNTS) &
       (tsse>=MIN_TSSE) & (frac_mito<=MAX_MITO) & (frac_dup<=MAX_DUP) )
if ("frip" in adata.obs) and (PEAKS_CLEAN and os.path.exists(PEAKS_CLEAN)):
    frip=np.asarray(adata.obs["frip"], dtype=float)
    mask = mask & (frip>=MIN_FRIP)
adata.obs["qc_pass"]=mask

idx=np.where(mask)[0]
obs_list=list(adata.obs_names)
pass_barcodes=[obs_list[i] for i in idx]

filtered_h5ad=os.path.join(OUTDIR, f"{PREFIX}.filtered.h5ad")
try:
    _=adata.subset(obs_indices=idx, inplace=False, out=filtered_h5ad, backend="hdf5")
    print("[INFO] filtered:", filtered_h5ad, "(subset/out)")
except Exception as e:
    print(f"[WARN] subset(out=...) failed ({e}); re-import whitelist fallback …")
    tmp=os.path.join(OUTDIR, f"{PREFIX}.filtered.tmp.h5ad")
    adata_f=snap.pp.import_fragments(
        FRAM=FRAG, chrom_sizes=genome_obj, file=tmp,
        sorted_by_barcode=False, whitelist=pass_barcodes,
        backend="hdf5", tempdir=BIG_TMP, chunk_size=CHUNK_SIZE, n_jobs=N_JOBS,
        min_num_fragments=MIN_FRAGS_IMPORT,
    )
    snap.metrics.tsse(adata_f, genome_obj)
    if PEAKS_CLEAN and os.path.exists(PEAKS_CLEAN):
        try:    snap.metrics.frip(adata_f, {"frip": PEAKS_CLEAN})
        except: 
            try: snap.metrics.frip(adata_f, PEAKS_CLEAN)
            except: pass
    os.replace(tmp, filtered_h5ad)
    if hasattr(adata_f, "close"):
        try: adata_f.close()
        except Exception: pass
    print("[INFO] filtered:", filtered_h5ad, "(re-import whitelist)")

pd.Series(pass_barcodes).to_csv(os.path.join(OUTDIR, f"{PREFIX}.filtered_barcodes.txt"), index=False, header=False)

# ---------------- scrublet doublet (tile/peak) ----------------
dbl_key=None
srcA=tileA if tileA is not None else (peakA if peakA is not None else None)
if srcA is None:
    try:
        print("[INFO] Building temporary 5kb tile matrix for doublet …")
        tmp_tile=os.path.join(OUTDIR, f"{PREFIX}.tile_5kb.doublet.zarr")
        srcA=snap.pp.make_tile_matrix(adata, file=tmp_tile, bin_size=5000)
    except Exception as e:
        print(f("[WARN] temp tile for doublet failed: {e}"))

if srcA is not None:
    try:
        snap.pp.scrublet(srcA, n_jobs=SCRUBLET_JOBS)
        s=pd.Series(np.asarray(srcA.obs["doublet_score"]), index=list(srcA.obs_names), name="doublet_score")
        adata.obs["doublet_score"]=s.reindex(adata.obs_names).astype(float).values
        adata.obs["doublet"]=np.asarray(adata.obs["doublet_score"], dtype=float)>float(DOUBLETS_THRESHOLD)
        dbl_key="doublet_score"
        print("[INFO] doublet_score written; threshold >", float(DOUBLETS_THRESHOLD))
    except Exception as e:
        print(f"[WARN] scrublet failed: {e}")

# ---------------- Unified PDF report ----------------
pdf_path=os.path.join(OUTDIR, f"{PREFIX}.qc_report.pdf") if WRITE_PDF else None
pdf=PdfPages(pdf_path) if WRITE_PDF else None
def savefig_or_inline():
    if pdf is not None: pdf.savefig(); plt.close()
    else: plt.show()

# Knee
barcode_knee_plot(adata.obs["n_fragment"].to_numpy(), "Barcode knee (n_fragment)"); savefig_or_inline()

# Fragment length (official or fallback)
ok_fsd=True
try:
    fsd=snap.metrics.frag_size_distr(adata); ok_fsd=(fsd is not None)
except Exception: ok_fsd=False
if ok_fsd:
    x=np.arange(len(fsd)); y=np.array(fsd, dtype=float); x=x[1:]; y=y[1:]
    plt.figure(); plt.plot(x,y); plt.yscale("log")
    plt.xlabel("fragment length (bp)"); plt.ylabel("count (log)")
    plt.title("Fragment length distribution"); plt.tight_layout(); savefig_or_inline()
    yp=y/(y.sum() if y.sum()>0 else 1.0)
    plt.figure(); plt.plot(x,yp)
    plt.xlabel("Size of Fragments (bp)"); plt.ylabel("Fragments (%)")
    plt.title("Fragment Size Distribution (percent)"); plt.tight_layout(); savefig_or_inline()
else:
    x,y,overflow=frag_len_histogram(FRAG, max_len=800, step=1)
    plt.figure(); plt.plot(x,y)
    plt.xlabel("Fragment size (bp)"); plt.ylabel("Fraction")
    plt.title(f"Fragment length distribution (fallback); overflow>{800}bp: {overflow:,}")
    plt.tight_layout(); savefig_or_inline()

# TSSe vs n_fragment
scatter_plot(np.asarray(adata.obs["n_fragment"]), np.asarray(adata.obs["tsse"]),
             "TSSe vs n_fragment", "n_fragment", "TSSe", logx=True); savefig_or_inline()

# Density panel
nf=np.asarray(adata.obs["n_fragment"], dtype=float); ts=np.asarray(adata.obs["tsse"], dtype=float)
log_nf=np.log10(np.clip(nf,1,None))
plt.figure(figsize=(6,6))
plt.hist2d(log_nf, ts, bins=200, cmap="viridis"); cb=plt.colorbar(); cb.set_label("density")
plt.axvline(np.log10(max(MIN_COUNTS,1)), ls="--", c="k"); plt.axhline(MIN_TSSE, ls="--", c="k")
med_nf=float(np.median(nf)); med_ts=float(np.median(ts)); n_pass=int(np.sum(np.asarray(adata.obs["qc_pass"], dtype=bool)))
txt=(f"{PREFIX}\nPass = {n_pass}\nMedian Frags = {med_nf:.1f}\nMedian TSSe = {med_ts:.4f}")
plt.text(log_nf.min()+0.05, ts.max()*0.95, txt, va="top")
plt.xlabel("Log10 (Unique Fragments)"); plt.ylabel("TSS Enrichment")
plt.title("TSSe vs Log10(Unique Fragments) — density"); plt.tight_layout(); savefig_or_inline()

# Histograms
if "frip" in adata.obs: hist_plot(np.asarray(adata.obs["frip"]), "Histogram: FRiP", "FRiP"); savefig_or_inline()
hist_plot(np.asarray(adata.obs["n_fragment"]), "Histogram: n_fragment", "n_fragment"); savefig_or_inline()
hist_plot(np.asarray(adata.obs["tsse"]),      "Histogram: TSSe",       "TSSe");       savefig_or_inline()
if "frac_mito" in adata.obs: hist_plot(np.asarray(adata.obs["frac_mito"]), "Histogram: frac_mito", "frac_mito"); savefig_or_inline()
if "frac_dup"  in adata.obs: hist_plot(np.asarray(adata.obs["frac_dup"]),  "Histogram: frac_dup",  "frac_dup");  savefig_or_inline()
if "n_peaks"   in adata.obs: hist_plot(np.asarray(adata.obs["n_peaks"]),   "Histogram: n_peaks per cell", "n_peaks"); savefig_or_inline()
if "n_tiles"   in adata.obs: hist_plot(np.asarray(adata.obs["n_tiles"]),   "Histogram: n_tiles per cell", "n_tiles"); savefig_or_inline()

# Correlations
if "frip" in adata.obs:
    scatter_plot(np.asarray(adata.obs["n_fragment"]), np.asarray(adata.obs["frip"]),
                 "FRiP vs n_fragment", "n_fragment", "FRiP", logx=True); savefig_or_inline()
    scatter_plot(np.asarray(adata.obs["frip"]), np.asarray(adata.obs["tsse"]),
                 "TSSe vs FRiP", "FRiP", "TSSe"); savefig_or_inline()
if "frac_mito" in adata.obs:
    scatter_plot(np.asarray(adata.obs["n_fragment"]), np.asarray(adata.obs["frac_mito"]),
                 "Mito fraction vs n_fragment", "n_fragment", "frac_mito", logx=True); savefig_or_inline()
if "frac_dup" in adata.obs:
    scatter_plot(np.asarray(adata.obs["n_fragment"]), np.asarray(adata.obs["frac_dup"]),
                 "Duplicate fraction vs n_fragment", "n_fragment", "frac_dup", logx=True); savefig_or_inline()

# Pass/fail overlay
try:
    x=np.asarray(adata.obs["n_fragment"], dtype=float); y=np.asarray(adata.obs["tsse"], dtype=float); keep=np.asarray(adata.obs["qc_pass"], dtype=bool)
    plt.figure()
    plt.scatter(x[~keep], y[~keep], s=3, alpha=0.35, label="fail")
    plt.scatter(x[keep],  y[keep],  s=3, alpha=0.65, label="pass")
    plt.xscale("log"); plt.xlabel("n_fragment"); plt.ylabel("TSSe"); plt.title("QC pass/fail (TSSe vs n_fragment)")
    plt.legend(); plt.tight_layout(); savefig_or_inline()
except Exception as e:
    print("[WARN] pass/fail plot failed:", e)

# Threshold overlays
for col, thr, label in [("tsse", MIN_TSSE, "TSSe ≥ MIN_TSSE"),
                        ("frac_mito", MAX_MITO, "frac_mito ≤ MAX_MITO"),
                        ("frac_dup",  MAX_DUP,  "frac_dup ≤ MAX_DUP")]:
    if col in adata.obs:
        v=np.asarray(adata.obs[col], dtype=float)
        plt.figure(); plt.hist(v[~np.isnan(v)], bins=60); plt.axvline(thr)
        plt.title(f"Threshold: {label}"); plt.xlabel(col); plt.ylabel("count"); plt.tight_layout(); savefig_or_inline()

if ("frip" in adata.obs) and (PEAKS is not None):
    v=np.asarray(adata.obs["frip"], dtype=float)
    plt.figure(); plt.hist(v[~np.isnan(v)], bins=60); plt.axvline(MIN_FRIP)
    plt.title("Threshold: FRiP ≥ MIN_FRIP"); plt.xlabel("frip"); plt.ylabel("count"); plt.tight_layout(); savefig_or_inline()

# Doublet figs
if "doublet_score" in adata.obs:
    v=np.asarray(adata.obs["doublet_score"], dtype=float); thr=float(DOUBLETS_THRESHOLD)
    plt.figure(); plt.hist(v[~np.isnan(v)], bins=60); plt.axvline(thr, ls="--", c="r")
    plt.xlabel("doublet_score"); plt.ylabel("count"); plt.title("Histogram: doublet_score")
    plt.tight_layout(); savefig_or_inline()

    plt.figure(); plt.scatter(np.asarray(adata.obs["n_fragment"]), v, s=3, alpha=0.5)
    plt.xscale("log"); plt.xlabel("n_fragment"); plt.ylabel("doublet_score")
    plt.title("doublet_score vs n_fragment"); plt.tight_layout(); savefig_or_inline()

    plt.figure(); plt.scatter(np.asarray(adata.obs["tsse"]), v, s=3, alpha=0.5)
    plt.xlabel("TSSe"); plt.ylabel("doublet_score"); plt.title("doublet_score vs TSSe")
    plt.tight_layout(); savefig_or_inline()

    plt.figure()
    kk=np.asarray(adata.obs.get("doublet", np.zeros_like(v, dtype=bool)), dtype=bool)
    plt.scatter(np.asarray(adata.obs["n_fragment"])[~kk], np.asarray(adata.obs["tsse"])[~kk], s=3, alpha=0.3, label="singlet")
    plt.scatter(np.asarray(adata.obs["n_fragment"])[ kk], np.asarray(adata.obs["tsse"])[ kk], s=5, alpha=0.6, label="doublet")
    plt.xscale("log"); plt.xlabel("n_fragment"); plt.ylabel("TSSe")
    plt.title("TSSe vs n_fragment (doublet overlay)"); plt.legend()
    plt.tight_layout(); savefig_or_inline()

if WRITE_PDF and pdf is not None:
    pdf.close(); print("[INFO] Unified QC PDF:", pdf_path)

# Close handles
for obj in [locals().get('adata', None), locals().get('peakA', None), locals().get('tileA', None)]:
    if hasattr(obj, "close"):
        try: obj.close()
        except Exception: pass
# =====================================================================================================


# In[2]:


# ===================== Multi-tech scATAC Benchmark (SnapATAC2, fragments-only, HDF5) =====================
# Uses your local mm10 annotation: gencode vM31 GTF/GZ
# Temp/cache/output directed to big scratch to avoid $HOME/TMP quota issues

# ---------- CONFIG ----------
BIG_TMP  = "/path/to/tmp"                           # temp/cache on big disk
OUTDIR   = "/path/to/benchmark_scATAC_hdf5"     # outputs on big disk
GENOME   = "mm10"                                                    # common build across datasets
GENE_ANNOT = "/path/to/gencode.vM31.basic.annotation.gtf.gz"  # your ref

# Uniform QC thresholds
MIN_COUNTS, MAX_COUNTS = 1000, 200000
MIN_TSSE,   MIN_FRIP   = 6.0,  0.20
MAX_MITO,   MAX_DUP    = 0.20, 0.80

# Import tuning (HDF5 only)
BACKEND          = "hdf5"
CHUNK_SIZE       = 1_000_000
N_JOBS_IMPORT    = 4
MIN_FRAGS_IMPORT = max(500, MIN_COUNTS)

# Doublets / extras
DO_NUCLEOSOME   = True
NUC_MAX_LINES   = None
SCRUBLET_JOBS   = 4
DOWNSAMPLE_PLOT = 10_000  # for TSSe vs log10(nFrags) panel

# Your 5 datasets (mm10)
DATASETS = {
    "10x_multiome": "/path/to/scATAC/brain/10xmultiome/M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_atac_fragments.tsv.gz",
    "10x_scatac"  : "/path/to/scATAC/brain/10xscatac/atac_v1_E18_brain_flash_5k_fragments.tsv.gz",
    "SeekGene"    : "/path/to/scATAC/brain/seekgene/atac_fragments.tsv.gz",
    "TXCI"        : "/path/to/scATAC/brain/txci-atac/GSM7852211_mm10.merged.fragments.tsv.gz",
    "Droplet"     : "/path/to/scATAC/brain/droplet/GSM3507342_Mouse1-Channel1.fragments.tsv.gz",
}

# ---------- ENV ----------
import os, sys, re, gzip, time, gc, warnings, logging
from pathlib import Path
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
os.environ.setdefault("HDF5_DISABLE_VERSION_CHECK", "2")
os.environ["TMPDIR"] = os.environ["TMP"] = os.environ["TEMP"] = BIG_TMP
Path(BIG_TMP).mkdir(parents=True, exist_ok=True)
Path(OUTDIR).mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", os.path.join(BIG_TMP, "mpl"))

# ---------- DEPS ----------
import numpy as np, pandas as pd, h5py
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
warnings.filterwarnings("ignore", message=".*_import_from_c.*")
try:
    import snapatac2 as snap
except Exception as e:
    print("ERROR: snapatac2 not installed. Conda: conda create -n snap2 -y python=3.10 snapatac2 h5py pandas matplotlib", file=sys.stderr); raise

# ---------- HELPERS ----------
def ensure_tbi_for_frag(frag_path: str):
    tbi = frag_path + ".tbi"
    tbi_gz = tbi + ".gz"
    if os.path.exists(tbi): return
    if os.path.exists(tbi_gz):
        print(f"[INFO] Decompressing {tbi_gz} -> {tbi}")
        with gzip.open(tbi_gz, "rb") as fi, open(tbi, "wb") as fo: fo.write(fi.read())
    else:
        print(f"[WARN] No .tbi alongside {frag_path} (SnapATAC2 may index at TMP).")

def close_all_handles():
    for v in list(globals().values()):
        if hasattr(v, "close"):
            try: v.close()
            except Exception: pass
    gc.collect()

def unique_store(base, name):
    ts=time.strftime("%Y%m%d-%H%M%S"); pid=os.getpid()
    return os.path.join(base, f"{name}.raw.{ts}.{pid}.h5ad")

def frag_len_histogram(frag_path, max_len=800, step=1):
    bins=np.zeros(max_len+1, dtype=np.int64)  # 0 collects overflow
    with gzip.open(frag_path,"rt") as fh:
        for ln in fh:
            if not ln or ln[0]=="#": continue
            p=ln.rstrip("\n").split("\t")
            if len(p)<3: continue
            try: L=int(p[2])-int(p[1])
            except: continue
            if 0 < L <= max_len: bins[L]+=1
            else: bins[0]+=1
    idx=np.arange(1,max_len+1,step)
    val=bins[1:max_len+1:step].astype(float)
    val=val/val.sum() if val.sum()>0 else val
    return idx, val, int(bins[0])

def compute_nucleosome_signal_by_streaming(frag_path, keep_barcodes=None,
                                           nfr_max=147, mono_min=147, mono_max=294, max_lines=None):
    nfr, mono = {}, {}
    i=0
    with gzip.open(frag_path, "rt") as fh:
        for line in fh:
            if max_lines is not None and i>=max_lines: break
            i+=1
            if not line or line[0]=="#": continue
            p=line.rstrip("\n").split("\t")
            if len(p)<4: continue
            try: L = int(p[2]) - int(p[1])
            except: continue
            bc=p[3]
            if keep_barcodes is not None and bc not in keep_barcodes: continue
            if L < nfr_max:               nfr[bc]  = nfr.get(bc,0)  + 1
            elif mono_min <= L < mono_max: mono[bc] = mono.get(bc,0) + 1
    ks = (keep_barcodes if keep_barcodes is not None else set(list(nfr.keys())+list(mono.keys())))
    return {bc: (mono.get(bc,0)/(nfr.get(bc,0) if nfr.get(bc,0)>0 else 1.0)) for bc in ks}

def summarize_series(x):
    x = x[~np.isnan(x)]
    if x.size==0: return dict(n=0, median=np.nan, q1=np.nan, q3=np.nan)
    q = np.quantile(x, [0.25, 0.5, 0.75])
    return dict(n=x.size, median=q[1], q1=q[0], q3=q[2])

def run_tsse_safe(ad, gene_anno_path, genome_obj):
    """Prefer local GTF/GFF; else fill NaN TSSe and continue."""
    try:
        if gene_anno_path and os.path.exists(gene_anno_path):
            snap.metrics.tsse(ad, gene_anno_path)  # local GTF/GFF (gz ok)
            return True
        raise FileNotFoundError(f"Local annotation not found: {gene_anno_path}")
    except Exception as e:
        print("[WARN] TSSe failed:", e)
        ad.obs["tsse"] = np.full(ad.n_obs, np.nan, dtype=float)
        return False

# ---------- GENOME ----------
genome_obj = snap.genome.mm10 if GENOME.lower()=="mm10" else getattr(snap.genome, GENOME.lower())
assert os.path.exists(GENE_ANNOT), f"Annotation not found: {GENE_ANNOT}"

# ---------- BENCHMARK ----------
summary_rows = []
pdf_path = os.path.join(OUTDIR, "benchmark_qc_report.pdf")
pdf = PdfPages(pdf_path)

for name, frag in DATASETS.items():
    print(f"\n=== {name} ===")
    assert os.path.exists(frag), f"{frag} not found"
    ensure_tbi_for_frag(frag)

    # Import (HDF5 only)
    close_all_handles()
    store = unique_store(OUTDIR, name)
    print("[INFO] import_fragments …", store)
    ad = snap.pp.import_fragments(
        frag,
        chrom_sizes = genome_obj,
        file        = store,
        backend     = "hdf5",
        sorted_by_barcode = False,
        whitelist   = None,
        tempdir     = BIG_TMP,
        chunk_size  = CHUNK_SIZE,
        n_jobs      = N_JOBS_IMPORT,
        min_num_fragments = MIN_FRAGS_IMPORT,
    )

    # TSSe using your local GENCODE vM31 annotation
    print("[INFO] TSSe …")
    _ok_tsse = run_tsse_safe(ad, GENE_ANNOT, genome_obj)

    # Tile matrix + doublets
    print("[INFO] Tile 5kb + Scrublet …")
    tile = snap.pp.make_tile_matrix(ad, file=os.path.join(OUTDIR, f"{name}.tile_5kb.h5ad"),
                                    bin_size=5000, backend="hdf5")
    try:
        snap.pp.scrublet(tile, n_jobs=SCRUBLET_JOBS)
        dbl = np.asarray(tile.obs["doublet_score"])
    except Exception as e:
        print("[WARN] scrublet failed:", e); dbl = np.full(ad.n_obs, np.nan)

    # Nucleosome signal (optional)
    if DO_NUCLEOSOME:
        nuc = compute_nucleosome_signal_by_streaming(frag, keep_barcodes=set(ad.obs_names), max_lines=NUC_MAX_LINES)
        nuc = pd.Series(nuc).reindex(ad.obs_names).astype(float).values
    else:
        nuc = np.full(ad.n_obs, np.nan)

    # Per-cell QC
    qc = pd.DataFrame(index=ad.obs_names)
    for k in ["n_fragment","frac_dup","frac_mito","tsse"]:
        if k in ad.obs: qc[k] = np.asarray(ad.obs[k])
    qc["doublet_score"]      = dbl
    qc["nucleosome_signal"]  = nuc
    qc.to_csv(os.path.join(OUTDIR, f"{name}.qc_cells.csv"))

    # Pass/fail
    nfr  = qc["n_fragment"].to_numpy(float)
    tsse = qc["tsse"].to_numpy(float)
    mito = qc["frac_mito"].to_numpy(float) if "frac_mito" in qc else np.full_like(tsse, np.nan)
    dup  = qc["frac_dup"].to_numpy(float)  if "frac_dup"  in qc else np.full_like(tsse, np.nan)
    mask = (nfr>=MIN_COUNTS) & (nfr<=MAX_COUNTS) & (tsse>=MIN_TSSE)
    if not np.isnan(mito).all(): mask &= (mito<=MAX_MITO)
    if not np.isnan(dup).all():  mask &= (dup <=MAX_DUP)
    qc["qc_pass"] = mask
    qc.to_csv(os.path.join(OUTDIR, f"{name}.qc_cells.with_pass.csv"))

    # Dataset summary
    row = {
        "dataset": name,
        "cells": int(qc.shape[0]),
        "pass_cells": int(mask.sum()),
        "pass_rate": float(mask.mean()),
    }
    for col in ["n_fragment","tsse","frac_mito","frac_dup","doublet_score","nucleosome_signal"]:
        if col in qc:
            s = summarize_series(qc[col].to_numpy(float))
            row.update({f"{col}_median": s["median"], f"{col}_q1": s["q1"], f"{col}_3q": s["q3"]})
    summary_rows.append(row)

    # ===== Plots into shared PDF =====
    # 1) Knee
    vals = np.sort(nfr)[::-1]; ranks=np.arange(1, len(vals)+1)
    plt.figure(); plt.plot(ranks, vals)
    plt.xscale("log"); plt.yscale("log"); plt.xlabel("barcode rank (log10)"); plt.ylabel("n_fragment (log10)")
    plt.title(f"{name} — Knee"); plt.tight_layout(); pdf.savefig(); plt.close()

    # 2) TSSe vs log10(nFrags) (downsample)
    ds = min(DOWNSAMPLE_PLOT, len(nfr))
    idx = np.random.choice(len(nfr), ds, replace=False) if ds<len(nfr) else np.arange(len(nfr))
    plt.figure(); plt.scatter(np.log10(np.clip(nfr[idx],1,None)), tsse[idx], s=3, alpha=0.5)
    plt.axvline(np.log10(MIN_COUNTS), ls="--"); plt.axhline(MIN_TSSE, ls="--")
    plt.xlabel("log10(n_fragment)"); plt.ylabel("TSSe"); plt.title(f"{name} — TSSe vs log10(nFrags)")
    plt.tight_layout(); pdf.savefig(); plt.close()

    # 3) Histograms
    def _hist(col, title, xlabel):
        if col in qc:
            x=qc[col].to_numpy(float)
            plt.figure(); plt.hist(x[~np.isnan(x)], bins=60)
            plt.title(f"{name} — {title}"); plt.xlabel(xlabel); plt.ylabel("count")
            plt.tight_layout(); pdf.savefig(); plt.close()
    _hist("n_fragment","n_fragment","n_fragment")
    _hist("tsse","TSSe","TSSe")
    _hist("frac_mito","frac_mito","frac_mito")
    _hist("frac_dup","frac_dup","frac_dup")
    _hist("doublet_score","doublet_score","doublet_score")
    _hist("nucleosome_signal","nucleosome_signal","nucleosome_signal")

# ---------- Save summary & cross-dataset plots ----------
summary = pd.DataFrame(summary_rows).sort_values("dataset")
summary_path = os.path.join(OUTDIR, "benchmark_summary.csv")
summary.to_csv(summary_path, index=False)
print("[INFO] Summary:", summary_path)

# Cross-dataset metric panels
def metric_panel(metric, ylabel):
    vals=[]; labels=[]
    for name,_ in DATASETS.items():
        fn=os.path.join(OUTDIR, f"{name}.qc_cells.csv")
        if not os.path.exists(fn): continue
        df=pd.read_csv(fn, index_col=0)
        if metric in df:
            x=df[metric].to_numpy(float); x=x[~np.isnan(x)]
            if x.size>0: vals.append(x); labels.append(name)
    if not vals: return
    plt.figure(figsize=(max(6,len(labels)*1.2),4))
    plt.boxplot(vals, showfliers=False)
    plt.xticks(range(1,len(labels)+1), labels, rotation=30, ha="right")
    plt.ylabel(ylabel); plt.title(f"{metric} — cross-dataset")
    plt.tight_layout(); pdf.savefig(); plt.close()

for m, y in [("n_fragment","n_fragment"),
             ("tsse","TSSe"),
             ("frac_mito","frac_mito"),
             ("frac_dup","frac_dup"),
             ("doublet_score","doublet_score")]:
    metric_panel(m, y)

pdf.close()
print("[INFO] Consolidated PDF:", pdf_path)

# Quick ranking previews
print("\n=== Ranking by TSSe & n_fragment (higher better) ===")
print(summary[["dataset","tsse_median","n_fragment_median","pass_rate"]]
      .sort_values(["tsse_median","n_fragment_median","pass_rate"], ascending=[False,False,False]))
print("\n=== Ranking by frac_mito/frac_dup/doublet_score (lower better) ===")
print(summary[["dataset","frac_mito_median","frac_dup_median","doublet_score_median"]]
      .sort_values(["frac_mito_median","frac_dup_median","doublet_score_median"]))


# In[3]:


# ===================== Multi-tech scATAC Benchmark (SnapATAC2, fragments-only, HDF5) =====================
# Uses your local mm10 annotation (GENCODE vM31 GTF/GZ)
# Temp/cache/output on big scratch; robust to API differences across SnapATAC2 versions
# Doublets: uses tile/bin matrix when available; otherwise skips gracefully

# ---------- CONFIG ----------
BIG_TMP  = "/path/to/tmp"                           # temp/cache on big disk
OUTDIR   = "/path/to/benchmark_scATAC_hdf5"     # outputs on big disk
GENOME   = "mm10"                                                    # common build across datasets
GENE_ANNOT = "/path/to/gencode.vM31.basic.annotation.gtf.gz"

# Uniform QC thresholds
MIN_COUNTS, MAX_COUNTS = 1000, 200000
MIN_TSSE,   MIN_FRIP   = 6.0,  0.20
MAX_MITO,   MAX_DUP    = 0.20, 0.80

# Import tuning (HDF5 only)
BACKEND          = "hdf5"
CHUNK_SIZE       = 1_000_000
N_JOBS_IMPORT    = 4
MIN_FRAGS_IMPORT = max(500, MIN_COUNTS)

# Doublets / extras
DO_NUCLEOSOME   = True
NUC_MAX_LINES   = None
SCRUBLET_JOBS   = 4
DOWNSAMPLE_PLOT = 10_000  # for TSSe vs log10(nFrags) panel

# Your 5 datasets (mm10)
DATASETS = {
    "10x_multiome": "/path/to/scATAC/brain/10xmultiome/M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_atac_fragments.tsv.gz",
    "10x_scatac"  : "/path/to/scATAC/brain/10xscatac/atac_v1_E18_brain_flash_5k_fragments.tsv.gz",
    "SeekGene"    : "/path/to/scATAC/brain/seekgene/atac_fragments.tsv.gz",
    "TXCI"        : "/path/to/scATAC/brain/txci-atac/GSM7852211_mm10.merged.fragments.tsv.gz",
    "Droplet"     : "/path/to/scATAC/brain/droplet/GSM3507342_Mouse1-Channel1.fragments.tsv.gz",
}

# ---------- ENV ----------
import os, sys, re, gzip, time, gc, warnings, logging
from pathlib import Path
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
os.environ.setdefault("HDF5_DISABLE_VERSION_CHECK", "2")
os.environ["TMPDIR"] = os.environ["TMP"] = os.environ["TEMP"] = BIG_TMP
Path(BIG_TMP).mkdir(parents=True, exist_ok=True)
Path(OUTDIR).mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", os.path.join(BIG_TMP, "mpl"))

# ---------- DEPS ----------
import numpy as np, pandas as pd, h5py
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
warnings.filterwarnings("ignore", message=".*_import_from_c.*")
try:
    import snapatac2 as snap
except Exception as e:
    print("ERROR: snapatac2 not installed. Conda: conda create -n snap2 -y python=3.10 snapatac2 h5py pandas matplotlib", file=sys.stderr); raise

# ---------- HELPERS ----------
def ensure_tbi_for_frag(frag_path: str):
    tbi = frag_path + ".tbi"
    tbi_gz = tbi + ".gz"
    if os.path.exists(tbi): return
    if os.path.exists(tbi_gz):
        print(f"[INFO] Decompressing {tbi_gz} -> {tbi}")
        with gzip.open(tbi_gz, "rb") as fi, open(tbi, "wb") as fo: fo.write(fi.read())
    else:
        print(f"[WARN] No .tbi alongside {frag_path} (SnapATAC2 may index at TMP).")

def close_all_handles():
    for v in list(globals().values()):
        if hasattr(v, "close"):
            try: v.close()
            except Exception: pass
    gc.collect()

def unique_store(base, name):
    ts=time.strftime("%Y%m%d-%H%M%S"); pid=os.getpid()
    return os.path.join(base, f"{name}.raw.{ts}.{pid}.h5ad")

def frag_len_histogram(frag_path, max_len=800, step=1):
    bins=np.zeros(max_len+1, dtype=np.int64)  # 0 collects overflow
    with gzip.open(frag_path,"rt") as fh:
        for ln in fh:
            if not ln or ln[0]=="#": continue
            p=ln.rstrip("\n").split("\t")
            if len(p)<3: continue
            try: L=int(p[2])-int(p[1])
            except: continue
            if 0 < L <= max_len: bins[L]+=1
            else: bins[0]+=1
    idx=np.arange(1,max_len+1,step)
    val=bins[1:max_len+1:step].astype(float)
    val=val/val.sum() if val.sum()>0 else val
    return idx, val, int(bins[0])

def compute_nucleosome_signal_by_streaming(frag_path, keep_barcodes=None,
                                           nfr_max=147, mono_min=147, mono_max=294, max_lines=None):
    nfr, mono = {}, {}
    i=0
    with gzip.open(frag_path, "rt") as fh:
        for line in fh:
            if max_lines is not None and i>=max_lines: break
            i+=1
            if not line or line[0]=="#": continue
            p=line.rstrip("\n").split("\t")
            if len(p)<4: continue
            try: L = int(p[2]) - int(p[1])
            except: continue
            bc=p[3]
            if keep_barcodes is not None and bc not in keep_barcodes: continue
            if L < nfr_max:               nfr[bc]  = nfr.get(bc,0)  + 1
            elif mono_min <= L < mono_max: mono[bc] = mono.get(bc,0) + 1
    ks = (keep_barcodes if keep_barcodes is not None else set(list(nfr.keys())+list(mono.keys())))
    return {bc: (mono.get(bc,0)/(nfr.get(bc,0) if nfr.get(bc,0)>0 else 1.0)) for bc in ks}

def summarize_series(x):
    x = x[~np.isnan(x)]
    if x.size==0: return dict(n=0, median=np.nan, q1=np.nan, q3=np.nan)
    q = np.quantile(x, [0.25, 0.5, 0.75])
    return dict(n=x.size, median=q[1], q1=q[0], q3=q[2])

def run_tsse_safe(ad, gene_anno_path, genome_obj):
    """Prefer local GTF/GFF; else fill NaN TSSe and continue (no cache downloads)."""
    try:
        if gene_anno_path and os.path.exists(gene_anno_path):
            snap.metrics.tsse(ad, gene_anno_path)  # local GTF/GFF (gz ok)
            return True
        raise FileNotFoundError(f"Local annotation not found: {gene_anno_path}")
    except Exception as e:
        print("[WARN] TSSe failed:", e)
        ad.obs["tsse"] = np.full(ad.n_obs, np.nan, dtype=float)
        return False

def build_tile_like_matrix(ad, out_path, bin_size=5000):
    """
    Return (tile_or_bin_adata, has_matrix) using whatever this SnapATAC2 exposes.
    Tries: pp.make_tile_matrix -> pp.make_bin_matrix -> fallback(None)
    """
    if hasattr(snap.pp, "make_tile_matrix"):
        return snap.pp.make_tile_matrix(ad, file=out_path, bin_size=bin_size, backend="hdf5"), True
    if hasattr(snap.pp, "make_bin_matrix"):
        return snap.pp.make_bin_matrix(ad, file=out_path, bin_size=bin_size, backend="hdf5"), True
    print("[WARN] SnapATAC2 has neither pp.make_tile_matrix nor pp.make_bin_matrix; skipping doublets.")
    return None, False

# ---------- GENOME ----------
genome_obj = snap.genome.mm10 if GENOME.lower()=="mm10" else getattr(snap.genome, GENOME.lower())
assert os.path.exists(GENE_ANNOT), f"Annotation not found: {GENE_ANNOT}"

# ---------- BENCHMARK ----------
summary_rows = []
pdf_path = os.path.join(OUTDIR, "benchmark_qc_report.pdf")
pdf = PdfPages(pdf_path)

for name, frag in DATASETS.items():
    print(f"\n=== {name} ===")
    assert os.path.exists(frag), f"{frag} not found"
    ensure_tbi_for_frag(frag)

    # Import (HDF5 only)
    close_all_handles()
    store = unique_store(OUTDIR, name)
    print("[INFO] import_fragments …", store)
    ad = snap.pp.import_fragments(
        frag,
        chrom_sizes = genome_obj,
        file        = store,
        backend     = "hdf5",
        sorted_by_barcode = False,
        whitelist   = None,
        tempdir     = BIG_TMP,
        chunk_size  = CHUNK_SIZE,
        n_jobs      = N_JOBS_IMPORT,
        min_num_fragments = MIN_FRAGS_IMPORT,
    )

    # TSSe using your local GENCODE vM31 annotation
    print("[INFO] TSSe …")
    _ok_tsse = run_tsse_safe(ad, GENE_ANNOT, genome_obj)

    # Tile/bin matrix + doublets (robust to API differences)
    print("[INFO] Tile 5kb + Scrublet …")
    tile_path = os.path.join(OUTDIR, f"{name}.tile_5kb.h5ad")
    tile, has_mat = build_tile_like_matrix(ad, tile_path, bin_size=5000)

    dbl = np.full(ad.n_obs, np.nan)
    if has_mat:
        try:
            snap.pp.scrublet(tile, n_jobs=SCRUBLET_JOBS)
            dbl = np.asarray(tile.obs["doublet_score"])
        except Exception as e:
            print("[WARN] scrublet failed:", e)

    # Nucleosome signal (optional)
    if DO_NUCLEOSOME:
        nuc = compute_nucleosome_signal_by_streaming(frag, keep_barcodes=set(ad.obs_names), max_lines=NUC_MAX_LINES)
        nuc = pd.Series(nuc).reindex(ad.obs_names).astype(float).values
    else:
        nuc = np.full(ad.n_obs, np.nan)

    # Per-cell QC
    qc = pd.DataFrame(index=ad.obs_names)
    for k in ["n_fragment","frac_dup","frac_mito","tsse"]:
        if k in ad.obs: qc[k] = np.asarray(ad.obs[k])
    qc["doublet_score"]      = dbl
    qc["nucleosome_signal"]  = nuc
    qc.to_csv(os.path.join(OUTDIR, f"{name}.qc_cells.csv"))

    # Pass/fail
    nfr  = qc["n_fragment"].to_numpy(float)
    tsse = qc["tsse"].to_numpy(float)
    mito = qc["frac_mito"].to_numpy(float) if "frac_mito" in qc else np.full_like(tsse, np.nan)
    dup  = qc["frac_dup"].to_numpy(float)  if "frac_dup"  in qc else np.full_like(tsse, np.nan)
    mask = (nfr>=MIN_COUNTS) & (nfr<=MAX_COUNTS) & (tsse>=MIN_TSSE)
    if not np.isnan(mito).all(): mask &= (mito<=MAX_MITO)
    if not np.isnan(dup).all():  mask &= (dup <=MAX_DUP)
    qc["qc_pass"] = mask
    qc.to_csv(os.path.join(OUTDIR, f"{name}.qc_cells.with_pass.csv"))

    # Dataset summary
    row = {
        "dataset": name,
        "cells": int(qc.shape[0]),
        "pass_cells": int(mask.sum()),
        "pass_rate": float(mask.mean()),
    }
    for col in ["n_fragment","tsse","frac_mito","frac_dup","doublet_score","nucleosome_signal"]:
        if col in qc:
            s = summarize_series(qc[col].to_numpy(float))
            row.update({f"{col}_median": s["median"], f"{col}_q1": s["q1"], f"{col}_q3": s["q3"]})
    summary_rows.append(row)

    # ===== Plots into shared PDF =====
    # 1) Knee
    vals = np.sort(nfr)[::-1]; ranks=np.arange(1, len(vals)+1)
    plt.figure(); plt.plot(ranks, vals)
    plt.xscale("log"); plt.yscale("log"); plt.xlabel("barcode rank (log10)"); plt.ylabel("n_fragment (log10)")
    plt.title(f"{name} — Knee"); plt.tight_layout(); pdf.savefig(); plt.close()

    # 2) TSSe vs log10(nFrags) (downsample)
    ds = min(DOWNSAMPLE_PLOT, len(nfr))
    idx = np.random.choice(len(nfr), ds, replace=False) if ds<len(nfr) else np.arange(len(nfr))
    plt.figure(); plt.scatter(np.log10(np.clip(nfr[idx],1,None)), tsse[idx], s=3, alpha=0.5)
    plt.axvline(np.log10(MIN_COUNTS), ls="--"); plt.axhline(MIN_TSSE, ls="--")
    plt.xlabel("log10(n_fragment)"); plt.ylabel("TSSe"); plt.title(f"{name} — TSSe vs log10(nFrags)")
    plt.tight_layout(); pdf.savefig(); plt.close()

    # 3) Histograms
    def _hist(col, title, xlabel):
        if col in qc:
            x=qc[col].to_numpy(float)
            plt.figure(); plt.hist(x[~np.isnan(x)], bins=60)
            plt.title(f"{name} — {title}"); plt.xlabel(xlabel); plt.ylabel("count")
            plt.tight_layout(); pdf.savefig(); plt.close()
    _hist("n_fragment","n_fragment","n_fragment")
    _hist("tsse","TSSe","TSSe")
    _hist("frac_mito","frac_mito","frac_mito")
    _hist("frac_dup","frac_dup","frac_dup")
    _hist("doublet_score","doublet_score","doublet_score")
    _hist("nucleosome_signal","nucleosome_signal","nucleosome_signal")

# ---------- Save summary & cross-dataset plots ----------
summary = pd.DataFrame(summary_rows).sort_values("dataset")
summary_path = os.path.join(OUTDIR, "benchmark_summary.csv")
summary.to_csv(summary_path, index=False)
print("[INFO] Summary:", summary_path)

# Cross-dataset metric panels
def metric_panel(metric, ylabel):
    vals=[]; labels=[]
    for name,_ in DATASETS.items():
        fn=os.path.join(OUTDIR, f"{name}.qc_cells.csv")
        if not os.path.exists(fn): continue
        df=pd.read_csv(fn, index_col=0)
        if metric in df:
            x=df[metric].to_numpy(float); x=x[~np.isnan(x)]
            if x.size>0: vals.append(x); labels.append(name)
    if not vals: return
    plt.figure(figsize=(max(6,len(labels)*1.2),4))
    plt.boxplot(vals, showfliers=False)
    plt.xticks(range(1,len(labels)+1), labels, rotation=30, ha="right")
    plt.ylabel(ylabel); plt.title(f"{metric} — cross-dataset")
    plt.tight_layout(); pdf.savefig(); plt.close()

for m, y in [("n_fragment","n_fragment"),
             ("tsse","TSSe"),
             ("frac_mito","frac_mito"),
             ("frac_dup","frac_dup"),
             ("doublet_score","doublet_score")]:
    metric_panel(m, y)

pdf.close()
print("[INFO] Consolidated PDF:", pdf_path)

# Quick ranking previews
print("\n=== Ranking by TSSe & n_fragment (higher better) ===")
print(summary[["dataset","tsse_median","n_fragment_median","pass_rate"]]
      .sort_values(["tsse_median","n_fragment_median","pass_rate"], ascending=[False,False,False]))
print("\n=== Ranking by frac_mito/frac_dup/doublet_score (lower better) ===")
print(summary[["dataset","frac_mito_median","frac_dup_median","doublet_score_median"]]
      .sort_values(["frac_mito_median","frac_dup_median","doublet_score_median"]))


# In[1]:


# ===================== Multi-tech scATAC Benchmark (SnapATAC2, fragments-only, HDF5) =====================
# Uses your local mm10 annotation (GENCODE vM31 GTF/GZ)
# Temp/cache/output on big scratch; robust to API differences across SnapATAC2 versions
# NOTE: Droplet dataset removed per request.

# ---------- CONFIG ----------
BIG_TMP  = "/path/to/tmp"                           # temp/cache on big disk
OUTDIR   = "/path/to/benchmark_scATAC_hdf5"     # outputs on big disk
GENOME   = "mm10"                                                    # common build across datasets
GENE_ANNOT = "/path/to/gencode.vM31.basic.annotation.gtf.gz"

# Uniform QC thresholds
MIN_COUNTS, MAX_COUNTS = 1000, 200000
MIN_TSSE,   MIN_FRIP   = 6.0,  0.20
MAX_MITO,   MAX_DUP    = 0.20, 0.80

# Import tuning (HDF5 only)
BACKEND          = "hdf5"
CHUNK_SIZE       = 1_000_000
N_JOBS_IMPORT    = 4
MIN_FRAGS_IMPORT = max(500, MIN_COUNTS)

# Doublets / extras
DO_NUCLEOSOME   = True
NUC_MAX_LINES   = None
SCRUBLET_JOBS   = 4
DOWNSAMPLE_PLOT = 10_000  # for TSSe vs log10(nFrags) panel

# Your 4 datasets (mm10) — Droplet removed
DATASETS = {
    "10x_multiome": "/path/to/scATAC/brain/10xmultiome/M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_atac_fragments.tsv.gz",
    "10x_scatac"  : "/path/to/scATAC/brain/10xscatac/atac_v1_E18_brain_flash_5k_fragments.tsv.gz",
    "SeekGene"    : "/path/to/scATAC/brain/seekgene/atac_fragments.tsv.gz",
    "TXCI"        : "/path/to/scATAC/brain/txci-atac/GSM7852211_mm10.merged.fragments.tsv.gz",
}

# ---------- ENV ----------
import os, sys, re, gzip, time, gc, warnings, logging
from pathlib import Path
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
os.environ.setdefault("HDF5_DISABLE_VERSION_CHECK", "2")
os.environ["TMPDIR"] = os.environ["TMP"] = os.environ["TEMP"] = BIG_TMP
Path(BIG_TMP).mkdir(parents=True, exist_ok=True)
Path(OUTDIR).mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", os.path.join(BIG_TMP, "mpl"))

# ---------- DEPS ----------
import numpy as np, pandas as pd, h5py
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
warnings.filterwarnings("ignore", message=".*_import_from_c.*")
try:
    import snapatac2 as snap
except Exception as e:
    print("ERROR: snapatac2 not installed. Conda: conda create -n snap2 -y python=3.10 snapatac2 h5py pandas matplotlib", file=sys.stderr); raise

# ---------- HELPERS ----------
def ensure_tbi_for_frag(frag_path: str):
    tbi = frag_path + ".tbi"
    tbi_gz = tbi + ".gz"
    if os.path.exists(tbi): return
    if os.path.exists(tbi_gz):
        print(f"[INFO] Decompressing {tbi_gz} -> {tbi}")
        with gzip.open(tbi_gz, "rb") as fi, open(tbi, "wb") as fo: fo.write(fi.read())
    else:
        print(f"[WARN] No .tbi alongside {frag_path} (SnapATAC2 may index at TMP).")

def close_all_handles():
    for v in list(globals().values()):
        if hasattr(v, "close"):
            try: v.close()
            except Exception: pass
    gc.collect()

def unique_store(base, name):
    ts=time.strftime("%Y%m%d-%H%M%S"); pid=os.getpid()
    return os.path.join(base, f"{name}.raw.{ts}.{pid}.h5ad")

def frag_len_histogram(frag_path, max_len=800, step=1):
    bins=np.zeros(max_len+1, dtype=np.int64)  # 0 collects overflow
    with gzip.open(frag_path,"rt") as fh:
        for ln in fh:
            if not ln or ln[0]=="#": continue
            p=ln.rstrip("\n").split("\t")
            if len(p)<3: continue
            try: L=int(p[2])-int(p[1])
            except: continue
            if 0 < L <= max_len: bins[L]+=1
            else: bins[0]+=1
    idx=np.arange(1,max_len+1,step)
    val=bins[1:max_len+1:step].astype(float)
    val=val/val.sum() if val.sum()>0 else val
    return idx, val, int(bins[0])

def compute_nucleosome_signal_by_streaming(frag_path, keep_barcodes=None,
                                           nfr_max=147, mono_min=147, mono_max=294, max_lines=None):
    nfr, mono = {}, {}
    i=0
    with gzip.open(frag_path, "rt") as fh:
        for line in fh:
            if max_lines is not None and i>=max_lines: break
            i+=1
            if not line or line[0]=="#": continue
            p=ln = line.rstrip("\n").split("\t")
            if len(p)<4: continue
            try: L = int(p[2]) - int(p[1])
            except: continue
            bc=p[3]
            if keep_barcodes is not None and bc not in keep_barcodes: continue
            if L < nfr_max:               nfr[bc]  = nfr.get(bc,0)  + 1
            elif mono_min <= L < mono_max: mono[bc] = mono.get(bc,0) + 1
    ks = (keep_barcodes if keep_barcodes is not None else set(list(nfr.keys())+list(mono.keys())))
    return {bc: (mono.get(bc,0)/(nfr.get(bc,0) if nfr.get(bc,0)>0 else 1.0)) for bc in ks}

def summarize_series(x):
    x = x[~np.isnan(x)]
    if x.size==0: return dict(n=0, median=np.nan, q1=np.nan, q3=np.nan)
    q = np.quantile(x, [0.25, 0.5, 0.75])
    return dict(n=x.size, median=q[1], q1=q[0], q3=q[2])

def run_tsse_safe(ad, gene_anno_path, genome_obj):
    """Prefer local GTF/GFF; else fill NaN TSSe and continue (no cache downloads)."""
    try:
        if gene_anno_path and os.path.exists(gene_anno_path):
            snap.metrics.tsse(ad, gene_anno_path)  # local GTF/GFF (gz ok)
            return True
        raise FileNotFoundError(f"Local annotation not found: {gene_anno_path}")
    except Exception as e:
        print("[WARN] TSSe failed:", e)
        ad.obs["tsse"] = np.full(ad.n_obs, np.nan, dtype=float)
        return False

def build_tile_like_matrix(ad, out_path, bin_size=5000):
    """
    Return (tile_or_bin_adata, has_matrix) using whatever this SnapATAC2 exposes.
    Tries: pp.make_tile_matrix -> pp.make_bin_matrix -> fallback(None)
    """
    if hasattr(snap.pp, "make_tile_matrix"):
        return snap.pp.make_tile_matrix(ad, file=out_path, bin_size=bin_size, backend="hdf5"), True
    if hasattr(snap.pp, "make_bin_matrix"):
        return snap.pp.make_bin_matrix(ad, file=out_path, bin_size=bin_size, backend="hdf5"), True
    print("[WARN] SnapATAC2 has neither pp.make_tile_matrix nor pp.make_bin_matrix; skipping doublets.")
    return None, False

# ---------- GENOME ----------
genome_obj = snap.genome.mm10 if GENOME.lower()=="mm10" else getattr(snap.genome, GENOME.lower())
assert os.path.exists(GENE_ANNOT), f"Annotation not found: {GENE_ANNOT}"

# ---------- BENCHMARK ----------
summary_rows = []
pdf_path = os.path.join(OUTDIR, "benchmark_qc_report.pdf")
pdf = PdfPages(pdf_path)

for name, frag in DATASETS.items():
    print(f"\n=== {name} ===")
    assert os.path.exists(frag), f"{frag} not found"
    ensure_tbi_for_frag(frag)

    # Import (HDF5 only)
    close_all_handles()
    store = unique_store(OUTDIR, name)
    print("[INFO] import_fragments …", store)
    ad = snap.pp.import_fragments(
        frag,
        chrom_sizes = genome_obj,
        file        = store,
        backend     = "hdf5",
        sorted_by_barcode = False,
        whitelist   = None,
        tempdir     = BIG_TMP,
        chunk_size  = CHUNK_SIZE,
        n_jobs      = N_JOBS_IMPORT,
        min_num_fragments = MIN_FRAGS_IMPORT,
    )

    # TSSe using your local GENCODE vM31 annotation
    print("[INFO] TSSe …")
    _ok_tsse = run_tsse_safe(ad, GENE_ANNOT, genome_obj)

    # Tile/bin matrix + doublets (robust to API differences)
    print("[INFO] Tile 5kb + Scrublet …")
    tile_path = os.path.join(OUTDIR, f"{name}.tile_5kb.h5ad")
    tile, has_mat = build_tile_like_matrix(ad, tile_path, bin_size=5000)

    dbl = np.full(ad.n_obs, np.nan)
    if has_mat:
        try:
            snap.pp.scrublet(tile, n_jobs=SCRUBLET_JOBS)
            dbl = np.asarray(tile.obs["doublet_score"])
        except Exception as e:
            print("[WARN] scrublet failed:", e)

    # Nucleosome signal (optional)
    if DO_NUCLEOSOME:
        nuc = compute_nucleosome_signal_by_streaming(frag, keep_barcodes=set(ad.obs_names), max_lines=NUC_MAX_LINES)
        nuc = pd.Series(nuc).reindex(ad.obs_names).astype(float).values
    else:
        nuc = np.full(ad.n_obs, np.nan)

    # Per-cell QC
    qc = pd.DataFrame(index=ad.obs_names)
    for k in ["n_fragment","frac_dup","frac_mito","tsse"]:
        if k in ad.obs: qc[k] = np.asarray(ad.obs[k])
    qc["doublet_score"]      = dbl
    qc["nucleosome_signal"]  = nuc
    qc.to_csv(os.path.join(OUTDIR, f"{name}.qc_cells.csv"))

    # Pass/fail
    nfr  = qc["n_fragment"].to_numpy(float)
    tsse = qc["tsse"].to_numpy(float)
    mito = qc["frac_mito"].to_numpy(float) if "frac_mito" in qc else np.full_like(tsse, np.nan)
    dup  = qc["frac_dup"].to_numpy(float)  if "frac_dup"  in qc else np.full_like(tsse, np.nan)
    mask = (nfr>=MIN_COUNTS) & (nfr<=MAX_COUNTS) & (tsse>=MIN_TSSE)
    if not np.isnan(mito).all(): mask &= (mito<=MAX_MITO)
    if not np.isnan(dup).all():  mask &= (dup <=MAX_DUP)
    qc["qc_pass"] = mask
    qc.to_csv(os.path.join(OUTDIR, f"{name}.qc_cells.with_pass.csv"))

    # Dataset summary
    row = {
        "dataset": name,
        "cells": int(qc.shape[0]),
        "pass_cells": int(mask.sum()),
        "pass_rate": float(mask.mean()),
    }
    for col in ["n_fragment","tsse","frac_mito","frac_dup","doublet_score","nucleosome_signal"]:
        if col in qc:
            s = summarize_series(qc[col].to_numpy(float))
            row.update({f"{col}_median": s["median"], f"{col}_q1": s["q1"], f"{col}_q3": s["q3"]})
    summary_rows.append(row)

    # ===== Plots into shared PDF =====
    # 1) Knee
    vals = np.sort(nfr)[::-1]; ranks=np.arange(1, len(vals)+1)
    plt.figure(); plt.plot(ranks, vals)
    plt.xscale("log"); plt.yscale("log"); plt.xlabel("barcode rank (log10)"); plt.ylabel("n_fragment (log10)")
    plt.title(f"{name} — Knee"); plt.tight_layout(); pdf.savefig(); plt.close()

    # 2) TSSe vs log10(nFrags) (downsample)
    ds = min(DOWNSAMPLE_PLOT, len(nfr))
    idx = np.random.choice(len(nfr), ds, replace=False) if ds<len(nfr) else np.arange(len(nfr))
    plt.figure(); plt.scatter(np.log10(np.clip(nfr[idx],1,None)), tsse[idx], s=3, alpha=0.5)
    plt.axvline(np.log10(MIN_COUNTS), ls="--"); plt.axhline(MIN_TSSE, ls="--")
    plt.xlabel("log10(n_fragment)"); plt.ylabel("TSSe"); plt.title(f"{name} — TSSe vs log10(nFrags)")
    plt.tight_layout(); pdf.savefig(); plt.close()

    # 3) Histograms
    def _hist(col, title, xlabel):
        if col in qc:
            x=qc[col].to_numpy(float)
            plt.figure(); plt.hist(x[~np.isnan(x)], bins=60)
            plt.title(f"{name} — {title}"); plt.xlabel(xlabel); plt.ylabel("count")
            plt.tight_layout(); pdf.savefig(); plt.close()
    _hist("n_fragment","n_fragment","n_fragment")
    _hist("tsse","TSSe","TSSe")
    _hist("frac_mito","frac_mito","frac_mito")
    _hist("frac_dup","frac_dup","frac_dup")
    _hist("doublet_score","doublet_score","doublet_score")
    _hist("nucleosome_signal","nucleosome_signal","nucleosome_signal")

# ---------- Save summary & cross-dataset plots ----------
summary = pd.DataFrame(summary_rows).sort_values("dataset")
summary_path = os.path.join(OUTDIR, "benchmark_summary.csv")
summary.to_csv(summary_path, index=False)
print("[INFO] Summary:", summary_path)

# Cross-dataset metric panels
def metric_panel(metric, ylabel):
    vals=[]; labels=[]
    for name,_ in DATASETS.items():
        fn=os.path.join(OUTDIR, f"{name}.qc_cells.csv")
        if not os.path.exists(fn): continue
        df=pd.read_csv(fn, index_col=0)
        if metric in df:
            x=df[metric].to_numpy(float); x=x[~np.isnan(x)]
            if x.size>0: vals.append(x); labels.append(name)
    if not vals: return
    plt.figure(figsize=(max(6,len(labels)*1.2),4))
    plt.boxplot(vals, showfliers=False)
    plt.xticks(range(1,len(labels)+1), labels, rotation=30, ha="right")
    plt.ylabel(ylabel); plt.title(f"{metric} — cross-dataset")
    plt.tight_layout(); pdf.savefig(); plt.close()

for m, y in [("n_fragment","n_fragment"),
             ("tsse","TSSe"),
             ("frac_mito","frac_mito"),
             ("frac_dup","frac_dup"),
             ("doublet_score","doublet_score")]:
    metric_panel(m, y)

pdf.close()
print("[INFO] Consolidated PDF:", pdf_path)

# Quick ranking previews
print("\n=== Ranking by TSSe & n_fragment (higher better) ===")
print(summary[["dataset","tsse_median","n_fragment_median","pass_rate"]]
      .sort_values(["tsse_median","n_fragment_median","pass_rate"], ascending=[False,False,False]))
print("\n=== Ranking by frac_mito/frac_dup/doublet_score (lower better) ===")
print(summary[["dataset","frac_mito_median","frac_dup_median","doublet_score_median"]]
      .sort_values(["frac_mito_median","frac_dup_median","doublet_score_median"]))


# In[2]:


# ===================== Multi-tech scATAC Benchmark + Scoring (SnapATAC2, fragments-only, HDF5) =====================
# Four technologies: 10x multiome, 10x scATAC, SeekGene, TXCI (Droplet removed)
# - Uniform thresholds: nFrags / TSSe / mito / dup
# - TSSe uses your local annotation (GENCODE vM31 .gtf.gz)
# - Doublets: if the current SnapATAC2 version has a tile/bin matrix API, run Scrublet, otherwise skip automatically
# - Outputs: per-dataset QC table, merged PDF report, summary CSV, and scoring and ranking plots for two weight sets (Neutral / SeekGene-lean)

# ---------- CONFIG ----------
BIG_TMP  = "/path/to/tmp"                           # large-disk TMP
OUTDIR   = "/path/to/benchmark_scATAC_hdf5"     # output directory
GENOME   = "mm10"
GENE_ANNOT = "/path/to/gencode.vM31.basic.annotation.gtf.gz"  # local annotation

# Uniform thresholds (adjustable)
MIN_COUNTS, MAX_COUNTS = 1000, 200000
MIN_TSSE,   MIN_FRIP   = 6.0,  0.20
MAX_MITO,   MAX_DUP    = 0.20, 0.80

# Import parameters (HDF5 only; your environment does not support zarr)
BACKEND          = "hdf5"
CHUNK_SIZE       = 1_000_000
N_JOBS_IMPORT    = 4
MIN_FRAGS_IMPORT = max(500, MIN_COUNTS)

# Additional analysis
DO_NUCLEOSOME   = True
NUC_MAX_LINES   = None
SCRUBLET_JOBS   = 4
DOWNSAMPLE_PLOT = 10_000  # number of downsampled points for TSSe vs log10(nFrags)

# Four datasets (Droplet removed)
DATASETS = {
    "10x_multiome": "/path/to/scATAC/brain/10xmultiome/M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_atac_fragments.tsv.gz",
    "10x_scatac"  : "/path/to/scATAC/brain/10xscatac/atac_v1_E18_brain_flash_5k_fragments.tsv.gz",
    "SeekGene"    : "/path/to/scATAC/brain/seekgene/atac_fragments.tsv.gz",
    "TXCI"        : "/path/to/scATAC/brain/txci-atac/GSM7852211_mm10.merged.fragments.tsv.gz",
}

# ---------- ENV ----------
import os, sys, re, gzip, time, gc, warnings, logging
from pathlib import Path
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
os.environ.setdefault("HDF5_DISABLE_VERSION_CHECK", "2")
os.environ["TMPDIR"] = os.environ["TMP"] = os.environ["TEMP"] = BIG_TMP
Path(BIG_TMP).mkdir(parents=True, exist_ok=True)
Path(OUTDIR).mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", os.path.join(BIG_TMP, "mpl"))

# ---------- DEPS ----------
import numpy as np, pandas as pd, h5py
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
warnings.filterwarnings("ignore", message=".*_import_from_c.*")
try:
    import snapatac2 as snap
except Exception as e:
    print("ERROR: snapatac2 not installed. Conda: conda create -n snap2 -y python=3.10 snapatac2 h5py pandas matplotlib", file=sys.stderr); raise

# ---------- HELPERS ----------
def ensure_tbi_for_frag(frag_path: str):
    """Ensure .tbi exists; if only .tbi.gz is present, decompress it."""
    tbi = frag_path + ".tbi"
    tbi_gz = tbi + ".gz"
    if os.path.exists(tbi): return
    if os.path.exists(tbi_gz):
        print(f"[INFO] Decompressing {tbi_gz} -> {tbi}")
        with gzip.open(tbi_gz, "rb") as fi, open(tbi, "wb") as fo: fo.write(fi.read())
    else:
        print(f"[WARN] No .tbi alongside {frag_path} (SnapATAC2 may index at TMP).")

def close_all_handles():
    """Close open backed handles as much as possible to avoid file locking."""
    for v in list(globals().values()):
        if hasattr(v, "close"):
            try: v.close()
            except Exception: pass
    gc.collect()

def unique_store(base, name):
    ts=time.strftime("%Y%m%d-%H%M%S"); pid=os.getpid()
    return os.path.join(base, f"{name}.raw.{ts}.{pid}.h5ad")

def frag_len_histogram(frag_path, max_len=800, step=1):
    """Fallback: stream-compute the fragment length distribution, returning (x, y_fraction, overflow_count)."""
    bins=np.zeros(max_len+1, dtype=np.int64)  # index 0 collects overflow
    with gzip.open(frag_path,"rt") as fh:
        for ln in fh:
            if not ln or ln[0]=="#": continue
            p=ln.rstrip("\n").split("\t")
            if len(p)<3: continue
            try: L=int(p[2])-int(p[1])
            except: continue
            if 0 < L <= max_len: bins[L]+=1
            else: bins[0]+=1
    idx=np.arange(1,max_len+1,step)
    val=bins[1:max_len+1:step].astype(float)
    val=val/val.sum() if val.sum()>0 else val
    return idx, val, int(bins[0])

def compute_nucleosome_signal_by_streaming(frag_path, keep_barcodes=None,
                                           nfr_max=147, mono_min=147, mono_max=294, max_lines=None):
    """nucleosome signal = mono/NFR, counting length bins per barcode."""
    nfr, mono = {}, {}
    i=0
    with gzip.open(frag_path, "rt") as fh:
        for line in fh:
            if max_lines is not None and i>=max_lines: break
            i+=1
            if not line or line[0]=="#": continue
            p=line.rstrip("\n").split("\t")
            if len(p)<4: continue
            try: L = int(p[2]) - int(p[1])
            except: continue
            bc=p[3]
            if keep_barcodes is not None and bc not in keep_barcodes: continue
            if L < nfr_max:               nfr[bc]  = nfr.get(bc,0)  + 1
            elif mono_min <= L < mono_max: mono[bc] = mono.get(bc,0) + 1
    ks = (keep_barcodes if keep_barcodes is not None else set(list(nfr.keys())+list(mono.keys())))
    return {bc: (mono.get(bc,0)/(nfr.get(bc,0) if nfr.get(bc,0)>0 else 1.0)) for bc in ks}

def summarize_series(x):
    x = x[~np.isnan(x)]
    if x.size==0: return dict(n=0, median=np.nan, q1=np.nan, q3=np.nan)
    q = np.quantile(x, [0.25, 0.5, 0.75])
    return dict(n=x.size, median=q[1], q1=q[0], q3=q[2])

def run_tsse_safe(ad, gene_anno_path, genome_obj):
    """Prefer local GTF/GFF; on failure write NaN to keep the pipeline running."""
    try:
        if gene_anno_path and os.path.exists(gene_anno_path):
            snap.metrics.tsse(ad, gene_anno_path)  # local annotation (gz is fine)
            return True
        raise FileNotFoundError(f"Local annotation not found: {gene_anno_path}")
    except Exception as e:
        print("[WARN] TSSe failed:", e)
        ad.obs["tsse"] = np.full(ad.n_obs, np.nan, dtype=float)
        return False

def build_tile_like_matrix(ad, out_path, bin_size=5000):
    """
    Compatible with different API versions:
      try pp.make_tile_matrix -> pp.make_bin_matrix -> fallback(None)
    returns (tile_or_bin_adata, has_matrix)
    """
    # Try with the backend parameter
    if hasattr(snap.pp, "make_tile_matrix"):
        try:
            return snap.pp.make_tile_matrix(ad, file=out_path, bin_size=bin_size, backend="hdf5"), True
        except TypeError:
            return snap.pp.make_tile_matrix(ad, file=out_path, bin_size=bin_size), True
    if hasattr(snap.pp, "make_bin_matrix"):
        try:
            return snap.pp.make_bin_matrix(ad, file=out_path, bin_size=bin_size, backend="hdf5"), True
        except TypeError:
            return snap.pp.make_bin_matrix(ad, file=out_path, bin_size=bin_size), True
    print("[WARN] SnapATAC2 has neither pp.make_tile_matrix nor pp.make_bin_matrix; skipping doublets.")
    return None, False

# ---------- GENOME ----------
genome_obj = snap.genome.mm10 if GENOME.lower()=="mm10" else getattr(snap.genome, GENOME.lower())
assert os.path.exists(GENE_ANNOT), f"Annotation not found: {GENE_ANNOT}"

# ---------- BENCHMARK LOOP ----------
summary_rows = []
pdf_path = os.path.join(OUTDIR, "benchmark_qc_report.pdf")
pdf = PdfPages(pdf_path)

for name, frag in DATASETS.items():
    print(f"\n=== {name} ===")
    assert os.path.exists(frag), f"{frag} not found"
    ensure_tbi_for_frag(frag)

    # Import
    close_all_handles()
    store = unique_store(OUTDIR, name)
    print("[INFO] import_fragments …", store)
    ad = snap.pp.import_fragments(
        frag,
        chrom_sizes = genome_obj,
        file        = store,
        backend     = "hdf5",
        sorted_by_barcode = False,
        whitelist   = None,
        tempdir     = BIG_TMP,
        chunk_size  = CHUNK_SIZE,
        n_jobs      = N_JOBS_IMPORT,
        min_num_fragments = MIN_FRAGS_IMPORT,
    )

    # TSSe
    print("[INFO] TSSe …")
    _ok_tsse = run_tsse_safe(ad, GENE_ANNOT, genome_obj)

    # Tile/bin + Scrublet (if available)
    print("[INFO] Tile 5kb + Scrublet …")
    tile_path = os.path.join(OUTDIR, f"{name}.tile_5kb.h5ad")
    tile, has_mat = build_tile_like_matrix(ad, tile_path, bin_size=5000)
    dbl = np.full(ad.n_obs, np.nan)
    if has_mat:
        try:
            snap.pp.scrublet(tile, n_jobs=SCRUBLET_JOBS)
            dbl = np.asarray(tile.obs["doublet_score"])
        except Exception as e:
            print("[WARN] scrublet failed:", e)

    # Nucleosome signal
    if DO_NUCLEOSOME:
        nuc = compute_nucleosome_signal_by_streaming(frag, keep_barcodes=set(ad.obs_names), max_lines=NUC_MAX_LINES)
        nuc = pd.Series(nuc).reindex(ad.obs_names).astype(float).values
    else:
        nuc = np.full(ad.n_obs, np.nan)

    # Per-cell QC & pass/fail
    qc = pd.DataFrame(index=ad.obs_names)
    for k in ["n_fragment","frac_dup","frac_mito","tsse"]:
        if k in ad.obs: qc[k] = np.asarray(ad.obs[k])
    qc["doublet_score"]      = dbl
    qc["nucleosome_signal"]  = nuc

    nfr  = qc["n_fragment"].to_numpy(float)
    tsse = qc["tsse"].to_numpy(float)
    mito = qc["frac_mito"].to_numpy(float) if "frac_mito" in qc else np.full_like(tsse, np.nan)
    dup  = qc["frac_dup"].to_numpy(float)  if "frac_dup"  in qc else np.full_like(tsse, np.nan)
    mask = (nfr>=MIN_COUNTS) & (nfr<=MAX_COUNTS) & (tsse>=MIN_TSSE)
    if not np.isnan(mito).all(): mask &= (mito<=MAX_MITO)
    if not np.isnan(dup).all():  mask &= (dup <=MAX_DUP)
    qc["qc_pass"] = mask

    qc.to_csv(os.path.join(OUTDIR, f"{name}.qc_cells.csv"))
    qc.to_csv(os.path.join(OUTDIR, f"{name}.qc_cells.with_pass.csv"))

    # Summary row
    row = {
        "dataset": name,
        "cells": int(qc.shape[0]),
        "pass_cells": int(mask.sum()),
        "pass_rate": float(mask.mean()),
        "n_fragment_median": float(np.nanmedian(qc["n_fragment"])),
        "tsse_median": float(np.nanmedian(qc["tsse"])),
        "frac_mito_median": float(np.nanmedian(qc["frac_mito"])) if "frac_mito" in qc else np.nan,
        "frac_dup_median" : float(np.nanmedian(qc["frac_dup"]))  if "frac_dup"  in qc else np.nan,
        "doublet_score_median": float(np.nanmedian(qc["doublet_score"])) if np.isfinite(dbl).any() else np.nan,
        "nucleosome_signal_median": float(np.nanmedian(qc["nucleosome_signal"])) if DO_NUCLEOSOME else np.nan,
    }
    summary_rows.append(row)

    # ===== Into PDF =====
    # 1) Knee
    vals = np.sort(nfr)[::-1]; ranks=np.arange(1, len(vals)+1)
    plt.figure(); plt.plot(ranks, vals)
    plt.xscale("log"); plt.yscale("log"); plt.xlabel("barcode rank (log10)"); plt.ylabel("n_fragment (log10)")
    plt.title(f"{name} — Knee"); plt.tight_layout(); pdf.savefig(); plt.close()

    # 2) Frag length: official or fallback
    ok_fsd=True
    try:
        fsd=snap.metrics.frag_size_distr(ad); ok_fsd=(fsd is not None)
    except Exception: ok_fsd=False
    if ok_fsd:
        x=np.arange(len(fsd)); y=np.array(fsd, dtype=float); x=x[1:]; y=y[1:]
        plt.figure(); plt.plot(x,y); plt.yscale("log")
        plt.xlabel("fragment length (bp)"); plt.ylabel("count (log)")
        plt.title(f"{name} — Fragment length distribution"); plt.tight_layout(); pdf.savefig(); plt.close()
        yp=y/(y.sum() if y.sum()>0 else 1.0)
        plt.figure(); plt.plot(x,yp)
        plt.xlabel("Size of Fragments (bp)"); plt.ylabel("Fragments (%)")
        plt.title(f"{name} — Fragment Size Distribution (percent)"); plt.tight_layout(); pdf.savefig(); plt.close()
    else:
        x,y,overflow=frag_len_histogram(frag, max_len=800, step=1)
        plt.figure(); plt.plot(x,y)
        plt.xlabel("Fragment size (bp)"); plt.ylabel("Fraction")
        plt.title(f"{name} — Frag length (fallback); overflow>{800}bp: {overflow:,}")
        plt.tight_layout(); pdf.savefig(); plt.close()

    # 3) TSSe vs log10(nFrags)
    ds = min(DOWNSAMPLE_PLOT, len(nfr))
    idx = np.random.choice(len(nfr), ds, replace=False) if ds<len(nfr) else np.arange(len(nfr))
    plt.figure(); plt.scatter(np.log10(np.clip(nfr[idx],1,None)), tsse[idx], s=3, alpha=0.5)
    plt.axvline(np.log10(MIN_COUNTS), ls="--"); plt.axhline(MIN_TSSE, ls="--")
    plt.xlabel("log10(n_fragment)"); plt.ylabel("TSSe"); plt.title(f"{name} — TSSe vs log10(nFrags)")
    plt.tight_layout(); pdf.savefig(); plt.close()

    # 4) Core histograms
    def _hist(col, title, xlabel):
        if col in qc:
            x=qc[col].to_numpy(float)
            plt.figure(); plt.hist(x[~np.isnan(x)], bins=60)
            plt.title(f"{name} — {title}"); plt.xlabel(xlabel); plt.ylabel("count")
            plt.tight_layout(); pdf.savefig(); plt.close()
    _hist("n_fragment","n_fragment","n_fragment")
    _hist("tsse","TSSe","TSSe")
    _hist("frac_mito","frac_mito","frac_mito")
    _hist("frac_dup","frac_dup","frac_dup")
    _hist("doublet_score","doublet_score","doublet_score")
    _hist("nucleosome_signal","nucleosome_signal","nucleosome_signal")

# ---------- Summary and cross-dataset plots ----------
summary = pd.DataFrame(summary_rows).sort_values("dataset")
summary_path = os.path.join(OUTDIR, "benchmark_summary.csv")
summary.to_csv(summary_path, index=False)
print("[INFO] Summary:", summary_path)

# Cross-dataset boxplot
def metric_panel(metric, ylabel):
    vals=[]; labels=[]
    for name,_ in DATASETS.items():
        fn=os.path.join(OUTDIR, f"{name}.qc_cells.csv")
        if not os.path.exists(fn): continue
        df=pd.read_csv(fn, index_col=0)
        if metric in df:
            x=df[metric].to_numpy(float); x=x[~np.isnan(x)]
            if x.size>0: vals.append(x); labels.append(name)
    if not vals: return
    plt.figure(figsize=(max(6,len(labels)*1.2),4))
    plt.boxplot(vals, showfliers=False)
    plt.xticks(range(1,len(labels)+1), labels, rotation=30, ha="right")
    plt.ylabel(ylabel); plt.title(f"{metric} — cross-dataset")
    plt.tight_layout(); pdf.savefig(); plt.close()

for m, y in [("n_fragment","n_fragment"),
             ("tsse","TSSe"),
             ("frac_mito","frac_mito"),
             ("frac_dup","frac_dup"),
             ("doublet_score","doublet_score"),
             ("nucleosome_signal","nucleosome_signal")]:
    metric_panel(m, y)

pdf.close()
print("[INFO] Consolidated PDF:", pdf_path)

# ---------- Scoring (Neutral / SeekGene-lean) ----------
def norm_higher_better(x):
    x = x.astype(float)
    lo, hi = np.nanmin(x), np.nanmax(x)
    if not np.isfinite(lo) or not np.isfinite(hi) or hi==lo:
        return np.full_like(x, np.nan, dtype=float)
    return 100*(x - lo)/(hi - lo)

def norm_lower_better(x):
    x = x.astype(float)
    lo, hi = np.nanmin(x), np.nanmax(x)
    if not np.isfinite(lo) or not np.isfinite(hi) or hi==lo:
        return np.full_like(x, np.nan, dtype=float)
    return 100*(hi - x)/(hi - lo)

# Read back each dataset's QC table and take the median
def load_medians(outdir, name):
    df = pd.read_csv(os.path.join(outdir, f"{name}.qc_cells.csv"), index_col=0)
    med = {}
    for col in ["n_fragment","tsse","frip","frac_mito","frac_dup","doublet_score","nucleosome_signal"]:
        med[col] = float(np.nanmedian(df[col])) if col in df.columns else np.nan
    # Pass rate
    dfp = pd.read_csv(os.path.join(outdir, f"{name}.qc_cells.with_pass.csv"), index_col=0)
    med["pass_rate"] = float(np.nanmean(dfp["qc_pass"])) if "qc_pass" in dfp.columns else np.nan
    med["log10_nfrag"] = np.log10(max(1.0, med["n_fragment"])) if np.isfinite(med["n_fragment"]) else np.nan
    return med

table = pd.DataFrame({n: load_medians(OUTDIR, n) for n in DATASETS}).T

scores = pd.DataFrame(index=table.index)
scores["S_tsse"]      = norm_higher_better(table["tsse"])
scores["S_frip"]      = norm_higher_better(table["frip"])
scores["S_nfrag"]     = norm_higher_better(table["log10_nfrag"])
scores["S_mito"]      = norm_lower_better(table["frac_mito"])
scores["S_dup"]       = norm_lower_better(table["frac_dup"])
scores["S_doublet"]   = norm_lower_better(table["doublet_score"])   # may be all NaN (if no tile/bin)
scores["S_passrate"]  = norm_higher_better(table["pass_rate"])

NEUTRAL = {
    "S_tsse":0.30, "S_frip":0.25, "S_nfrag":0.20,
    "S_mito":0.10, "S_dup":0.10, "S_doublet":0.03, "S_passrate":0.02
}
SEEKGENE_LEAN = {
    "S_tsse":0.35, "S_frip":0.30, "S_nfrag":0.20,
    "S_mito":0.07, "S_dup":0.05, "S_doublet":0.02, "S_passrate":0.01
}

def weighted_total(S, W):
    w = pd.Series(W, index=S.columns, dtype=float)
    def row_total(r):
        valid = ~r.isna()
        if not valid.any(): return np.nan
        ww = w[valid]; ww = ww/ww.sum()
        return float((r[valid]*ww).sum())
    return S.apply(row_total, axis=1)

scores["Total_neutral"]  = weighted_total(scores, NEUTRAL)
scores["Total_seekgene"] = weighted_total(scores, SEEKGENE_LEAN)

# Output the scoring table and rankings
scores_out = os.path.join(OUTDIR, "benchmark_scoring.csv")
pd.concat([table, scores], axis=1).to_csv(scores_out)
print("[INFO] Scoring table ->", scores_out)

rank_neutral  = scores["Total_neutral"].sort_values(ascending=False)
rank_seekgene = scores["Total_seekgene"].sort_values(ascending=False)

def barplot(series, title, png):
    plt.figure(figsize=(6,3.6))
    series.plot(kind="barh", color="#446ecf")
    plt.gca().invert_yaxis()
    plt.xlabel("Score (0–100)"); plt.title(title); plt.tight_layout()
    plt.savefig(png, dpi=150); plt.close()

barplot(rank_neutral,  "Benchmark Ranking — Neutral",      os.path.join(OUTDIR, "ranking_neutral.png"))
barplot(rank_seekgene, "Benchmark Ranking — SeekGene-lean", os.path.join(OUTDIR, "ranking_seekgene.png"))

print("\n=== Ranking (Neutral) ===")
print(rank_neutral)
print("\n=== Ranking (SeekGene-lean) ===")
print(rank_seekgene)
print("[INFO] Plots:",
      os.path.join(OUTDIR, "ranking_neutral.png"), ",",
      os.path.join(OUTDIR, "ranking_seekgene.png"))


# In[5]:


os.chdir("/path/to/benchmark_scATAC_hdf5/")


# In[6]:


# ===== Make barplots from benchmark_summary.csv & benchmark_scoring.csv, save to a multi-page PDF =====
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# ---- File paths (if not in the current directory, change to absolute paths) ----
SUMMARY_CSV = "benchmark_summary.csv"
SCORING_CSV = "benchmark_scoring.csv"
OUT_PDF     = "benchmark_bars.pdf"

# ---- Read ----
summary = pd.read_csv(SUMMARY_CSV)
scoring = pd.read_csv(SCORING_CSV, index_col=0)

# Keep dataset order consistent
order = list(summary["dataset"])
scoring = scoring.reindex(order)

# Utility functions
def _to_pct(x):
    """Convert a 0-1 ratio to a percentage value (keep NaN if NaN)"""
    x = x.astype(float)
    return x * 100.0

def barplot_series(ax, series, title, ylabel, is_percent=False, rotate=30):
    s = series.astype(float)
    if is_percent:
        s = _to_pct(s)
        ylabel = ylabel if ylabel else "Percent (%)"
    s = s.replace([np.inf, -np.inf], np.nan)
    s.plot(kind="bar", ax=ax)
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.set_xlabel("")
    ax.set_xticklabels(s.index, rotation=rotate, ha="right")
    # Annotate values on the bars
    for p in ax.patches:
        val = p.get_height()
        if np.isfinite(val):
            ax.annotate(f"{val:.2f}", (p.get_x()+p.get_width()/2, val),
                        xytext=(0, 3), textcoords="offset points", ha="center", va="bottom", fontsize=8)
    ax.margins(x=0.02)

# ---- Metrics to plot (by page) ----
# From summary: median/pass rate
pages = []

# 1) n_fragment_median (also include log10 median for clarity)
if "n_fragment_median" in summary.columns:
    s_nf = pd.Series(summary["n_fragment_median"].values, index=summary["dataset"], name="n_fragment_median")
    pages.append(("Median Unique Fragments", s_nf, "Fragments (median)", False))
    # log10
    s_nf_log = np.log10(np.clip(s_nf, 1, None))
    pages.append(("Median log10(Unique Fragments)", s_nf_log, "log10(Fragments)", False))

# 2) tsse_median
if "tsse_median" in summary.columns:
    s_tsse = pd.Series(summary["tsse_median"].values, index=summary["dataset"], name="tsse_median")
    pages.append(("Median TSSe", s_tsse, "TSSe (median)", False))

# 3) frac_mito_median / frac_dup_median (displayed as percentages)
if "frac_mito_median" in summary.columns:
    s_mito = pd.Series(summary["frac_mito_median"].values, index=summary["dataset"], name="frac_mito_median")
    pages.append(("Median Mito Fraction", s_mito, "Mito (%)", True))
if "frac_dup_median" in summary.columns:
    s_dup = pd.Series(summary["frac_dup_median"].values, index=summary["dataset"], name="frac_dup_median")
    pages.append(("Median Duplicate Fraction", s_dup, "Duplicates (%)", True))

# 4) pass_rate (percentage)
if "pass_rate" in summary.columns:
    s_pass = pd.Series(summary["pass_rate"].values, index=summary["dataset"], name="pass_rate")
    pages.append(("Pass Rate", s_pass, "Pass (%)", True))

# 5) doublet_score_median / nucleosome_signal_median (if present)
if "doublet_score_median" in summary.columns:
    s_dbl = pd.Series(summary["doublet_score_median"].values, index=summary["dataset"], name="doublet_score_median")
    if not s_dbl.isna().all():
        pages.append(("Median Doublet Score", s_dbl, "Doublet score (median)", False))
if "nucleosome_signal_median" in summary.columns:
    s_nuc = pd.Series(summary["nucleosome_signal_median"].values, index=summary["dataset"], name="nucleosome_signal_median")
    if not s_nuc.isna().all():
        pages.append(("Median Nucleosome Signal (mono/NFR)", s_nuc, "Nucleosome signal", False))

# From scoring: total score and individual subscores
if "Total_neutral" in scoring.columns:
    s_neut = scoring["Total_neutral"]
    pages.append(("Total Score — Neutral", s_neut, "Score (0–100)", False))
if "Total_seekgene" in scoring.columns:
    s_seek = scoring["Total_seekgene"]
    pages.append(("Total Score — SeekGene-lean", s_seek, "Score (0–100)", False))

# Optional: individual subscores (if you want to see the composition)
for col, title in [
    ("S_tsse",     "Subscore: TSSe"),
    ("S_frip",     "Subscore: FRiP"),
    ("S_nfrag",    "Subscore: log10(nFrags)"),
    ("S_mito",     "Subscore: Mito (lower better)"),
    ("S_dup",      "Subscore: Duplicate (lower better)"),
    ("S_doublet",  "Subscore: Doublet (lower better)"),
    ("S_passrate", "Subscore: Pass rate"),
]:
    if col in scoring.columns and not scoring[col].isna().all():
        pages.append((title, scoring[col], "Score (0–100)", False))

# ---- Generate PDF ----
with PdfPages(OUT_PDF) as pdf:
    # Cover page: simple overview
    fig = plt.figure(figsize=(8.0, 4.8))
    plt.axis("off")
    txt = f"scATAC Benchmark Barplots\n\nSummary file: {os.path.abspath(SUMMARY_CSV)}\nScoring file: {os.path.abspath(SCORING_CSV)}\nDatasets: {', '.join(order)}"
    plt.text(0.02, 0.8, txt, fontsize=12, va="top")
    pdf.savefig(fig); plt.close(fig)

    # One metric per page
    for title, series, ylabel, is_percent in pages:
        fig, ax = plt.subplots(figsize=(8.0, 4.8))
        barplot_series(ax, series.reindex(order), title, ylabel, is_percent=is_percent)
        pdf.savefig(fig); plt.close(fig)

print(f"[INFO] Saved barplots PDF -> {os.path.abspath(OUT_PDF)}")

# Also print both rankings (to screen)
def print_ranking(series, name):
    s = series.sort_values(ascending=False)
    print(f"\n=== Ranking: {name} ===")
    print(s)

if "Total_neutral" in scoring.columns:
    print_ranking(scoring["Total_neutral"], "Total_neutral")
if "Total_seekgene" in scoring.columns:
    print_ranking(scoring["Total_seekgene"], "Total_seekgene")


# In[9]:


# ===================== scATAC Benchmark — Extras (Triplet FragLen / FRiP with Union Peaks / Pseudobulk Repro) =====================
# Dependencies: snapatac2, pandas, numpy, matplotlib; (optional) bedtools/macs2 (if you build a union peak set via CLI)
# Prerequisite: you have run the "Multi-tech scATAC Benchmark + Scoring" main script and produced *.qc_cells.csv and raw .h5ad under OUTDIR

# ---------- CONFIG ----------
BIG_TMP  = "/path/to/tmp"
OUTDIR   = "/path/to/benchmark_scATAC_hdf5"
GENOME   = "mm10"
# Union peak set (strongly recommended to provide a union BED; if not, see the CLI comments below on building one with macs2+bedtools)
UNION_PEAKS = "/path/to/benchmark_scATAC_hdf5/union_peaks.bed"   # <- change to your path

# Datasets (same as the main script; Droplet removed)
DATASETS = {
    "10x_multiome": "/path/to/scATAC/brain/10xmultiome/M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_atac_fragments.tsv.gz",
    "10x_scatac"  : "/path/to/scATAC/brain/10xscatac/atac_v1_E18_brain_flash_5k_fragments.tsv.gz",
    "SeekGene"    : "/path/to/scATAC/brain/seekgene/atac_fragments.tsv.gz",
    "TXCI"        : "/path/to/scATAC/brain/txci-atac/GSM7852211_mm10.merged.fragments.tsv.gz",
}
FRIP_THRESH = 0.20             # FRiP threshold (uniform-threshold filtering)
BIN_SIZE    = 5000             # pseudobulk tile/bin size (if API available)

# ---------- Imports ----------
import os, re, gzip, glob, gc, warnings
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

warnings.filterwarnings("ignore", message=".*_import_from_c.*")
import snapatac2 as snap

# ---------- Utilities ----------
def latest_h5ad_for(name, outdir):
    patt = os.path.join(outdir, f"{name}.raw.*.h5ad")
    files = sorted(glob.glob(patt))
    return files[-1] if files else None

def frag_len_histogram(frag_path, max_len=800, step=1):
    bins=np.zeros(max_len+1, dtype=np.int64)  # 0 collects overflow
    with gzip.open(frag_path,"rt") as fh:
        for ln in fh:
            if not ln or ln[0]=="#": continue
            p=ln.rstrip("\n").split("\t")
            if len(p)<3: continue
            try: L=int(p[2])-int(p[1])
            except: continue
            if 0 < L <= max_len: bins[L]+=1
            else: bins[0]+=1
    idx=np.arange(1,max_len+1,step)
    val=bins[1:max_len+1:step].astype(float)
    tot=val.sum(); frac=val/(tot if tot>0 else 1.0)
    return idx, frac, int(bins[0])

def plot_triplet_length_panels(datasets, out_pdf):
    with PdfPages(out_pdf) as pdf:
        for name, frag in datasets.items():
            x, frac, overflow = frag_len_histogram(frag, max_len=800, step=1)
            # NFR/mono/di bins
            nfr_m = (x<147); mono_m=(x>=147)&(x<294); di_m=(x>=294)&(x<441)
            f_nfr, f_mono, f_di = frac[nfr_m].sum(), frac[mono_m].sum(), frac[di_m].sum()
            # Plot
            plt.figure(figsize=(7.5,4.2))
            plt.plot(x, frac, lw=1.2)
            for thr,lab in [(147,"147bp"),(294,"294bp"),(441,"441bp")]:
                plt.axvline(thr, ls="--", c="grey", lw=0.8)
                plt.text(thr+5, max(frac)*0.9, lab, fontsize=8, rotation=90, va="top", color="grey")
            plt.title(f"{name} — Fragment Length Triplet\nNFR={f_nfr:.2%}  Mono={f_mono:.2%}  Di={f_di:.2%}  (overflow>{800}bp:{overflow:,})")
            plt.xlabel("Fragment size (bp)"); plt.ylabel("Fraction")
            plt.tight_layout(); pdf.savefig(); plt.close()

# ---------- 1) Triplet FragLen (NFR/Mono/Di peak shapes & proportions) ----------
triplet_pdf = os.path.join(OUTDIR, "fraglen_triplet_panels.pdf")
plot_triplet_length_panels(DATASETS, triplet_pdf)
print("[INFO] FragLen triplet panels ->", triplet_pdf)

# ---------- 2) FRiP with Union Peaks (per-cell FRiP + uniform threshold) ----------
if not (UNION_PEAKS and os.path.exists(UNION_PEAKS)):
    print("[WARN] UNION_PEAKS not set/found; skip FRiP. You can provide a union peaks.bed and rerun this section.")
else:
    for name, frag in DATASETS.items():
        h5 = latest_h5ad_for(name, OUTDIR)
        assert h5, f"No h5ad found for {name} under {OUTDIR}"
        print(f"[INFO] FRiP — {name} ({h5})")
        ad = snap.read(h5) if hasattr(snap, "read") else snap.pp.load_file(h5)  # compatible with different API versions
        try:
            snap.metrics.frip(ad, {"frip": UNION_PEAKS})
        except TypeError:
            snap.metrics.frip(ad, UNION_PEAKS)
        # Write back to the QC table (append FRiP column) and update the pass mask by uniform threshold
        qc_path = os.path.join(OUTDIR, f"{name}.qc_cells.csv")
        df = pd.read_csv(qc_path, index_col=0)
        frip = np.asarray(ad.obs["frip"], dtype=float)
        df["frip"] = frip
        # Pass-rate file (if it exists)
        pass_path = os.path.join(OUTDIR, f"{name}.qc_cells.with_pass.csv")
        if os.path.exists(pass_path):
            dpf = pd.read_csv(pass_path, index_col=0)
        else:
            dpf = df.copy()
        # Recompute pass: reuse previous thresholds + FRiP threshold
        nfr  = df["n_fragment"].to_numpy(float)
        tsse = df["tsse"].to_numpy(float)
        mito = df["frac_mito"].to_numpy(float) if "frac_mito" in df else np.full_like(tsse, np.nan)
        dup  = df["frac_dup"].to_numpy(float)  if "frac_dup"  in df else np.full_like(tsse, np.nan)
        mask = (nfr>=MIN_COUNTS) & (nfr<=MAX_COUNTS) & (tsse>=MIN_TSSE)
        if not np.isnan(mito).all(): mask &= (mito<=MAX_MITO)
        if not np.isnan(dup).all():  mask &= (dup<=MAX_DUP)
        mask &= (frip>=FRIP_THRESH)
        dpf["qc_pass"] = mask
        # Save
        df.to_csv(qc_path)
        dpf.to_csv(pass_path)
        # Emit a FRiP histogram page (alongside the main report)
        plt.figure(figsize=(6.5,4))
        v = df["frip"].to_numpy(float)
        plt.hist(v[~np.isnan(v)], bins=60)
        plt.axvline(FRIP_THRESH, ls="--", c="red", label=f"FRiP ≥ {FRIP_THRESH}")
        plt.title(f"{name} — FRiP histogram")
        plt.xlabel("FRiP"); plt.ylabel("count"); plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(OUTDIR, f"{name}.frip_hist.png"), dpi=150)
        plt.close()
    print("[INFO] FRiP per-cell computed & qc_pass updated (threshold=%.2f)." % FRIP_THRESH)

# ---------- 3) Reproducibility (pseudobulk Spearman heatmap; based on the union peak set) ----------
# Note: if the current SnapATAC2 version supports make_peak_matrix/make_bin_matrix (similar to the tile/bin compatibility above), we build a peak matrix per dataset,
#       then summarize cells into pseudobulk (sum across cells), and build a Spearman correlation matrix across all datasets.
def build_peak_matrix_compat(ad, peak_bed, out_path):
    """Peak matrix construction compatible with different APIs; returns (AnnData_or_None, success)"""
    if hasattr(snap.pp, "make_peak_matrix"):
        try:
            return snap.pp.make_peak_matrix(ad, file=out_path, peak_file=peak_bed), True
        except TypeError:
            return snap.pp.make_peak_matrix(ad, file=out_path, peaks=peak_bed), True
        except Exception:
            return None, False
    elif hasattr(snap.pp, "make_bin_matrix"):
        # Older versions may have no peak matrix, only bin; in that case skip reproducibility assessment
        return None, False
    return None, False

if UNION_PEAKS and os.path.exists(UNION_PEAKS):
    pseudo = {}   # name -> 1D counts array over union peaks
    peak_n = None
    for name,_ in DATASETS.items():
        h5 = latest_h5ad_for(name, OUTDIR)
        assert h5, f"No h5ad found for {name}"
        print(f"[INFO] Peak matrix (union) — {name}")
        ad = snap.read(h5) if hasattr(snap, "read") else snap.pp.load_file(h5)
        peak_path = os.path.join(OUTDIR, f"{name}.union_peak_matrix.h5ad")
        pm, ok = build_peak_matrix_compat(ad, UNION_PEAKS, peak_path)
        if not ok or pm is None:
            print(f"[WARN] make_peak_matrix unavailable/failed for {name}; skip reproducibility.")
            pseudo = {}
            break
        X = pm.X
        if hasattr(X, "getnnz"):
            c = np.array(X.sum(axis=0)).ravel()   # pseudobulk (sum over peaks)
        else:
            c = np.asarray(X).sum(axis=0).ravel()
        pseudo[name] = c
        peak_n = len(c)
        # Cleanup
        if hasattr(pm, "close"):
            try: pm.close()
            except Exception: pass
        del pm, X; gc.collect()

    if pseudo:
        df = pd.DataFrame(pseudo)    # rows=peaks, cols=samples
        # Spearman correlation matrix
        corr = df.corr(method="spearman")
        # Output heatmap
        heat_pdf = os.path.join(OUTDIR, "pseudobulk_repro_spearman.pdf")
        with PdfPages(heat_pdf) as pdfh:
            plt.figure(figsize=(5.2,4.8))
            im = plt.imshow(corr.values, cmap="viridis", vmin=0, vmax=1)
            plt.colorbar(im, fraction=0.046, pad=0.04, label="Spearman ρ")
            plt.xticks(range(len(corr.columns)), corr.columns, rotation=30, ha="right")
            plt.yticks(range(len(corr.index)), corr.index)
            plt.title("Pseudobulk reproducibility (Spearman)")
            plt.tight_layout(); pdfh.savefig(); plt.close()
        corr.to_csv(os.path.join(OUTDIR, "pseudobulk_spearman_matrix.csv"))
        print("[INFO] Pseudobulk Spearman heatmap ->", heat_pdf)
    else:
        print("[WARN] Repro heatmap skipped (no peak-matrix API).")
else:
    print("[WARN] UNION_PEAKS not found; reproducibility skipped.")

# ---------- Scoring (based on the existing scoring file; auto-updates if FRiP was added) ----------
scoring_csv = os.path.join(OUTDIR, "benchmark_scoring.csv")
if os.path.exists(scoring_csv):
    scoring = pd.read_csv(scoring_csv, index_col=0)
    # If FRiP/pass-rate updates were newly added, recompute the SeekGene-lean total score (same for Neutral)
    def norm_higher_better(x):
        x = x.astype(float); 
        if x.isna().all(): return x
        lo, hi = np.nanmin(x), np.nanmax(x); 
        return 100*(x-lo)/(hi-lo) if hi>lo else np.full_like(x, np.nan, dtype=float)
    def norm_lower_better(x):
        x = x.astype(float); 
        if x.isna().all(): return x
        lo, hi = np.nanmin(x), np.nanmax(x); 
        return 100*(hi-x)/(hi-lo) if hi>lo else np.full_like(x, np.nan, dtype=float)

    # Read medians from the summary/scoring files
    summary = pd.read_csv(os.path.join(OUTDIR, "benchmark_summary.csv"))
    order = list(summary["dataset"])
    # Recompute medians from the latest qc_cells.csv (especially FRiP and pass_rate)
    def load_medians(name):
        df = pd.read_csv(os.path.join(OUTDIR, f"{name}.qc_cells.csv"), index_col=0)
        med = {}
        for col in ["n_fragment","tsse","frip","frac_mito","frac_dup","doublet_score","nucleosome_signal"]:
            med[col]=float(np.nanmedian(df[col])) if col in df.columns else np.nan
        dpf = pd.read_csv(os.path.join(OUTDIR, f"{name}.qc_cells.with_pass.csv"), index_col=0)
        med["pass_rate"]=float(np.nanmean(dpf["qc_pass"])) if "qc_pass" in dpf.columns else np.nan
        med["log10_nfrag"]=np.log10(max(1.0, med["n_fragment"])) if np.isfinite(med["n_fragment"]) else np.nan
        return med
    table = pd.DataFrame({n: load_medians(n) for n in order}).T

    S = pd.DataFrame(index=table.index)
    S["S_tsse"]     = norm_higher_better(table["tsse"])
    S["S_frip"]     = norm_higher_better(table["frip"])
    S["S_nfrag"]    = norm_higher_better(table["log10_nfrag"])
    S["S_mito"]     = norm_lower_better(table["frac_mito"])
    S["S_dup"]      = norm_lower_better(table["frac_dup"])
    S["S_doublet"]  = norm_lower_better(table["doublet_score"])
    S["S_passrate"] = norm_higher_better(table["pass_rate"])

    NEUTRAL = {"S_tsse":0.30, "S_frip":0.25, "S_nfrag":0.20,
               "S_mito":0.10, "S_dup":0.10, "S_doublet":0.03, "S_passrate":0.02}
    SEEKGENE_LEAN = {"S_tsse":0.35, "S_frip":0.30, "S_nfrag":0.20,
                     "S_mito":0.07, "S_dup":0.05, "S_doublet":0.02, "S_passrate":0.01}

    def weighted_total(S, W):
        w = pd.Series(W, index=S.columns, dtype=float)
        def row_total(r):
            valid = ~r.isna()
            if not valid.any(): return np.nan
            ww = w[valid]; ww = ww/ww.sum()
            return float((r[valid]*ww).sum())
        return S.apply(row_total, axis=1)

    S["Total_neutral"]  = weighted_total(S, NEUTRAL)
    S["Total_seekgene"] = weighted_total(S, SEEKGENE_LEAN)

    out = pd.concat([table, S], axis=1)
    out.to_csv(scoring_csv)
    print("[INFO] Updated scoring ->", scoring_csv)

    # Two ranking plots
    def barplot(series, title, png):
        plt.figure(figsize=(6,3.6))
        series.sort_values(ascending=False).plot(kind="barh", color="#446ecf")
        plt.gca().invert_yaxis()
        plt.xlabel("Score (0–100)"); plt.title(title); plt.tight_layout()
        plt.savefig(png, dpi=150); plt.close()
    barplot(S["Total_neutral"],  "Benchmark Ranking — Neutral",
            os.path.join(OUTDIR, "ranking_neutral.updated.png"))
    barplot(S["Total_seekgene"], "Benchmark Ranking — SeekGene-lean",
            os.path.join(OUTDIR, "ranking_seekgene.updated.png"))
    print("[INFO] Rankings updated:",
          os.path.join(OUTDIR, "ranking_neutral.updated.png"), ",",
          os.path.join(OUTDIR, "ranking_seekgene.updated.png"))

# ---------- Notes ----------
# If you need to build a "union peak set":
#  A) If you already have multiple sets of peaks (.bed or narrowPeak), merge them with bedtools:
#     cat peaks1.bed peaks2.bed ... | sort -k1,1 -k2,2n | bedtools merge > union_peaks.bed
#  B) If you have no peaks, you can do pseudobulk + macs2 peak calling per dataset (may be large and slow):
#     zcat sample.fragments.tsv.gz | cut -f1-3 > sample.pseudo.bed
#     macs2 callpeak -t sample.pseudo.bed -f BED -g mm -n sample -B --nomodel --shift -100 --extsize 200 -q 0.01
#     then merge all sample_peaks.narrowPeak (same as A)
#
#  IDR (two replicates):
#    After macs2 peak calling on the two replicates of the same technology to obtain narrowPeak, use the idr tool:
#    idr --samples rep1.narrowPeak rep2.narrowPeak --peak-list union_peaks.bed --plot \
#        --output-file idr.txt --log-output-file idr.log
#    Peaks with IDR<0.05 in the result are considered stably reproducible; you can compute the proportion / plot a curve.


# In[2]:


# ===================== Multi-tech scATAC Benchmark + Extras + Scoring (SeekGene-lean option) =====================
# Four datasets: 10x multiome / 10x scATAC / SeekGene / TXCI (Droplet removed)
# End-to-end: fragments import -> TSSe (local GENCODE vM31) -> (if supported) tile/bin+Scrublet -> fragment-length triple peaks
#        (optional) union peak set FRiP -> (optional) pseudobulk reproducibility -> summary/scoring (FAIR & DISPLAY) -> ranking plots & PDF
# Note: this script is compatible with different SnapATAC2 API versions (make_tile_matrix/make_bin_matrix/peak_matrix may not exist)

# ---------- CONFIG ----------
BIG_TMP  = "/path/to/tmp"                           # large-disk TMP
OUTDIR   = "/path/to/benchmark_scATAC_hdf5"     # output directory
GENOME   = "mm10"
GENE_ANNOT = "/path/to/gencode.vM31.basic.annotation.gtf.gz"

# Four datasets (Droplet removed)
DATASETS = {
    "10x_multiome": "/path/to/scATAC/brain/10xmultiome/M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_atac_fragments.tsv.gz",
    "10x_scatac"  : "/path/to/scATAC/brain/10xscatac/atac_v1_E18_brain_flash_5k_fragments.tsv.gz",
    "SeekGene"    : "/path/to/scATAC/brain/seekgene/atac_fragments.tsv.gz",
    "TXCI"        : "/path/to/scATAC/brain/txci-atac/GSM7852211_mm10.merged.fragments.tsv.gz",
}

# (optional) union peak set (if absent, set to ""/None; this script auto-skips FRiP & reproducibility)
UNION_PEAKS = ""  # e.g.: "/path/to/benchmark_scATAC_hdf5/union_peaks.bed"

# Uniform thresholds (closer to the scale of your current data; can be fine-tuned later based on results)
MIN_COUNTS, MAX_COUNTS = 1000, 200000
MIN_TSSE_FOR_PASS      = 1.0
FRIP_THRESH            = 0.20
MAX_MITO               = 0.20

# Import parameters (HDF5 only; your environment does not support zarr)
BACKEND          = "hdf5"
CHUNK_SIZE       = 1_000_000
N_JOBS_IMPORT    = 4
MIN_FRAGS_IMPORT = max(500, MIN_COUNTS)

# Additional analysis
DO_NUCLEOSOME   = True
NUC_MAX_LINES   = None
SCRUBLET_JOBS   = 4
DOWNSAMPLE_PLOT = 10_000

# ---------- ENV ----------
import os, sys, re, gzip, time, gc, warnings, glob
from pathlib import Path
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
os.environ.setdefault("HDF5_DISABLE_VERSION_CHECK", "2")
os.environ["TMPDIR"] = os.environ["TMP"] = os.environ["TEMP"] = BIG_TMP
Path(BIG_TMP).mkdir(parents=True, exist_ok=True)
Path(OUTDIR).mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", os.path.join(BIG_TMP, "mpl"))

# ---------- DEPS ----------
import numpy as np, pandas as pd, h5py
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
warnings.filterwarnings("ignore", message=".*_import_from_c.*")
try:
    import snapatac2 as snap
except Exception as e:
    print("ERROR: snapatac2 not installed. Conda: conda create -n snap2 -y python=3.10 snapatac2 h5py pandas matplotlib", file=sys.stderr); raise

# ---------- HELPERS ----------
def ensure_tbi_for_frag(frag_path: str):
    tbi = frag_path + ".tbi"
    tbi_gz = tbi + ".gz"
    if os.path.exists(tbi): return
    if os.path.exists(tbi_gz):
        print(f"[INFO] Decompressing {tbi_gz} -> {tbi}")
        with gzip.open(tbi_gz, "rb") as fi, open(tbi, "wb") as fo: fo.write(fi.read())
    else:
        print(f"[WARN] No .tbi alongside {frag_path} (SnapATAC2 may index at TMP).")

def close_all_handles():
    for v in list(globals().values()):
        if hasattr(v, "close"):
            try: v.close()
            except Exception: pass
    gc.collect()

def unique_store(base, name):
    ts=time.strftime("%Y%m%d-%H%M%S"); pid=os.getpid()
    return os.path.join(base, f"{name}.raw.{ts}.{pid}.h5ad")

def frag_len_histogram(frag_path, max_len=800, step=1):
    bins=np.zeros(max_len+1, dtype=np.int64)  # 0 collects overflow
    with gzip.open(frag_path,"rt") as fh:
        for ln in fh:
            if not ln or ln[0]=="#": continue
            p=ln.rstrip("\n").split("\t")
            if len(p)<3: continue
            try: L=int(p[2])-int(p[1])
            except: continue
            if 0 < L <= max_len: bins[L]+=1
            else: bins[0]+=1
    idx=np.arange(1,max_len+1,step)
    val=bins[1:max_len+1:step].astype(float)
    tot=val.sum(); frac=val/(tot if tot>0 else 1.0)
    return idx, frac, int(bins[0])

def compute_nucleosome_signal_by_streaming(frag_path, keep_barcodes=None,
                                           nfr_max=147, mono_min=147, mono_max=294, max_lines=None):
    nfr, mono = {}, {}
    i=0
    with gzip.open(frag_path, "rt") as fh:
        for line in fh:
            if max_lines is not None and i>=max_lines: break
            i+=1
            if not line or line[0]=="#": continue
            p=line.rstrip("\n").split("\t")
            if len(p)<4: continue
            try: L = int(p[2]) - int(p[1])
            except: continue
            bc=p[3]
            if keep_barcodes is not None and bc not in keep_barcodes: continue
            if L < nfr_max:               nfr[bc]  = nfr.get(bc,0)  + 1
            elif mono_min <= L < mono_max: mono[bc] = mono.get(bc,0) + 1
    ks = (keep_barcodes if keep_barcodes is not None else set(list(nfr.keys())+list(mono.keys())))
    return {bc: (mono.get(bc,0)/(nfr.get(bc,0) if nfr.get(bc,0)>0 else 1.0)) for bc in ks}

def summarize_series(x):
    x = x[~np.isnan(x)]
    if x.size==0: return dict(n=0, median=np.nan, q1=np.nan, q3=np.nan)
    q = np.quantile(x, [0.25, 0.5, 0.75])
    return dict(n=x.size, median=q[1], q1=q[0], q3=q[2])

def latest_h5ad_for(name, outdir):
    patt = os.path.join(outdir, f"{name}.raw.*.h5ad")
    files = sorted(glob.glob(patt))
    return files[-1] if files else None

def run_tsse_safe(ad, gene_anno_path, genome_obj):
    """Prefer local GTF/GFF; else write NaN TSSe and continue (avoid cache-download failures)"""
    try:
        if gene_anno_path and os.path.exists(gene_anno_path):
            snap.metrics.tsse(ad, gene_anno_path)
            return True
        raise FileNotFoundError(f"Local annotation not found: {gene_anno_path}")
    except Exception as e:
        print("[WARN] TSSe failed:", e)
        ad.obs["tsse"] = np.full(ad.n_obs, np.nan, dtype=float)
        return False

def build_tile_like_matrix(ad, out_path, bin_size=5000):
    """Compatible with different APIs: pp.make_tile_matrix -> pp.make_bin_matrix -> None"""
    if hasattr(snap.pp, "make_tile_matrix"):
        try:    return snap.pp.make_tile_matrix(ad, file=out_path, bin_size=bin_size, backend="hdf5"), True
        except TypeError: return snap.pp.make_tile_matrix(ad, file=out_path, bin_size=bin_size), True
    if hasattr(snap.pp, "make_bin_matrix"):
        try:    return snap.pp.make_bin_matrix(ad, file=out_path, bin_size=bin_size, backend="hdf5"), True
        except TypeError: return snap.pp.make_bin_matrix(ad, file=out_path, bin_size=bin_size), True
    print("[WARN] No tile/bin matrix API; skipping doublets.")
    return None, False

def build_peak_matrix_compat(ad, peak_bed, out_path):
    """Peak matrix (returns None, False if the version does not support it)"""
    if hasattr(snap.pp, "make_peak_matrix"):
        try:
            return snap.pp.make_peak_matrix(ad, file=out_path, peak_file=peak_bed), True
        except TypeError:
            return snap.pp.make_peak_matrix(ad, file=out_path, peaks=peak_bed), True
        except Exception:
            return None, False
    return None, False

# ---------- BENCHMARK LOOP ----------
summary_rows = []
pdf_path = os.path.join(OUTDIR, "benchmark_full_report.pdf")
pdf = PdfPages(pdf_path)

genome_obj = snap.genome.mm10 if GENOME.lower()=="mm10" else getattr(snap.genome, GENOME.lower())
assert os.path.exists(GENE_ANNOT), f"Annotation not found: {GENE_ANNOT}"

for name, frag in DATASETS.items():
    print(f"\n=== {name} ===")
    assert os.path.exists(frag), f"{frag} not found"
    ensure_tbi_for_frag(frag)

    # Import
    close_all_handles()
    store = unique_store(OUTDIR, name)
    print("[INFO] import_fragments …", store)
    ad = snap.pp.import_fragments(
        frag,
        chrom_sizes = genome_obj,
        file        = store,
        backend     = "hdf5",
        sorted_by_barcode = False,
        whitelist   = None,
        tempdir     = BIG_TMP,
        chunk_size  = CHUNK_SIZE,
        n_jobs      = N_JOBS_IMPORT,
        min_num_fragments = MIN_FRAGS_IMPORT,
    )

    # TSSe
    print("[INFO] TSSe …")
    _ = run_tsse_safe(ad, GENE_ANNOT, genome_obj)

    # Tile/bin + Scrublet (if available)
    print("[INFO] Tile 5kb + Scrublet …")
    tile_path = os.path.join(OUTDIR, f"{name}.tile_5kb.h5ad")
    tile, has_mat = build_tile_like_matrix(ad, tile_path, bin_size=5000)
    dbl = np.full(ad.n_obs, np.nan)
    if has_mat:
        try:
            snap.pp.scrublet(tile, n_jobs=SCRUBLET_JOBS)
            dbl = np.asarray(tile.obs["doublet_score"])
        except Exception as e:
            print("[WARN] scrublet failed:", e)

    # Nucleosome signal
    if DO_NUCLEOSOME:
        nuc = compute_nucleosome_signal_by_streaming(frag, keep_barcodes=set(ad.obs_names), max_lines=NUC_MAX_LINES)
        nuc = pd.Series(nuc).reindex(ad.obs_names).astype(float).values
    else:
        nuc = np.full(ad.n_obs, np.nan)

    # Per-cell QC & pass/fail (use 1.0 for the TSSe threshold; dup not included in the pass decision)
    qc = pd.DataFrame(index=ad.obs_names)
    for k in ["n_fragment","frac_dup","frac_mito","tsse"]:
        if k in ad.obs: qc[k] = np.asarray(ad.obs[k])
    qc["doublet_score"]      = dbl
    qc["nucleosome_signal"]  = nuc

    nfr  = qc["n_fragment"].to_numpy(float) if "n_fragment" in qc else np.full(ad.n_obs, np.nan)
    tsse = qc["tsse"].to_numpy(float)       if "tsse"       in qc else np.full(ad.n_obs, np.nan)
    mito = qc["frac_mito"].to_numpy(float)  if "frac_mito"  in qc else np.full(ad.n_obs, np.nan)

    mask = (nfr>=MIN_COUNTS) & (nfr<=MAX_COUNTS) & (tsse>=MIN_TSSE_FOR_PASS)
    if not np.isnan(mito).all(): mask &= (mito<=MAX_MITO)
    qc["qc_pass"] = mask

    qc_path  = os.path.join(OUTDIR, f"{name}.qc_cells.csv")
    pass_path= os.path.join(OUTDIR, f"{name}.qc_cells.with_pass.csv")
    qc.to_csv(qc_path); qc.to_csv(pass_path)

    # Summary row
    row = {
        "dataset": name,
        "cells": int(qc.shape[0]),
        "pass_cells": int(mask.sum()),
        "pass_rate": float(mask.mean()),
        "n_fragment_median": float(np.nanmedian(qc["n_fragment"])) if "n_fragment" in qc else np.nan,
        "tsse_median": float(np.nanmedian(qc["tsse"])) if "tsse" in qc else np.nan,
        "frac_mito_median": float(np.nanmedian(qc["frac_mito"])) if "frac_mito" in qc else np.nan,
        "frac_dup_median" : float(np.nanmedian(qc["frac_dup"]))  if "frac_dup"  in qc else np.nan,
        "doublet_score_median": float(np.nanmedian(qc["doublet_score"])) if "doublet_score" in qc else np.nan,
        "nucleosome_signal_median": float(np.nanmedian(qc["nucleosome_signal"])) if "nucleosome_signal" in qc else np.nan,
    }
    summary_rows.append(row)

    # ===== PDF: knee + length distribution (count/percent/triple peaks) + TSSe vs log10(nFrags) + histogram =====
    # 1) Knee
    vals = np.sort(nfr)[::-1]; ranks=np.arange(1, len(vals)+1)
    plt.figure(); plt.plot(ranks, vals)
    plt.xscale("log"); plt.yscale("log"); plt.xlabel("barcode rank (log10)"); plt.ylabel("n_fragment (log10)")
    plt.title(f"{name} — Knee"); plt.tight_layout(); pdf.savefig(); plt.close()

    # 2) Frag length
    ok_fsd=True
    try:
        fsd=snap.metrics.frag_size_distr(ad); ok_fsd=(fsd is not None)
    except Exception: ok_fsd=False
    if ok_fsd:
        x=np.arange(len(fsd)); y=np.array(fsd, dtype=float); x=x[1:]; y=y[1:]
        y_pct=y/(y.sum() if y.sum()>0 else 1.0)
        # counts
        plt.figure(); plt.plot(x,y); plt.yscale("log")
        plt.xlabel("fragment length (bp)"); plt.ylabel("count (log)")
        plt.title(f"{name} — Fragment length distribution"); plt.tight_layout(); pdf.savefig(); plt.close()
        # percent + triplet bands
        nfr_m=(x<147); mono_m=(x>=147)&(x<294); di_m=(x>=294)&(x<441)
        f_nfr, f_mono, f_di = y_pct[nfr_m].sum(), y_pct[mono_m].sum(), y_pct[di_m].sum()
        plt.figure(); plt.plot(x,y_pct, lw=1.2)
        for thr,lab in [(147,"147bp"),(294,"294bp"),(441,"441bp")]:
            plt.axvline(thr, ls="--", c="grey", lw=0.8)
            plt.text(thr+5, max(y_pct)*0.9, lab, fontsize=8, rotation=90, va="top", color="grey")
        plt.xlabel("Size of Fragments (bp)"); plt.ylabel("Fragments (%)")
        plt.title(f"{name} — Frag Size (percent)  NFR={f_nfr:.2%}  Mono={f_mono:.2%}  Di={f_di:.2%}")
        plt.tight_layout(); pdf.savefig(); plt.close()
    else:
        x,y_pct,overflow=frag_len_histogram(frag, max_len=800, step=1)
        nfr_m=(x<147); mono_m=(x>=147)&(x<294); di_m=(x>=294)&(x<441)
        f_nfr, f_mono, f_di = y_pct[nfr_m].sum(), y_pct[mono_m].sum(), y_pct[di_m].sum()
        plt.figure(); plt.plot(x,y_pct, lw=1.2)
        for thr,lab in [(147,"147bp"),(294,"294bp"),(441,"441bp")]:
            plt.axvline(thr, ls="--", c="grey", lw=0.8)
            plt.text(thr+5, max(y_pct)*0.9, lab, fontsize=8, rotation=90, va="top", color="grey")
        plt.xlabel("Size of Fragments (bp)"); plt.ylabel("Fragments (%)")
        plt.title(f"{name} — Frag Size (fallback)  NFR={f_nfr:.2%}  Mono={f_mono:.2%}  Di={f_di:.2%}  overflow>{800}bp:{overflow:,}")
        plt.tight_layout(); pdf.savefig(); plt.close()

    # 3) TSSe vs log10(nFrags)
    ds = min(DOWNSAMPLE_PLOT, len(nfr))
    idx = np.random.choice(len(nfr), ds, replace=False) if ds<len(nfr) else np.arange(len(nfr))
    plt.figure(); plt.scatter(np.log10(np.clip(nfr[idx],1,None)), tsse[idx], s=3, alpha=0.5)
    plt.axvline(np.log10(MIN_COUNTS), ls="--"); plt.axhline(MIN_TSSE_FOR_PASS, ls="--")
    plt.xlabel("log10(n_fragment)"); plt.ylabel("TSSe"); plt.title(f"{name} — TSSe vs log10(nFrags)")
    plt.tight_layout(); pdf.savefig(); plt.close()

    # 4) Histograms
    def _hist(col, title, xlabel):
        if col in qc:
            x=qc[col].to_numpy(float)
            plt.figure(); plt.hist(x[~np.isnan(x)], bins=60)
            plt.title(f"{name} — {title}"); plt.xlabel(xlabel); plt.ylabel("count")
            plt.tight_layout(); pdf.savefig(); plt.close()
    _hist("n_fragment","n_fragment","n_fragment")
    _hist("tsse","TSSe","TSSe")
    _hist("frac_mito","frac_mito","frac_mito")
    _hist("frac_dup","frac_dup","frac_dup")
    _hist("doublet_score","doublet_score","doublet_score")
    _hist("nucleosome_signal","nucleosome_signal","nucleosome_signal")

pdf.close()
print("[INFO] Full PDF ->", pdf_path)

# ---------- FRiP with union peaks (optional; only when UNION_PEAKS exists) ----------
has_union = bool(UNION_PEAKS) and os.path.exists(UNION_PEAKS)
if has_union:
    for name in DATASETS:
        h5 = latest_h5ad_for(name, OUTDIR)
        if not h5: continue
        print(f"[INFO] FRiP — {name}")
        ad = snap.read(h5) if hasattr(snap, "read") else snap.pp.load_file(h5)
        try:    snap.metrics.frip(ad, {"frip": UNION_PEAKS})
        except TypeError: snap.metrics.frip(ad, UNION_PEAKS)

        # Update FRiP & pass mask
        qc_path  = os.path.join(OUTDIR, f"{name}.qc_cells.csv")
        pass_path= os.path.join(OUTDIR, f"{name}.qc_cells.with_pass.csv")
        df  = pd.read_csv(qc_path, index_col=0)
        dpf = pd.read_csv(pass_path, index_col=0) if os.path.exists(pass_path) else df.copy()
        frip= np.asarray(ad.obs["frip"], dtype=float)
        df["frip"] = frip

        nfr  = df["n_fragment"].to_numpy(float) if "n_fragment" in df else np.full(len(df), np.nan)
        tsse = df["tsse"].to_numpy(float)       if "tsse"       in df else np.full(len(df), np.nan)
        mito = df["frac_mito"].to_numpy(float)  if "frac_mito"  in df else np.full(len(df), np.nan)

        mask = (nfr>=MIN_COUNTS) & (nfr<=MAX_COUNTS) & (tsse>=MIN_TSSE_FOR_PASS)
        if not np.isnan(mito).all(): mask &= (mito<=MAX_MITO)
        mask &= (frip>=FRIP_THRESH)
        dpf["qc_pass"] = mask

        df.to_csv(qc_path); dpf.to_csv(pass_path)

        # FRiP histogram (separate PNG)
        plt.figure(figsize=(6.5,4))
        v = df["frip"].to_numpy(float)
        plt.hist(v[~np.isnan(v)], bins=60)
        plt.axvline(FRIP_THRESH, ls="--", c="red", label=f"FRiP ≥ {FRIP_THRESH}")
        plt.title(f"{name} — FRiP histogram"); plt.xlabel("FRiP"); plt.ylabel("count"); plt.legend()
        plt.tight_layout(); plt.savefig(os.path.join(OUTDIR, f"{name}.frip_hist.png"), dpi=150); plt.close()
    print("[INFO] FRiP computed & qc_pass updated.")

# ---------- Pseudobulk reproducibility (optional; requires peak matrix API & UNION_PEAKS) ----------
if has_union:
    pseudo = {}
    for name in DATASETS:
        h5 = latest_h5ad_for(name, OUTDIR)
        if not h5: continue
        ad = snap.read(h5) if hasattr(snap, "read") else snap.pp.load_file(h5)
        peak_path = os.path.join(OUTDIR, f"{name}.union_peak_matrix.h5ad")
        pm, ok = build_peak_matrix_compat(ad, UNION_PEAKS, peak_path)
        if not ok or pm is None:
            print(f"[WARN] make_peak_matrix unavailable/failed for {name}; reproducibility skipped.")
            pseudo = {}
            break
        X = pm.X
        c = np.array(X.sum(axis=0)).ravel() if hasattr(X, "getnnz") else np.asarray(X).sum(axis=0).ravel()
        pseudo[name]=c
        if hasattr(pm, "close"):
            try: pm.close()
            except Exception: pass
        del pm, X; gc.collect()

    if pseudo:
        dfp = pd.DataFrame(pseudo)    # rows=peaks, cols=samples
        corr = dfp.corr(method="spearman")
        heat_pdf = os.path.join(OUTDIR, "pseudobulk_repro_spearman.pdf")
        with PdfPages(heat_pdf) as pdfh:
            plt.figure(figsize=(5.2,4.8))
            im = plt.imshow(corr.values, cmap="viridis", vmin=0, vmax=1)
            plt.colorbar(im, fraction=0.046, pad=0.04, label="Spearman ρ")
            plt.xticks(range(len(corr.columns)), corr.columns, rotation=30, ha="right")
            plt.yticks(range(len(corr.index)), corr.index)
            plt.title("Pseudobulk reproducibility (Spearman)")
            plt.tight_layout(); pdfh.savefig(); plt.close()
        corr.to_csv(os.path.join(OUTDIR, "pseudobulk_spearman_matrix.csv"))
        print("[INFO] Pseudobulk heatmap ->", heat_pdf)

# ---------- Summary CSV ----------
rows=[]
for name in DATASETS:
    qcf = os.path.join(OUTDIR, f"{name}.qc_cells.csv")
    qpf = os.path.join(OUTDIR, f"{name}.qc_cells.with_pass.csv")
    if not os.path.exists(qcf): continue
    df = pd.read_csv(qcf, index_col=0)
    dpf= pd.read_csv(qpf, index_col=0) if os.path.exists(qpf) else df
    def med(col): 
        return float(np.nanmedian(df[col])) if col in df.columns and not np.isnan(df[col]).all() else np.nan
    rows.append({
        "dataset": name,
        "cells": int(df.shape[0]),
        "pass_cells": int(np.nan_to_num(dpf["qc_pass"]).sum()) if "qc_pass" in dpf else 0,
        "pass_rate": float(np.nanmean(dpf["qc_pass"])) if "qc_pass" in dpf else np.nan,
        "n_fragment_median": med("n_fragment"),
        "tsse_median": med("tsse"),
        "frac_mito_median": med("frac_mito"),
        "frac_dup_median" : med("frac_dup"),
        "doublet_score_median": med("doublet_score"),
        "nucleosome_signal_median": med("nucleosome_signal"),
    })
summary = pd.DataFrame(rows).sort_values("dataset")
summary.to_csv(os.path.join(OUTDIR, "benchmark_summary.csv"), index=False)
print("[INFO] Summary CSV ->", os.path.join(OUTDIR, "benchmark_summary.csv"))

# ---------- Scoring（FAIR & DISPLAY/SeekGene-lean） ----------
def safe_nanmedian(a):
    arr = np.asarray(a, dtype=float)
    return np.nan if np.isnan(arr).all() else float(np.nanmedian(arr))

def norm_hi(x):
    x = x.astype(float); 
    if x.isna().all(): return x
    lo,hi=np.nanmin(x),np.nanmax(x)
    return 100*(x-lo)/(hi-lo) if hi>lo else np.full_like(x, np.nan, dtype=float)

def norm_lo(x):
    x = x.astype(float); 
    if x.isna().all(): return x
    lo,hi=np.nanmin(x),np.nanmax(x)
    return 100*(hi-x)/(hi-lo) if hi>lo else np.full_like(x, np.nan, dtype=float)

def med_table(names):
    rows=[]
    for n in names:
        qcf = os.path.join(OUTDIR, f"{n}.qc_cells.csv")
        qpf = os.path.join(OUTDIR, f"{n}.qc_cells.with_pass.csv")
        if not os.path.exists(qcf): continue
        df  = pd.read_csv(qcf, index_col=0)
        dpf = pd.read_csv(qpf, index_col=0) if os.path.exists(qpf) else df
        med = {
            "n_fragment"        : safe_nanmedian(df["n_fragment"]) if "n_fragment" in df else np.nan,
            "tsse"              : safe_nanmedian(df["tsse"])       if "tsse"       in df else np.nan,
            "frip"              : safe_nanmedian(df["frip"])       if "frip"       in df else np.nan,
            "frac_mito"         : safe_nanmedian(df["frac_mito"])  if "frac_mito"  in df else np.nan,
            "frac_dup"          : safe_nanmedian(df["frac_dup"])   if "frac_dup"   in df else np.nan,
            "doublet_score"     : safe_nanmedian(df["doublet_score"]) if "doublet_score" in df else np.nan,
            "nucleosome_signal" : safe_nanmedian(df["nucleosome_signal"]) if "nucleosome_signal" in df else np.nan,
            "pass_rate"         : float(np.nanmean(dpf["qc_pass"])) if "qc_pass" in dpf else np.nan,
        }
        med["log10_nfrag"] = np.log10(med["n_fragment"]) if np.isfinite(med["n_fragment"]) and med["n_fragment"]>0 else np.nan
        rows.append((n, med))
    return pd.DataFrame(dict(rows)).T

tab = med_table(DATASETS.keys())

# FAIR: objective (drop dup/doublet, downweight pass; if FRiP is all NaN it is auto-redistributed to the remaining items)
FAIR = {"S_tsse":0.40, "S_frip":0.35, "S_nfrag":0.20, "S_mito":0.05, "S_passrate":0.0, "S_doublet":0.0, "S_dup":0.0}
# DISPLAY: SeekGene-lean (emphasizes FRiP + TSSe + complexity more)
DISPLAY = {"S_tsse":0.45, "S_frip":0.40, "S_nfrag":0.15, "S_mito":0.0, "S_passrate":0.0, "S_doublet":0.0, "S_dup":0.0}

def build_subscores(table):
    S = pd.DataFrame(index=table.index)
    S["S_tsse"]     = norm_hi(table["tsse"])
    S["S_frip"]     = norm_hi(table["frip"])
    S["S_nfrag"]    = norm_hi(table["log10_nfrag"])
    S["S_mito"]     = norm_lo(table["frac_mito"])
    S["S_dup"]      = norm_lo(table["frac_dup"])
    S["S_doublet"]  = norm_lo(table["doublet_score"])
    S["S_passrate"] = norm_hi(table["pass_rate"])
    return S

def row_weighted_total(row, weights: dict):
    vals = row.reindex(weights.keys())
    valid = vals.notna()
    if not valid.any(): return np.nan
    w = pd.Series(weights)[valid]
    w = w / w.sum()
    return float((vals[valid] * w).sum())

S = build_subscores(tab)
Sf = pd.DataFrame(index=S.index)
Sd = pd.DataFrame(index=S.index)
Sf["Total_fair"]    = S.apply(lambda r: row_weighted_total(r, FAIR), axis=1)
Sd["Total_display"] = S.apply(lambda r: row_weighted_total(r, DISPLAY), axis=1)

# Write back to scoring CSV
scoring_csv = os.path.join(OUTDIR, "benchmark_scoring.csv")
pd.concat([tab, S, Sf, Sd], axis=1).to_csv(scoring_csv)
print("[INFO] Scoring CSV ->", scoring_csv)

# Ranking plots
def barplot(series, title, png):
    s = series.sort_values(ascending=False)
    plt.figure(figsize=(6,3.6))
    s.plot(kind="barh")
    plt.gca().invert_yaxis()
    plt.xlabel("Score (0–100)"); plt.title(title); plt.tight_layout()
    plt.savefig(png, dpi=150); plt.close()

barplot(Sf["Total_fair"],    "Benchmark Ranking — FAIR",    os.path.join(OUTDIR, "ranking_fair.png"))
barplot(Sd["Total_display"], "Benchmark Ranking — DISPLAY (SeekGene-lean)", os.path.join(OUTDIR, "ranking_display.png"))
print("[INFO] Rankings ->", os.path.join(OUTDIR, "ranking_fair.png"), ",", os.path.join(OUTDIR, "ranking_display.png"))

# Also save a summary bar-chart PDF (from summary/scoring)
summary = pd.read_csv(os.path.join(OUTDIR, "benchmark_summary.csv"))
order = list(summary["dataset"])
sc = pd.read_csv(scoring_csv, index_col=0).reindex(order)

bars_pdf = os.path.join(OUTDIR, "benchmark_bars.pdf")
with PdfPages(bars_pdf) as bpdf:
    # Cover page
    fig = plt.figure(figsize=(8,4))
    plt.axis("off"); plt.text(0.02,0.9,"scATAC Benchmark — Barplots", fontsize=14)
    plt.text(0.02,0.7,f"Summary: {os.path.join(OUTDIR,'benchmark_summary.csv')}\nScoring: {scoring_csv}", fontsize=10)
    bpdf.savefig(fig); plt.close(fig)

    def page(series, title, ylabel=None, pct=False):
        s = pd.Series(series.values, index=order).astype(float)
        if pct: s = s*100.0
        fig, ax = plt.subplots(figsize=(8,4))
        s.plot(kind="bar", ax=ax)
        ax.set_title(title); ax.set_ylabel(ylabel or ("Percent (%)" if pct else "")); ax.set_xlabel("")
        ax.set_xticklabels(order, rotation=25, ha="right")
        for p in ax.patches:
            v=p.get_height()
            if np.isfinite(v): ax.annotate(f"{v:.2f}", (p.get_x()+p.get_width()/2, v), xytext=(0,3),
                                           textcoords="offset points", ha="center", va="bottom", fontsize=8)
        fig.tight_layout(); bpdf.savefig(fig); plt.close(fig)

    # summary metrics
    page(summary["n_fragment_median"], "Median Unique Fragments", "Fragments")
    page(np.log10(np.clip(summary["n_fragment_median"],1,None)), "Median log10(Unique Fragments)", "log10(Fragments)")
    page(summary["tsse_median"], "Median TSSe", "TSSe")
    page(summary["frac_mito_median"], "Median Mito Fraction", "Mito (%)", pct=True)
    page(summary["frac_dup_median"], "Median Duplicate Fraction", "Dup (%)", pct=True)
    page(summary["pass_rate"], "Pass Rate", "Pass (%)", pct=True)

    # scoring total
    if "Total_fair" in sc.columns:     page(sc["Total_fair"], "Total Score — FAIR", "Score (0–100)")
    if "Total_display" in sc.columns:  page(sc["Total_display"], "Total Score — DISPLAY", "Score (0–100)")

print("[INFO] Barplots PDF ->", bars_pdf)


# In[1]:


# ===================== Multi-tech scATAC Benchmark + Extras + Scoring (SeekGene-lean option) =====================
# Four datasets: 10x multiome / 10x scATAC / SeekGene / TXCI (Droplet removed)
# End-to-end: fragments import -> TSSe (local GENCODE vM31) -> (if supported) tile/bin+Scrublet -> fragment-length triple peaks
#        (optional) union peak set FRiP -> (optional) pseudobulk reproducibility -> summary/scoring (FAIR & DISPLAY) -> ranking plots & PDF
# Note: this script is compatible with different SnapATAC2 API versions (make_tile_matrix/make_bin_matrix/peak_matrix may not exist)

# ---------- CONFIG ----------
BIG_TMP  = "/path/to/tmp"                           # large-disk TMP
OUTDIR   = "/path/to/benchmark_scATAC_hdf5"     # output directory
GENOME   = "mm10"
GENE_ANNOT = "/path/to/gencode.vM31.basic.annotation.gtf.gz"

# Four datasets (Droplet removed)
DATASETS = {
    "10x_multiome": "/path/to/scATAC/brain/10xmultiome/M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_atac_fragments.tsv.gz",
    "10x_scatac"  : "/path/to/scATAC/brain/10xscatac/atac_v1_E18_brain_flash_5k_fragments.tsv.gz",
    "SeekGene"    : "/path/to/scATAC/brain/seekgene/atac_fragments.tsv.gz",
    "TXCI"        : "/path/to/scATAC/brain/txci-atac/GSM7852211_mm10.merged.fragments.tsv.gz",
}

# (optional) union peak set (if absent, set to ""/None; this script auto-skips FRiP & reproducibility)
UNION_PEAKS = ""  # e.g.: "/path/to/benchmark_scATAC_hdf5/union_peaks.bed"

# Uniform thresholds (closer to the scale of your current data; can be fine-tuned later based on results)
MIN_COUNTS, MAX_COUNTS = 1000, 200000
MIN_TSSE_FOR_PASS      = 1.0
FRIP_THRESH            = 0.20
MAX_MITO               = 0.20

# Import parameters (HDF5 only; your environment does not support zarr)
BACKEND          = "hdf5"
CHUNK_SIZE       = 1_000_000
N_JOBS_IMPORT    = 4
MIN_FRAGS_IMPORT = max(500, MIN_COUNTS)

# Additional analysis
DO_NUCLEOSOME   = True
NUC_MAX_LINES   = None
SCRUBLET_JOBS   = 4
DOWNSAMPLE_PLOT = 10_000

# ---------- ENV ----------
import os, sys, re, gzip, time, gc, warnings, glob
from pathlib import Path
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
os.environ.setdefault("HDF5_DISABLE_VERSION_CHECK", "2")
os.environ["TMPDIR"] = os.environ["TMP"] = os.environ["TEMP"] = BIG_TMP
Path(BIG_TMP).mkdir(parents=True, exist_ok=True)
Path(OUTDIR).mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", os.path.join(BIG_TMP, "mpl"))

# ---------- DEPS ----------
import numpy as np, pandas as pd, h5py
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
warnings.filterwarnings("ignore", message=".*_import_from_c.*")
try:
    import snapatac2 as snap
except Exception as e:
    print("ERROR: snapatac2 not installed. Conda: conda create -n snap2 -y python=3.10 snapatac2 h5py pandas matplotlib", file=sys.stderr); raise

# Publication-quality PDF embedding
plt.rcParams.update({
    "pdf.fonttype": 42,  # Illustrator-friendly
    "ps.fonttype": 42,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
})

# ---------- HELPERS ----------
def ensure_tbi_for_frag(frag_path: str):
    tbi = frag_path + ".tbi"
    tbi_gz = tbi + ".gz"
    if os.path.exists(tbi): return
    if os.path.exists(tbi_gz):
        print(f"[INFO] Decompressing {tbi_gz} -> {tbi}")
        with gzip.open(tbi_gz, "rb") as fi, open(tbi, "wb") as fo: fo.write(fi.read())
    else:
        print(f"[WARN] No .tbi alongside {frag_path} (SnapATAC2 may index at TMP).")

def close_all_handles():
    for v in list(globals().values()):
        if hasattr(v, "close"):
            try: v.close()
            except Exception: pass
    gc.collect()

def unique_store(base, name):
    ts=time.strftime("%Y%m%d-%H%M%S"); pid=os.getpid()
    return os.path.join(base, f"{name}.raw.{ts}.{pid}.h5ad")

def frag_len_histogram(frag_path, max_len=800, step=1):
    bins=np.zeros(max_len+1, dtype=np.int64)  # 0 collects overflow
    with gzip.open(frag_path,"rt") as fh:
        for ln in fh:
            if not ln or ln[0]=="#": continue
            p=ln.rstrip("\n").split("\t")
            if len(p)<3: continue
            try: L=int(p[2])-int(p[1])
            except: continue
            if 0 < L <= max_len: bins[L]+=1
            else: bins[0]+=1
    idx=np.arange(1,max_len+1,step)
    val=bins[1:max_len+1:step].astype(float)
    tot=val.sum(); frac=val/(tot if tot>0 else 1.0)
    return idx, frac, int(bins[0])

def compute_nucleosome_signal_by_streaming(frag_path, keep_barcodes=None,
                                           nfr_max=147, mono_min=147, mono_max=294, max_lines=None):
    nfr, mono = {}, {}
    i=0
    with gzip.open(frag_path, "rt") as fh:
        for line in fh:
            if max_lines is not None and i>=max_lines: break
            i+=1
            if not line or line[0]=="#": continue
            p=line.rstrip("\n").split("\t")
            if len(p)<4: continue
            try: L = int(p[2]) - int(p[1])
            except: continue
            bc=p[3]
            if keep_barcodes is not None and bc not in keep_barcodes: continue
            if L < nfr_max:               nfr[bc]  = nfr.get(bc,0)  + 1
            elif mono_min <= L < mono_max: mono[bc] = mono.get(bc,0) + 1
    ks = (keep_barcodes if keep_barcodes is not None else set(list(nfr.keys())+list(mono.keys())))
    return {bc: (mono.get(bc,0)/(nfr.get(bc,0) if nfr.get(bc,0)>0 else 1.0)) for bc in ks}

def summarize_series(x):
    x = x[~np.isnan(x)]
    if x.size==0: return dict(n=0, median=np.nan, q1=np.nan, q3=np.nan)
    q = np.quantile(x, [0.25, 0.5, 0.75])
    return dict(n=x.size, median=q[1], q1=q[0], q3=q[2])

def latest_h5ad_for(name, outdir):
    patt = os.path.join(outdir, f"{name}.raw.*.h5ad")
    files = sorted(glob.glob(patt))
    return files[-1] if files else None

def run_tsse_safe(ad, gene_anno_path, genome_obj):
    """Prefer local GTF/GFF; else write NaN TSSe and continue"""
    try:
        if gene_anno_path and os.path.exists(gene_anno_path):
            snap.metrics.tsse(ad, gene_anno_path)
            return True
        raise FileNotFoundError(f"Local annotation not found: {gene_anno_path}")
    except Exception as e:
        print("[WARN] TSSe failed:", e)
        ad.obs["tsse"] = np.full(ad.n_obs, np.nan, dtype=float)
        return False

def build_tile_like_matrix(ad, out_path, bin_size=5000):
    """Compatible with different APIs: pp.make_tile_matrix -> pp.make_bin_matrix -> None"""
    if hasattr(snap.pp, "make_tile_matrix"):
        try:    return snap.pp.make_tile_matrix(ad, file=out_path, bin_size=bin_size, backend="hdf5"), True
        except TypeError: return snap.pp.make_tile_matrix(ad, file=out_path, bin_size=bin_size), True
    if hasattr(snap.pp, "make_bin_matrix"):
        try:    return snap.pp.make_bin_matrix(ad, file=out_path, bin_size=bin_size, backend="hdf5"), True
        except TypeError: return snap.pp.make_bin_matrix(ad, file=out_path, bin_size=bin_size), True
    print("[WARN] No tile/bin matrix API; skipping doublets.")
    return None, False

def build_peak_matrix_compat(ad, peak_bed, out_path):
    """Peak matrix (returns None, False if the version does not support it)"""
    if hasattr(snap.pp, "make_peak_matrix"):
        try:
            return snap.pp.make_peak_matrix(ad, file=out_path, peak_file=peak_bed), True
        except TypeError:
            return snap.pp.make_peak_matrix(ad, file=out_path, peaks=peak_bed), True
        except Exception:
            return None, False
    return None, False

# ---------- BENCHMARK LOOP ----------
summary_rows = []
pdf_path = os.path.join(OUTDIR, "benchmark_full_report.pdf")
pdf = PdfPages(pdf_path)

genome_obj = snap.genome.mm10 if GENOME.lower()=="mm10" else getattr(snap.genome, GENOME.lower())
assert os.path.exists(GENE_ANNOT), f"Annotation not found: {GENE_ANNOT}"

for name, frag in DATASETS.items():
    print(f"\n=== {name} ===")
    assert os.path.exists(frag), f"{frag} not found"
    ensure_tbi_for_frag(frag)

    # Import
    close_all_handles()
    store = unique_store(OUTDIR, name)
    print("[INFO] import_fragments …", store)
    ad = snap.pp.import_fragments(
        frag,
        chrom_sizes = genome_obj,
        file        = store,
        backend     = "hdf5",
        sorted_by_barcode = False,
        whitelist   = None,
        tempdir     = BIG_TMP,
        chunk_size  = CHUNK_SIZE,
        n_jobs      = N_JOBS_IMPORT,
        min_num_fragments = MIN_FRAGS_IMPORT,
    )

    # TSSe
    print("[INFO] TSSe …")
    _ = run_tsse_safe(ad, GENE_ANNOT, genome_obj)

    # Tile/bin + Scrublet (if available)
    print("[INFO] Tile 5kb + Scrublet …")
    tile_path = os.path.join(OUTDIR, f"{name}.tile_5kb.h5ad")
    tile, has_mat = build_tile_like_matrix(ad, tile_path, bin_size=5000)
    dbl = np.full(ad.n_obs, np.nan)
    if has_mat:
        try:
            snap.pp.scrublet(tile, n_jobs=SCRUBLET_JOBS)
            dbl = np.asarray(tile.obs["doublet_score"])
        except Exception as e:
            print("[WARN] scrublet failed:", e)

    # Nucleosome signal
    if DO_NUCLEOSOME:
        nuc = compute_nucleosome_signal_by_streaming(frag, keep_barcodes=set(ad.obs_names), max_lines=NUC_MAX_LINES)
        nuc = pd.Series(nuc).reindex(ad.obs_names).astype(float).values
    else:
        nuc = np.full(ad.n_obs, np.nan)

    # Per-cell QC & pass/fail
    qc = pd.DataFrame(index=ad.obs_names)
    for k in ["n_fragment","frac_dup","frac_mito","tsse"]:
        if k in ad.obs: qc[k] = np.asarray(ad.obs[k])
    qc["doublet_score"]      = dbl
    qc["nucleosome_signal"]  = nuc

    nfr  = qc["n_fragment"].to_numpy(float) if "n_fragment" in qc else np.full(ad.n_obs, np.nan)
    tsse = qc["tsse"].to_numpy(float)       if "tsse"       in qc else np.full(ad.n_obs, np.nan)
    mito = qc["frac_mito"].to_numpy(float)  if "frac_mito"  in qc else np.full(ad.n_obs, np.nan)

    mask = (nfr>=MIN_COUNTS) & (nfr<=MAX_COUNTS) & (tsse>=MIN_TSSE_FOR_PASS)
    if not np.isnan(mito).all(): mask &= (mito<=MAX_MITO)
    qc["qc_pass"] = mask

    qc_path  = os.path.join(OUTDIR, f"{name}.qc_cells.csv")
    pass_path= os.path.join(OUTDIR, f"{name}.qc_cells.with_pass.csv")
    qc.to_csv(qc_path); qc.to_csv(pass_path)

    # Summary row
    row = {
        "dataset": name,
        "cells": int(qc.shape[0]),
        "pass_cells": int(mask.sum()),
        "pass_rate": float(mask.mean()),
        "n_fragment_median": float(np.nanmedian(qc["n_fragment"])) if "n_fragment" in qc else np.nan,
        "tsse_median": float(np.nanmedian(qc["tsse"])) if "tsse" in qc else np.nan,
        "frac_mito_median": float(np.nanmedian(qc["frac_mito"])) if "frac_mito" in qc else np.nan,
        "frac_dup_median" : float(np.nanmedian(qc["frac_dup"]))  if "frac_dup"  in qc else np.nan,
        "doublet_score_median": float(np.nanmedian(qc["doublet_score"])) if "doublet_score" in qc else np.nan,
        "nucleosome_signal_median": float(np.nanmedian(qc["nucleosome_signal"])) if "nucleosome_signal" in qc else np.nan,
    }
    summary_rows.append(row)

    # ===== PDF: knee + length distribution + TSSe vs log10(nFrags) + histogram =====
    # 1) Knee
    vals = np.sort(nfr)[::-1]; ranks=np.arange(1, len(vals)+1)
    plt.figure(); plt.plot(ranks, vals)
    plt.xscale("log"); plt.yscale("log"); plt.xlabel("barcode rank (log10)"); plt.ylabel("n_fragment (log10)")
    plt.title(f"{name} — Knee"); plt.tight_layout(); pdf.savefig(); plt.close()

    # 2) Frag length
    ok_fsd=True
    try:
        fsd=snap.metrics.frag_size_distr(ad); ok_fsd=(fsd is not None)
    except Exception: ok_fsd=False
    if ok_fsd:
        x=np.arange(len(fsd)); y=np.array(fsd, dtype=float); x=x[1:]; y=y[1:]
        y_pct=y/(y.sum() if y.sum()>0 else 1.0)
        # counts
        plt.figure(); plt.plot(x,y); plt.yscale("log")
        plt.xlabel("fragment length (bp)"); plt.ylabel("count (log)")
        plt.title(f"{name} — Fragment length distribution"); plt.tight_layout(); pdf.savefig(); plt.close()
        # percent + triplet bands
        nfr_m=(x<147); mono_m=(x>=147)&(x<294); di_m=(x>=294)&(x<441)
        f_nfr, f_mono, f_di = y_pct[nfr_m].sum(), y_pct[mono_m].sum(), y_pct[di_m].sum()
        plt.figure(); plt.plot(x,y_pct, lw=1.2)
        for thr,lab in [(147,"147bp"),(294,"294bp"),(441,"441bp")]:
            plt.axvline(thr, ls="--", c="grey", lw=0.8)
            plt.text(thr+5, max(y_pct)*0.9, lab, fontsize=8, rotation=90, va="top", color="grey")
        plt.xlabel("Size of Fragments (bp)"); plt.ylabel("Fragments (%)")
        plt.title(f"{name} — Frag Size (percent)  NFR={f_nfr:.2%}  Mono={f_mono:.2%}  Di={f_di:.2%}")
        plt.tight_layout(); pdf.savefig(); plt.close()
    else:
        x,y_pct,overflow=frag_len_histogram(frag, max_len=800, step=1)
        nfr_m=(x<147); mono_m=(x>=147)&(x<294); di_m=(x>=294)&(x<441)
        f_nfr, f_mono, f_di = y_pct[nfr_m].sum(), y_pct[mono_m].sum(), y_pct[di_m].sum()
        plt.figure(); plt.plot(x,y_pct, lw=1.2)
        for thr,lab in [(147,"147bp"),(294,"294bp"),(441,"441bp")]:
            plt.axvline(thr, ls="--", c="grey", lw=0.8)
            plt.text(thr+5, max(y_pct)*0.9, lab, fontsize=8, rotation=90, va="top", color="grey")
        plt.xlabel("Size of Fragments (bp)"); plt.ylabel("Fragments (%)")
        plt.title(f"{name} — Frag Size (fallback)  NFR={f_nfr:.2%}  Mono={f_mono:.2%}  Di={f_di:.2%}  overflow>{800}bp:{overflow:,}")
        plt.tight_layout(); pdf.savefig(); plt.close()

    # 3) TSSe vs log10(nFrags)
    ds = min(DOWNSAMPLE_PLOT, len(nfr))
    idx = np.random.choice(len(nfr), ds, replace=False) if ds<len(nfr) else np.arange(len(nfr))
    plt.figure(); plt.scatter(np.log10(np.clip(nfr[idx],1,None)), tsse[idx], s=3, alpha=0.5)
    plt.axvline(np.log10(MIN_COUNTS), ls="--"); plt.axhline(MIN_TSSE_FOR_PASS, ls="--")
    plt.xlabel("log10(n_fragment)"); plt.ylabel("TSSe"); plt.title(f"{name} — TSSe vs log10(nFrags)")
    plt.tight_layout(); pdf.savefig(); plt.close()

    # 4) Histograms
    def _hist(col, title, xlabel):
        if col in qc:
            x=qc[col].to_numpy(float)
            plt.figure(); plt.hist(x[~np.isnan(x)], bins=60)
            plt.title(f"{name} — {title}"); plt.xlabel(xlabel); plt.ylabel("count")
            plt.tight_layout(); pdf.savefig(); plt.close()
    _hist("n_fragment","n_fragment","n_fragment")
    _hist("tsse","TSSe","TSSe")
    _hist("frac_mito","frac_mito","frac_mito")
    _hist("frac_dup","frac_dup","frac_dup")
    _hist("doublet_score","doublet_score","doublet_score")
    _hist("nucleosome_signal","nucleosome_signal","nucleosome_signal")

pdf.close()
print("[INFO] Full PDF ->", pdf_path)

# ---------- FRiP with union peaks (optional; only when UNION_PEAKS exists) ----------
has_union = bool(UNION_PEAKS) and os.path.exists(UNION_PEAKS)
if has_union:
    for name in DATASETS:
        h5 = latest_h5ad_for(name, OUTDIR)
        if not h5: continue
        print(f"[INFO] FRiP — {name}")
        ad = snap.read(h5) if hasattr(snap, "read") else snap.pp.load_file(h5)
        try:    snap.metrics.frip(ad, {"frip": UNION_PEAKS})
        except TypeError: snap.metrics.frip(ad, UNION_PEAKS)

        # Update FRiP & pass mask
        qc_path  = os.path.join(OUTDIR, f"{name}.qc_cells.csv")
        pass_path= os.path.join(OUTDIR, f"{name}.qc_cells.with_pass.csv")
        df  = pd.read_csv(qc_path, index_col=0)
        dpf = pd.read_csv(pass_path, index_col=0) if os.path.exists(pass_path) else df.copy()
        frip= np.asarray(ad.obs["frip"], dtype=float)
        df["frip"] = frip

        nfr  = df["n_fragment"].to_numpy(float) if "n_fragment" in df else np.full(len(df), np.nan)
        tsse = df["tsse"].to_numpy(float)       if "tsse"       in df else np.full(len(df), np.nan)
        mito = df["frac_mito"].to_numpy(float)  if "frac_mito"  in df else np.full(len(df), np.nan)

        mask = (nfr>=MIN_COUNTS) & (nfr<=MAX_COUNTS) & (tsse>=MIN_TSSE_FOR_PASS)
        if not np.isnan(mito).all(): mask &= (mito<=MAX_MITO)
        mask &= (frip>=FRIP_THRESH)
        dpf["qc_pass"] = mask

        df.to_csv(qc_path); dpf.to_csv(pass_path)

        # FRiP histogram (separate PNG)
        plt.figure(figsize=(6.5,4))
        v = df["frip"].to_numpy(float)
        plt.hist(v[~np.isnan(v)], bins=60)
        plt.axvline(FRIP_THRESH, ls="--", c="red", label=f"FRiP ≥ {FRIP_THRESH}")
        plt.title(f"{name} — FRiP histogram"); plt.xlabel("FRiP"); plt.ylabel("count"); plt.legend()
        plt.tight_layout(); plt.savefig(os.path.join(OUTDIR, f"{name}.frip_hist.png"), dpi=150); plt.close()
    print("[INFO] FRiP computed & qc_pass updated.")

# ---------- Pseudobulk reproducibility (optional; requires peak matrix API & UNION_PEAKS) ----------
if has_union:
    pseudo = {}
    for name in DATASETS:
        h5 = latest_h5ad_for(name, OUTDIR)
        if not h5: continue
        ad = snap.read(h5) if hasattr(snap, "read") else snap.pp.load_file(h5)
        peak_path = os.path.join(OUTDIR, f"{name}.union_peak_matrix.h5ad")
        pm, ok = build_peak_matrix_compat(ad, UNION_PEAKS, peak_path)
        if not ok or pm is None:
            print(f"[WARN] make_peak_matrix unavailable/failed for {name}; reproducibility skipped.")
            pseudo = {}
            break
        X = pm.X
        c = np.array(X.sum(axis=0)).ravel() if hasattr(X, "getnnz") else np.asarray(X).sum(axis=0).ravel()
        pseudo[name]=c
        if hasattr(pm, "close"):
            try: pm.close()
            except Exception: pass
        del pm, X; gc.collect()

    if pseudo:
        dfp = pd.DataFrame(pseudo)    # rows=peaks, cols=samples
        corr = dfp.corr(method="spearman")
        heat_pdf = os.path.join(OUTDIR, "pseudobulk_repro_spearman.pdf")
        with PdfPages(heat_pdf) as pdfh:
            plt.figure(figsize=(5.2,4.8))
            im = plt.imshow(corr.values, cmap="viridis", vmin=0, vmax=1)
            plt.colorbar(im, fraction=0.046, pad=0.04, label="Spearman ρ")
            plt.xticks(range(len(corr.columns)), corr.columns, rotation=30, ha="right")
            plt.yticks(range(len(corr.index)), corr.index)
            plt.title("Pseudobulk reproducibility (Spearman)")
            plt.tight_layout(); pdfh.savefig(); plt.close()
        corr.to_csv(os.path.join(OUTDIR, "pseudobulk_spearman_matrix.csv"))
        print("[INFO] Pseudobulk heatmap ->", heat_pdf)

# ---------- Summary CSV ----------
rows=[]
for name in DATASETS:
    qcf = os.path.join(OUTDIR, f"{name}.qc_cells.csv")
    qpf = os.path.join(OUTDIR, f"{name}.qc_cells.with_pass.csv")
    if not os.path.exists(qcf): continue
    df = pd.read_csv(qcf, index_col=0)
    dpf= pd.read_csv(qpf, index_col=0) if os.path.exists(qpf) else df
    def med(col): 
        return float(np.nanmedian(df[col])) if col in df.columns and not np.isnan(df[col]).all() else np.nan
    rows.append({
        "dataset": name,
        "cells": int(df.shape[0]),
        "pass_cells": int(np.nan_to_num(dpf["qc_pass"]).sum()) if "qc_pass" in dpf else 0,
        "pass_rate": float(np.nanmean(dpf["qc_pass"])) if "qc_pass" in dpf else np.nan,
        "n_fragment_median": med("n_fragment"),
        "tsse_median": med("tsse"),
        "frac_mito_median": med("frac_mito"),
        "frac_dup_median" : med("frac_dup"),
        "doublet_score_median": med("doublet_score"),
        "nucleosome_signal_median": med("nucleosome_signal"),
    })
summary = pd.DataFrame(rows).sort_values("dataset")
summary.to_csv(os.path.join(OUTDIR, "benchmark_summary.csv"), index=False)
print("[INFO] Summary CSV ->", os.path.join(OUTDIR, "benchmark_summary.csv"))

# ---------- Scoring（FAIR & DISPLAY/SeekGene-lean） ----------
def safe_nanmedian(a):
    arr = np.asarray(a, dtype=float)
    return np.nan if np.isnan(arr).all() else float(np.nanmedian(arr))

def norm_hi(x):
    x = x.astype(float)
    if x.isna().all(): return x
    lo,hi=np.nanmin(x),np.nanmax(x)
    return 100*(x-lo)/(hi-lo) if hi>lo else np.full_like(x, np.nan, dtype=float)

def norm_lo(x):
    x = x.astype(float)
    if x.isna().all(): return x
    lo,hi=np.nanmin(x),np.nanmax(x)
    return 100*(hi-x)/(hi-lo) if hi>lo else np.full_like(x, np.nan, dtype=float)

def med_table(names):
    rows=[]
    for n in names:
        qcf = os.path.join(OUTDIR, f"{n}.qc_cells.csv")
        qpf = os.path.join(OUTDIR, f"{n}.qc_cells.with_pass.csv")
        if not os.path.exists(qcf): continue
        df  = pd.read_csv(qcf, index_col=0)
        dpf = pd.read_csv(qpf, index_col=0) if os.path.exists(qpf) else df
        med = {
            "n_fragment"        : safe_nanmedian(df["n_fragment"]) if "n_fragment" in df else np.nan,
            "tsse"              : safe_nanmedian(df["tsse"])       if "tsse"       in df else np.nan,
            "frip"              : safe_nanmedian(df["frip"])       if "frip"       in df else np.nan,
            "frac_mito"         : safe_nanmedian(df["frac_mito"])  if "frac_mito"  in df else np.nan,
            "frac_dup"          : safe_nanmedian(df["frac_dup"])   if "frac_dup"   in df else np.nan,
            "doublet_score"     : safe_nanmedian(df["doublet_score"]) if "doublet_score" in df else np.nan,
            "nucleosome_signal" : safe_nanmedian(df["nucleosome_signal"]) if "nucleosome_signal" in df else np.nan,
            "pass_rate"         : float(np.nanmean(dpf["qc_pass"])) if "qc_pass" in dpf else np.nan,
        }
        med["log10_nfrag"] = np.log10(med["n_fragment"]) if np.isfinite(med["n_fragment"]) and med["n_fragment"]>0 else np.nan
        rows.append((n, med))
    return pd.DataFrame(dict(rows)).T

tab = med_table(DATASETS.keys())

# FAIR: objective (drop dup/doublet, downweight pass; if FRiP is all NaN it is auto-redistributed to the remaining items)
FAIR = {"S_tsse":0.40, "S_frip":0.35, "S_nfrag":0.20, "S_mito":0.05, "S_passrate":0.0, "S_doublet":0.0, "S_dup":0.0}
# DISPLAY: SeekGene-lean (emphasizes FRiP + TSSe + complexity more)
DISPLAY = {"S_tsse":0.45, "S_frip":0.40, "S_nfrag":0.15, "S_mito":0.0, "S_passrate":0.0, "S_doublet":0.0, "S_dup":0.0}

def build_subscores(table):
    S = pd.DataFrame(index=table.index)
    S["S_tsse"]     = norm_hi(table["tsse"])
    S["S_frip"]     = norm_hi(table["frip"])
    S["S_nfrag"]    = norm_hi(table["log10_nfrag"])
    S["S_mito"]     = norm_lo(table["frac_mito"])
    S["S_dup"]      = norm_lo(table["frac_dup"])
    S["S_doublet"]  = norm_lo(table["doublet_score"])
    S["S_passrate"] = norm_hi(table["pass_rate"])
    return S

def row_weighted_total(row, weights: dict):
    vals = row.reindex(weights.keys())
    valid = vals.notna()
    if not valid.any(): return np.nan
    w = pd.Series(weights)[valid]
    w = w / w.sum()
    return float((vals[valid] * w).sum())

S = build_subscores(tab)
Sf = pd.DataFrame(index=S.index)
Sd = pd.DataFrame(index=S.index)
Sf["Total_fair"]    = S.apply(lambda r: row_weighted_total(r, FAIR), axis=1)
Sd["Total_display"] = S.apply(lambda r: row_weighted_total(r, DISPLAY), axis=1)

# Write back to scoring CSV
scoring_csv = os.path.join(OUTDIR, "benchmark_scoring.csv")
pd.concat([tab, S, Sf, Sd], axis=1).to_csv(scoring_csv)
print("[INFO] Scoring CSV ->", scoring_csv)

# Ranking plots (original PNG version)
def barplot(series, title, png):
    s = series.sort_values(ascending=False)
    plt.figure(figsize=(6,3.6))
    s.plot(kind="barh")
    plt.gca().invert_yaxis()
    plt.xlabel("Score (0–100)"); plt.title(title); plt.tight_layout()
    plt.savefig(png, dpi=150); plt.close()

barplot(Sf["Total_fair"],    "Benchmark Ranking — FAIR",    os.path.join(OUTDIR, "ranking_fair.png"))
barplot(Sd["Total_display"], "Benchmark Ranking — DISPLAY (SeekGene-lean)", os.path.join(OUTDIR, "ranking_display.png"))
print("[INFO] Rankings ->", os.path.join(OUTDIR, "ranking_fair.png"), ",", os.path.join(OUTDIR, "ranking_display.png"))

# Also save a summary bar-chart PDF (from summary/scoring)
summary = pd.read_csv(os.path.join(OUTDIR, "benchmark_summary.csv"))
order = list(summary["dataset"])
sc = pd.read_csv(scoring_csv, index_col=0).reindex(order)

bars_pdf = os.path.join(OUTDIR, "benchmark_bars.pdf")
with PdfPages(bars_pdf) as bpdf:
    # Cover page
    fig = plt.figure(figsize=(8,4))
    plt.axis("off"); plt.text(0.02,0.9,"scATAC Benchmark — Barplots", fontsize=14)
    plt.text(0.02,0.7,f"Summary: {os.path.join(OUTDIR,'benchmark_summary.csv')}\nScoring: {scoring_csv}", fontsize=10)
    bpdf.savefig(fig); plt.close(fig)

    def page(series, title, ylabel=None, pct=False):
        s = pd.Series(series.values, index=order).astype(float)
        if pct: s = s*100.0
        fig, ax = plt.subplots(figsize=(8,4))
        s.plot(kind="bar", ax=ax)
        ax.set_title(title); ax.set_ylabel(ylabel or ("Percent (%)" if pct else "")); ax.set_xlabel("")
        ax.set_xticklabels(order, rotation=25, ha="right")
        for p in ax.patches:
            v=p.get_height()
            if np.isfinite(v): ax.annotate(f"{v:.2f}", (p.get_x()+p.get_width()/2, v), xytext=(0,3),
                                           textcoords="offset points", ha="center", va="bottom", fontsize=8)
        fig.tight_layout(); bpdf.savefig(fig); plt.close(fig)

    # summary metrics
    page(summary["n_fragment_median"], "Median Unique Fragments", "Fragments")
    page(np.log10(np.clip(summary["n_fragment_median"],1,None)), "Median log10(Unique Fragments)", "log10(Fragments)")
    page(summary["tsse_median"], "Median TSSe", "TSSe")
    page(summary["frac_mito_median"], "Median Mito Fraction", "Mito (%)", pct=True)
    page(summary["frac_dup_median"], "Median Duplicate Fraction", "Dup (%)", pct=True)
    page(summary["pass_rate"], "Pass Rate", "Pass (%)", pct=True)

    # scoring total
    if "Total_fair" in sc.columns:     page(sc["Total_fair"], "Total Score — FAIR", "Score (0–100)")
    if "Total_display" in sc.columns:  page(sc["Total_display"], "Total Score — DISPLAY", "Score (0–100)")

print("[INFO] Barplots PDF ->", bars_pdf)

# ======================= Scoring Visualizations: multi-page report + standalone PDF =======================
# Palette (by your four datasets)
SCATAC_PALETTE = {
    "10x_multiome": "#86c7b4",
    "10x_scatac":   "#fdbf6f",
    "SeekGene":     "#d8a0c0",
    "TXCI":         "#005AC8",
}

# Subscore colors (stacked breakdown page)
METRIC_COLORS = {
    "S_tsse": "#1f77b4",
    "S_frip": "#ff7f0e",
    "S_nfrag": "#2ca02c",
    "S_mito": "#d62728",
    "S_dup": "#9467bd",
    "S_doublet": "#8c564b",
    "S_passrate": "#7f7f7f",
}

def _safe_cols(df, cols):
    return [c for c in cols if c in df.columns]

def _order_for(series: pd.Series):
    try:
        return series.sort_values(ascending=False).index.tolist()
    except Exception:
        return series.index.tolist()

def _weighted_contrib(row, weights: dict):
    """Return the weighted contribution of each metric (normalize weights over non-NA items)."""
    vals = row.reindex(weights.keys())
    ok = vals.notna()
    if not ok.any():
        return pd.Series(index=weights.keys(), dtype=float)
    w = pd.Series(weights).reindex(vals.index)[ok]
    w = w / w.sum()
    return (vals[ok] * w)

def make_scoring_report(tab, S, Sf, Sd, FAIR, DISPLAY, outdir=OUTDIR,
                        palette=SCATAC_PALETTE, metric_colors=METRIC_COLORS):
    # Sort
    order_fair    = _order_for(Sf["Total_fair"]) if "Total_fair" in Sf else list(S.index)
    order_display = _order_for(Sd["Total_display"]) if "Total_display" in Sd else list(S.index)
    order_union   = list(dict.fromkeys(order_display + order_fair))  # merged order, display shown first

    report_pdf = os.path.join(outdir, "benchmark_scoring_report.pdf")
    os.makedirs(outdir, exist_ok=True)
    with PdfPages(report_pdf) as rep:

        # ---- Cover page ----
        fig = plt.figure(figsize=(8.5, 5)); plt.axis("off")
        txt = "scATAC Benchmark — Scoring Report\n\n"
        if "Total_display" in Sd:
            top = Sd["Total_display"].sort_values(ascending=False).head(3)
            txt += "Top (DISPLAY): " + ", ".join([f"{k}={v:.1f}" for k,v in top.items()]) + "\n"
        if "Total_fair" in Sf:
            top = Sf["Total_fair"].sort_values(ascending=False).head(3)
            txt += "Top (FAIR): " + ", ".join([f"{k}={v:.1f}" for k,v in top.items()]) + "\n"
        plt.text(0.02, 0.85, txt, fontsize=14, va="top")
        rep.savefig(fig); plt.close(fig)

        # ---- Subscores heatmap ----
        cols = _safe_cols(S, ["S_tsse","S_frip","S_nfrag","S_mito","S_passrate","S_doublet","S_dup"])
        if len(cols) > 0:
            M = S.loc[order_union, cols].astype(float)
            fig, ax = plt.subplots(figsize=(max(6, 0.8*len(cols)+2), max(3, 0.45*len(order_union)+1)))
            im = ax.imshow(M.fillna(0).values, cmap="viridis", vmin=0, vmax=100, aspect="auto")
            # Mask NaN with a white box
            for r in range(M.shape[0]):
                for c in range(M.shape[1]):
                    if pd.isna(M.iat[r,c]):
                        ax.add_patch(plt.Rectangle((c-0.5, r-0.5), 1, 1, color="white", zorder=2))
            plt.colorbar(im, fraction=0.046, pad=0.04, label="Subscore (0–100)")
            ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, rotation=30, ha="right")
            ax.set_yticks(range(len(order_union))); ax.set_yticklabels(order_union)
            ax.set_title("Subscores heatmap")
            plt.tight_layout(); rep.savefig(fig); fig.savefig(os.path.join(outdir, "subscores_heatmap.pdf"), format="pdf"); plt.close(fig)

        # ---- Lollipop ranking plot (FAIR, DISPLAY) ----
        def lollipop(series, title, savepath):
            s = series.dropna().sort_values(ascending=True)  # low -> high, highest on the right
            fig, ax = plt.subplots(figsize=(6.6, 3.8))
            y = np.arange(len(s))
            cols = [palette.get(k, "#999999") for k in s.index]
            ax.hlines(y, xmin=0, xmax=s.values, color=cols, linewidth=2)
            ax.plot(s.values, y, "o", color="black", markersize=4)
            ax.set_yticks(y); ax.set_yticklabels(s.index)
            ax.set_xlabel("Score (0–100)"); ax.set_title(title)
            ax.grid(True, axis='x', linestyle=':', linewidth=0.8, alpha=0.6)
            plt.tight_layout()
            rep.savefig(fig); fig.savefig(savepath, format="pdf"); plt.close(fig)

        if "Total_fair" in Sf:
            lollipop(Sf["Total_fair"], "Ranking — FAIR", os.path.join(outdir, "ranking_fair.pdf"))
        if "Total_display" in Sd:
            lollipop(Sd["Total_display"], "Ranking — DISPLAY (SeekGene-lean)", os.path.join(outdir, "ranking_display.pdf"))

        # ---- TSSe vs FRiP scatter (marker size ~ log10 fragments) ----
        fig, ax = plt.subplots(figsize=(6.2, 4.8))
        x = tab["tsse"].astype(float)
        y = tab["frip"].astype(float)
        s = tab["log10_nfrag"].astype(float)
        if s.notna().any():
            s_norm = (s - np.nanmin(s)) / max(1e-6, (np.nanmax(s) - np.nanmin(s)))
        else:
            s_norm = s*0+0.5
        sizes = 40 + 160 * s_norm.fillna(0.3)
        for name in order_union:
            if name not in x.index: continue
            if pd.isna(x.loc[name]) or pd.isna(y.loc[name]): continue
            ax.scatter(x.loc[name], y.loc[name], s=float(sizes.loc[name]),
                       color=palette.get(name, "#999999"), alpha=0.8, edgecolor="white", linewidth=0.6, zorder=3)
            ax.text(x.loc[name], y.loc[name], f"  {name}", va="center", fontsize=8)
        ax.set_xlabel("Median TSSe")
        ax.set_ylabel("Median FRiP")
        ax.set_title("TSSe vs FRiP (marker size ~ log10 unique fragments)")
        ax.grid(True, linestyle=":", linewidth=0.8, alpha=0.6)
        plt.tight_layout(); rep.savefig(fig); fig.savefig(os.path.join(outdir, "tsse_vs_frip.pdf"), format="pdf"); plt.close(fig)

        # ---- Weighted contribution stacked plot (DISPLAY & FAIR) ----
        def stacked_contrib(weights, title, fname):
            metrics = [k for k in weights.keys() if k in S.columns]
            if not metrics: return
            C = []
            for ds in order_union:
                contrib = _weighted_contrib(S.loc[ds, :], weights).reindex(metrics)
                C.append(contrib.values)
            C = np.array(C).T  # (n_metrics, n_ds)
            fig, ax = plt.subplots(figsize=(max(6, 0.6*len(order_union)+1.8), 4.0))
            bottoms = np.zeros(C.shape[1])
            x = np.arange(C.shape[1])
            for i, m in enumerate(metrics):
                ax.bar(x, C[i], bottom=bottoms, color=metric_colors.get(m, "#888888"), edgecolor="none", label=m)
                bottoms += C[i]
            # Dataset outline (by sample palette)
            for i, ds in enumerate(order_union):
                ax.bar(x[i], bottoms[i], fill=False, edgecolor=palette.get(ds, "#666666"), linewidth=1.0)
            ax.set_xticks(x); ax.set_xticklabels(order_union, rotation=25, ha="right")
            ax.set_ylabel("Weighted contribution to total (0–100)")
            ax.set_title(title)
            ax.legend(frameon=False, ncol=min(3, len(metrics)))
            ax.grid(True, axis='y', linestyle=':', linewidth=0.8, alpha=0.6)
            plt.tight_layout(); rep.savefig(fig); fig.savefig(fname, format="pdf"); plt.close(fig)

        if "Total_display" in Sd:
            stacked_contrib(DISPLAY, "Contribution breakdown — DISPLAY (SeekGene-lean)",
                            os.path.join(outdir, "contrib_display.pdf"))
        if "Total_fair" in Sf:
            stacked_contrib(FAIR, "Contribution breakdown — FAIR",
                            os.path.join(outdir, "contrib_fair.pdf"))

        # ---- Radar plot (one page per dataset) ----
        cols = _safe_cols(S, ["S_tsse","S_frip","S_nfrag","S_mito","S_passrate","S_doublet","S_dup"])
        if len(cols) >= 3:
            theta = np.linspace(0, 2*np.pi, len(cols)+1)  # close loop
            for ds in order_union:
                vals = S.loc[ds, cols].astype(float).values
                fig = plt.figure(figsize=(4.8, 4.8))
                ax = plt.subplot(111, polar=True)
                v = np.append(vals, vals[0])  # close loop
                ax.plot(theta, v, color=palette.get(ds, "#999999"), linewidth=2)
                ax.fill(theta, v, color=palette.get(ds, "#999999"), alpha=0.20)
                ax.set_xticks(theta[:-1]); ax.set_xticklabels(cols, fontsize=8)
                ax.set_yticks([20,40,60,80]); ax.set_yticklabels(["20","40","60","80"])
                ax.set_ylim(0,100)
                ax.set_title(f"Radar — {ds}", va="bottom")
                plt.tight_layout(); rep.savefig(fig); plt.close(fig)

    print("[INFO] Scoring report PDF ->", report_pdf)

# === Generate scoring report and standalone PDF ===
make_scoring_report(tab, S, Sf, Sd, FAIR, DISPLAY, outdir=OUTDIR, palette=SCATAC_PALETTE)


# In[ ]:


# ===================== Multi-tech scATAC Benchmark + Extras + Scoring (SeekGene-lean option) =====================
# Four datasets: 10x multiome / 10x scATAC / SeekGene / TXCI (Droplet removed)
# End-to-end: fragments import -> TSSe (local GENCODE vM31) -> (if supported) tile/bin+Scrublet -> fragment-length triple peaks
#        (optional) union peak set FRiP -> (optional) pseudobulk reproducibility -> summary/scoring (FAIR & DISPLAY) -> ranking plots & PDF
# Note: compatible with different SnapATAC2 API versions (make_tile_matrix/make_bin_matrix/peak_matrix may not exist)

# ---------- CONFIG ----------
BIG_TMP  = "/path/to/tmp"
OUTDIR   = "/path/to/benchmark_scATAC_hdf5"
GENOME   = "mm10"
GENE_ANNOT = "/path/to/gencode.vM31.basic.annotation.gtf.gz"

DATASETS = {
    "10x_multiome": "/path/to/scATAC/brain/10xmultiome/M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_atac_fragments.tsv.gz",
    "10x_scatac"  : "/path/to/scATAC/brain/10xscatac/atac_v1_E18_brain_flash_5k_fragments.tsv.gz",
    "SeekGene"    : "/path/to/scATAC/brain/seekgene/atac_fragments.tsv.gz",
    "TXCI"        : "/path/to/scATAC/brain/txci-atac/GSM7852211_mm10.merged.fragments.tsv.gz",
}

# (optional) uniform peak set
UNION_PEAKS = ""  # e.g.: "/path/to/benchmark_scATAC_hdf5/union_peaks.bed"

# QC thresholds
MIN_COUNTS, MAX_COUNTS = 1000, 200000
MIN_TSSE_FOR_PASS      = 1.0
FRIP_THRESH            = 0.20
MAX_MITO               = 0.20

# Import parameters
BACKEND          = "hdf5"
CHUNK_SIZE       = 1_000_000
N_JOBS_IMPORT    = 4
MIN_FRAGS_IMPORT = max(500, MIN_COUNTS)

# Additional analysis
DO_NUCLEOSOME   = True
NUC_MAX_LINES   = None
SCRUBLET_JOBS   = 4
DOWNSAMPLE_PLOT = 10_000

# ---------- ENV ----------
import os, sys, re, gzip, time, gc, warnings, glob
from pathlib import Path
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
os.environ.setdefault("HDF5_DISABLE_VERSION_CHECK", "2")
os.environ["TMPDIR"] = os.environ["TMP"] = os.environ["TEMP"] = BIG_TMP
Path(BIG_TMP).mkdir(parents=True, exist_ok=True)
Path(OUTDIR).mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", os.path.join(BIG_TMP, "mpl"))

# ---------- DEPS ----------
import numpy as np, pandas as pd, h5py
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
warnings.filterwarnings("ignore", message=".*_import_from_c.*")
try:
    import snapatac2 as snap
except Exception as e:
    print("ERROR: install snapatac2 env: conda create -n snap2 -y python=3.10 snapatac2 h5py pandas matplotlib", file=sys.stderr); raise

# Publication-quality PDF
plt.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})

# ===== Fixed color scheme (as requested) =====
SCATAC_PALETTE = {
    "TXCI":         "#d4de9c",  # 1
    "SeekGene":     "#d8a0c0",  # 2
    "10x_multiome": "#86c7b4",  # 3
    "10x_scatac":   "#9cd2ed",  # 4
}
def c_for(name: str) -> str:
    return SCATAC_PALETTE.get(name, "#999999")

# ---------- HELPERS ----------
def ensure_tbi_for_frag(frag_path: str):
    tbi = frag_path + ".tbi"; tbi_gz = tbi + ".gz"
    if os.path.exists(tbi): return
    if os.path.exists(tbi_gz):
        print(f"[INFO] Decompress {tbi_gz} -> {tbi}")
        with gzip.open(tbi_gz, "rb") as fi, open(tbi, "wb") as fo: fo.write(fi.read())
    else:
        print(f"[WARN] No .tbi for {frag_path} (SnapATAC2 may index at TMP).")

def close_all_handles():
    for v in list(globals().values()):
        if hasattr(v, "close"):
            try: v.close()
            except Exception: pass
    gc.collect()

def unique_store(base, name):
    ts=time.strftime("%Y%m%d-%H%M%S"); pid=os.getpid()
    return os.path.join(base, f"{name}.raw.{ts}.{pid}.h5ad")

def frag_len_histogram(frag_path, max_len=800, step=1):
    bins=np.zeros(max_len+1, dtype=np.int64)
    with gzip.open(frag_path,"rt") as fh:
        for ln in fh:
            if not ln or ln[0]=="#": continue
            p=ln.rstrip("\n").split("\t")
            if len(p)<3: continue
            try: L=int(p[2])-int(p[1])
            except: continue
            if 0 < L <= max_len: bins[L]+=1
            else: bins[0]+=1
    idx=np.arange(1,max_len+1,step)
    val=bins[1:max_len+1:step].astype(float)
    tot=val.sum(); frac=val/(tot if tot>0 else 1.0)
    return idx, frac, int(bins[0])

def compute_nucleosome_signal_by_streaming(frag_path, keep_barcodes=None,
                                           nfr_max=147, mono_min=147, mono_max=294, max_lines=None):
    nfr, mono = {}, {}
    i=0
    with gzip.open(frag_path, "rt") as fh:
        for line in fh:
            if max_lines is not None and i>=max_lines: break
            i+=1
            if not line or line[0]=="#": continue
            p=line.rstrip("\n").split("\t")
            if len(p)<4: continue
            try: L = int(p[2]) - int(p[1])
            except: continue
            bc=p[3]
            if keep_barcodes is not None and bc not in keep_barcodes: continue
            if L < nfr_max:               nfr[bc]  = nfr.get(bc,0)  + 1
            elif mono_min <= L < mono_max: mono[bc] = mono.get(bc,0) + 1
    ks = (keep_barcodes if keep_barcodes is not None else set(list(nfr.keys())+list(mono.keys())))
    return {bc: (mono.get(bc,0)/(nfr.get(bc,0) if nfr.get(bc,0)>0 else 1.0)) for bc in ks}

def latest_h5ad_for(name, outdir):
    patt = os.path.join(outdir, f"{name}.raw.*.h5ad")
    files = sorted(glob.glob(patt))
    return files[-1] if files else None

def run_tsse_safe(ad, gene_anno_path, genome_obj):
    try:
        if gene_anno_path and os.path.exists(gene_anno_path):
            snap.metrics.tsse(ad, gene_anno_path); return True
        raise FileNotFoundError(f"Local annotation not found: {gene_anno_path}")
    except Exception as e:
        print("[WARN] TSSe failed:", e); ad.obs["tsse"] = np.full(ad.n_obs, np.nan, dtype=float); return False

def build_tile_like_matrix(ad, out_path, bin_size=5000):
    if hasattr(snap.pp, "make_tile_matrix"):
        try:    return snap.pp.make_tile_matrix(ad, file=out_path, bin_size=bin_size, backend="hdf5"), True
        except TypeError: return snap.pp.make_tile_matrix(ad, file=out_path, bin_size=bin_size), True
    if hasattr(snap.pp, "make_bin_matrix"):
        try:    return snap.pp.make_bin_matrix(ad, file=out_path, bin_size=bin_size, backend="hdf5"), True
        except TypeError: return snap.pp.make_bin_matrix(ad, file=out_path, bin_size=bin_size), True
    print("[WARN] No tile/bin matrix API; skipping doublets."); return None, False

def build_peak_matrix_compat(ad, peak_bed, out_path):
    if hasattr(snap.pp, "make_peak_matrix"):
        try:    return snap.pp.make_peak_matrix(ad, file=out_path, peak_file=peak_bed), True
        except TypeError: return snap.pp.make_peak_matrix(ad, file=out_path, peaks=peak_bed), True
        except Exception: return None, False
    return None, False

# ---------- BENCHMARK LOOP ----------
summary_rows = []
pdf_path = os.path.join(OUTDIR, "benchmark_full_report.pdf")
pdf = PdfPages(pdf_path)

genome_obj = snap.genome.mm10 if GENOME.lower()=="mm10" else getattr(snap.genome, GENOME.lower())
assert os.path.exists(GENE_ANNOT), f"Annotation not found: {GENE_ANNOT}"

for name, frag in DATASETS.items():
    print(f"\n=== {name} ===")
    assert os.path.exists(frag), f"{frag} not found"
    ensure_tbi_for_frag(frag)

    # Import
    close_all_handles()
    store = unique_store(OUTDIR, name)
    print("[INFO] import_fragments …", store)
    ad = snap.pp.import_fragments(
        frag,
        chrom_sizes = genome_obj,
        file        = store,
        backend     = "hdf5",
        sorted_by_barcode = False,
        whitelist   = None,
        tempdir     = BIG_TMP,
        chunk_size  = CHUNK_SIZE,
        n_jobs      = N_JOBS_IMPORT,
        min_num_fragments = MIN_FRAGS_IMPORT,
    )

    # TSSe
    print("[INFO] TSSe …")
    _ = run_tsse_safe(ad, GENE_ANNOT, genome_obj)

    # Tile/bin + Scrublet
    print("[INFO] Tile 5kb + Scrublet …")
    tile_path = os.path.join(OUTDIR, f"{name}.tile_5kb.h5ad")
    tile, has_mat = build_tile_like_matrix(ad, tile_path, bin_size=5000)
    dbl = np.full(ad.n_obs, np.nan)
    if has_mat:
        try:
            snap.pp.scrublet(tile, n_jobs=SCRUBLET_JOBS)
            dbl = np.asarray(tile.obs["doublet_score"])
        except Exception as e:
            print("[WARN] scrublet failed:", e)

    # Nucleosome signal
    if DO_NUCLEOSOME:
        nuc = compute_nucleosome_signal_by_streaming(frag, keep_barcodes=set(ad.obs_names), max_lines=NUC_MAX_LINES)
        nuc = pd.Series(nuc).reindex(ad.obs_names).astype(float).values
    else:
        nuc = np.full(ad.n_obs, np.nan)

    # Per-cell QC
    qc = pd.DataFrame(index=ad.obs_names)
    for k in ["n_fragment","frac_dup","frac_mito","tsse"]:
        if k in ad.obs: qc[k] = np.asarray(ad.obs[k])
    qc["doublet_score"]      = dbl
    qc["nucleosome_signal"]  = nuc

    nfr  = qc["n_fragment"].to_numpy(float) if "n_fragment" in qc else np.full(ad.n_obs, np.nan)
    tsse = qc["tsse"].to_numpy(float)       if "tsse"       in qc else np.full(ad.n_obs, np.nan)
    mito = qc["frac_mito"].to_numpy(float)  if "frac_mito"  in qc else np.full(ad.n_obs, np.nan)

    mask = (nfr>=MIN_COUNTS) & (nfr<=MAX_COUNTS) & (tsse>=MIN_TSSE_FOR_PASS)
    if not np.isnan(mito).all(): mask &= (mito<=MAX_MITO)
    qc["qc_pass"] = mask

    qc_path  = os.path.join(OUTDIR, f"{name}.qc_cells.csv")
    pass_path= os.path.join(OUTDIR, f"{name}.qc_cells.with_pass.csv")
    qc.to_csv(qc_path); qc.to_csv(pass_path)

    # Summary row
    summary_rows.append({
        "dataset": name,
        "cells": int(qc.shape[0]),
        "pass_cells": int(mask.sum()),
        "pass_rate": float(mask.mean()),
        "n_fragment_median": float(np.nanmedian(qc["n_fragment"])) if "n_fragment" in qc else np.nan,
        "tsse_median": float(np.nanmedian(qc["tsse"])) if "tsse" in qc else np.nan,
        "frac_mito_median": float(np.nanmedian(qc["frac_mito"])) if "frac_mito" in qc else np.nan,
        "frac_dup_median" : float(np.nanmedian(qc["frac_dup"]))  if "frac_dup"  in qc else np.nan,
        "doublet_score_median": float(np.nanmedian(qc["doublet_score"])) if "doublet_score" in qc else np.nan,
        "nucleosome_signal_median": float(np.nanmedian(qc["nucleosome_signal"])) if "nucleosome_signal" in qc else np.nan,
    })

    # ===== PDF (all colored) =====
    # 1) Knee
    vals = np.sort(nfr)[::-1]; ranks=np.arange(1, len(vals)+1)
    plt.figure(); plt.plot(ranks, vals, color=c_for(name))
    plt.xscale("log"); plt.yscale("log")
    plt.xlabel("barcode rank (log10)"); plt.ylabel("n_fragment (log10)")
    plt.title(f"{name} — Knee"); plt.tight_layout(); pdf.savefig(); plt.close()

    # 2) Frag length
    ok_fsd=True
    try:
        fsd=snap.metrics.frag_size_distr(ad); ok_fsd=(fsd is not None)
    except Exception: ok_fsd=False
    if ok_fsd:
        x=np.arange(len(fsd)); y=np.array(fsd, dtype=float); x=x[1:]; y=y[1:]
        y_pct=y/(y.sum() if y.sum()>0 else 1.0)
        # counts
        plt.figure(); plt.plot(x,y, color=c_for(name))
        plt.yscale("log"); plt.xlabel("fragment length (bp)"); plt.ylabel("count (log)")
        plt.title(f"{name} — Fragment length distribution"); plt.tight_layout(); pdf.savefig(); plt.close()
        # percent + triplet bands
        nfr_m=(x<147); mono_m=(x>=147)&(x<294); di_m=(x>=294)&(x<441)
        f_nfr, f_mono, f_di = y_pct[nfr_m].sum(), y_pct[mono_m].sum(), y_pct[di_m].sum()
        plt.figure(); plt.plot(x,y_pct, lw=1.2, color=c_for(name))
        for thr,lab in [(147,"147bp"),(294,"294bp"),(441,"441bp")]:
            plt.axvline(thr, ls="--", c="grey", lw=0.8)
            plt.text(thr+5, max(y_pct)*0.9, lab, fontsize=8, rotation=90, va="top", color="grey")
        plt.xlabel("Size of Fragments (bp)"); plt.ylabel("Fragments (%)")
        plt.title(f"{name} — Frag Size (%)  NFR={f_nfr:.2%}  Mono={f_mono:.2%}  Di={f_di:.2%}")
        plt.tight_layout(); pdf.savefig(); plt.close()
    else:
        x,y_pct,overflow=frag_len_histogram(frag, max_len=800, step=1)
        nfr_m=(x<147); mono_m=(x>=147)&(x<294); di_m=(x>=294)&(x<441)
        f_nfr, f_mono, f_di = y_pct[nfr_m].sum(), y_pct[mono_m].sum(), y_pct[di_m].sum()
        plt.figure(); plt.plot(x,y_pct, lw=1.2, color=c_for(name))
        for thr,lab in [(147,"147bp"),(294,"294bp"),(441,"441bp")]:
            plt.axvline(thr, ls="--", c="grey", lw=0.8)
            plt.text(thr+5, max(y_pct)*0.9, lab, fontsize=8, rotation=90, va="top", color="grey")
        plt.xlabel("Size of Fragments (bp)"); plt.ylabel("Fragments (%)")
        plt.title(f"{name} — Frag Size (fallback)  NFR={f_nfr:.2%}  Mono={f_mono:.2%}  Di={f_di:.2%}  overflow>{800}bp:{overflow:,}")
        plt.tight_layout(); pdf.savefig(); plt.close()

    # 3) TSSe vs log10(nFrags)
    ds = min(DOWNSAMPLE_PLOT, len(nfr))
    idx = np.random.choice(len(nfr), ds, replace=False) if ds<len(nfr) else np.arange(len(nfr))
    plt.figure(); plt.scatter(np.log10(np.clip(nfr[idx],1,None)), tsse[idx], s=6, alpha=0.5, color=c_for(name))
    plt.axvline(np.log10(MIN_COUNTS), ls="--"); plt.axhline(MIN_TSSE_FOR_PASS, ls="--")
    plt.xlabel("log10(n_fragment)"); plt.ylabel("TSSe"); plt.title(f"{name} — TSSe vs log10(nFrags)")
    plt.tight_layout(); pdf.savefig(); plt.close()

    # 4) histograms (all colored)
    def _hist(col, title, xlabel):
        if col in qc:
            x=qc[col].to_numpy(float)
            plt.figure(); plt.hist(x[~np.isnan(x)], bins=60, color=c_for(name))
            plt.title(f"{name} — {title}"); plt.xlabel(xlabel); plt.ylabel("count")
            plt.tight_layout(); pdf.savefig(); plt.close()
    _hist("n_fragment","n_fragment","n_fragment")
    _hist("tsse","TSSe","TSSe")
    _hist("frac_mito","frac_mito","frac_mito")
    _hist("frac_dup","frac_dup","frac_dup")
    _hist("doublet_score","doublet_score","doublet_score")
    _hist("nucleosome_signal","nucleosome_signal","nucleosome_signal")

pdf.close()
print("[INFO] Full PDF ->", pdf_path)

# ---------- FRiP (optional; requires UNION_PEAKS) ----------
has_union = bool(UNION_PEAKS) and os.path.exists(UNION_PEAKS)
if has_union:
    for name in DATASETS:
        h5 = latest_h5ad_for(name, OUTDIR)
        if not h5: continue
        print(f"[INFO] FRiP — {name}")
        ad = snap.read(h5) if hasattr(snap, "read") else snap.pp.load_file(h5)
        try:    snap.metrics.frip(ad, {"frip": UNION_PEAKS})
        except TypeError: snap.metrics.frip(ad, UNION_PEAKS)

        qc_path  = os.path.join(OUTDIR, f"{name}.qc_cells.csv")
        pass_path= os.path.join(OUTDIR, f"{name}.qc_cells.with_pass.csv")
        df  = pd.read_csv(qc_path, index_col=0)
        dpf = pd.read_csv(pass_path, index_col=0) if os.path.exists(pass_path) else df.copy()
        frip= np.asarray(ad.obs["frip"], dtype=float)
        df["frip"] = frip

        nfr  = df["n_fragment"].to_numpy(float) if "n_fragment" in df else np.full(len(df), np.nan)
        tsse = df["tsse"].to_numpy(float)       if "tsse"       in df else np.full(len(df), np.nan)
        mito = df["frac_mito"].to_numpy(float)  if "frac_mito"  in df else np.full(len(df), np.nan)

        mask = (nfr>=MIN_COUNTS) & (nfr<=MAX_COUNTS) & (tsse>=MIN_TSSE_FOR_PASS)
        if not np.isnan(mito).all(): mask &= (mito<=MAX_MITO)
        mask &= (frip>=FRIP_THRESH)
        dpf["qc_pass"] = mask

        df.to_csv(qc_path); dpf.to_csv(pass_path)

        # FRiP histogram (colored)
        plt.figure(figsize=(6.5,4))
        v = df["frip"].to_numpy(float)
        plt.hist(v[~np.isnan(v)], bins=60, color=c_for(name))
        plt.axvline(FRIP_THRESH, ls="--", c="red", label=f"FRiP ≥ {FRIP_THRESH}")
        plt.title(f"{name} — FRiP histogram"); plt.xlabel("FRiP"); plt.ylabel("count"); plt.legend()
        plt.tight_layout(); plt.savefig(os.path.join(OUTDIR, f"{name}.frip_hist.png"), dpi=150); plt.close()
    print("[INFO] FRiP computed & qc_pass updated.")

# ---------- Summary CSV ----------
rows=[]
for name in DATASETS:
    qcf = os.path.join(OUTDIR, f"{name}.qc_cells.csv")
    qpf = os.path.join(OUTDIR, f"{name}.qc_cells.with_pass.csv")
    if not os.path.exists(qcf): continue
    df = pd.read_csv(qcf, index_col=0)
    dpf= pd.read_csv(qpf, index_col=0) if os.path.exists(qpf) else df
    def med(col):
        return float(np.nanmedian(df[col])) if col in df.columns and not np.isnan(df[col]).all() else np.nan
    rows.append({
        "dataset": name,
        "cells": int(df.shape[0]),
        "pass_cells": int(np.nan_to_num(dpf["qc_pass"]).sum()) if "qc_pass" in dpf else 0,
        "pass_rate": float(np.nanmean(dpf["qc_pass"])) if "qc_pass" in dpf else np.nan,
        "n_fragment_median": med("n_fragment"),
        "tsse_median": med("tsse"),
        "frac_mito_median": med("frac_mito"),
        "frac_dup_median" : med("frac_dup"),
        "doublet_score_median": med("doublet_score"),
        "nucleosome_signal_median": med("nucleosome_signal"),
    })
summary = pd.DataFrame(rows).sort_values("dataset")
summary.to_csv(os.path.join(OUTDIR, "benchmark_summary.csv"), index=False)
print("[INFO] Summary CSV ->", os.path.join(OUTDIR, "benchmark_summary.csv"))

# ---------- Scoring（FAIR & DISPLAY/SeekGene-lean） ----------
def safe_nanmedian(a):
    arr = np.asarray(a, dtype=float)
    return np.nan if np.isnan(arr).all() else float(np.nanmedian(arr))

def norm_hi(x):
    x = x.astype(float)
    if x.isna().all(): return x
    lo,hi=np.nanmin(x),np.nanmax(x)
    return 100*(x-lo)/(hi-lo) if hi>lo else np.full_like(x, np.nan, dtype=float)

def norm_lo(x):
    x = x.astype(float)
    if x.isna().all(): return x
    lo,hi=np.nanmin(x),np.nanmax(x)
    return 100*(hi-x)/(hi-lo) if hi>lo else np.full_like(x, np.nan, dtype=float)

def med_table(names):
    rows=[]
    for n in names:
        qcf = os.path.join(OUTDIR, f"{n}.qc_cells.csv")
        qpf = os.path.join(OUTDIR, f"{n}.qc_cells.with_pass.csv")
        if not os.path.exists(qcf): continue
        df  = pd.read_csv(qcf, index_col=0)
        dpf = pd.read_csv(qpf, index_col=0) if os.path.exists(qpf) else df
        med = {
            "n_fragment"        : safe_nanmedian(df["n_fragment"]) if "n_fragment" in df else np.nan,
            "tsse"              : safe_nanmedian(df["tsse"])       if "tsse"       in df else np.nan,
            "frip"              : safe_nanmedian(df["frip"])       if "frip"       in df else np.nan,
            "frac_mito"         : safe_nanmedian(df["frac_mito"])  if "frac_mito"  in df else np.nan,
            "frac_dup"          : safe_nanmedian(df["frac_dup"])   if "frac_dup"   in df else np.nan,
            "doublet_score"     : safe_nanmedian(df["doublet_score"]) if "doublet_score" in df else np.nan,
            "nucleosome_signal" : safe_nanmedian(df["nucleosome_signal"]) if "nucleosome_signal" in df else np.nan,
            "pass_rate"         : float(np.nanmean(dpf["qc_pass"])) if "qc_pass" in dpf else np.nan,
        }
        med["log10_nfrag"] = np.log10(med["n_fragment"]) if np.isfinite(med["n_fragment"]) and med["n_fragment"]>0 else np.nan
        rows.append((n, med))
    return pd.DataFrame(dict(rows)).T

tab = med_table(DATASETS.keys())

FAIR    = {"S_tsse":0.40, "S_frip":0.35, "S_nfrag":0.20, "S_mito":0.05, "S_passrate":0.0, "S_doublet":0.0, "S_dup":0.0}
DISPLAY = {"S_tsse":0.45, "S_frip":0.40, "S_nfrag":0.15, "S_mito":0.0, "S_passrate":0.0, "S_doublet":0.0, "S_dup":0.0}

def build_subscores(table):
    S = pd.DataFrame(index=table.index)
    S["S_tsse"]     = norm_hi(table["tsse"])
    S["S_frip"]     = norm_hi(table["frip"])
    S["S_nfrag"]    = norm_hi(table["log10_nfrag"])
    S["S_mito"]     = norm_lo(table["frac_mito"])
    S["S_dup"]      = norm_lo(table["frac_dup"])
    S["S_doublet"]  = norm_lo(table["doublet_score"])
    S["S_passrate"] = norm_hi(table["pass_rate"])
    return S

def row_weighted_total(row, weights: dict):
    vals = row.reindex(weights.keys()); valid = vals.notna()
    if not valid.any(): return np.nan
    w = pd.Series(weights)[valid]; w = w / w.sum()
    return float((vals[valid] * w).sum())

S  = build_subscores(tab)
Sf = pd.DataFrame(index=S.index); Sd = pd.DataFrame(index=S.index)
Sf["Total_fair"]    = S.apply(lambda r: row_weighted_total(r, FAIR), axis=1)
Sd["Total_display"] = S.apply(lambda r: row_weighted_total(r, DISPLAY), axis=1)

# Ranking plots (PNG) -- also colored by dataset
def barplot(series, title, png):
    s = series.sort_values(ascending=False)
    fig, ax = plt.subplots(figsize=(6,3.6))
    y = np.arange(len(s))
    cols = [c_for(k) for k in s.index]
    ax.barh(y, s.values, color=cols, edgecolor="none")
    ax.set_yticks(y); ax.set_yticklabels(s.index)
    ax.invert_yaxis()
    ax.set_xlabel("Score (0–100)"); ax.set_title(title)
    ax.grid(True, axis="x", linestyle=":", linewidth=0.8, alpha=0.6)
    fig.tight_layout(); fig.savefig(png, dpi=150); plt.close(fig)

barplot(Sf["Total_fair"],    "Benchmark Ranking — FAIR",    os.path.join(OUTDIR, "ranking_fair.png"))
barplot(Sd["Total_display"], "Benchmark Ranking — DISPLAY (SeekGene-lean)", os.path.join(OUTDIR, "ranking_display.png"))
print("[INFO] Rankings ->", os.path.join(OUTDIR, "ranking_fair.png"), ",", os.path.join(OUTDIR, "ranking_display.png"))

# Bars PDF (colored bars)
summary = pd.read_csv(os.path.join(OUTDIR, "benchmark_summary.csv"))
order = list(summary["dataset"])
scoring_csv = os.path.join(OUTDIR, "benchmark_scoring.csv")
sc = pd.read_csv(scoring_csv, index_col=0).reindex(order)

bars_pdf = os.path.join(OUTDIR, "benchmark_bars.pdf")
with PdfPages(bars_pdf) as bpdf:
    fig = plt.figure(figsize=(8,4))
    plt.axis("off"); plt.text(0.02,0.9,"scATAC Benchmark — Barplots", fontsize=14)
    plt.text(0.02,0.7,f"Summary: {os.path.join(OUTDIR,'benchmark_summary.csv')}\nScoring: {scoring_csv}", fontsize=10)
    bpdf.savefig(fig); plt.close(fig)

    def page(series, title, ylabel=None, pct=False):
        s = pd.Series(series.values, index=order).astype(float)
        if pct: s = s*100.0
        fig, ax = plt.subplots(figsize=(8,4))
        x = np.arange(len(s))
        cols = [c_for(k) for k in s.index]
        ax.bar(x, s.values, color=cols, edgecolor="none")
        ax.set_title(title); ax.set_ylabel(ylabel or ("Percent (%)" if pct else "")); ax.set_xlabel("")
        ax.set_xticks(x); ax.set_xticklabels(order, rotation=25, ha="right")
        for i, v in enumerate(s.values):
            if np.isfinite(v):
                ax.annotate(f"{v:.2f}", (x[i], v), xytext=(0,3), textcoords="offset points",
                            ha="center", va="bottom", fontsize=8)
        fig.tight_layout(); bpdf.savefig(fig); plt.close(fig)

    page(summary["n_fragment_median"], "Median Unique Fragments", "Fragments")
    page(np.log10(np.clip(summary["n_fragment_median"],1,None)), "Median log10(Unique Fragments)", "log10(Fragments)")
    page(summary["tsse_median"], "Median TSSe", "TSSe")
    page(summary["frac_mito_median"], "Median Mito Fraction", "Mito (%)", pct=True)
    page(summary["frac_dup_median"], "Median Duplicate Fraction", "Dup (%)", pct=True)
    page(summary["pass_rate"], "Pass Rate", "Pass (%)", pct=True)

    if "Total_fair" in sc.columns:     page(sc["Total_fair"], "Total Score — FAIR", "Score (0–100)")
    if "Total_display" in sc.columns:  page(sc["Total_display"], "Total Score — DISPLAY", "Score (0–100)")

print("[INFO] Barplots PDF ->", bars_pdf)

# ======================= Scoring Visualizations: multi-page report + standalone PDF =======================
# This reuses the SCATAC_PALETTE defined above
METRIC_COLORS = {
    "S_tsse": "#1f77b4", "S_frip": "#ff7f0e", "S_nfrag": "#2ca02c",
    "S_mito": "#d62728", "S_dup": "#9467bd", "S_doublet": "#8c564b", "S_passrate": "#7f7f7f",
}
def _safe_cols(df, cols): return [c for c in cols if c in df.columns]
def _order_for(series: pd.Series):
    try: return series.sort_values(ascending=False).index.tolist()
    except Exception: return series.index.tolist()
def _weighted_contrib(row, weights: dict):
    vals = row.reindex(weights.keys()); ok = vals.notna()
    if not ok.any(): return pd.Series(index=weights.keys(), dtype=float)
    w = pd.Series(weights).reindex(vals.index)[ok]; w = w / w.sum()
    return (vals[ok] * w)

def make_scoring_report(tab, S, Sf, Sd, FAIR, DISPLAY, outdir=OUTDIR,
                        palette=SCATAC_PALETTE, metric_colors=METRIC_COLORS):
    order_fair    = _order_for(Sf["Total_fair"]) if "Total_fair" in Sf else list(S.index)
    order_display = _order_for(Sd["Total_display"]) if "Total_display" in Sd else list(S.index)
    order_union   = list(dict.fromkeys(order_display + order_fair))

    report_pdf = os.path.join(outdir, "benchmark_scoring_report.pdf")
    os.makedirs(outdir, exist_ok=True)
    with PdfPages(report_pdf) as rep:
        # Cover page
        fig = plt.figure(figsize=(8.5, 5)); plt.axis("off")
        txt = "scATAC Benchmark — Scoring Report\n\n"
        if "Total_display" in Sd:
            top = Sd["Total_display"].sort_values(ascending=False).head(3)
            txt += "Top (DISPLAY): " + ", ".join([f"{k}={v:.1f}" for k,v in top.items()]) + "\n"
        if "Total_fair" in Sf:
            top = Sf["Total_fair"].sort_values(ascending=False).head(3)
            txt += "Top (FAIR): " + ", ".join([f"{k}={v:.1f}" for k,v in top.items()]) + "\n"
        plt.text(0.02, 0.85, txt, fontsize=14, va="top"); rep.savefig(fig); plt.close(fig)

        # Subscore heatmap
        cols = _safe_cols(S, ["S_tsse","S_frip","S_nfrag","S_mito","S_passrate","S_doublet","S_dup"])
        if len(cols) > 0:
            M = S.loc[order_union, cols].astype(float)
            fig, ax = plt.subplots(figsize=(max(6, 0.8*len(cols)+2), max(3, 0.45*len(order_union)+1)))
            im = ax.imshow(M.fillna(0).values, cmap="viridis", vmin=0, vmax=100, aspect="auto")
            # NaN mask
            for r in range(M.shape[0]):
                for c in range(M.shape[1]):
                    if pd.isna(M.iat[r,c]): ax.add_patch(plt.Rectangle((c-0.5, r-0.5), 1, 1, color="white", zorder=2))
            plt.colorbar(im, fraction=0.046, pad=0.04, label="Subscore (0–100)")
            ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, rotation=30, ha="right")
            ax.set_yticks(range(len(order_union))); ax.set_yticklabels(order_union)
            ax.set_title("Subscores heatmap")
            plt.tight_layout(); rep.savefig(fig); fig.savefig(os.path.join(outdir, "subscores_heatmap.pdf"), format="pdf"); plt.close(fig)

        # Lollipop ranking plot (FAIR/DISPLAY), lines colored by dataset
        def lollipop(series, title, savepath):
            s = series.dropna().sort_values(ascending=True)
            fig, ax = plt.subplots(figsize=(6.6, 3.8))
            y = np.arange(len(s)); cols = [palette.get(k, "#999999") for k in s.index]
            ax.hlines(y, xmin=0, xmax=s.values, color=cols, linewidth=2)
            ax.plot(s.values, y, "o", color="black", markersize=4)
            ax.set_yticks(y); ax.set_yticklabels(s.index)
            ax.set_xlabel("Score (0–100)"); ax.set_title(title)
            ax.grid(True, axis='x', linestyle=':', linewidth=0.8, alpha=0.6)
            plt.tight_layout(); rep.savefig(fig); fig.savefig(savepath, format="pdf"); plt.close(fig)

        if "Total_fair" in Sf:    lollipop(Sf["Total_fair"],    "Ranking — FAIR",    os.path.join(outdir, "ranking_fair.pdf"))
        if "Total_display" in Sd: lollipop(Sd["Total_display"], "Ranking — DISPLAY (SeekGene-lean)", os.path.join(outdir, "ranking_display.pdf"))

        # TSSe vs FRiP (marker color = dataset color, marker size ~ log10 frags)
        fig, ax = plt.subplots(figsize=(6.2, 4.8))
        x = tab["tsse"].astype(float); y = tab["frip"].astype(float); s = tab["log10_nfrag"].astype(float)
        if s.notna().any(): s_norm = (s - np.nanmin(s)) / max(1e-6, (np.nanmax(s) - np.nanmin(s)))
        else:               s_norm = s*0+0.5
        sizes = 40 + 160 * s_norm.fillna(0.3)
        for ds in order_union:
            if ds not in x.index: continue
            if pd.isna(x.loc[ds]) or pd.isna(y.loc[ds]): continue
            ax.scatter(x.loc[ds], y.loc[ds], s=float(sizes.loc[ds]),
                       color=palette.get(ds, "#999999"), alpha=0.8, edgecolor="white", linewidth=0.6, zorder=3)
            ax.text(x.loc[ds], y.loc[ds], f"  {ds}", va="center", fontsize=8)
        ax.set_xlabel("Median TSSe"); ax.set_ylabel("Median FRiP")
        ax.set_title("TSSe vs FRiP (marker size ~ log10 unique fragments)")
        ax.grid(True, linestyle=":", linewidth=0.8, alpha=0.6)
        plt.tight_layout(); rep.savefig(fig); fig.savefig(os.path.join(outdir, "tsse_vs_frip.pdf"), format="pdf"); plt.close(fig)

        # Weighted contribution stacked plot (border colored by dataset)
        def stacked_contrib(weights, title, fname):
            metrics = [k for k in weights.keys() if k in S.columns]
            if not metrics: return
            C = []
            for ds in order_union:
                contrib = _weighted_contrib(S.loc[ds, :], weights).reindex(metrics)
                C.append(contrib.values)
            C = np.array(C).T
            fig, ax = plt.subplots(figsize=(max(6, 0.6*len(order_union)+1.8), 4.0))
            bottoms = np.zeros(C.shape[1]); x = np.arange(C.shape[1])
            for i, m in enumerate(metrics):
                ax.bar(x, C[i], bottom=bottoms, color=metric_colors.get(m, "#888888"), edgecolor="none", label=m)
                bottoms += C[i]
            for i, ds in enumerate(order_union):
                ax.bar(x[i], bottoms[i], fill=False, edgecolor=palette.get(ds, "#666666"), linewidth=1.0)
            ax.set_xticks(x); ax.set_xticklabels(order_union, rotation=25, ha="right")
            ax.set_ylabel("Weighted contribution to total (0–100)")
            ax.set_title(title)
            ax.legend(frameon=False, ncol=min(3, len(metrics)))
            ax.grid(True, axis='y', linestyle=':', linewidth=0.8, alpha=0.6)
            plt.tight_layout(); rep.savefig(fig); fig.savefig(fname, format="pdf"); plt.close(fig)

        if "Total_display" in Sd: stacked_contrib(DISPLAY, "Contribution breakdown — DISPLAY (SeekGene-lean)", os.path.join(outdir, "contrib_display.pdf"))
        if "Total_fair" in Sf:    stacked_contrib(FAIR,    "Contribution breakdown — FAIR",                    os.path.join(outdir, "contrib_fair.pdf"))

        # Radar plot (one page per dataset)
        cols = _safe_cols(S, ["S_tsse","S_frip","S_nfrag","S_mito","S_passrate","S_doublet","S_dup"])
        if len(cols) >= 3:
            theta = np.linspace(0, 2*np.pi, len(cols)+1)
            for ds in order_union:
                vals = S.loc[ds, cols].astype(float).values
                fig = plt.figure(figsize=(4.8, 4.8)); ax = plt.subplot(111, polar=True)
                v = np.append(vals, vals[0])
                ax.plot(theta, v, color=palette.get(ds, "#999999"), linewidth=2)
                ax.fill(theta, v, color=palette.get(ds, "#999999"), alpha=0.20)
                ax.set_xticks(theta[:-1]); ax.set_xticklabels(cols, fontsize=8)
                ax.set_yticks([20,40,60,80]); ax.set_yticklabels(["20","40","60","80"])
                ax.set_ylim(0,100); ax.set_title(f"Radar — {ds}", va="bottom")
                plt.tight_layout(); rep.savefig(fig); plt.close(fig)

    print("[INFO] Scoring report PDF ->", report_pdf)

# Generate scoring report
make_scoring_report(tab, S, Sf, Sd, FAIR, DISPLAY, outdir=OUTDIR, palette=SCATAC_PALETTE)


# In[ ]:



