#!/usr/bin/env python3
"""Apply strict skin filtering to AlphaGenome results and export per-SNV boxplot stats."""

from __future__ import annotations

import argparse
import logging
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


LOGGER = logging.getLogger(__name__)

STRICT_SKIN_TERMS = (
    "skin",
    "suprapubic",
    "lower leg",
    "epiderm",
    "dermis",
    "keratinocyte",
    "melanocyte",
    "foreskin",
)

ESSENTIAL_COLUMNS = [
    "input_snv",
    "gene_name",
    "track_name",
    "biosample_name",
    "gtex_tissue",
    "raw_score",
    "quantile_score",
    "chrom",
    "pos",
    "ref",
    "alt",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Re-filter AlphaGenome RNA/ATAC results with a stricter skin-only rule "
            "and export per-SNV summary statistics suitable for boxplot rendering."
        )
    )
    parser.add_argument("--rnaseq-all", required=True, type=Path, help="AlphaGenome RNA all-scores CSV.")
    parser.add_argument("--atac-all", required=True, type=Path, help="AlphaGenome ATAC all-scores CSV.")
    parser.add_argument("--output-dir", required=True, type=Path, help="Output directory.")
    parser.add_argument(
        "--chunksize",
        type=int,
        default=200_000,
        help="CSV chunk size used for streaming large RNA files.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"],
        help="Logging verbosity.",
    )
    return parser.parse_args()


def configure_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def strict_skin_mask(df: pd.DataFrame) -> pd.Series:
    biosample = df.get("biosample_name", pd.Series("", index=df.index)).fillna("").astype(str).str.lower()
    gtex = df.get("gtex_tissue", pd.Series("", index=df.index)).fillna("").astype(str).str.lower()
    combined = biosample + " || " + gtex
    mask = pd.Series(False, index=df.index)
    for term in STRICT_SKIN_TERMS:
        mask = mask | combined.str.contains(term, regex=False)
    return mask


def append_filtered_rows(df: pd.DataFrame, output_csv: Path, first_write: bool) -> None:
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    write_df = df.copy()
    missing_cols = [col for col in ESSENTIAL_COLUMNS if col not in write_df.columns]
    for col in missing_cols:
        write_df[col] = np.nan
    write_df = write_df.loc[:, ESSENTIAL_COLUMNS]
    write_df.to_csv(output_csv, mode="w" if first_write else "a", header=first_write, index=False)


def update_stats(stats: dict[str, list[float]], snvs: pd.Series, values: pd.Series) -> None:
    for snv, value in zip(snvs.astype(str), values.astype(float)):
        if np.isfinite(value):
            stats[snv].append(float(value))


def summarize_value_lists(value_map: dict[str, list[float]], prefix: str) -> pd.DataFrame:
    rows = []
    for snv, values in value_map.items():
        arr = np.asarray(values, dtype=float)
        if arr.size == 0:
            continue
        rows.append(
            {
                "input_snv": snv,
                f"{prefix}_n": int(arr.size),
                f"{prefix}_mean": float(arr.mean()),
                f"{prefix}_std": float(arr.std(ddof=0)),
                f"{prefix}_min": float(arr.min()),
                f"{prefix}_q1": float(np.quantile(arr, 0.25)),
                f"{prefix}_median": float(np.quantile(arr, 0.50)),
                f"{prefix}_q3": float(np.quantile(arr, 0.75)),
                f"{prefix}_max": float(arr.max()),
            }
        )
    return pd.DataFrame(rows)


def process_score_file(
    csv_path: Path,
    filtered_output_csv: Path,
    prefix: str,
    chunksize: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    LOGGER.info("Strict-skin filtering %s: %s", prefix, csv_path)
    abs_quantile_map: dict[str, list[float]] = defaultdict(list)
    raw_abs_map: dict[str, list[float]] = defaultdict(list)
    track_counter: dict[str, int] = defaultdict(int)
    biosample_counter: dict[str, int] = defaultdict(int)
    first_write = True
    kept_rows = 0

    for chunk_idx, chunk in enumerate(pd.read_csv(csv_path, chunksize=chunksize), start=1):
        mask = strict_skin_mask(chunk)
        keep = chunk.loc[mask].copy()
        if keep.empty:
            continue

        kept_rows += len(keep)
        keep["abs_quantile_score"] = pd.to_numeric(keep["quantile_score"], errors="coerce").abs()
        keep["abs_raw_score"] = pd.to_numeric(keep["raw_score"], errors="coerce").abs()
        update_stats(abs_quantile_map, keep["input_snv"], keep["abs_quantile_score"])
        update_stats(raw_abs_map, keep["input_snv"], keep["abs_raw_score"])

        for name, count in keep["track_name"].fillna("").astype(str).value_counts().items():
            if name:
                track_counter[name] += int(count)
        for name, count in keep["biosample_name"].fillna("").astype(str).value_counts().items():
            if name:
                biosample_counter[name] += int(count)

        append_filtered_rows(keep, filtered_output_csv, first_write=first_write)
        first_write = False

        if chunk_idx % 20 == 0:
            LOGGER.info("%s chunks processed: %d | kept rows so far: %d", prefix, chunk_idx, kept_rows)

    abs_df = summarize_value_lists(abs_quantile_map, prefix=f"{prefix}_abs_quantile")
    raw_df = summarize_value_lists(raw_abs_map, prefix=f"{prefix}_abs_raw")
    summary_df = abs_df.merge(raw_df, on="input_snv", how="outer")

    metadata_rows = []
    for track_name, count in sorted(track_counter.items(), key=lambda item: item[1], reverse=True):
        metadata_rows.append({"type": "track_name", "name": track_name, "rows": count})
    for biosample_name, count in sorted(biosample_counter.items(), key=lambda item: item[1], reverse=True):
        metadata_rows.append({"type": "biosample_name", "name": biosample_name, "rows": count})
    metadata_df = pd.DataFrame(metadata_rows)

    LOGGER.info("%s strict-skin kept rows: %d", prefix, kept_rows)
    return summary_df, metadata_df


def main() -> None:
    args = parse_args()
    configure_logging(args.log_level)

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    rna_filtered_csv = output_dir / "alphagenome_rnaseq_strict_skin_scores.csv"
    atac_filtered_csv = output_dir / "alphagenome_atac_strict_skin_scores.csv"

    rna_summary, rna_meta = process_score_file(
        csv_path=args.rnaseq_all,
        filtered_output_csv=rna_filtered_csv,
        prefix="rna",
        chunksize=args.chunksize,
    )
    atac_summary, atac_meta = process_score_file(
        csv_path=args.atac_all,
        filtered_output_csv=atac_filtered_csv,
        prefix="atac",
        chunksize=args.chunksize,
    )

    merged = pd.merge(rna_summary, atac_summary, on="input_snv", how="outer").sort_values("input_snv")
    merged.to_csv(output_dir / "alphagenome_strict_skin_boxplot_stats.csv", index=False)
    rna_meta.to_csv(output_dir / "alphagenome_rnaseq_strict_skin_metadata.csv", index=False)
    atac_meta.to_csv(output_dir / "alphagenome_atac_strict_skin_metadata.csv", index=False)

    LOGGER.info("Wrote strict-skin boxplot stats to %s", output_dir / "alphagenome_strict_skin_boxplot_stats.csv")


if __name__ == "__main__":
    main()
