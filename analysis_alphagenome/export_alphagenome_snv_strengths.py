#!/usr/bin/env python3
"""Export compact AlphaGenome SNV strength summary."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd


LOGGER = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Read the merged VESPA-AlphaGenome summary table and export a compact "
            "3-column CSV with input_snv, RNA strength, and ATAC strength."
        )
    )
    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Path to vespa_alphagenome_merged_summary.csv.",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output CSV path.",
    )
    parser.add_argument(
        "--rna-column",
        default="rna_max_abs_quantile",
        help="Column used as AlphaGenome RNA strength.",
    )
    parser.add_argument(
        "--atac-column",
        default="atac_max_abs_quantile",
        help="Column used as AlphaGenome ATAC strength.",
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


def main() -> None:
    args = parse_args()
    configure_logging(args.log_level)

    LOGGER.info("Loading merged summary: %s", args.input)
    df = pd.read_csv(args.input)

    required = ["input_snv", args.rna_column, args.atac_column]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    out = df.loc[:, required].copy()
    out = out.rename(
        columns={
            args.rna_column: "alphagenome_rna_strength",
            args.atac_column: "alphagenome_atac_strength",
        }
    )
    out = out.sort_values("input_snv").reset_index(drop=True)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    LOGGER.info("Wrote %d rows to %s", len(out), args.output)


if __name__ == "__main__":
    main()
