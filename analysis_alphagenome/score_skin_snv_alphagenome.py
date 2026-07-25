#!/usr/bin/env python3
"""Score SNV effects in skin-related AlphaGenome tracks."""

from __future__ import annotations

import argparse
import csv
import logging
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


LOGGER = logging.getLogger(__name__)
SNV_PATTERN = re.compile(r"^(chr[^:]+|[^:]+):(\d+)_([ACGTN]+)>([ACGTN]+)$", re.IGNORECASE)
DEFAULT_SKIN_TERMS = ("skin", "suprapubic", "epidermis", "dermis", "keratinocyte", "fibroblast")
SCORER_ALIASES = {
    "rna_seq": "RNA_SEQ",
    "cage": "CAGE",
    "procap": "PROCAP",
    "atac": "ATAC",
    "dnase": "DNASE",
    "chip_histone": "CHIP_HISTONE",
    "chip_tf": "CHIP_TF",
    "polyadenylation": "POLYADENYLATION",
    "splice_sites": "SPLICE_SITES",
    "splice_site_usage": "SPLICE_SITE_USAGE",
    "splice_junctions": "SPLICE_JUNCTIONS",
}


@dataclass(frozen=True)
class SNV:
    variant_id: str
    chrom: str
    pos: int
    ref: str
    alt: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Read an SNV list and use AlphaGenome to score variant effects, "
            "then export all scores plus a skin-track subset."
        )
    )
    parser.add_argument(
        "--snv-list",
        required=True,
        type=Path,
        help=(
            "Input SNV table. Supported forms: columns variant_id/CHROM/POS/REF/ALT, "
            "columns chrom/pos/ref/alt, or a single SNV column like chr1:123_A>G."
        ),
    )
    parser.add_argument("--output-all", required=True, type=Path, help="Output CSV for all AlphaGenome scores.")
    parser.add_argument(
        "--output-skin",
        required=True,
        type=Path,
        help="Output CSV filtered to skin-related AlphaGenome metadata rows.",
    )
    parser.add_argument(
        "--api-key",
        default=None,
        help="AlphaGenome API key. If omitted, ALPHAGENOME_API_KEY is used.",
    )
    parser.add_argument(
        "--sequence-length",
        default="1MB",
        choices=["16KB", "100KB", "500KB", "1MB"],
        help="Sequence context centered on each SNV.",
    )
    parser.add_argument(
        "--organism",
        default="human",
        choices=["human", "mouse"],
        help="Genome/model organism. Human uses hg38 coordinates; mouse uses mm10 coordinates.",
    )
    parser.add_argument(
        "--api-timeout",
        type=float,
        default=120.0,
        help="AlphaGenome API timeout in seconds.",
    )
    parser.add_argument(
        "--skin-term",
        action="append",
        dest="skin_terms",
        default=None,
        help="Case-insensitive term used to keep skin rows. Can be supplied multiple times.",
    )
    parser.add_argument(
        "--scorer",
        action="append",
        dest="scorers",
        default=None,
        choices=[
            "rna_seq",
            "cage",
            "procap",
            "atac",
            "dnase",
            "chip_histone",
            "chip_tf",
            "polyadenylation",
            "splice_sites",
            "splice_site_usage",
            "splice_junctions",
        ],
        help="AlphaGenome scorer to run. Can be supplied multiple times. Defaults to all recommended scorers.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=100,
        help="Number of SNVs per AlphaGenome batch.",
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=5,
        help="Parallel workers used by AlphaGenome score_variants.",
    )
    parser.add_argument(
        "--plot-dir",
        type=Path,
        default=None,
        help="Optional output directory for summary figures.",
    )
    parser.add_argument(
        "--failed-snvs",
        type=Path,
        default=None,
        help="Optional CSV path for SNVs that still fail after single-SNV retry.",
    )
    parser.add_argument(
        "--top-n-plot",
        type=int,
        default=30,
        help="Number of top SNVs shown in bar plots.",
    )
    parser.add_argument(
        "--no-progress",
        action="store_true",
        help="Disable AlphaGenome progress bars.",
    )
    parser.add_argument(
        "--mock-client",
        action="store_true",
        help="Use deterministic mock scoring for local smoke tests without alphagenome or an API key.",
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
    logging.getLogger("fontTools").setLevel(logging.WARNING)
    logging.getLogger("matplotlib").setLevel(logging.WARNING)


def normalize_header(name: str) -> str:
    return name.strip().lower().replace("#", "").replace(" ", "").replace("-", "_")


def sniff_delimiter(path: Path) -> str:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        sample = handle.read(4096)
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters="\t,;")
        return dialect.delimiter
    except csv.Error:
        return "\t" if "\t" in sample else ","


def parse_snv_text(text: str) -> tuple[str, int, str, str]:
    match = SNV_PATTERN.match(text.strip())
    if match is None:
        raise ValueError(f"Unsupported SNV format: {text}")
    chrom, pos, ref, alt = match.groups()
    return chrom, int(pos), ref.upper(), alt.upper()


def alphagenome_variant_id(snv: SNV) -> str:
    return f"{snv.chrom}:{snv.pos}:{snv.ref}>{snv.alt}"


def first_present(row: dict[str, str], names: Sequence[str]) -> str | None:
    for name in names:
        value = row.get(name)
        if value is not None and str(value).strip():
            return str(value).strip()
    return None


def read_snvs(path: Path) -> list[SNV]:
    delimiter = sniff_delimiter(path)
    snvs: list[SNV] = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError(f"Input file has no header: {path}")

        original_names = list(reader.fieldnames)
        header_is_snv_list = len(original_names) == 1
        if header_is_snv_list:
            try:
                chrom, pos, ref, alt = parse_snv_text(original_names[0])
                snvs.append(SNV(variant_id=original_names[0], chrom=chrom, pos=pos, ref=ref, alt=alt))
            except ValueError:
                header_is_snv_list = False

        normalized_names = [normalize_header(name) for name in reader.fieldnames]
        reader.fieldnames = normalized_names

        for row_number, row in enumerate(reader, start=2):
            row = {normalize_header(key): value for key, value in row.items() if key is not None}
            if header_is_snv_list:
                variant_id = first_present(row, (normalized_names[0],))
                if variant_id:
                    chrom, pos, ref, alt = parse_snv_text(variant_id)
                    snvs.append(SNV(variant_id=variant_id, chrom=chrom, pos=pos, ref=ref, alt=alt))
                continue

            chrom = first_present(row, ("chrom", "chr", "chromosome", "chromosome_name"))
            pos = first_present(row, ("pos", "position", "start"))
            ref = first_present(row, ("ref", "reference", "reference_bases", "reference_allele"))
            alt = first_present(row, ("alt", "alternate", "alternate_bases", "alternate_allele"))
            variant_id = first_present(row, ("variant_id", "id", "snv", "variant", "name"))

            if chrom and pos and ref and alt:
                snv = SNV(
                    variant_id=variant_id or f"{chrom}:{int(pos)}_{ref.upper()}>{alt.upper()}",
                    chrom=chrom,
                    pos=int(pos),
                    ref=ref.upper(),
                    alt=alt.upper(),
                )
            elif variant_id:
                parsed_chrom, parsed_pos, parsed_ref, parsed_alt = parse_snv_text(variant_id)
                snv = SNV(
                    variant_id=variant_id,
                    chrom=parsed_chrom,
                    pos=parsed_pos,
                    ref=parsed_ref,
                    alt=parsed_alt,
                )
            else:
                raise ValueError(
                    f"Row {row_number} must contain CHROM/POS/REF/ALT columns or a parseable SNV/variant_id column."
                )

            if len(snv.ref) != 1 or len(snv.alt) != 1:
                raise ValueError(f"Row {row_number} is not an SNV: {snv.variant_id}")
            if snv.ref == snv.alt:
                raise ValueError(f"Row {row_number} has identical REF and ALT: {snv.variant_id}")
            snvs.append(snv)

    if not snvs:
        raise ValueError(f"No SNVs found in {path}")
    return snvs


def load_alphagenome_model(api_key: str, timeout: float):
    try:
        from alphagenome.models import dna_client
    except ImportError as exc:
        raise RuntimeError("alphagenome is not installed. Install it before running real API scoring.") from exc
    return dna_client.create(api_key, timeout=timeout)


def score_snvs_with_alphagenome(
    snvs: Sequence[SNV],
    api_key: str,
    sequence_length: str,
    organism_name: str,
    scorer_names: Sequence[str] | None,
    api_timeout: float,
    batch_size: int,
    max_workers: int,
    progress_bar: bool,
    failed_snvs_path: Path | None,
) -> list[dict[str, object]]:
    try:
        from alphagenome.data import genome
        from alphagenome.models import dna_client, variant_scorers
    except ImportError as exc:
        raise RuntimeError("alphagenome is not installed. Install it before running real API scoring.") from exc

    model = load_alphagenome_model(api_key, timeout=api_timeout)
    sequence_key = f"SEQUENCE_LENGTH_{sequence_length}"
    sequence_size = dna_client.SUPPORTED_SEQUENCE_LENGTHS[sequence_key]
    organism_map = {
        "human": dna_client.Organism.HOMO_SAPIENS,
        "mouse": dna_client.Organism.MUS_MUSCULUS,
    }
    organism = organism_map[organism_name]

    all_scorers = variant_scorers.RECOMMENDED_VARIANT_SCORERS
    if scorer_names:
        scorers = [all_scorers[SCORER_ALIASES[name]] for name in scorer_names]
    else:
        scorers = list(all_scorers.values())
    selected_scorers = [
        scorer
        for scorer in scorers
        if (organism.value in variant_scorers.SUPPORTED_ORGANISMS[scorer.base_variant_scorer])
        and not (
            scorer.requested_output == dna_client.OutputType.PROCAP
            and organism == dna_client.Organism.MUS_MUSCULUS
        )
    ]

    records: list[dict[str, object]] = []
    failed_records: list[dict[str, object]] = []
    total = len(snvs)

    def append_tidy_scores(score_groups: Sequence[Sequence[object]], scored_snvs: Sequence[SNV], label: str) -> None:
        nonempty_scores = [[adata for adata in variant_scores if getattr(adata, "n_obs", 0) > 0] for variant_scores in score_groups]
        nonempty_scores = [variant_scores for variant_scores in nonempty_scores if variant_scores]
        if not nonempty_scores:
            LOGGER.warning("AlphaGenome returned no score rows for %s.", label)
            return

        tidy_df = variant_scorers.tidy_scores(nonempty_scores)
        if tidy_df is None or tidy_df.empty:
            LOGGER.warning("AlphaGenome returned no tidy score rows for %s.", label)
            return
        records.extend(attach_input_metadata(tidy_df.to_dict("records"), scored_snvs))

    for start in range(0, total, batch_size):
        chunk = list(snvs[start : start + batch_size])
        LOGGER.info("Scoring SNV batch %d-%d of %d", start + 1, start + len(chunk), total)
        variants = [
            genome.Variant(
                chromosome=snv.chrom,
                position=snv.pos,
                reference_bases=snv.ref,
                alternate_bases=snv.alt,
                name=snv.variant_id,
            )
            for snv in chunk
        ]
        intervals = [variant.reference_interval.resize(sequence_size) for variant in variants]
        try:
            batch_scores = model.score_variants(
                intervals=intervals,
                variants=variants,
                variant_scorers=selected_scorers,
                organism=organism,
                progress_bar=progress_bar,
                max_workers=max_workers,
            )
            append_tidy_scores(batch_scores, chunk, f"batch {start + 1}-{start + len(chunk)}")
        except Exception as exc:
            LOGGER.warning(
                "Batch %d-%d failed with %s. Retrying SNVs one by one.",
                start + 1,
                start + len(chunk),
                exc,
            )
            for snv, variant, interval in zip(chunk, variants, intervals):
                try:
                    single_scores = model.score_variant(
                        interval=interval,
                        variant=variant,
                        variant_scorers=selected_scorers,
                        organism=organism,
                    )
                    append_tidy_scores([single_scores], [snv], snv.variant_id)
                except Exception as single_exc:
                    LOGGER.error("SNV failed after single retry: %s | %s", snv.variant_id, single_exc)
                    failed_records.append(
                        {
                            "variant_id": snv.variant_id,
                            "chrom": snv.chrom,
                            "pos": snv.pos,
                            "ref": snv.ref,
                            "alt": snv.alt,
                            "error": str(single_exc),
                        }
                    )

    if failed_snvs_path is not None:
        write_records(
            failed_snvs_path,
            failed_records,
            fieldnames=["variant_id", "chrom", "pos", "ref", "alt", "error"],
        )
        if failed_records:
            LOGGER.warning("Wrote %d failed SNVs to %s", len(failed_records), failed_snvs_path)

    if not records:
        LOGGER.warning("AlphaGenome returned no score rows for the selected variants/scorers.")
        return []
    return records


def score_snvs_with_mock(snvs: Sequence[SNV]) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for index, snv in enumerate(snvs, start=1):
        signed_effect = round((index * 0.031) * (1 if snv.alt in {"A", "G"} else -1), 6)
        records.append(
            {
                "variant_id": snv.variant_id,
                "chrom": snv.chrom,
                "pos": snv.pos,
                "ref": snv.ref,
                "alt": snv.alt,
                "score_name": "mock_center_mask_diff_mean",
                "requested_output": "RNA_SEQ",
                "biosample_name": "skin of body, suprapubic",
                "ontology_curie": "UBERON:0036149",
                "value": signed_effect,
            }
        )
        records.append(
            {
                "variant_id": snv.variant_id,
                "chrom": snv.chrom,
                "pos": snv.pos,
                "ref": snv.ref,
                "alt": snv.alt,
                "score_name": "mock_center_mask_diff_mean",
                "requested_output": "RNA_SEQ",
                "biosample_name": "whole blood",
                "ontology_curie": "UBERON:0000178",
                "value": round(-signed_effect / 2, 6),
            }
        )
    return records


def attach_input_metadata(records: list[dict[str, object]], snvs: Sequence[SNV]) -> list[dict[str, object]]:
    annotated: list[dict[str, object]] = []
    snv_by_index = list(snvs)
    snv_by_alpha_id = {alphagenome_variant_id(snv): snv for snv in snv_by_index}
    snv_by_input_id = {snv.variant_id: snv for snv in snv_by_index}
    for record_index, record in enumerate(records):
        out_record = dict(record)
        if "variant_id" in out_record:
            record_variant_id = str(out_record["variant_id"])
            snv = snv_by_alpha_id.get(record_variant_id) or snv_by_input_id.get(record_variant_id)
        else:
            snv = snv_by_index[min(record_index, len(snv_by_index) - 1)]
        if snv is not None:
            out_record.setdefault("alphagenome_variant_id", alphagenome_variant_id(snv))
            out_record.setdefault("input_snv", snv.variant_id)
            out_record.setdefault("chrom", snv.chrom)
            out_record.setdefault("pos", snv.pos)
            out_record.setdefault("ref", snv.ref)
            out_record.setdefault("alt", snv.alt)
        annotated.append(out_record)
    return annotated


def row_matches_skin(row: dict[str, object], skin_terms: Iterable[str]) -> bool:
    searchable = " ".join("" if value is None else str(value) for value in row.values()).lower()
    return any(term.lower() in searchable for term in skin_terms)


def write_records(path: Path, records: Sequence[dict[str, object]], fieldnames: Sequence[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    output_fieldnames: list[str] = list(fieldnames or [])
    for record in records:
        for key in record:
            if key not in output_fieldnames:
                output_fieldnames.append(key)
    if not output_fieldnames:
        output_fieldnames = ["variant_id", "chrom", "pos", "ref", "alt", "score_name", "raw_score"]

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_fieldnames, extrasaction="ignore")
        writer.writeheader()
        for record in records:
            writer.writerow(record)


def plot_summary(records: Sequence[dict[str, object]], plot_dir: Path, top_n: int) -> None:
    if not records:
        LOGGER.warning("Skip plotting because there are no score records.")
        return

    import matplotlib.pyplot as plt
    import pandas as pd
    import seaborn as sns

    plot_dir.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(records)
    if "raw_score" not in df.columns:
        LOGGER.warning("Skip plotting because raw_score is not present.")
        return
    df["raw_score"] = pd.to_numeric(df["raw_score"], errors="coerce")
    df["abs_raw_score"] = df["raw_score"].abs()
    if "quantile_score" in df.columns:
        df["quantile_score"] = pd.to_numeric(df["quantile_score"], errors="coerce")

    sns.set_theme(style="whitegrid", context="paper", font="Arial")
    plt.rcParams.update(
        {
            "axes.edgecolor": "#222222",
            "axes.linewidth": 0.8,
            "axes.labelcolor": "#222222",
            "xtick.color": "#222222",
            "ytick.color": "#222222",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig, ax = plt.subplots(figsize=(4.2, 3.0))
    sns.histplot(df["raw_score"].dropna(), bins=40, color="#3B6EA8", edgecolor="white", linewidth=0.3, ax=ax)
    ax.set_xlabel("AlphaGenome raw score")
    ax.set_ylabel("Track-SNV count")
    sns.despine(ax=ax)
    fig.tight_layout()
    fig.savefig(plot_dir / "skin_raw_score_distribution.pdf")
    fig.savefig(plot_dir / "skin_raw_score_distribution.png", dpi=300)
    plt.close(fig)

    group_col = "input_snv" if "input_snv" in df.columns else "variant_id"
    top = (
        df.dropna(subset=["abs_raw_score"])
        .groupby(group_col, as_index=False)["abs_raw_score"]
        .max()
        .sort_values("abs_raw_score", ascending=False)
        .head(top_n)
    )
    if not top.empty:
        fig_height = max(3.2, 0.16 * len(top) + 1.0)
        fig, ax = plt.subplots(figsize=(5.2, fig_height))
        sns.barplot(data=top, y=group_col, x="abs_raw_score", color="#2E7D6B", ax=ax)
        ax.set_xlabel("Max absolute raw score")
        ax.set_ylabel("SNV")
        sns.despine(ax=ax)
        fig.tight_layout()
        fig.savefig(plot_dir / "top_skin_snv_abs_raw_score.pdf")
        fig.savefig(plot_dir / "top_skin_snv_abs_raw_score.png", dpi=300)
        plt.close(fig)

    if "quantile_score" in df.columns:
        fig, ax = plt.subplots(figsize=(3.4, 3.2))
        sns.scatterplot(
            data=df,
            x="raw_score",
            y="quantile_score",
            hue="biosample_name" if "biosample_name" in df.columns and df["biosample_name"].nunique() <= 8 else None,
            s=18,
            linewidth=0,
            alpha=0.8,
            palette="Set2",
            ax=ax,
        )
        ax.axvline(0, color="#777777", lw=0.7, ls="--")
        ax.axhline(0, color="#777777", lw=0.7, ls="--")
        ax.set_xlabel("Raw score")
        ax.set_ylabel("Quantile score")
        if ax.get_legend() is not None:
            ax.legend(frameon=False, fontsize=7, title=None)
        sns.despine(ax=ax)
        fig.tight_layout()
        fig.savefig(plot_dir / "skin_raw_vs_quantile_score.pdf")
        fig.savefig(plot_dir / "skin_raw_vs_quantile_score.png", dpi=300)
        plt.close(fig)


def main() -> None:
    args = parse_args()
    configure_logging(args.log_level)

    snvs = read_snvs(args.snv_list)
    LOGGER.info("Loaded %d SNVs from %s", len(snvs), args.snv_list)

    if args.mock_client:
        records = score_snvs_with_mock(snvs)
    else:
        api_key = args.api_key or os.environ.get("ALPHAGENOME_API_KEY")
        if not api_key:
            raise RuntimeError("Missing AlphaGenome API key. Pass --api-key or set ALPHAGENOME_API_KEY.")
        records = score_snvs_with_alphagenome(
            snvs=snvs,
            api_key=api_key,
            sequence_length=args.sequence_length,
            organism_name=args.organism,
            scorer_names=args.scorers,
            api_timeout=args.api_timeout,
            batch_size=args.batch_size,
            max_workers=args.max_workers,
            progress_bar=not args.no_progress,
            failed_snvs_path=args.failed_snvs
            or args.output_skin.with_name(f"{args.output_skin.stem}_failed_snvs.csv"),
        )

    skin_terms = tuple(args.skin_terms or DEFAULT_SKIN_TERMS)
    skin_records = [record for record in records if row_matches_skin(record, skin_terms)]
    write_records(args.output_all, records)
    write_records(args.output_skin, skin_records, fieldnames=list(records[0].keys()) if records else None)
    if args.plot_dir is not None:
        plot_summary(skin_records, args.plot_dir, args.top_n_plot)

    LOGGER.info("Wrote %d total score rows to %s", len(records), args.output_all)
    LOGGER.info("Wrote %d skin-related score rows to %s", len(skin_records), args.output_skin)
    if not skin_records:
        LOGGER.warning("No rows matched skin terms: %s", ", ".join(skin_terms))


if __name__ == "__main__":
    main()
