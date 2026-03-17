#!/usr/bin/env python
"""Rank cross-species H2A expression differences across shared genes and tissues."""

from __future__ import annotations

import argparse
import itertools
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import pandas as pd

from gene_compare_common import (
    DEFAULT_DETAIL_INDEX,
    DEFAULT_HEATMAP_DIR,
    DEFAULT_PROCESSED_DIR,
    DEFAULT_RANKING_OUT_DIR,
    DEFAULT_SHARED_INDEX,
    GENERIC_TISSUES,
    build_gene_long_dataframe,
    build_shared_gene_summary,
    collect_gene_rows,
    load_detail_index,
)


DEFAULT_MIN_GENE_SPECIES_COUNT = 4
MIN_CONSERVED_TISSUES = 10
ALL_SHARED_SHORTLIST_LIMIT = 10
VARIANT_SHORTLIST_LIMIT = 10
CONSERVED_SHORTLIST_LIMIT = 3


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Rank strong cross-species H2A expression differences across shared "
            "canonical gene names and produce manuscript-ready shortlists."
        )
    )
    parser.add_argument(
        "--shared-index",
        default=str(DEFAULT_SHARED_INDEX),
        help="Shared gene summary CSV with species_count per canonical gene name.",
    )
    parser.add_argument(
        "--detail-index",
        default=str(DEFAULT_DETAIL_INDEX),
        help="Detail index CSV built by summarize_shared_h2a_gene_names_across_species.py.",
    )
    parser.add_argument(
        "--heatmap-dir",
        default=str(DEFAULT_HEATMAP_DIR),
        help="Directory with per-species heatmap subdirectories.",
    )
    parser.add_argument(
        "--processed-dir",
        default=str(DEFAULT_PROCESSED_DIR),
        help="Directory with per-species processed *_present_gold.tsv files.",
    )
    parser.add_argument(
        "--out-dir",
        default=str(DEFAULT_RANKING_OUT_DIR),
        help="Directory for ranking tables and manuscript shortlist.",
    )
    parser.add_argument(
        "--min-species-per-tissue",
        default=4,
        type=int,
        help="Keep only tissues represented in at least this many species for candidate selection.",
    )
    parser.add_argument(
        "--candidate-quantile",
        default=0.90,
        type=float,
        help="Global range quantile threshold for conservative candidates.",
    )
    parser.add_argument(
        "--high-confidence-quantile",
        default=0.95,
        type=float,
        help="Global range quantile threshold for high-confidence candidates.",
    )
    parser.add_argument(
        "--scope",
        default="two-tier",
        choices=["two-tier", "all-shared", "variants-only"],
        help="Shortlist scope for manuscript-oriented output.",
    )
    return parser.parse_args()


def load_shared_index(shared_index_path: Path, detail_df: pd.DataFrame, detail_index_path: Path) -> pd.DataFrame:
    if shared_index_path.exists():
        shared_df = pd.read_csv(shared_index_path, dtype=str)
    else:
        shared_df = build_shared_gene_summary(detail_df, detail_index_path)
    if shared_df.empty:
        return shared_df

    for col in ["gene_name", "canonical_gene_name"]:
        if col not in shared_df.columns:
            shared_df[col] = ""
        shared_df[col] = shared_df[col].fillna("").astype(str).str.strip()
    if "species_count" not in shared_df.columns:
        shared_df["species_count"] = 0
    shared_df["species_count"] = pd.to_numeric(shared_df["species_count"], errors="coerce").fillna(0)
    return shared_df


def choose_target_genes(shared_df: pd.DataFrame, min_gene_species_count: int) -> List[str]:
    filtered = shared_df[shared_df["species_count"] >= int(min_gene_species_count)].copy()
    if filtered.empty:
        return []
    filtered = filtered.sort_values(
        by=["species_count", "canonical_gene_name"],
        ascending=[False, True],
    )
    return filtered["canonical_gene_name"].tolist()


def summarize_gene_class(gene_rows: pd.DataFrame) -> str:
    classes = sorted(
        {
            value
            for value in gene_rows["class"].fillna("").astype(str).str.strip().tolist()
            if value
        }
    )
    if not classes:
        return ""
    return ",".join(classes)


def build_species_level_long(detail_df: pd.DataFrame, target_genes: Iterable[str]) -> pd.DataFrame:
    expr_cache: Dict[str, pd.DataFrame] = {}
    frames: List[pd.DataFrame] = []
    for gene_name in target_genes:
        gene_rows, _ = collect_gene_rows(gene_name, detail_df)
        gene_class = summarize_gene_class(gene_rows)
        gene_long_df = build_gene_long_dataframe(
            gene_name,
            gene_rows,
            expr_cache=expr_cache,
            aggregate="mean",
        )
        species_level = (
            gene_long_df[
                [
                    "gene_name",
                    "species_dir",
                    "species_name",
                    "tissue",
                    "agg_expression_score",
                ]
            ]
            .drop_duplicates(subset=["gene_name", "species_dir", "tissue"], keep="first")
            .rename(columns={"agg_expression_score": "expression_score"})
        )
        species_level["gene_class"] = gene_class
        species_level["gene_species_count"] = int(gene_rows["species_dir"].nunique())
        frames.append(species_level)

    if not frames:
        return pd.DataFrame(
            columns=[
                "gene_name",
                "species_dir",
                "species_name",
                "tissue",
                "expression_score",
                "gene_class",
                "gene_species_count",
            ]
        )

    long_df = pd.concat(frames, ignore_index=True)
    long_df["expression_score"] = pd.to_numeric(long_df["expression_score"], errors="coerce")
    long_df = long_df[long_df["expression_score"].notna()].copy()
    long_df["gene_species_count"] = pd.to_numeric(
        long_df["gene_species_count"], errors="coerce"
    ).fillna(0).astype(int)
    long_df = long_df[~long_df["tissue"].isin(GENERIC_TISSUES)].copy()
    return long_df


def build_gene_tissue_summary(species_level_long_df: pd.DataFrame) -> pd.DataFrame:
    comparable_df = species_level_long_df.copy()
    if comparable_df.empty:
        return pd.DataFrame()

    summary_df = (
        comparable_df.groupby(["gene_name", "tissue"], as_index=False)
        .agg(
            gene_class=("gene_class", "first"),
            gene_species_count=("gene_species_count", "first"),
            species_n=("species_dir", "nunique"),
            mean_score=("expression_score", "mean"),
            median_score=("expression_score", "median"),
            std_score=("expression_score", "std"),
            min_score=("expression_score", "min"),
            max_score=("expression_score", "max"),
        )
    )
    summary_df = summary_df[summary_df["species_n"] >= 2].copy()
    summary_df["range"] = summary_df["max_score"] - summary_df["min_score"]

    idx_max = comparable_df.groupby(["gene_name", "tissue"])["expression_score"].idxmax()
    idx_min = comparable_df.groupby(["gene_name", "tissue"])["expression_score"].idxmin()

    max_df = comparable_df.loc[
        idx_max,
        ["gene_name", "tissue", "species_dir", "species_name", "expression_score"],
    ].rename(
        columns={
            "species_dir": "max_species",
            "species_name": "max_species_name",
            "expression_score": "_max_score_check",
        }
    )
    min_df = comparable_df.loc[
        idx_min,
        ["gene_name", "tissue", "species_dir", "species_name", "expression_score"],
    ].rename(
        columns={
            "species_dir": "min_species",
            "species_name": "min_species_name",
            "expression_score": "_min_score_check",
        }
    )

    summary_df = summary_df.merge(max_df, on=["gene_name", "tissue"], how="left")
    summary_df = summary_df.merge(min_df, on=["gene_name", "tissue"], how="left")
    summary_df["std_score"] = summary_df["std_score"].fillna(0.0)
    return summary_df.sort_values(
        by=["range", "species_n", "gene_name", "tissue"],
        ascending=[False, False, True, True],
    ).reset_index(drop=True)


def assign_candidate_flags(
    summary_df: pd.DataFrame,
    min_species_per_tissue: int,
    candidate_quantile: float,
    high_confidence_quantile: float,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, float, float]:
    conservative_base = summary_df[summary_df["species_n"] >= int(min_species_per_tissue)].copy()
    if conservative_base.empty:
        raise RuntimeError("No gene/tissue combinations met the conservative species filter.")

    candidate_threshold = float(conservative_base["range"].quantile(candidate_quantile))
    high_conf_threshold = float(conservative_base["range"].quantile(high_confidence_quantile))

    summary_df = summary_df.copy()
    summary_df["passes_conservative_species_filter"] = summary_df["species_n"] >= int(
        min_species_per_tissue
    )
    summary_df["is_conservative_candidate"] = summary_df["passes_conservative_species_filter"] & (
        summary_df["range"] >= candidate_threshold
    )
    summary_df["is_high_confidence_candidate"] = summary_df[
        "passes_conservative_species_filter"
    ] & (summary_df["range"] >= high_conf_threshold)
    summary_df["candidate_level"] = "background"
    summary_df.loc[summary_df["is_conservative_candidate"], "candidate_level"] = "conservative"
    summary_df.loc[
        summary_df["is_high_confidence_candidate"], "candidate_level"
    ] = "high_confidence"

    conservative_df = summary_df[summary_df["is_conservative_candidate"]].copy()
    high_conf_df = summary_df[summary_df["is_high_confidence_candidate"]].copy()
    conservative_df = conservative_df.sort_values(
        by=["range", "species_n", "gene_name", "tissue"],
        ascending=[False, False, True, True],
    ).reset_index(drop=True)
    high_conf_df = high_conf_df.sort_values(
        by=["range", "species_n", "gene_name", "tissue"],
        ascending=[False, False, True, True],
    ).reset_index(drop=True)
    return summary_df, conservative_df, high_conf_df, candidate_threshold, high_conf_threshold


def build_pairwise_contrasts(
    species_level_long_df: pd.DataFrame,
    candidate_df: pd.DataFrame,
) -> pd.DataFrame:
    empty_cols = [
        "gene_name",
        "gene_class",
        "tissue",
        "candidate_level",
        "species_n",
        "combo_range",
        "max_species",
        "max_species_name",
        "min_species",
        "min_species_name",
        "species_high",
        "species_high_name",
        "high_score",
        "species_low",
        "species_low_name",
        "low_score",
        "abs_diff",
    ]
    if candidate_df.empty:
        return pd.DataFrame(columns=empty_cols)

    candidate_lookup = {
        (row["gene_name"], row["tissue"]): row for row in candidate_df.to_dict(orient="records")
    }
    rows: List[dict] = []
    for (gene_name, tissue), combo_df in species_level_long_df.groupby(["gene_name", "tissue"]):
        candidate_row = candidate_lookup.get((gene_name, tissue))
        if candidate_row is None:
            continue

        records = combo_df[
            ["species_dir", "species_name", "expression_score"]
        ].to_dict(orient="records")
        for left, right in itertools.combinations(records, 2):
            if left["expression_score"] >= right["expression_score"]:
                high = left
                low = right
            else:
                high = right
                low = left
            rows.append(
                {
                    "gene_name": gene_name,
                    "gene_class": candidate_row["gene_class"],
                    "tissue": tissue,
                    "candidate_level": candidate_row["candidate_level"],
                    "species_n": int(candidate_row["species_n"]),
                    "combo_range": float(candidate_row["range"]),
                    "max_species": candidate_row["max_species"],
                    "max_species_name": candidate_row["max_species_name"],
                    "min_species": candidate_row["min_species"],
                    "min_species_name": candidate_row["min_species_name"],
                    "species_high": high["species_dir"],
                    "species_high_name": high["species_name"],
                    "high_score": float(high["expression_score"]),
                    "species_low": low["species_dir"],
                    "species_low_name": low["species_name"],
                    "low_score": float(low["expression_score"]),
                    "abs_diff": float(high["expression_score"] - low["expression_score"]),
                }
            )

    if not rows:
        return pd.DataFrame(columns=empty_cols)

    return pd.DataFrame(rows).sort_values(
        by=["abs_diff", "combo_range", "gene_name", "tissue"],
        ascending=[False, False, True, True],
    ).reset_index(drop=True)


def build_gene_variability_summary(summary_df: pd.DataFrame, min_species_per_tissue: int) -> pd.DataFrame:
    conservative_summary = summary_df[summary_df["species_n"] >= int(min_species_per_tissue)].copy()
    if conservative_summary.empty:
        return pd.DataFrame()

    gene_df = (
        conservative_summary.groupby("gene_name", as_index=False)
        .agg(
            gene_class=("gene_class", "first"),
            gene_species_count=("gene_species_count", "first"),
            tissues_compared=("tissue", "count"),
            max_range=("range", "max"),
            median_range=("range", "median"),
            mean_range=("range", "mean"),
        )
    )
    return gene_df.sort_values(
        by=["median_range", "max_range", "gene_name"],
        ascending=[True, True, True],
    ).reset_index(drop=True)


def take_shortlist_rows(
    high_conf_df: pd.DataFrame,
    conservative_df: pd.DataFrame,
    limit: int,
    gene_class: Optional[str] = None,
) -> pd.DataFrame:
    frames: List[pd.DataFrame] = []
    if gene_class is None:
        high_subset = high_conf_df
        cons_subset = conservative_df
    else:
        high_subset = high_conf_df[high_conf_df["gene_class"].eq(gene_class)].copy()
        cons_subset = conservative_df[conservative_df["gene_class"].eq(gene_class)].copy()
    frames.extend([high_subset, cons_subset])

    picked_keys = set()
    rows: List[dict] = []
    for frame in frames:
        for record in frame.to_dict(orient="records"):
            key = (record["gene_name"], record["tissue"])
            if key in picked_keys:
                continue
            picked_keys.add(key)
            rows.append(record)
            if len(rows) >= limit:
                return pd.DataFrame(rows)
    return pd.DataFrame(rows)


def take_ranked_summary_rows(
    summary_df: pd.DataFrame,
    limit: int,
    min_species_per_tissue: int,
    gene_class: Optional[str] = None,
) -> pd.DataFrame:
    ranked_df = summary_df[summary_df["species_n"] >= int(min_species_per_tissue)].copy()
    if gene_class is not None:
        ranked_df = ranked_df[ranked_df["gene_class"].eq(gene_class)].copy()
    ranked_df = ranked_df.sort_values(
        by=["range", "species_n", "gene_name", "tissue"],
        ascending=[False, False, True, True],
    )
    return ranked_df.head(limit).reset_index(drop=True)


def format_candidate_sentence(row: pd.Series) -> str:
    return (
        f'In {row["tissue"]}, {row["gene_name"]} expression was higher in '
        f'{row["max_species_name"]} ({row["max_score"]:.2f}) than in '
        f'{row["min_species_name"]} ({row["min_score"]:.2f}), yielding a '
        f'cross-species range of {row["range"]:.2f} across {int(row["species_n"])} species.'
    )


def write_candidate_section(lines: List[str], title: str, section_df: pd.DataFrame) -> None:
    lines.append(f"## {title}")
    if section_df.empty:
        lines.append("")
        lines.append("- No candidates met the current filtering rules.")
        lines.append("")
        return

    lines.append("")
    for row in section_df.to_dict(orient="records"):
        label = (
            f'- `{row["gene_name"]}` in `{row["tissue"]}`: '
            f'{row["min_species"]} {row["min_score"]:.2f} vs '
            f'{row["max_species"]} {row["max_score"]:.2f} '
            f'(range {row["range"]:.2f}, {int(row["species_n"])} species, '
            f'{row["candidate_level"]}).'
        )
        lines.append(label)
        lines.append(f'  Template: "{format_candidate_sentence(pd.Series(row))}"')
    lines.append("")


def write_conserved_section(lines: List[str], gene_variability_df: pd.DataFrame) -> None:
    lines.append("## Conserved genes")
    lines.append("")
    conserved_df = gene_variability_df[gene_variability_df["tissues_compared"] >= MIN_CONSERVED_TISSUES]
    conserved_df = conserved_df.head(CONSERVED_SHORTLIST_LIMIT)
    if conserved_df.empty:
        lines.append("- No genes met the minimum tissue-count rule for the conserved block.")
        lines.append("")
        return

    for row in conserved_df.to_dict(orient="records"):
        sentence = (
            f'{row["gene_name"]} showed comparatively stable expression across species '
            f'(median range {row["median_range"]:.2f}, max range {row["max_range"]:.2f}, '
            f'{int(row["tissues_compared"])} tissues).'
        )
        lines.append(
            f'- `{row["gene_name"]}` ({row["gene_class"] or "unclassified"}): '
            f'median range {row["median_range"]:.2f}, max range {row["max_range"]:.2f}, '
            f'{int(row["tissues_compared"])} tissues.'
        )
        lines.append(f'  Template: "{sentence}"')
    lines.append("")


def build_manuscript_shortlist(
    *,
    out_path: Path,
    scope: str,
    summary_df: pd.DataFrame,
    conservative_df: pd.DataFrame,
    high_conf_df: pd.DataFrame,
    gene_variability_df: pd.DataFrame,
    candidate_threshold: float,
    high_conf_threshold: float,
    min_species_per_tissue: int,
    target_gene_count: int,
) -> None:
    lines: List[str] = []
    lines.append("# Cross-species H2A manuscript shortlist")
    lines.append("")
    lines.append(
        f"- Genes analyzed: {target_gene_count} shared genes with species_count >= {DEFAULT_MIN_GENE_SPECIES_COUNT}."
    )
    lines.append(
        f"- Tissue filter for candidate selection: species_n >= {min_species_per_tissue}."
    )
    lines.append(
        f"- Generic tissues removed: {', '.join(sorted(GENERIC_TISSUES))}."
    )
    lines.append(
        f"- Conservative threshold: range >= global p90 = {candidate_threshold:.2f}."
    )
    lines.append(
        f"- High-confidence threshold: range >= global p95 = {high_conf_threshold:.2f}."
    )
    lines.append(
        f"- Candidate counts: {len(conservative_df)} conservative, {len(high_conf_df)} high-confidence."
    )
    lines.append("")

    if scope in {"two-tier", "all-shared"}:
        all_shared_df = take_shortlist_rows(
            high_conf_df,
            conservative_df,
            limit=ALL_SHARED_SHORTLIST_LIMIT,
        )
        write_candidate_section(lines, "All shared genes", all_shared_df)

    if scope in {"two-tier", "variants-only"}:
        variant_df = take_ranked_summary_rows(
            summary_df,
            limit=VARIANT_SHORTLIST_LIMIT,
            min_species_per_tissue=min_species_per_tissue,
            gene_class="variant",
        )
        write_candidate_section(lines, "Variant genes", variant_df)

    write_conserved_section(lines, gene_variability_df)
    out_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()

    if not 0 < args.candidate_quantile < 1:
        raise ValueError("--candidate-quantile must be between 0 and 1")
    if not 0 < args.high_confidence_quantile < 1:
        raise ValueError("--high-confidence-quantile must be between 0 and 1")
    if args.candidate_quantile > args.high_confidence_quantile:
        raise ValueError("--candidate-quantile must be <= --high-confidence-quantile")

    shared_index_path = Path(args.shared_index)
    detail_index_path = Path(args.detail_index)
    heatmap_dir = Path(args.heatmap_dir)
    processed_dir = Path(args.processed_dir)
    out_dir = Path(args.out_dir)
    tables_dir = out_dir / "tables"
    reports_dir = out_dir / "reports"

    if not heatmap_dir.exists():
        raise FileNotFoundError(heatmap_dir)
    if not processed_dir.exists():
        raise FileNotFoundError(processed_dir)

    tables_dir.mkdir(parents=True, exist_ok=True)
    reports_dir.mkdir(parents=True, exist_ok=True)

    detail_df = load_detail_index(
        detail_index_path,
        heatmap_dir=heatmap_dir,
        processed_dir=processed_dir,
    )
    shared_df = load_shared_index(shared_index_path, detail_df, detail_index_path)
    target_genes = choose_target_genes(shared_df, DEFAULT_MIN_GENE_SPECIES_COUNT)
    if not target_genes:
        raise RuntimeError("No shared genes met the conservative gene-count filter.")

    species_level_long_df = build_species_level_long(detail_df, target_genes)
    summary_df = build_gene_tissue_summary(species_level_long_df)
    summary_df, conservative_df, high_conf_df, candidate_threshold, high_conf_threshold = (
        assign_candidate_flags(
            summary_df,
            min_species_per_tissue=args.min_species_per_tissue,
            candidate_quantile=args.candidate_quantile,
            high_confidence_quantile=args.high_confidence_quantile,
        )
    )

    candidate_union_df = pd.concat([conservative_df, high_conf_df], ignore_index=True)
    candidate_union_df = candidate_union_df.drop_duplicates(
        subset=["gene_name", "tissue"], keep="last"
    ).reset_index(drop=True)
    pairwise_df = build_pairwise_contrasts(species_level_long_df, candidate_union_df)
    gene_variability_df = build_gene_variability_summary(
        summary_df,
        min_species_per_tissue=args.min_species_per_tissue,
    )

    out_summary_csv = tables_dir / "gene_tissue_summary.csv"
    out_conservative_csv = tables_dir / "conservative_candidates.csv"
    out_high_conf_csv = tables_dir / "high_confidence_candidates.csv"
    out_pairwise_csv = tables_dir / "pairwise_contrasts.csv"
    out_shortlist_md = reports_dir / "manuscript_shortlist.md"

    summary_df.to_csv(out_summary_csv, index=False, encoding="utf-8")
    conservative_df.to_csv(out_conservative_csv, index=False, encoding="utf-8")
    high_conf_df.to_csv(out_high_conf_csv, index=False, encoding="utf-8")
    pairwise_df.to_csv(out_pairwise_csv, index=False, encoding="utf-8")
    build_manuscript_shortlist(
        out_path=out_shortlist_md,
        scope=args.scope,
        summary_df=summary_df,
        conservative_df=conservative_df,
        high_conf_df=high_conf_df,
        gene_variability_df=gene_variability_df,
        candidate_threshold=candidate_threshold,
        high_conf_threshold=high_conf_threshold,
        min_species_per_tissue=args.min_species_per_tissue,
        target_gene_count=len(target_genes),
    )

    print(f"Saved summary to {out_summary_csv}")
    print(f"Saved conservative candidates to {out_conservative_csv}")
    print(f"Saved high-confidence candidates to {out_high_conf_csv}")
    print(f"Saved pairwise contrasts to {out_pairwise_csv}")
    print(f"Saved manuscript shortlist to {out_shortlist_md}")


if __name__ == "__main__":
    main()
