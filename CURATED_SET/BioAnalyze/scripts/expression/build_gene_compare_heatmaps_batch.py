#!/usr/bin/env python
"""Build cross-species gene_compare heatmaps for multiple shared H2A genes."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List

import pandas as pd

from build_gene_compare_heatmap import run_gene_compare_heatmap
from gene_compare_common import (
    DEFAULT_DETAIL_INDEX,
    DEFAULT_GENE_COMPARE_DATA_ROOT,
    DEFAULT_GENE_COMPARE_FIG_ROOT,
    DEFAULT_HEATMAP_DIR,
    DEFAULT_PROCESSED_DIR,
    DEFAULT_SHARED_INDEX,
    build_shared_gene_summary,
    load_detail_index,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build gene_compare heatmaps for shared H2A genes using the same "
            "logic as the single-gene CLI."
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
        "--out-fig-root",
        default=str(DEFAULT_GENE_COMPARE_FIG_ROOT),
        help="Root directory for gene_compare heatmap figures.",
    )
    parser.add_argument(
        "--out-data-root",
        default=str(DEFAULT_GENE_COMPARE_DATA_ROOT),
        help="Root directory for gene_compare output tables.",
    )
    parser.add_argument(
        "--min-species-count",
        default=4,
        type=int,
        help="Only build heatmaps for genes present in at least this many species.",
    )
    parser.add_argument(
        "--genes",
        nargs="*",
        default=[],
        help="Optional gene names to run explicitly. Supports spaces or comma-separated values.",
    )
    parser.add_argument(
        "--tissue-mode",
        default="union",
        choices=["union"],
        help="Tissue set mode for the heatmap.",
    )
    parser.add_argument(
        "--include-absent-species",
        action="store_true",
        help="Include all known species as empty columns when a gene is absent.",
    )
    parser.add_argument(
        "--aggregate",
        default="mean",
        choices=["mean"],
        help="Aggregation for multiple rows per tissue/species.",
    )
    return parser.parse_args()


def parse_gene_values(raw_values: List[str]) -> List[str]:
    genes: List[str] = []
    for value in raw_values:
        for chunk in value.split(","):
            gene = chunk.strip()
            if gene:
                genes.append(gene)
    seen = set()
    ordered: List[str] = []
    for gene in genes:
        norm = gene.casefold()
        if norm in seen:
            continue
        seen.add(norm)
        ordered.append(gene)
    return ordered


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


def choose_genes(
    shared_df: pd.DataFrame,
    explicit_genes: List[str],
    min_species_count: int,
) -> List[str]:
    if explicit_genes:
        gene_lookup = {
            row["canonical_gene_name"].casefold(): row["canonical_gene_name"]
            for row in shared_df.to_dict(orient="records")
            if row.get("canonical_gene_name", "")
        }
        gene_lookup.update(
            {
                row["gene_name"].casefold(): row["canonical_gene_name"] or row["gene_name"]
                for row in shared_df.to_dict(orient="records")
                if row.get("gene_name", "")
            }
        )
        chosen: List[str] = []
        missing: List[str] = []
        for gene in explicit_genes:
            canonical = gene_lookup.get(gene.casefold())
            if canonical is None:
                missing.append(gene)
                continue
            if canonical not in chosen:
                chosen.append(canonical)
        if missing:
            raise RuntimeError(f"Gene(s) not found in shared index: {', '.join(missing)}")
        return chosen

    filtered = shared_df[shared_df["species_count"] >= int(min_species_count)].copy()
    if filtered.empty:
        return []
    filtered = filtered.sort_values(
        by=["species_count", "canonical_gene_name"],
        ascending=[False, True],
    )
    return filtered["canonical_gene_name"].tolist()


def main() -> None:
    args = parse_args()

    shared_index_path = Path(args.shared_index)
    detail_index_path = Path(args.detail_index)
    heatmap_dir = Path(args.heatmap_dir)
    processed_dir = Path(args.processed_dir)
    out_fig_root = Path(args.out_fig_root)
    out_data_root = Path(args.out_data_root)

    if not heatmap_dir.exists():
        raise FileNotFoundError(heatmap_dir)
    if not processed_dir.exists():
        raise FileNotFoundError(processed_dir)

    expr_cache = {}
    detail_df = load_detail_index(
        detail_index_path,
        heatmap_dir=heatmap_dir,
        processed_dir=processed_dir,
        expr_cache=expr_cache,
    )
    shared_df = load_shared_index(shared_index_path, detail_df, detail_index_path)
    target_genes = choose_genes(
        shared_df,
        explicit_genes=parse_gene_values(args.genes),
        min_species_count=args.min_species_count,
    )

    if not target_genes:
        raise RuntimeError("No genes matched the batch selection criteria.")

    for gene_name in target_genes:
        outputs = run_gene_compare_heatmap(
            gene_name=gene_name,
            heatmap_dir=heatmap_dir,
            processed_dir=processed_dir,
            out_fig_root=out_fig_root,
            out_data_root=out_data_root,
            detail_index_path=detail_index_path,
            tissue_mode=args.tissue_mode,
            include_absent_species=args.include_absent_species,
            aggregate=args.aggregate,
            expr_cache=expr_cache,
            detail_df=detail_df,
        )
        print(f"[ok] {gene_name} -> {outputs['out_png']}")

    print(f"Built {len(target_genes)} gene_compare heatmap(s).")


if __name__ == "__main__":
    main()
