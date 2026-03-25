#!/usr/bin/env python
"""Build a cross-species heatmap for one canonical gene name."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from gene_compare_common import (
    DEFAULT_DETAIL_INDEX,
    DEFAULT_GENE_COMPARE_DATA_ROOT,
    DEFAULT_GENE_COMPARE_FIG_ROOT,
    DEFAULT_HEATMAP_DIR,
    DEFAULT_PROCESSED_DIR,
    build_gene_long_dataframe,
    build_matrix,
    collect_gene_rows,
    ensure_mean_aggregate,
    infer_species_name,
    load_detail_index,
    path_str,
    safe_slug,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a cross-species heatmap for one canonical gene name using "
            "the normalized H2A cell files and the canonical variant maps."
        )
    )
    parser.add_argument(
        "--gene-name",
        required=True,
        help="Canonical gene name to compare across species, e.g. H2AJ.",
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
        "--detail-index",
        default=str(DEFAULT_DETAIL_INDEX),
        help="Detail index CSV built by summarize_shared_h2a_gene_names_across_species.py.",
    )
    parser.add_argument(
        "--tissue-mode",
        default="union",
        choices=["union"],
        help="Tissue set mode for the heatmap.",
    )
    parser.add_argument(
        "--include-absent-species",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Include all known species as empty columns when the gene is absent.",
    )
    parser.add_argument(
        "--aggregate",
        default="mean",
        choices=["mean"],
        help="Aggregation for multiple rows per tissue/species.",
    )
    return parser.parse_args()


def plot_heatmap(matrix_df: pd.DataFrame, gene_name: str, out_png: Path, out_svg: Path) -> None:
    if matrix_df.empty:
        raise RuntimeError(f"Heatmap matrix is empty for gene '{gene_name}'.")

    plot_df = np.log10(matrix_df + 1)
    fig_w = max(9, plot_df.shape[1] * 1.2)
    fig_h = max(12, plot_df.shape[0] * 0.22)

    sns.set(style="whitegrid")
    plt.figure(figsize=(fig_w, fig_h))
    ax = sns.heatmap(
        plot_df,
        cmap="viridis",
        mask=plot_df.isna(),
        linewidths=0.1,
        cbar_kws={"label": "log10(Expression score + 1)"},
        xticklabels=1,
        yticklabels=1,
    )
    ax.set_title(f"{gene_name} Expression Across Species (normalized cells)")
    ax.set_xlabel("Species")
    ax.set_ylabel("Anatomical entity name")
    plt.xticks(rotation=45, ha="right", fontsize=9)
    plt.yticks(rotation=0, fontsize=7)
    plt.tight_layout()

    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_svg, bbox_inches="tight")
    plt.close()


def write_metadata(
    *,
    gene_name: str,
    gene_rows: pd.DataFrame,
    long_df: pd.DataFrame,
    matrix_df: pd.DataFrame,
    species_order: List[str],
    tissue_mode: str,
    aggregate: str,
    include_absent_species: bool,
    detail_index_path: Path,
    metadata_path: Path,
) -> None:
    species_entries: List[dict] = []
    for species_dir in species_order:
        species_df = gene_rows[gene_rows["species_dir"].eq(species_dir)].copy()
        if species_df.empty:
            species_entries.append(
                {
                    "species_dir": species_dir,
                    "species_name": infer_species_name(species_dir),
                    "ensembl_gene_ids": [],
                    "present_gold_path": "",
                    "map_path": "",
                    "row_count": 0,
                    "tissue_count": 0,
                }
            )
            continue

        expr_path = species_df["present_gold_path"].iloc[0]
        map_path = species_df["map_path"].iloc[0]
        species_entries.append(
            {
                "species_dir": species_dir,
                "species_name": species_df["species_name"].iloc[0],
                "ensembl_gene_ids": sorted(
                    species_df["ensembl_gene_id"].dropna().astype(str).unique().tolist()
                ),
                "present_gold_path": expr_path,
                "map_path": map_path,
                "row_count": int(long_df[long_df["species_dir"].eq(species_dir)].shape[0]),
                "tissue_count": int(
                    long_df[long_df["species_dir"].eq(species_dir)]["tissue"].nunique()
                ),
            }
        )

    metadata = {
        "gene_name": gene_name,
        "tissue_mode": tissue_mode,
        "aggregate": aggregate,
        "include_absent_species": include_absent_species,
        "detail_index_path": path_str(detail_index_path),
        "species_count": len(species_order),
        "tissue_count": int(matrix_df.shape[0]),
        "species_order": species_order,
        "species": species_entries,
    }

    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")


def run_gene_compare_heatmap(
    *,
    gene_name: str,
    heatmap_dir: Path,
    processed_dir: Path,
    out_fig_root: Path,
    out_data_root: Path,
    detail_index_path: Path,
    tissue_mode: str = "union",
    include_absent_species: bool = False,
    aggregate: str = "mean",
    expr_cache: Optional[Dict[str, pd.DataFrame]] = None,
    detail_df: Optional[pd.DataFrame] = None,
) -> Dict[str, Path]:
    ensure_mean_aggregate(aggregate)
    if tissue_mode != "union":
        raise ValueError(f"Unsupported tissue_mode: {tissue_mode}")
    if not gene_name.strip():
        raise RuntimeError("Gene name must not be empty.")
    if not heatmap_dir.exists():
        raise FileNotFoundError(heatmap_dir)
    if not processed_dir.exists():
        raise FileNotFoundError(processed_dir)

    if detail_df is None:
        detail_df = load_detail_index(
            detail_index_path,
            heatmap_dir,
            processed_dir,
            expr_cache=expr_cache,
        )

    gene_rows, species_order = collect_gene_rows(
        gene_name,
        detail_df,
        include_absent_species=include_absent_species,
        heatmap_dir=heatmap_dir,
    )
    canonical_gene_name = gene_rows["canonical_gene_name"].iloc[0]
    gene_slug = safe_slug(canonical_gene_name)

    long_df = build_gene_long_dataframe(
        canonical_gene_name,
        gene_rows,
        expr_cache=expr_cache,
        aggregate=aggregate,
    )
    matrix_df = build_matrix(long_df, species_order, aggregate=aggregate)

    fig_dir = out_fig_root / gene_slug
    data_dir = out_data_root / gene_slug
    out_png = fig_dir / f"{gene_slug}_gene_compare_heatmap.png"
    out_svg = fig_dir / f"{gene_slug}_gene_compare_heatmap.svg"
    out_long_csv = data_dir / f"{gene_slug}_gene_compare_long.csv"
    out_matrix_csv = data_dir / f"{gene_slug}_gene_compare_matrix.csv"
    out_metadata_json = data_dir / f"{gene_slug}_gene_compare_metadata.json"

    fig_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)

    long_df.to_csv(out_long_csv, index=False, encoding="utf-8")
    matrix_df.to_csv(out_matrix_csv, encoding="utf-8")
    plot_heatmap(matrix_df, canonical_gene_name, out_png, out_svg)
    write_metadata(
        gene_name=canonical_gene_name,
        gene_rows=gene_rows,
        long_df=long_df,
        matrix_df=matrix_df,
        species_order=species_order,
        tissue_mode=tissue_mode,
        aggregate=aggregate,
        include_absent_species=include_absent_species,
        detail_index_path=detail_index_path,
        metadata_path=out_metadata_json,
    )

    return {
        "out_long_csv": out_long_csv,
        "out_matrix_csv": out_matrix_csv,
        "out_png": out_png,
        "out_svg": out_svg,
        "out_metadata_json": out_metadata_json,
    }


def main() -> None:
    args = parse_args()
    outputs = run_gene_compare_heatmap(
        gene_name=args.gene_name.strip(),
        heatmap_dir=Path(args.heatmap_dir),
        processed_dir=Path(args.processed_dir),
        out_fig_root=Path(args.out_fig_root),
        out_data_root=Path(args.out_data_root),
        detail_index_path=Path(args.detail_index),
        tissue_mode=args.tissue_mode,
        include_absent_species=args.include_absent_species,
        aggregate=args.aggregate,
    )

    print(f"Saved long table to {outputs['out_long_csv']}")
    print(f"Saved matrix table to {outputs['out_matrix_csv']}")
    print(f"Saved heatmap to {outputs['out_png']}")
    print(f"Saved heatmap to {outputs['out_svg']}")
    print(f"Saved metadata to {outputs['out_metadata_json']}")


if __name__ == "__main__":
    main()
