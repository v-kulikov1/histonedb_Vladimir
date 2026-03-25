#!/usr/bin/env python
"""Summarize shared H2A gene names across species and build a reusable index."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from gene_compare_common import (
    DEFAULT_HEATMAP_DIR,
    DEFAULT_INDEX_DIR,
    DEFAULT_PROCESSED_DIR,
    DEFAULT_SHARED_STATS_DIR,
    DETAIL_INDEX_FILENAME,
    build_detail_index,
    build_shared_gene_summary,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Count in how many species each H2A Gene name appears in the "
            "species-level normalized H2A TSV files used for heatmaps and build "
            "a reusable detail index for gene_compare workflows."
        )
    )
    parser.add_argument(
        "--heatmap-dir",
        default=str(DEFAULT_HEATMAP_DIR),
        help="Directory with per-species heatmap subdirectories.",
    )
    parser.add_argument(
        "--processed-dir",
        default=str(DEFAULT_PROCESSED_DIR),
        help="Directory with per-species normalized H2A TSV files.",
    )
    parser.add_argument(
        "--out-dir",
        default=str(DEFAULT_SHARED_STATS_DIR),
        help="Directory for summary table and plot.",
    )
    parser.add_argument(
        "--out-index-dir",
        default=str(DEFAULT_INDEX_DIR),
        help="Directory for the detailed reusable gene index.",
    )
    return parser.parse_args()


def plot_summary(summary_df: pd.DataFrame, out_png: Path, out_svg: Path) -> None:
    if summary_df.empty:
        raise RuntimeError("No gene names were shared across more than one species.")

    plot_df = summary_df.sort_values(
        by=["species_count", "gene_name"], ascending=[True, False]
    ).reset_index(drop=True)

    fig_h = max(6, len(plot_df) * 0.28)
    fig_w = max(10, min(18, plot_df["species_count"].max() * 1.2))

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    bars = ax.barh(plot_df["gene_name"], plot_df["species_count"], color="#3C6E71")

    ax.set_title("Shared H2A Gene Names Across Species")
    ax.set_xlabel("Number of species with gene in normalized H2A TSV")
    ax.set_ylabel("Gene name")
    ax.set_xlim(0, max(int(plot_df["species_count"].max()) + 1, 2))

    for bar, value in zip(bars, plot_df["species_count"]):
        ax.text(
            bar.get_width() + 0.05,
            bar.get_y() + bar.get_height() / 2,
            str(int(value)),
            va="center",
            fontsize=9,
        )

    plt.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()

    heatmap_dir = Path(args.heatmap_dir)
    processed_dir = Path(args.processed_dir)
    out_dir = Path(args.out_dir)
    out_index_dir = Path(args.out_index_dir)

    if not heatmap_dir.exists():
        raise FileNotFoundError(heatmap_dir)
    if not processed_dir.exists():
        raise FileNotFoundError(processed_dir)

    out_dir.mkdir(parents=True, exist_ok=True)
    out_index_dir.mkdir(parents=True, exist_ok=True)

    detail_index_path = out_index_dir / DETAIL_INDEX_FILENAME
    expr_cache = {}
    detail_df = build_detail_index(heatmap_dir, processed_dir, expr_cache=expr_cache)
    detail_df.to_csv(detail_index_path, index=False, encoding="utf-8")

    summary_df = build_shared_gene_summary(detail_df, detail_index_path)

    out_csv = out_dir / "shared_h2a_gene_names_across_species.csv"
    out_png = out_dir / "shared_h2a_gene_names_across_species.png"
    out_svg = out_dir / "shared_h2a_gene_names_across_species.svg"

    summary_df.to_csv(out_csv, index=False, encoding="utf-8")
    plot_summary(summary_df, out_png, out_svg)

    print(f"Saved detail index with {len(detail_df)} rows to {detail_index_path}")
    print(f"Saved {len(summary_df)} shared gene names to {out_csv}")
    print(f"Saved plot to {out_png}")
    print(f"Saved plot to {out_svg}")


if __name__ == "__main__":
    main()
