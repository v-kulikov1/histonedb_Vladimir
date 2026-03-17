#!/usr/bin/env python
"""Plot overview and focused panels for strong cross-species H2A candidates."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Tuple

import matplotlib.pyplot as plt
import pandas as pd

from gene_compare_common import (
    DEFAULT_GENE_COMPARE_DATA_ROOT,
    DEFAULT_RANKING_PLOTS_DIR,
    DEFAULT_RANKING_TABLES_DIR,
    safe_slug,
)


DEFAULT_OUT_DIR = DEFAULT_RANKING_PLOTS_DIR
DEFAULT_TOP_OVERALL = 2
DEFAULT_TOP_VARIANT = 4
DEFAULT_MIN_SPECIES = 4
DEFAULT_QUANTILE = 0.95
DEFAULT_PANELS_PER_PAGE = 6
DEFAULT_PREFERRED_VARIANTS = [
    ("H2AX", "adipose tissue"),
    ("MACROH2A2", "zone of skin"),
    ("H2AJ", "testis"),
    ("H2AJ", "pituitary gland"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create an overview scatter plot and focused species-level panels "
            "for strong cross-species H2A candidates."
        )
    )
    parser.add_argument(
        "--ranking-dir",
        default=str(DEFAULT_RANKING_TABLES_DIR),
        help="Directory with gene_tissue_summary.csv and candidate tables.",
    )
    parser.add_argument(
        "--gene-compare-data-root",
        default=str(DEFAULT_GENE_COMPARE_DATA_ROOT),
        help="Root directory with per-gene gene_compare long CSV files.",
    )
    parser.add_argument(
        "--out-dir",
        default=str(DEFAULT_OUT_DIR),
        help="Directory for output figures.",
    )
    parser.add_argument(
        "--top-overall",
        default=DEFAULT_TOP_OVERALL,
        type=int,
        help="Number of top overall candidates to include when --candidates is not used.",
    )
    parser.add_argument(
        "--top-variant",
        default=DEFAULT_TOP_VARIANT,
        type=int,
        help="Number of top variant candidates to include when --candidates is not used.",
    )
    parser.add_argument(
        "--candidates",
        nargs="*",
        default=[],
        help="Optional explicit candidates in GENE::tissue format.",
    )
    parser.add_argument(
        "--include-class-p95",
        action="store_true",
        help="Also build paginated class-specific p95 panels for clustered and variant genes.",
    )
    parser.add_argument(
        "--min-species",
        default=DEFAULT_MIN_SPECIES,
        type=int,
        help="Minimum species_n for candidate selection in overview and class-specific p95 panels.",
    )
    parser.add_argument(
        "--quantile",
        default=DEFAULT_QUANTILE,
        type=float,
        help="Quantile threshold for class-specific candidate pages when --include-class-p95 is used.",
    )
    parser.add_argument(
        "--panels-per-page",
        default=DEFAULT_PANELS_PER_PAGE,
        type=int,
        help="Maximum number of candidate panels per page for class-specific outputs.",
    )
    return parser.parse_args()


def parse_candidate_values(values: List[str]) -> List[Tuple[str, str]]:
    pairs: List[Tuple[str, str]] = []
    for value in values:
        if "::" not in value:
            raise ValueError(f"Invalid candidate format: {value}. Use GENE::tissue")
        gene_name, tissue = value.split("::", 1)
        gene_name = gene_name.strip()
        tissue = tissue.strip()
        if not gene_name or not tissue:
            raise ValueError(f"Invalid candidate format: {value}. Use GENE::tissue")
        pairs.append((gene_name, tissue))
    return pairs


def choose_default_candidates(
    summary_df: pd.DataFrame,
    high_conf_df: pd.DataFrame,
    top_overall: int,
    top_variant: int,
) -> pd.DataFrame:
    overall_df = high_conf_df.sort_values(
        by=["range", "species_n", "gene_name", "tissue"],
        ascending=[False, False, True, True],
    ).head(top_overall)
    variant_df = (
        summary_df[(summary_df["gene_class"] == "variant") & (summary_df["species_n"] >= 4)]
        .sort_values(
            by=["range", "species_n", "gene_name", "tissue"],
            ascending=[False, False, True, True],
        )
    )

    preferred_frames: List[pd.DataFrame] = []
    for gene_name, tissue in DEFAULT_PREFERRED_VARIANTS:
        preferred = variant_df[
            variant_df["gene_name"].eq(gene_name) & variant_df["tissue"].eq(tissue)
        ].copy()
        if not preferred.empty:
            preferred_frames.append(preferred.head(1))
    preferred_df = pd.concat(preferred_frames, ignore_index=True) if preferred_frames else pd.DataFrame()

    ranked_variant_df = variant_df.copy()
    chosen_variant_df = pd.concat(
        [preferred_df, ranked_variant_df],
        ignore_index=True,
    ).drop_duplicates(subset=["gene_name", "tissue"], keep="first")
    chosen_variant_df = chosen_variant_df.head(top_variant)

    chosen_df = pd.concat([overall_df, chosen_variant_df], ignore_index=True)
    chosen_df = chosen_df.drop_duplicates(subset=["gene_name", "tissue"], keep="first")
    return chosen_df.reset_index(drop=True)


def select_class_quantile_candidates(
    summary_df: pd.DataFrame,
    gene_class: str,
    min_species: int,
    quantile: float,
) -> pd.DataFrame:
    class_df = summary_df[
        (summary_df["gene_class"] == gene_class) & (summary_df["species_n"] >= int(min_species))
    ].copy()
    if class_df.empty:
        return class_df
    threshold = float(class_df["range"].quantile(quantile))
    class_df = class_df[class_df["range"] >= threshold].copy()
    return class_df.sort_values(
        by=["range", "species_n", "gene_name", "tissue"],
        ascending=[False, False, True, True],
    ).reset_index(drop=True)


def load_candidate_scores(
    gene_compare_data_root: Path,
    gene_name: str,
    tissue: str,
) -> pd.DataFrame:
    gene_slug = safe_slug(gene_name)
    long_path = gene_compare_data_root / gene_slug / f"{gene_slug}_gene_compare_long.csv"
    if not long_path.exists():
        raise FileNotFoundError(long_path)

    long_df = pd.read_csv(long_path)
    filtered = long_df[long_df["tissue"].eq(tissue)].copy()
    if filtered.empty:
        raise RuntimeError(f"No rows found for {gene_name} / {tissue} in {long_path}")

    plot_df = (
        filtered.groupby(["species_dir", "species_name"], as_index=False)["agg_expression_score"]
        .mean()
        .rename(columns={"agg_expression_score": "expression_score"})
        .sort_values(by=["expression_score", "species_name"], ascending=[True, True])
        .reset_index(drop=True)
    )
    return plot_df


def plot_overview(
    summary_df: pd.DataFrame,
    highlight_df: pd.DataFrame,
    out_png: Path,
    out_svg: Path,
    min_species: int = DEFAULT_MIN_SPECIES,
) -> None:
    plot_df = summary_df[summary_df["species_n"] >= int(min_species)].copy()
    color_map = {"clustered": "#5B7C8D", "variant": "#C06C3E"}

    fig, ax = plt.subplots(figsize=(11, 7))
    for gene_class, class_df in plot_df.groupby("gene_class", dropna=False):
        ax.scatter(
            class_df["species_n"],
            class_df["range"],
            s=class_df["species_n"] * 28,
            alpha=0.65,
            color=color_map.get(gene_class, "#888888"),
            label=gene_class or "unclassified",
            edgecolor="white",
            linewidth=0.6,
        )

    for row in highlight_df.to_dict(orient="records"):
        ax.scatter(
            row["species_n"],
            row["range"],
            s=row["species_n"] * 42,
            color="#B22222",
            edgecolor="black",
            linewidth=0.8,
            zorder=5,
        )
        ax.annotate(
            f'{row["gene_name"]}\n{row["tissue"]}',
            (row["species_n"], row["range"]),
            xytext=(6, 4),
            textcoords="offset points",
            fontsize=8,
        )

    ax.set_title("Cross-species H2A candidates by range and species coverage")
    ax.set_xlabel("Number of species with tissue-level score")
    ax.set_ylabel("Expression score range across species")
    ax.legend(title="Gene class")
    ax.grid(alpha=0.2)
    plt.tight_layout()

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    plt.close(fig)


def plot_candidate_panels(
    candidate_df: pd.DataFrame,
    gene_compare_data_root: Path,
    out_png: Path,
    out_svg: Path,
) -> None:
    if candidate_df.empty:
        raise RuntimeError("No candidates available for focused panels.")

    n_panels = len(candidate_df)
    n_cols = 2
    n_rows = (n_panels + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(14, max(4 * n_rows, 6)))
    axes_flat = axes.flatten() if hasattr(axes, "flatten") else [axes]
    color_map = {"clustered": "#5B7C8D", "variant": "#C06C3E"}

    for ax, row in zip(axes_flat, candidate_df.to_dict(orient="records")):
        plot_df = load_candidate_scores(
            gene_compare_data_root,
            row["gene_name"],
            row["tissue"],
        )
        bars = ax.barh(
            plot_df["species_name"],
            plot_df["expression_score"],
            color=color_map.get(row["gene_class"], "#888888"),
            alpha=0.9,
        )
        ax.set_title(
            f'{row["gene_name"]} | {row["tissue"]}\n'
            f'range {row["range"]:.2f}, {int(row["species_n"])} species'
        )
        ax.set_xlabel("Expression score")
        ax.set_ylabel("")
        ax.set_xlim(0, max(100, float(plot_df["expression_score"].max()) + 5))

        for bar, value in zip(bars, plot_df["expression_score"]):
            ax.text(
                value + 0.8,
                bar.get_y() + bar.get_height() / 2,
                f"{value:.2f}",
                va="center",
                fontsize=8,
            )

    for ax in axes_flat[n_panels:]:
        ax.axis("off")

    plt.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    plt.close(fig)


def write_paginated_panels(
    candidate_df: pd.DataFrame,
    gene_compare_data_root: Path,
    out_dir: Path,
    stem_prefix: str,
    panels_per_page: int,
) -> List[Path]:
    if candidate_df.empty:
        return []

    written_paths: List[Path] = []
    total = len(candidate_df)
    for page_start in range(0, total, panels_per_page):
        page_idx = page_start // panels_per_page + 1
        page_df = candidate_df.iloc[page_start : page_start + panels_per_page].reset_index(drop=True)
        out_png = out_dir / f"{stem_prefix}_page{page_idx}.png"
        out_svg = out_dir / f"{stem_prefix}_page{page_idx}.svg"
        plot_candidate_panels(page_df, gene_compare_data_root, out_png, out_svg)
        written_paths.extend([out_png, out_svg])
    return written_paths


def main() -> None:
    args = parse_args()

    if not 0 < args.quantile < 1:
        raise ValueError("--quantile must be between 0 and 1")
    if args.panels_per_page < 1:
        raise ValueError("--panels-per-page must be >= 1")

    ranking_dir = Path(args.ranking_dir)
    gene_compare_data_root = Path(args.gene_compare_data_root)
    out_dir = Path(args.out_dir)

    summary_path = ranking_dir / "gene_tissue_summary.csv"
    high_conf_path = ranking_dir / "high_confidence_candidates.csv"
    if not summary_path.exists():
        raise FileNotFoundError(summary_path)
    if not high_conf_path.exists():
        raise FileNotFoundError(high_conf_path)

    summary_df = pd.read_csv(summary_path)
    high_conf_df = pd.read_csv(high_conf_path)

    explicit_pairs = parse_candidate_values(args.candidates)
    if explicit_pairs:
        key_df = pd.DataFrame(explicit_pairs, columns=["gene_name", "tissue"])
        candidate_df = key_df.merge(summary_df, on=["gene_name", "tissue"], how="left")
        if candidate_df.isna().any(axis=None):
            missing = candidate_df[candidate_df["range"].isna()][["gene_name", "tissue"]]
            raise RuntimeError(
                "Missing candidate(s) in summary: "
                + ", ".join(f"{r.gene_name}::{r.tissue}" for r in missing.itertuples())
            )
    else:
        candidate_df = choose_default_candidates(
            summary_df,
            high_conf_df,
            top_overall=args.top_overall,
            top_variant=args.top_variant,
        )

    out_overview_png = out_dir / "candidate_range_overview.png"
    out_overview_svg = out_dir / "candidate_range_overview.svg"
    out_panels_png = out_dir / "candidate_focus_panels.png"
    out_panels_svg = out_dir / "candidate_focus_panels.svg"

    plot_overview(
        summary_df,
        candidate_df,
        out_overview_png,
        out_overview_svg,
        min_species=args.min_species,
    )
    plot_candidate_panels(candidate_df, gene_compare_data_root, out_panels_png, out_panels_svg)

    print(f"Saved overview plot to {out_overview_png}")
    print(f"Saved overview plot to {out_overview_svg}")
    print(f"Saved candidate panels to {out_panels_png}")
    print(f"Saved candidate panels to {out_panels_svg}")

    if args.include_class_p95:
        for gene_class in ["clustered", "variant"]:
            class_df = select_class_quantile_candidates(
                summary_df,
                gene_class=gene_class,
                min_species=args.min_species,
                quantile=args.quantile,
            )
            written = write_paginated_panels(
                class_df,
                gene_compare_data_root=gene_compare_data_root,
                out_dir=out_dir,
                stem_prefix=f"{gene_class}_p{int(args.quantile * 100):02d}_panels",
                panels_per_page=args.panels_per_page,
            )
            if written:
                print(
                    f"Saved {len(class_df)} {gene_class} candidate panel(s) across "
                    f"{len(written) // 2} page(s) for p{int(args.quantile * 100):02d}."
                )
            else:
                print(f"No {gene_class} candidates met the class-specific p{int(args.quantile * 100):02d} threshold.")


if __name__ == "__main__":
    main()
