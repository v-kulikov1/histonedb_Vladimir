#!/usr/bin/env python
"""Build a human H2A boxplot from the heatmap-aligned normalized tissue table."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Sequence

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd

from gene_compare_common import GENERIC_TISSUES, normalize_map_df
from normalized_expression_common import (
    build_species_heatmap_display_index,
    load_processed_expression_cells,
    sort_gene_labels,
)


DEFAULT_EXPR_TSV = Path(
    r"CURATED_SET/BioAnalyze/data/processed/homo_sapiens/homo_sapiens_expr_advanced_H2A_present_gold.tsv"
)
DEFAULT_MAP_TSV = Path(
    r"CURATED_SET/BioAnalyze/data/processed/homo_sapiens/h2a_hs_canonical_variant_map.tsv"
)
DEFAULT_OUT_DIR = Path(r"CURATED_SET/BioAnalyze/figures/boxplot/human")
DEFAULT_GENERIC_TISSUES = sorted(GENERIC_TISSUES)
PRESENTATION_GENE_ORDER = [
    "H2AC1",
    "H2AC4",
    "H2AC6",
    "H2AC7",
    "H2AC8",
    "H2AC11",
    "H2AC12",
    "H2AC13",
    "H2AC14",
    "H2AC15",
    "H2AC16",
    "H2AC17",
    "H2AC18",
    "H2AC19",
    "H2AC20",
    "H2AC21",
    "H2AC25",
    "H2AB1",
    "H2AB2",
    "H2AB3",
    "H2AJ",
    "H2AL1Q",
    "H2AL3",
    "H2AP",
    "H2AX",
    "H2AZ1",
    "H2AZ2",
    "MACROH2A1",
    "MACROH2A2",
]
CLASS_COLORS = {
    "clustered": "#6ea6c8",
    "variant": "#d7b46a",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a human H2A boxplot from the normalized Gene x tissue table used "
            "for the species heatmaps."
        )
    )
    parser.add_argument(
        "--expr-tsv",
        default=str(DEFAULT_EXPR_TSV),
        help="Processed human H2A present_gold TSV with one normalized row per Gene x tissue.",
    )
    parser.add_argument(
        "--map-tsv",
        default=str(DEFAULT_MAP_TSV),
        help="Canonical/variant map TSV for the human H2A heatmap rows.",
    )
    parser.add_argument(
        "--out-dir",
        default=str(DEFAULT_OUT_DIR),
        help="Output directory for boxplot figures and tables.",
    )
    parser.add_argument(
        "--exclude-tissue",
        action="append",
        default=[],
        help=(
            "Additional tissue name to exclude. May be passed multiple times. "
            "Default exclusions are the generic tissues used by BioAnalyze."
        ),
    )
    return parser.parse_args()


def add_presentation_numbers(display_df: pd.DataFrame) -> pd.DataFrame:
    gene_names = display_df["gene_name"].astype(str).tolist()
    missing = [gene for gene in PRESENTATION_GENE_ORDER if gene not in gene_names]
    extras = sorted(set(gene_names) - set(PRESENTATION_GENE_ORDER))
    if missing or extras:
        details: List[str] = []
        if missing:
            details.append(f"missing={','.join(missing)}")
        if extras:
            details.append(f"extras={','.join(extras)}")
        raise RuntimeError(
            "Presentation order does not match the discovered human H2A genes: "
            + "; ".join(details)
        )

    presentation_map = {
        gene_name: idx + 1 for idx, gene_name in enumerate(PRESENTATION_GENE_ORDER)
    }
    display_df = display_df.copy()
    display_df["presentation_number"] = (
        display_df["gene_name"].map(presentation_map).astype(int)
    )
    return display_df


def build_display_index(map_tsv: Path, expr_tsv: Path) -> pd.DataFrame:
    map_df = normalize_map_df(map_tsv)
    expr_df = load_processed_expression_cells(expr_tsv)
    display_df = build_species_heatmap_display_index(
        map_df,
        expr_df,
        preferred_id_col="hgnc_id",
    )
    if display_df.empty:
        raise RuntimeError("No displayable human H2A genes were found in the provided files.")

    label_map = dict(
        zip(
            display_df["ensembl_gene_id"].astype(str),
            display_df["label"].astype(str),
        )
    )
    ordered_labels = sort_gene_labels(
        label_map,
        display_df["ensembl_gene_id"].astype(str).tolist(),
    )
    order_map = {label: idx + 1 for idx, label in enumerate(ordered_labels)}
    display_df["display_order"] = display_df["label"].map(order_map)
    display_df = display_df.sort_values(
        by=["display_order", "gene_name", "ensembl_gene_id"],
        ascending=[True, True, True],
        kind="stable",
    ).reset_index(drop=True)
    display_df = add_presentation_numbers(display_df)
    return display_df[
        [
            "ensembl_gene_id",
            "gene_name",
            "class",
            "label",
            "display_order",
            "presentation_number",
        ]
    ].copy()


def build_source_table(
    expr_tsv: Path,
    display_df: pd.DataFrame,
    excluded_tissues: List[str],
) -> pd.DataFrame:
    expr_df = load_processed_expression_cells(expr_tsv).copy()
    expr_df = expr_df.rename(
        columns={
            "Gene ID": "ensembl_gene_id",
            "Gene name": "raw_gene_name",
            "Anatomical entity name": "tissue",
        }
    )
    expr_df["ensembl_gene_id"] = expr_df["ensembl_gene_id"].fillna("").astype(str).str.strip()
    expr_df["tissue"] = expr_df["tissue"].fillna("").astype(str).str.strip()
    expr_df = expr_df[
        expr_df["ensembl_gene_id"].isin(display_df["ensembl_gene_id"])
        & expr_df["tissue"].ne("")
        & ~expr_df["tissue"].isin(excluded_tissues)
    ].copy()

    if expr_df.empty:
        raise RuntimeError("All rows were filtered out before building the human H2A boxplot.")

    source_df = expr_df.merge(display_df, on="ensembl_gene_id", how="inner")
    source_df["gene_name"] = source_df["gene_name"].fillna("").astype(str).str.strip()
    source_df["label"] = source_df["label"].fillna("").astype(str).str.strip()
    source_df["class"] = source_df["class"].fillna("").astype(str).str.strip()
    source_df["display_order"] = pd.to_numeric(
        source_df["display_order"], errors="coerce"
    ).fillna(0).astype(int)
    source_df["presentation_number"] = pd.to_numeric(
        source_df["presentation_number"], errors="coerce"
    ).fillna(0).astype(int)
    source_df["cell_mean_score"] = pd.to_numeric(
        source_df["cell_mean_score"], errors="coerce"
    ).astype(float)
    source_df["cell_std_score"] = pd.to_numeric(
        source_df["cell_std_score"], errors="coerce"
    ).fillna(0.0).astype(float)
    source_df["cell_n"] = pd.to_numeric(source_df["cell_n"], errors="coerce").fillna(0).astype(int)
    source_df["cell_status"] = source_df["cell_status"].fillna("").astype(str).str.strip()

    source_df = source_df[
        [
            "display_order",
            "presentation_number",
            "ensembl_gene_id",
            "gene_name",
            "label",
            "class",
            "tissue",
            "cell_mean_score",
            "cell_std_score",
            "cell_n",
            "cell_status",
        ]
    ].sort_values(
        by=["display_order", "gene_name", "tissue"],
        ascending=[True, True, True],
        kind="stable",
    ).reset_index(drop=True)

    missing_ids = sorted(set(display_df["ensembl_gene_id"]) - set(source_df["ensembl_gene_id"]))
    if missing_ids:
        raise RuntimeError(
            "Some displayable human H2A genes have no non-generic tissues after filtering: "
            + ", ".join(missing_ids)
        )
    return source_df


def compute_boxplot_stats(values: pd.Series) -> Dict[str, float]:
    arr = np.sort(values.astype(float).to_numpy())
    if arr.size == 0:
        raise RuntimeError("Cannot compute boxplot statistics for an empty value set.")

    q1 = float(np.percentile(arr, 25))
    median = float(np.percentile(arr, 50))
    q3 = float(np.percentile(arr, 75))
    iqr = q3 - q1
    lower_fence = q1 - 1.5 * iqr
    upper_fence = q3 + 1.5 * iqr

    inside = arr[(arr >= lower_fence) & (arr <= upper_fence)]
    if inside.size == 0:
        lower_whisker = float(arr.min())
        upper_whisker = float(arr.max())
    else:
        lower_whisker = float(inside.min())
        upper_whisker = float(inside.max())

    outlier_count = int(((arr < lower_fence) | (arr > upper_fence)).sum())
    std_value = float(np.std(arr, ddof=1)) if arr.size > 1 else 0.0

    return {
        "n_tissues": int(arr.size),
        "mean": float(arr.mean()),
        "std": std_value,
        "min": float(arr.min()),
        "q1": q1,
        "median": median,
        "q3": q3,
        "max": float(arr.max()),
        "iqr": float(iqr),
        "lower_whisker": lower_whisker,
        "upper_whisker": upper_whisker,
        "outlier_count": outlier_count,
    }


def build_stats_table(source_df: pd.DataFrame, display_df: pd.DataFrame) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for display_row in display_df.to_dict(orient="records"):
        gene_values = source_df[
            source_df["ensembl_gene_id"].eq(display_row["ensembl_gene_id"])
        ]["cell_mean_score"]
        stats = compute_boxplot_stats(gene_values)
        rows.append(
            {
                "display_order": int(display_row["display_order"]),
                "presentation_number": int(display_row["presentation_number"]),
                "ensembl_gene_id": display_row["ensembl_gene_id"],
                "gene_name": display_row["gene_name"],
                "label": display_row["label"],
                "class": display_row["class"],
                **stats,
            }
        )

    return pd.DataFrame(rows).sort_values(
        by=["display_order", "gene_name"],
        ascending=[True, True],
        kind="stable",
    ).reset_index(drop=True)


def build_presentation_stats(stats_df: pd.DataFrame) -> pd.DataFrame:
    return stats_df.sort_values(
        by=["presentation_number", "gene_name"],
        ascending=[True, True],
        kind="stable",
    ).reset_index(drop=True)


def plot_boxplot(
    source_df: pd.DataFrame,
    stats_df: pd.DataFrame,
    out_png: Path,
    out_svg: Path,
    *,
    tick_labels: Sequence[str],
    y_label: str,
    x_label_scale: float = 1.0,
    x_tick_scale: float = 1.0,
    y_tick_scale: float = 1.0,
    x_grid_color: str = "#d9dee2",
    x_grid_width: float = 0.8,
    y_grid_color: str = "#ffffff",
    y_grid_width: float = 0.0,
    layout_rect: tuple[float, float, float, float] = (0.0, 0.0, 1.0, 1.0),
    title: str = "Human H2A expression across tissues",
) -> None:
    gene_order = stats_df["gene_name"].astype(str).tolist()
    gene_class_map = stats_df.set_index("gene_name")["class"].to_dict()
    data = [
        source_df.loc[source_df["gene_name"].eq(gene_name), "cell_mean_score"].astype(float).to_numpy()
        for gene_name in gene_order
    ]

    fig_w = 14
    fig_h = max(8, len(gene_order) * 0.42)

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    box = ax.boxplot(
        data,
        vert=False,
        tick_labels=list(tick_labels),
        patch_artist=True,
        widths=0.68,
        showfliers=True,
        flierprops={
            "marker": "o",
            "markerfacecolor": "#4c566a",
            "markeredgecolor": "#4c566a",
            "markersize": 3,
            "alpha": 0.55,
        },
        medianprops={"color": "#111111", "linewidth": 1.6},
        whiskerprops={"color": "#444444", "linewidth": 1.1},
        capprops={"color": "#444444", "linewidth": 1.1},
        boxprops={"edgecolor": "#444444", "linewidth": 1.1},
    )

    for patch, gene_name in zip(box["boxes"], gene_order):
        patch.set_facecolor(CLASS_COLORS.get(gene_class_map.get(gene_name, ""), "#b7c0c8"))
        patch.set_alpha(0.95)

    ax.invert_yaxis()
    ax.set_xlabel("Expression score")
    base_x_label_size = float(ax.xaxis.label.get_size())
    ax.xaxis.label.set_size(base_x_label_size * float(x_label_scale))
    ax.set_ylabel(y_label)
    ax.set_title(title)
    x_tick_size = float(ax.xaxis.get_ticklabels()[0].get_size()) if ax.xaxis.get_ticklabels() else 10.0
    y_tick_size = float(ax.yaxis.get_ticklabels()[0].get_size()) if ax.yaxis.get_ticklabels() else 10.0
    ax.tick_params(axis="x", labelsize=x_tick_size * float(x_tick_scale))
    ax.tick_params(axis="y", labelsize=y_tick_size * float(y_tick_scale))
    ax.grid(axis="x", color=x_grid_color, linewidth=x_grid_width)
    ax.grid(axis="y", color=y_grid_color, linewidth=y_grid_width)
    ax.set_axisbelow(True)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=layout_rect)
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    plt.close(fig)


def plot_number_legend(stats_df: pd.DataFrame, out_png: Path) -> None:
    fig_h = max(8, len(stats_df) * 0.38)
    fig, ax = plt.subplots(figsize=(8.2, fig_h))
    ax.axis("off")

    n_rows = len(stats_df)
    top = 0.965
    bottom = 0.035
    step = (top - bottom) / max(n_rows - 1, 1)
    x_number = 0.12
    x_swatch = 0.17
    swatch_size = 0.026
    x_gene = 0.23

    for idx, row in enumerate(stats_df.to_dict(orient="records"), start=0):
        y = top - idx * step
        color = CLASS_COLORS.get(str(row["class"]), "#333333")
        ax.text(
            x_number,
            y,
            f"{int(row['presentation_number'])}",
            transform=ax.transAxes,
            va="center",
            ha="right",
            fontsize=17,
            fontweight="bold",
            color="#111111",
        )
        ax.add_patch(
            Rectangle(
                (x_swatch, y - swatch_size / 2),
                swatch_size,
                swatch_size,
                transform=ax.transAxes,
                facecolor=color,
                edgecolor="#444444",
                linewidth=0.8,
            )
        )
        ax.text(
            x_gene,
            y,
            f"{row['gene_name']}",
            transform=ax.transAxes,
            va="center",
            ha="left",
            fontsize=17,
            fontweight="normal",
            color="#111111",
        )

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=(0.02, 0.01, 0.995, 0.995))
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


def run(expr_tsv: Path, map_tsv: Path, out_dir: Path, excluded_tissues: List[str]) -> Dict[str, Path]:
    display_df = build_display_index(map_tsv, expr_tsv)
    source_df = build_source_table(expr_tsv, display_df, excluded_tissues)
    stats_df = build_stats_table(source_df, display_df)
    presentation_stats_df = build_presentation_stats(stats_df)

    out_dir.mkdir(parents=True, exist_ok=True)
    out_png = out_dir / "h2a_human_boxplot.png"
    out_svg = out_dir / "h2a_human_boxplot.svg"
    out_source_tsv = out_dir / "h2a_human_boxplot_source_values.tsv"
    out_stats_tsv = out_dir / "h2a_human_boxplot_stats.tsv"
    out_presentation_png = out_dir / "h2a_human_boxplot_presentation.png"
    out_presentation_svg = out_dir / "h2a_human_boxplot_presentation.svg"
    out_legend_png = out_dir / "h2a_human_boxplot_number_legend.png"

    source_df.to_csv(out_source_tsv, sep="\t", index=False)
    stats_df.to_csv(out_stats_tsv, sep="\t", index=False)
    plot_boxplot(
        source_df,
        stats_df,
        out_png,
        out_svg,
        tick_labels=stats_df["gene_name"].astype(str).tolist(),
        y_label="Gene name",
    )
    plot_boxplot(
        source_df,
        presentation_stats_df,
        out_presentation_png,
        out_presentation_svg,
        tick_labels=presentation_stats_df["presentation_number"].astype(str).tolist(),
        y_label="",
        x_label_scale=2.0,
        x_tick_scale=2.0,
        y_tick_scale=2.0,
        x_grid_color="#bec8d2",
        x_grid_width=1.1,
        y_grid_color="#e3e9ee",
        y_grid_width=0.7,
        layout_rect=(0.08, 0.08, 0.995, 0.99),
        title="",
    )
    plot_number_legend(presentation_stats_df, out_legend_png)

    return {
        "out_png": out_png,
        "out_svg": out_svg,
        "out_source_tsv": out_source_tsv,
        "out_stats_tsv": out_stats_tsv,
        "out_presentation_png": out_presentation_png,
        "out_presentation_svg": out_presentation_svg,
        "out_legend_png": out_legend_png,
    }


def main() -> None:
    args = parse_args()
    excluded_tissues = sorted(
        {
            tissue.strip()
            for tissue in (DEFAULT_GENERIC_TISSUES + list(args.exclude_tissue))
            if tissue.strip()
        }
    )
    outputs = run(
        expr_tsv=Path(args.expr_tsv),
        map_tsv=Path(args.map_tsv),
        out_dir=Path(args.out_dir),
        excluded_tissues=excluded_tissues,
    )

    print(f"Excluded tissues: {', '.join(excluded_tissues)}")
    print(f"Saved boxplot PNG to {outputs['out_png']}")
    print(f"Saved boxplot SVG to {outputs['out_svg']}")
    print(f"Saved source table to {outputs['out_source_tsv']}")
    print(f"Saved stats table to {outputs['out_stats_tsv']}")
    print(f"Saved presentation PNG to {outputs['out_presentation_png']}")
    print(f"Saved presentation SVG to {outputs['out_presentation_svg']}")
    print(f"Saved legend PNG to {outputs['out_legend_png']}")


if __name__ == "__main__":
    main()
