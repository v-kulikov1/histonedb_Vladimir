#!/usr/bin/env python
"""Build presentation-friendly human/pan H2A heatmaps with shared legends."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colorbar import ColorbarBase
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
import seaborn as sns

BIOANALYZE_SCRIPTS_ROOT = Path(__file__).resolve().parents[1]
if str(BIOANALYZE_SCRIPTS_ROOT) not in sys.path:
    sys.path.insert(0, str(BIOANALYZE_SCRIPTS_ROOT))

from bioanalyze_paths import get_bioanalyze_data_root, get_bioanalyze_figures_root
from gene_compare_common import normalize_map_df
from normalized_expression_common import load_processed_expression_cells


DEFAULT_ALIGNED_EXPR = (
    get_bioanalyze_data_root() / "processed" / "intersections" / "homo_sapiens_expr_advanced_H2A_present_gold_intersection.tsv"
)
DEFAULT_PAN_EXPR = get_bioanalyze_data_root() / "processed" / "pan_troglodytes" / "pan_troglodytes_expr_advanced_H2A_present_gold.tsv"
DEFAULT_PAN_MAP = get_bioanalyze_data_root() / "processed" / "pan_troglodytes" / "pan_troglodytes_h2a_canonical_variant_map.tsv"
DEFAULT_OUT_DIR = get_bioanalyze_figures_root() / "heatmaps" / "presentation_human_pan"
CLASS_COLORS = {
    "clustered": "#6ea6c8",
    "variant": "#d7b46a",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build presentation-friendly shared H2A heatmaps for Homo sapiens "
            "and Pan troglodytes with separate legends."
        )
    )
    parser.add_argument(
        "--aligned-expr",
        default=str(DEFAULT_ALIGNED_EXPR),
        help="Processed aligned/intersection human TSV used by hs_aligned_all.",
    )
    parser.add_argument(
        "--pan-expr",
        default=str(DEFAULT_PAN_EXPR),
        help="Processed Pan troglodytes TSV used by h2a_pan_troglodytes_all.",
    )
    parser.add_argument(
        "--pan-map",
        default=str(DEFAULT_PAN_MAP),
        help="Pan troglodytes canonical/variant map TSV for canonical gene names.",
    )
    parser.add_argument(
        "--out-dir",
        default=str(DEFAULT_OUT_DIR),
        help="Output directory for presentation heatmaps, legends, and mapping TSVs.",
    )
    return parser.parse_args()


def _clean_text_columns(df: pd.DataFrame, columns: Sequence[str]) -> pd.DataFrame:
    for col in columns:
        if col not in df.columns:
            df[col] = ""
        df[col] = df[col].fillna("").astype(str).str.strip()
    return df


def load_aligned_expression(expr_tsv: Path) -> pd.DataFrame:
    df = pd.read_csv(expr_tsv, sep="\t", low_memory=True)
    df = _clean_text_columns(
        df,
        ["Gene ID", "Gene name", "Anatomical entity name", "merged_gene_name", "class"],
    )
    for col in ["cell_mean_score", "Expression score"]:
        if col not in df.columns:
            df[col] = np.nan
        df[col] = pd.to_numeric(df[col], errors="coerce")

    if "cell_mean_score" not in df.columns or df["cell_mean_score"].isna().all():
        df["cell_mean_score"] = df["Expression score"]
    df = df[
        df["Gene ID"].ne("")
        & df["Anatomical entity name"].ne("")
        & df["merged_gene_name"].ne("")
        & df["class"].ne("")
        & df["cell_mean_score"].notna()
    ].copy()
    return df


def build_pan_gene_metadata(
    pan_expr_df: pd.DataFrame,
    pan_map_df: pd.DataFrame,
) -> pd.DataFrame:
    pan_expr_df = pan_expr_df.copy()
    pan_expr_df = _clean_text_columns(
        pan_expr_df,
        ["Gene ID", "Gene name", "Anatomical entity name"],
    )
    pan_gene_ids = (
        pan_expr_df["Gene ID"].dropna().astype(str).str.strip().loc[lambda s: s.ne("")].unique().tolist()
    )
    meta_df = pan_map_df.copy()
    meta_df = _clean_text_columns(
        meta_df,
        ["ensembl_gene_id", "gene_name", "class", "label"],
    )
    meta_df = meta_df[meta_df["ensembl_gene_id"].isin(pan_gene_ids)].copy()
    meta_df = meta_df.sort_values(
        by=["gene_name", "label", "ensembl_gene_id"],
        ascending=[True, True, True],
        kind="stable",
    )
    meta_df = meta_df.drop_duplicates(subset=["gene_name"], keep="first")
    return meta_df[["ensembl_gene_id", "gene_name", "class"]].copy()


def build_shared_metadata(
    aligned_df: pd.DataFrame,
    pan_expr_df: pd.DataFrame,
    pan_map_df: pd.DataFrame,
) -> Tuple[pd.DataFrame, List[str], List[str]]:
    aligned_meta = (
        aligned_df[["Gene ID", "merged_gene_name", "class"]]
        .drop_duplicates()
        .rename(
            columns={
                "Gene ID": "aligned_gene_id",
                "merged_gene_name": "gene_name",
                "class": "class",
            }
        )
    )
    aligned_meta = _clean_text_columns(aligned_meta, ["aligned_gene_id", "gene_name", "class"])
    aligned_meta = aligned_meta[aligned_meta["gene_name"].ne("")].copy()
    aligned_meta = aligned_meta.sort_values(
        by=["gene_name", "aligned_gene_id"],
        ascending=[True, True],
        kind="stable",
    ).drop_duplicates(subset=["gene_name"], keep="first")

    pan_meta = build_pan_gene_metadata(pan_expr_df, pan_map_df).rename(
        columns={"ensembl_gene_id": "pan_gene_id"}
    )
    pan_meta = _clean_text_columns(pan_meta, ["pan_gene_id", "gene_name", "class"])

    shared_gene_names = sorted(
        set(aligned_meta["gene_name"].tolist()) & set(pan_meta["gene_name"].tolist())
    )
    if not shared_gene_names:
        raise RuntimeError("No shared gene names were found between aligned human and pan inputs.")

    shared_meta = aligned_meta[aligned_meta["gene_name"].isin(shared_gene_names)].copy()
    shared_meta = shared_meta.merge(
        pan_meta[["pan_gene_id", "gene_name"]],
        on="gene_name",
        how="left",
    )
    shared_meta = shared_meta.sort_values(
        by=["gene_name", "aligned_gene_id", "pan_gene_id"],
        ascending=[True, True, True],
        kind="stable",
    ).reset_index(drop=True)
    shared_meta["gene_number"] = np.arange(1, len(shared_meta) + 1, dtype=int)

    aligned_tissues = sorted(
        aligned_df["Anatomical entity name"]
        .dropna()
        .astype(str)
        .str.strip()
        .loc[lambda s: s.ne("")]
        .unique()
        .tolist()
    )
    pan_tissues = set(
        pan_expr_df["Anatomical entity name"]
        .dropna()
        .astype(str)
        .str.strip()
        .loc[lambda s: s.ne("")]
        .unique()
        .tolist()
    )
    shared_tissues = [tissue for tissue in aligned_tissues if tissue in pan_tissues]
    if not shared_tissues:
        raise RuntimeError("No shared tissues remained after intersecting aligned and pan inputs.")
    if len(shared_tissues) > 26:
        raise RuntimeError(
            "Too many shared tissues for single-letter labels; extend the labeling scheme first."
        )

    shared_meta = shared_meta[
        shared_meta["pan_gene_id"].fillna("").astype(str).str.strip().ne("")
    ].reset_index(drop=True)
    if shared_meta.empty:
        raise RuntimeError("Shared gene metadata is empty after joining aligned and pan IDs.")

    tissue_letters = [chr(ord("A") + idx) for idx in range(len(shared_tissues))]
    return shared_meta, shared_tissues, tissue_letters


def build_matrix(
    expr_df: pd.DataFrame,
    *,
    gene_id_col: str,
    gene_name_col: str,
    shared_gene_names: Sequence[str],
    tissue_names: Sequence[str],
) -> pd.DataFrame:
    df = expr_df.copy()
    df = _clean_text_columns(df, [gene_id_col, gene_name_col, "Anatomical entity name"])
    if "cell_mean_score" not in df.columns:
        if "Expression score" not in df.columns:
            raise RuntimeError("Input expression table is missing both cell_mean_score and Expression score.")
        df["cell_mean_score"] = pd.to_numeric(df["Expression score"], errors="coerce")
    else:
        df["cell_mean_score"] = pd.to_numeric(df["cell_mean_score"], errors="coerce")

    df = df[
        df[gene_name_col].isin(shared_gene_names)
        & df["Anatomical entity name"].isin(tissue_names)
        & df["cell_mean_score"].notna()
    ].copy()
    if df.empty:
        raise RuntimeError("Filtered expression table is empty while building the presentation matrix.")

    matrix = df.pivot_table(
        index=gene_name_col,
        columns="Anatomical entity name",
        values="cell_mean_score",
        aggfunc="mean",
    )
    matrix = matrix.reindex(index=list(shared_gene_names), columns=list(tissue_names))
    if matrix.shape != (len(shared_gene_names), len(tissue_names)):
        raise RuntimeError(
            f"Unexpected matrix shape {matrix.shape}; expected ({len(shared_gene_names)}, {len(tissue_names)})."
        )
    if matrix.notna().sum().sum() == 0:
        raise RuntimeError("Presentation matrix contains no non-missing values.")
    return matrix


def build_display_maps(
    shared_meta: pd.DataFrame,
    tissue_names: Sequence[str],
    tissue_letters: Sequence[str],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    gene_map = shared_meta[
        ["gene_number", "gene_name", "class", "aligned_gene_id", "pan_gene_id"]
    ].copy()
    gene_map = gene_map.sort_values("gene_number", kind="stable").reset_index(drop=True)

    tissue_map = pd.DataFrame(
        {
            "tissue_letter": list(tissue_letters),
            "tissue_name": list(tissue_names),
        }
    )
    tissue_map["axis"] = "OX"
    gene_map["axis"] = "OY"
    return gene_map, tissue_map


def save_map_tables(gene_map: pd.DataFrame, tissue_map: pd.DataFrame, out_dir: Path) -> None:
    gene_map = gene_map[
        ["axis", "gene_number", "gene_name", "class", "aligned_gene_id", "pan_gene_id"]
    ].copy()
    tissue_map = tissue_map[["axis", "tissue_letter", "tissue_name"]].copy()
    gene_map.to_csv(out_dir / "human_pan_gene_number_map.tsv", sep="\t", index=False)
    tissue_map.to_csv(out_dir / "human_pan_tissue_letter_map.tsv", sep="\t", index=False)


def plot_heatmap(
    matrix: pd.DataFrame,
    *,
    species_title: str,
    out_png: Path,
    out_svg: Path,
    gene_numbers: Sequence[int],
    tissue_letters: Sequence[str],
) -> None:
    sns.set_theme(style="white")
    mat_log = np.log10(matrix + 1.0)

    fig_w = max(16.0, len(tissue_letters) * 0.62)
    fig_h = max(10.0, len(gene_numbers) * 0.62)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    sns.heatmap(
        mat_log,
        cmap="viridis",
        vmin=0.0,
        vmax=2.0,
        mask=mat_log.isna(),
        cbar=False,
        linewidths=0.35,
        linecolor="#f0f3f6",
        xticklabels=list(tissue_letters),
        yticklabels=[str(num) for num in gene_numbers],
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title(species_title, fontsize=24, fontweight="bold", pad=16)
    ax.tick_params(axis="x", labelrotation=0, labelsize=16, length=0, pad=8)
    ax.tick_params(axis="y", labelrotation=0, labelsize=16, length=0, pad=8)
    for spine in ax.spines.values():
        spine.set_visible(False)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=(0.02, 0.02, 0.995, 0.96))
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    plt.close(fig)


def plot_gene_legend(gene_map: pd.DataFrame, out_png: Path, out_svg: Path) -> None:
    fig_h = max(7.0, len(gene_map) * 0.45 + 1.0)
    fig, ax = plt.subplots(figsize=(8.8, fig_h))
    ax.axis("off")

    ax.text(
        0.05,
        0.985,
        "OY: Gene index",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=20,
        fontweight="bold",
        color="#111111",
    )

    n_rows = len(gene_map)
    top = 0.91
    bottom = 0.05
    step = (top - bottom) / max(n_rows - 1, 1)
    x_number = 0.14
    x_swatch = 0.19
    swatch_size = 0.03
    x_gene = 0.255

    for idx, row in enumerate(gene_map.to_dict(orient="records")):
        y = top - idx * step
        color = CLASS_COLORS.get(str(row["class"]), "#333333")
        ax.text(
            x_number,
            y,
            str(int(row["gene_number"])),
            transform=ax.transAxes,
            va="center",
            ha="right",
            fontsize=18,
            fontweight="bold",
            color="#111111",
        )
        ax.add_patch(
            Rectangle(
                (x_swatch, y - swatch_size / 2.0),
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
            str(row["gene_name"]),
            transform=ax.transAxes,
            va="center",
            ha="left",
            fontsize=18,
            color="#111111",
        )

    fig.tight_layout(rect=(0.02, 0.01, 0.995, 0.995))
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    plt.close(fig)


def plot_tissue_legend(tissue_map: pd.DataFrame, out_png: Path, out_svg: Path) -> None:
    entries = tissue_map.to_dict(orient="records")
    if not entries:
        raise RuntimeError("Tissue map is empty; cannot build tissue legend.")

    fig_h = max(8.0, len(entries) * 0.50 + 1.2)
    fig_w = 10.5
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.axis("off")

    ax.text(
        0.05,
        0.985,
        "OX: Anatomical entity",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=24,
        fontweight="bold",
        color="#111111",
    )

    top = 0.92
    bottom = 0.04
    step = (top - bottom) / max(len(entries) - 1, 1)
    x_letter = 0.13
    x_name = 0.21

    for idx, row in enumerate(entries):
        y = top - idx * step
        ax.text(
            x_letter,
            y,
            str(row["tissue_letter"]),
            transform=ax.transAxes,
            va="center",
            ha="right",
            fontsize=18,
            fontweight="bold",
            color="#111111",
        )
        ax.text(
            x_name,
            y,
            str(row["tissue_name"]),
            transform=ax.transAxes,
            va="center",
            ha="left",
            fontsize=18,
            color="#111111",
        )

    fig.tight_layout(rect=(0.02, 0.01, 0.995, 0.995))
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    plt.close(fig)


def plot_color_legend(out_png: Path, out_svg: Path) -> None:
    fig = plt.figure(figsize=(8.0, 2.2))
    ax = fig.add_axes([0.08, 0.40, 0.84, 0.24])
    norm = mpl.colors.Normalize(vmin=0.0, vmax=2.0)
    ColorbarBase(
        ax,
        cmap=plt.get_cmap("viridis"),
        norm=norm,
        orientation="horizontal",
        ticks=[0.0, 0.5, 1.0, 1.5, 2.0],
    )
    ax.set_title(
        "Color scale: log10(Expression score + 1)",
        fontsize=18,
        fontweight="bold",
        pad=12,
    )
    ax.tick_params(labelsize=14)
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    aligned_expr_path = Path(args.aligned_expr)
    pan_expr_path = Path(args.pan_expr)
    pan_map_path = Path(args.pan_map)
    out_dir = Path(args.out_dir)

    for path in [aligned_expr_path, pan_expr_path, pan_map_path]:
        if not path.exists():
            raise FileNotFoundError(path)

    aligned_df = load_aligned_expression(aligned_expr_path)
    pan_expr_df = load_processed_expression_cells(pan_expr_path)
    pan_map_df = normalize_map_df(pan_map_path)

    shared_meta, tissue_names, tissue_letters = build_shared_metadata(
        aligned_df,
        pan_expr_df,
        pan_map_df,
    )
    shared_gene_names = shared_meta["gene_name"].astype(str).tolist()
    gene_numbers = shared_meta["gene_number"].astype(int).tolist()

    aligned_matrix = build_matrix(
        aligned_df,
        gene_id_col="Gene ID",
        gene_name_col="merged_gene_name",
        shared_gene_names=shared_gene_names,
        tissue_names=tissue_names,
    )
    pan_named_df = pan_expr_df.merge(
        shared_meta[["gene_name", "pan_gene_id"]],
        left_on="Gene ID",
        right_on="pan_gene_id",
        how="left",
    )
    pan_matrix = build_matrix(
        pan_named_df,
        gene_id_col="Gene ID",
        gene_name_col="gene_name",
        shared_gene_names=shared_gene_names,
        tissue_names=tissue_names,
    )

    if aligned_matrix.shape != pan_matrix.shape:
        raise RuntimeError(
            f"Presentation matrices do not align: human={aligned_matrix.shape}, pan={pan_matrix.shape}"
        )

    gene_map, tissue_map = build_display_maps(shared_meta, tissue_names, tissue_letters)
    out_dir.mkdir(parents=True, exist_ok=True)
    save_map_tables(gene_map, tissue_map, out_dir)

    plot_heatmap(
        aligned_matrix,
        species_title="Homo sapiens",
        out_png=out_dir / "hs_aligned_all_presentation.png",
        out_svg=out_dir / "hs_aligned_all_presentation.svg",
        gene_numbers=gene_numbers,
        tissue_letters=tissue_letters,
    )
    plot_heatmap(
        pan_matrix,
        species_title="Pan troglodytes",
        out_png=out_dir / "h2a_pan_troglodytes_all_presentation.png",
        out_svg=out_dir / "h2a_pan_troglodytes_all_presentation.svg",
        gene_numbers=gene_numbers,
        tissue_letters=tissue_letters,
    )
    plot_gene_legend(
        gene_map,
        out_dir / "human_pan_gene_number_legend.png",
        out_dir / "human_pan_gene_number_legend.svg",
    )
    plot_tissue_legend(
        tissue_map,
        out_dir / "human_pan_tissue_letter_legend.png",
        out_dir / "human_pan_tissue_letter_legend.svg",
    )
    plot_color_legend(
        out_dir / "human_pan_color_scale_legend.png",
        out_dir / "human_pan_color_scale_legend.svg",
    )

    print(f"HUMAN_MATRIX_SHAPE={aligned_matrix.shape[0]}x{aligned_matrix.shape[1]}")
    print(f"PAN_MATRIX_SHAPE={pan_matrix.shape[0]}x{pan_matrix.shape[1]}")
    print(f"SHARED_GENES={len(shared_gene_names)}")
    print(f"SHARED_TISSUES={len(tissue_names)}")
    print(f"REMOVED_FROM_PAN={sorted(set(pan_expr_df['Anatomical entity name']) - set(tissue_names))}")
    print(f"OUT_DIR={out_dir}")


if __name__ == "__main__":
    main()
