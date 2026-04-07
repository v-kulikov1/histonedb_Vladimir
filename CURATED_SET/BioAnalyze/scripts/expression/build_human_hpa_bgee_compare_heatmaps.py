#!/usr/bin/env python
"""Build presentation-friendly human cH2A comparison heatmaps for HPA nTPM vs Bgee."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Sequence

import matplotlib as mpl

mpl.use("Agg")

from matplotlib.colorbar import ColorbarBase
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

BIOANALYZE_SCRIPTS_ROOT = Path(__file__).resolve().parents[1]
if str(BIOANALYZE_SCRIPTS_ROOT) not in sys.path:
    sys.path.insert(0, str(BIOANALYZE_SCRIPTS_ROOT))

from bioanalyze_paths import get_bioanalyze_data_root, get_bioanalyze_figures_root
from normalized_expression_common import load_processed_expression_cells


DEFAULT_HPA_TSV = get_bioanalyze_data_root() / "expression_nTPM" / "human" / "h2a_human_gene_ntpm_cells.tsv"
DEFAULT_BGEE_TSV = (
    get_bioanalyze_data_root() / "processed" / "homo_sapiens" / "Homo_sapiens_expr_advanced_H2A_present_gold.tsv"
)
DEFAULT_MAP_TSV = get_bioanalyze_data_root() / "processed" / "homo_sapiens" / "h2a_hs_canonical_variant_map.tsv"
DEFAULT_OUT_FIG_DIR = get_bioanalyze_figures_root() / "heatmaps" / "compare_nTPM_bgee"
DEFAULT_OUT_DATA_DIR = get_bioanalyze_data_root() / "processed" / "intersections" / "human_hpa_bgee"

GENE_ORDER = [
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
]
EXACT_TISSUES = [
    "adipose tissue",
    "adrenal gland",
    "amygdala",
    "cerebellum",
    "cerebral cortex",
    "colon",
    "endometrium",
    "esophagus",
    "fallopian tube",
    "hypothalamus",
    "kidney",
    "liver",
    "lung",
    "nucleus accumbens",
    "ovary",
    "pancreas",
    "pituitary gland",
    "putamen",
    "retina",
    "small intestine",
    "spinal cord",
    "spleen",
    "stomach",
    "substantia nigra",
    "testis",
    "thyroid gland",
    "urinary bladder",
    "vagina",
]
SAFE_TISSUE_SYNONYMS = {
    "breast": "mammary gland",
    "caudate": "caudate nucleus",
    "cervix": "uterine cervix",
    "heart muscle": "myocardium",
    "prostate": "prostate gland",
    "salivary gland": "saliva-secreting gland",
    "skeletal muscle": "skeletal muscle tissue",
}
EXCLUDED_HPA_TISSUES = ["hippocampus", "skin"]
UNDERLINE_COMBINING_MARK = "\u0332"
CLASS_COLORS = {
    "clustered": "#6ea6c8",
    "variant": "#d7b46a",
}
ZERO_OUTLINE_COLOR = "#d62728"
ZERO_OUTLINE_WIDTH = 0.75


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build presentation-friendly human cH2A comparison heatmaps with "
            "compact axes, separate legends, and source-native labels."
        )
    )
    parser.add_argument(
        "--hpa-tsv",
        default=str(DEFAULT_HPA_TSV),
        help="HPA-derived human cH2A gene x tissue nTPM cells TSV.",
    )
    parser.add_argument(
        "--bgee-tsv",
        default=str(DEFAULT_BGEE_TSV),
        help="Processed human Bgee H2A present-gold TSV.",
    )
    parser.add_argument(
        "--map-tsv",
        default=str(DEFAULT_MAP_TSV),
        help="Human canonical Bgee map TSV used as the canonical gene-name source.",
    )
    parser.add_argument(
        "--out-fig-dir",
        default=str(DEFAULT_OUT_FIG_DIR),
        help="Output directory for comparison heatmaps, legends, and colorbars.",
    )
    parser.add_argument(
        "--out-data-dir",
        default=str(DEFAULT_OUT_DATA_DIR),
        help="Output directory for aligned comparison tables and metadata.",
    )
    return parser.parse_args()


def clean_text_columns(df: pd.DataFrame, columns: Sequence[str]) -> pd.DataFrame:
    for column in columns:
        if column not in df.columns:
            df[column] = ""
        df[column] = df[column].fillna("").astype(str).str.strip()
    return df


def require_columns(df: pd.DataFrame, columns: Sequence[str], *, label: str) -> None:
    missing = [column for column in columns if column not in df.columns]
    if missing:
        raise RuntimeError(f"{label} is missing required columns: {', '.join(missing)}")


def underline_label(text: str) -> str:
    chars: List[str] = []
    for char in str(text):
        if char in {" ", "\n"}:
            chars.append(char)
        else:
            chars.append(f"{char}{UNDERLINE_COMBINING_MARK}")
    return "".join(chars)


def index_to_letter_code(index_zero_based: int) -> str:
    index = index_zero_based + 1
    code = ""
    while index > 0:
        index, remainder = divmod(index - 1, 26)
        code = chr(ord("A") + remainder) + code
    return code


def load_hpa_df(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path, sep="\t", low_memory=True)
    require_columns(df, ["gene_name", "tissue", "mean_nTPM"], label="HPA TSV")
    df = clean_text_columns(df, ["gene_name", "tissue"])
    df["mean_nTPM"] = pd.to_numeric(df["mean_nTPM"], errors="coerce")
    df["nTPM"] = df["mean_nTPM"]
    df = df[df["gene_name"].ne("") & df["tissue"].ne("")].copy()
    return df


def load_bgee_df(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    df = load_processed_expression_cells(path).copy()
    require_columns(
        df,
        ["Gene ID", "Gene name", "Anatomical entity name", "cell_mean_score"],
        label="Bgee TSV",
    )
    df = clean_text_columns(df, ["Gene ID", "Gene name", "Anatomical entity name"])
    df["cell_mean_score"] = pd.to_numeric(df["cell_mean_score"], errors="coerce")
    df = df[df["Gene ID"].ne("") & df["Anatomical entity name"].ne("")].copy()
    return df


def load_map_df(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path, sep="\t", dtype=str)
    require_columns(
        df,
        ["ensembl_gene_id", "gene_name", "class", "label"],
        label="Canonical map TSV",
    )
    df = clean_text_columns(df, ["ensembl_gene_id", "gene_name", "class", "label"])
    df = df[df["ensembl_gene_id"].ne("") & df["gene_name"].ne("")].copy()
    return df


def build_gene_mapping(map_df: pd.DataFrame, bgee_df: pd.DataFrame) -> pd.DataFrame:
    map_subset = (
        map_df[map_df["gene_name"].isin(GENE_ORDER)][["ensembl_gene_id", "gene_name", "class"]]
        .sort_values(["gene_name", "ensembl_gene_id"], ascending=[True, True], kind="stable")
        .drop_duplicates(subset=["gene_name"], keep="first")
        .copy()
    )
    missing_genes = [gene_name for gene_name in GENE_ORDER if gene_name not in set(map_subset["gene_name"])]
    if missing_genes:
        raise RuntimeError(
            "Canonical map is missing required cH2A genes: " + ", ".join(missing_genes)
        )

    bgee_label_df = (
        bgee_df[["Gene ID", "Gene name"]]
        .assign(_len=lambda x: x["Gene name"].str.len())
        .sort_values(["Gene ID", "_len", "Gene name"], ascending=[True, False, True], kind="stable")
        .drop_duplicates(subset=["Gene ID"], keep="first")
        .rename(columns={"Gene ID": "ensembl_gene_id", "Gene name": "bgee_gene_label"})
        .drop(columns=["_len"])
    )

    gene_mapping = map_subset.merge(bgee_label_df, on="ensembl_gene_id", how="left")
    gene_mapping["canonical_gene_name"] = gene_mapping["gene_name"]
    gene_mapping["hpa_gene_label"] = gene_mapping["canonical_gene_name"]
    gene_mapping["bgee_gene_label"] = gene_mapping["bgee_gene_label"].fillna("").astype(str).str.strip()
    empty_bgee = gene_mapping["bgee_gene_label"].eq("")
    gene_mapping.loc[empty_bgee, "bgee_gene_label"] = gene_mapping.loc[
        empty_bgee, "canonical_gene_name"
    ]
    gene_mapping["label_match_type"] = np.where(
        gene_mapping["bgee_gene_label"].eq(gene_mapping["canonical_gene_name"]),
        "exact",
        "alias",
    )
    gene_mapping["gene_order"] = gene_mapping["canonical_gene_name"].map(
        {gene_name: idx for idx, gene_name in enumerate(GENE_ORDER, start=1)}
    )
    gene_mapping["gene_number"] = gene_mapping["gene_order"].astype(int)
    gene_mapping["legend_label"] = np.where(
        gene_mapping["label_match_type"].eq("alias"),
        gene_mapping["canonical_gene_name"] + " [Bgee: " + gene_mapping["bgee_gene_label"] + "]",
        gene_mapping["canonical_gene_name"],
    )
    gene_mapping = gene_mapping.sort_values("gene_order", kind="stable").reset_index(drop=True)
    return gene_mapping[
        [
            "gene_order",
            "gene_number",
            "canonical_gene_name",
            "ensembl_gene_id",
            "class",
            "hpa_gene_label",
            "bgee_gene_label",
            "label_match_type",
            "legend_label",
        ]
    ].copy()


def build_tissue_mapping(hpa_df: pd.DataFrame, bgee_df: pd.DataFrame) -> pd.DataFrame:
    hpa_tissues = set(hpa_df["tissue"].astype(str))
    bgee_tissues = set(bgee_df["Anatomical entity name"].astype(str))

    rows: List[dict] = []
    for tissue in EXACT_TISSUES:
        if tissue not in hpa_tissues:
            raise RuntimeError(f"Exact HPA tissue is missing from the HPA input: {tissue}")
        if tissue not in bgee_tissues:
            raise RuntimeError(f"Exact Bgee tissue is missing from the Bgee input: {tissue}")
        rows.append(
            {
                "tissue_order": len(rows) + 1,
                "compare_tissue_key": tissue,
                "hpa_tissue_label": tissue,
                "bgee_tissue_label": tissue,
                "match_type": "exact",
                "is_underlined": False,
            }
        )

    for hpa_label, bgee_label in SAFE_TISSUE_SYNONYMS.items():
        if hpa_label not in hpa_tissues:
            raise RuntimeError(f"Safe-synonym HPA tissue is missing from the HPA input: {hpa_label}")
        if bgee_label not in bgee_tissues:
            raise RuntimeError(f"Safe-synonym Bgee tissue is missing from the Bgee input: {bgee_label}")
        rows.append(
            {
                "tissue_order": len(rows) + 1,
                "compare_tissue_key": bgee_label,
                "hpa_tissue_label": hpa_label,
                "bgee_tissue_label": bgee_label,
                "match_type": "safe_synonym",
                "is_underlined": True,
            }
        )

    tissue_mapping = pd.DataFrame(rows)
    tissue_mapping["tissue_letter"] = [
        index_to_letter_code(idx) for idx in range(len(tissue_mapping))
    ]
    remaining_hpa = sorted(hpa_tissues - set(tissue_mapping["hpa_tissue_label"]))
    expected_remaining = sorted(EXCLUDED_HPA_TISSUES)
    if remaining_hpa != expected_remaining:
        raise RuntimeError(
            "Unexpected HPA tissues remain outside the fixed compare mapping. "
            f"Expected {expected_remaining}, found {remaining_hpa}"
        )
    return tissue_mapping


def build_hpa_aligned_long(
    hpa_df: pd.DataFrame,
    gene_mapping: pd.DataFrame,
    tissue_mapping: pd.DataFrame,
) -> pd.DataFrame:
    observed = hpa_df.rename(
        columns={
            "gene_name": "canonical_gene_name",
            "tissue": "hpa_tissue_label",
        }
    ).copy()
    observed = observed.merge(
        tissue_mapping[
            [
                "compare_tissue_key",
                "tissue_letter",
                "hpa_tissue_label",
                "match_type",
                "is_underlined",
                "tissue_order",
                "bgee_tissue_label",
            ]
        ],
        on="hpa_tissue_label",
        how="inner",
    )
    observed = observed[observed["canonical_gene_name"].isin(GENE_ORDER)].copy()
    observed["nTPM"] = pd.to_numeric(observed["nTPM"], errors="coerce")
    observed = observed[
        [
            "canonical_gene_name",
            "compare_tissue_key",
            "tissue_letter",
            "hpa_tissue_label",
            "bgee_tissue_label",
            "match_type",
            "is_underlined",
            "nTPM",
        ]
    ].copy()

    grid = gene_mapping[
        [
            "gene_order",
            "gene_number",
            "canonical_gene_name",
            "ensembl_gene_id",
            "hpa_gene_label",
        ]
    ].merge(
        tissue_mapping[
            [
                "tissue_order",
                "tissue_letter",
                "compare_tissue_key",
                "hpa_tissue_label",
                "bgee_tissue_label",
                "match_type",
                "is_underlined",
            ]
        ],
        how="cross",
    )
    aligned = grid.merge(
        observed,
        on=[
            "canonical_gene_name",
            "compare_tissue_key",
            "tissue_letter",
            "hpa_tissue_label",
            "bgee_tissue_label",
            "match_type",
            "is_underlined",
        ],
        how="left",
    )
    aligned = aligned.sort_values(
        by=["gene_order", "tissue_order"],
        ascending=[True, True],
        kind="stable",
    ).reset_index(drop=True)
    return aligned[
        [
            "gene_order",
            "gene_number",
            "tissue_order",
            "tissue_letter",
            "canonical_gene_name",
            "ensembl_gene_id",
            "hpa_gene_label",
            "compare_tissue_key",
            "hpa_tissue_label",
            "bgee_tissue_label",
            "match_type",
            "is_underlined",
            "nTPM",
        ]
    ].copy()


def build_bgee_aligned_long(
    bgee_df: pd.DataFrame,
    gene_mapping: pd.DataFrame,
    tissue_mapping: pd.DataFrame,
) -> pd.DataFrame:
    observed = bgee_df.merge(
        gene_mapping[
            [
                "canonical_gene_name",
                "ensembl_gene_id",
                "bgee_gene_label",
            ]
        ],
        left_on="Gene ID",
        right_on="ensembl_gene_id",
        how="inner",
    )
    observed = observed.rename(
        columns={
            "Anatomical entity name": "bgee_tissue_label",
        }
    ).copy()
    observed = observed.merge(
        tissue_mapping[
            [
                "compare_tissue_key",
                "tissue_letter",
                "hpa_tissue_label",
                "bgee_tissue_label",
                "match_type",
                "is_underlined",
                "tissue_order",
            ]
        ],
        on="bgee_tissue_label",
        how="inner",
    )
    observed["cell_mean_score"] = pd.to_numeric(observed["cell_mean_score"], errors="coerce")
    observed = observed[
        [
            "canonical_gene_name",
            "ensembl_gene_id",
            "bgee_gene_label",
            "compare_tissue_key",
            "tissue_letter",
            "hpa_tissue_label",
            "bgee_tissue_label",
            "match_type",
            "is_underlined",
            "cell_mean_score",
        ]
    ].copy()

    grid = gene_mapping[
        [
            "gene_order",
            "gene_number",
            "canonical_gene_name",
            "ensembl_gene_id",
            "bgee_gene_label",
        ]
    ].merge(
        tissue_mapping[
            [
                "tissue_order",
                "tissue_letter",
                "compare_tissue_key",
                "hpa_tissue_label",
                "bgee_tissue_label",
                "match_type",
                "is_underlined",
            ]
        ],
        how="cross",
    )
    aligned = grid.merge(
        observed,
        on=[
            "canonical_gene_name",
            "ensembl_gene_id",
            "bgee_gene_label",
            "compare_tissue_key",
            "tissue_letter",
            "hpa_tissue_label",
            "bgee_tissue_label",
            "match_type",
            "is_underlined",
        ],
        how="left",
    )
    aligned = aligned.sort_values(
        by=["gene_order", "tissue_order"],
        ascending=[True, True],
        kind="stable",
    ).reset_index(drop=True)
    return aligned[
        [
            "gene_order",
            "gene_number",
            "tissue_order",
            "tissue_letter",
            "canonical_gene_name",
            "ensembl_gene_id",
            "bgee_gene_label",
            "compare_tissue_key",
            "hpa_tissue_label",
            "bgee_tissue_label",
            "match_type",
            "is_underlined",
            "cell_mean_score",
        ]
    ].copy()


def build_matrix(
    aligned_df: pd.DataFrame,
    *,
    value_col: str,
    gene_col: str = "canonical_gene_name",
    tissue_col: str = "compare_tissue_key",
) -> pd.DataFrame:
    matrix = aligned_df.pivot_table(
        index=gene_col,
        columns=tissue_col,
        values=value_col,
        aggfunc="mean",
    )
    matrix = matrix.reindex(index=GENE_ORDER, columns=EXACT_TISSUES + list(SAFE_TISSUE_SYNONYMS.values()))
    return matrix


def validate_matrix_shape(matrix: pd.DataFrame, *, name: str) -> None:
    expected_shape = (len(GENE_ORDER), len(EXACT_TISSUES) + len(SAFE_TISSUE_SYNONYMS))
    if matrix.shape != expected_shape:
        raise RuntimeError(f"{name} matrix has shape {matrix.shape}; expected {expected_shape}")
    empty_rows = matrix.index[matrix.notna().sum(axis=1).eq(0)].tolist()
    if empty_rows:
        raise RuntimeError(f"{name} matrix has empty gene rows: {empty_rows}")
    empty_cols = matrix.columns[matrix.notna().sum(axis=0).eq(0)].tolist()
    if empty_cols:
        raise RuntimeError(f"{name} matrix has empty tissue columns: {empty_cols}")


def build_gene_legend_map(gene_mapping: pd.DataFrame) -> pd.DataFrame:
    gene_map = gene_mapping.sort_values("gene_order", kind="stable").copy()
    gene_map["axis"] = "OY"
    return gene_map[
        [
            "axis",
            "gene_number",
            "canonical_gene_name",
            "hpa_gene_label",
            "bgee_gene_label",
            "legend_label",
            "class",
            "label_match_type",
            "ensembl_gene_id",
        ]
    ].copy()


def build_tissue_legend_map(
    tissue_mapping: pd.DataFrame,
    *,
    source_col: str,
    label_col_name: str,
) -> pd.DataFrame:
    tissue_map = tissue_mapping.sort_values("tissue_order", kind="stable").copy()
    tissue_map["axis"] = "OX"
    return tissue_map[
        [
            "axis",
            "tissue_letter",
            "compare_tissue_key",
            source_col,
            "match_type",
            "is_underlined",
        ]
    ].rename(columns={source_col: label_col_name}).copy()


def build_xticklabels(tissue_mapping: pd.DataFrame) -> List[str]:
    labels: List[str] = []
    for row in tissue_mapping.sort_values("tissue_order", kind="stable").to_dict(orient="records"):
        label = str(row["tissue_letter"])
        if bool(row["is_underlined"]):
            label = underline_label(label)
        labels.append(label)
    return labels


def build_yticklabels(gene_mapping: pd.DataFrame) -> List[str]:
    ordered = gene_mapping.sort_values("gene_order", kind="stable")
    return ordered["gene_number"].astype(int).astype(str).tolist()


def save_table(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def style_axis(ax: plt.Axes) -> None:
    for spine in ax.spines.values():
        spine.set_visible(False)


def plot_heatmap(
    matrix: pd.DataFrame,
    *,
    xticklabels: Sequence[str],
    yticklabels: Sequence[str],
    title: str,
    out_png: Path,
    out_svg: Path,
    legend_note: str,
    vmax: float,
    transform_mode: str = "log10p1",
) -> None:
    raw_matrix = matrix.astype(float).copy()
    if transform_mode == "log10p1":
        plot_df = np.log10(raw_matrix + 1.0)
    elif transform_mode == "none":
        plot_df = raw_matrix.copy()
    else:
        raise RuntimeError(f"Unsupported transform_mode: {transform_mode}")

    zero_mask = raw_matrix.notna() & raw_matrix.eq(0.0)

    fig_w = max(16.0, plot_df.shape[1] * 0.38 + 2.5)
    fig_h = max(9.5, plot_df.shape[0] * 0.52 + 1.6)

    sns.set_theme(style="white")
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    sns.heatmap(
        plot_df,
        cmap="viridis",
        mask=plot_df.isna(),
        linewidths=0.35,
        linecolor="#eef2f5",
        vmin=0.0,
        vmax=vmax,
        cbar=False,
        xticklabels=list(xticklabels),
        yticklabels=list(yticklabels),
        ax=ax,
    )
    for row_idx, col_idx in np.argwhere(zero_mask.to_numpy(dtype=bool)):
        ax.add_patch(
            Rectangle(
                (col_idx, row_idx),
                1.0,
                1.0,
                fill=False,
                edgecolor=ZERO_OUTLINE_COLOR,
                linewidth=ZERO_OUTLINE_WIDTH,
                joinstyle="miter",
                capstyle="butt",
            )
        )
    ax.set_title(title, fontsize=20, fontweight="bold", pad=12)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", labelrotation=0, labelsize=13, length=0, pad=6)
    ax.tick_params(axis="y", labelrotation=0, labelsize=13, length=0, pad=6)
    style_axis(ax)
    fig.text(
        0.995,
        0.012,
        legend_note,
        ha="right",
        va="bottom",
        fontsize=10,
        color="#333333",
    )
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=(0.01, 0.03, 0.995, 0.985))
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    plt.close(fig)


def plot_gene_legend(gene_map: pd.DataFrame, out_png: Path, out_svg: Path) -> None:
    fig_h = max(7.0, len(gene_map) * 0.48 + 1.25)
    fig, ax = plt.subplots(figsize=(11.8, fig_h))
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
    if gene_map["label_match_type"].eq("alias").any():
        ax.text(
            0.05,
            0.948,
            "Alias shown inline when the Bgee label differs.",
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=11,
            color="#4a4a4a",
        )

    top = 0.90
    bottom = 0.04
    step = (top - bottom) / max(len(gene_map) - 1, 1)
    x_number = 0.12
    x_swatch = 0.17
    swatch_size = 0.03
    x_gene = 0.235

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
            fontsize=17,
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
            str(row["legend_label"]),
            transform=ax.transAxes,
            va="center",
            ha="left",
            fontsize=16,
            color="#111111",
        )

    fig.tight_layout(rect=(0.02, 0.01, 0.995, 0.995))
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    plt.close(fig)


def plot_tissue_legend(
    tissue_map: pd.DataFrame,
    *,
    title: str,
    out_png: Path,
    out_svg: Path,
) -> None:
    entries = tissue_map.to_dict(orient="records")
    if not entries:
        raise RuntimeError("Tissue map is empty; cannot build tissue legend.")

    fig_h = max(10.5, len(entries) * 0.48 + 1.6)
    fig, ax = plt.subplots(figsize=(11.2, fig_h))
    ax.axis("off")

    ax.text(
        0.05,
        0.985,
        title,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=22,
        fontweight="bold",
        color="#111111",
    )
    ax.text(
        0.05,
        0.948,
        "Underlined labels = safe synonym mapping",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=11,
        color="#4a4a4a",
    )

    top = 0.91
    bottom = 0.04
    step = (top - bottom) / max(len(entries) - 1, 1)
    x_letter = 0.13
    x_name = 0.22
    label_col = [col for col in tissue_map.columns if col.endswith("_tissue_label")][0]

    for idx, row in enumerate(entries):
        y = top - idx * step
        label = str(row[label_col])
        if bool(row["is_underlined"]):
            label = underline_label(label)
        letter = str(row["tissue_letter"])
        if bool(row["is_underlined"]):
            letter = underline_label(letter)
        ax.text(
            x_letter,
            y,
            letter,
            transform=ax.transAxes,
            va="center",
            ha="right",
            fontsize=17,
            fontweight="bold",
            color="#111111",
        )
        ax.text(
            x_name,
            y,
            label,
            transform=ax.transAxes,
            va="center",
            ha="left",
            fontsize=16,
            color="#111111",
        )

    fig.tight_layout(rect=(0.02, 0.01, 0.995, 0.995))
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    plt.close(fig)


def plot_colorbar(
    *,
    title: str,
    vmax: float,
    out_png: Path,
    out_svg: Path,
) -> None:
    fig = plt.figure(figsize=(8.2, 2.25))
    ax = fig.add_axes([0.08, 0.40, 0.84, 0.24])
    ticks = np.linspace(0.0, vmax, num=5) if vmax > 0 else np.array([0.0])
    norm = mpl.colors.Normalize(vmin=0.0, vmax=vmax if vmax > 0 else 1.0)
    ColorbarBase(
        ax,
        cmap=plt.get_cmap("viridis"),
        norm=norm,
        orientation="horizontal",
        ticks=ticks,
    )
    ax.set_title(title, fontsize=16, fontweight="bold", pad=10)
    ax.tick_params(labelsize=12)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    plt.close(fig)


def build_metadata(
    *,
    hpa_tsv: Path,
    bgee_tsv: Path,
    map_tsv: Path,
    out_fig_dir: Path,
    out_data_dir: Path,
    gene_mapping: pd.DataFrame,
    tissue_mapping: pd.DataFrame,
    hpa_matrix: pd.DataFrame,
    bgee_matrix: pd.DataFrame,
) -> Dict[str, object]:
    alias_rows = gene_mapping[gene_mapping["label_match_type"].eq("alias")].copy()
    alias_notes = [
        {
            "canonical_gene_name": row["canonical_gene_name"],
            "ensembl_gene_id": row["ensembl_gene_id"],
            "bgee_gene_label": row["bgee_gene_label"],
        }
        for row in alias_rows.to_dict(orient="records")
    ]

    metadata: Dict[str, object] = {
        "hpa_tsv": hpa_tsv.as_posix(),
        "bgee_tsv": bgee_tsv.as_posix(),
        "map_tsv": map_tsv.as_posix(),
        "out_fig_dir": out_fig_dir.as_posix(),
        "out_data_dir": out_data_dir.as_posix(),
        "gene_order": list(GENE_ORDER),
        "gene_numbers": gene_mapping[["gene_number", "canonical_gene_name"]].to_dict(orient="records"),
        "exact_tissues": list(EXACT_TISSUES),
        "safe_tissue_synonyms": dict(SAFE_TISSUE_SYNONYMS),
        "excluded_hpa_tissues": list(EXCLUDED_HPA_TISSUES),
        "tissue_letters": tissue_mapping[
            [
                "tissue_letter",
                "compare_tissue_key",
                "hpa_tissue_label",
                "bgee_tissue_label",
                "match_type",
                "is_underlined",
            ]
        ].to_dict(orient="records"),
        "gene_count": int(len(gene_mapping)),
        "tissue_count": int(len(tissue_mapping)),
        "hpa_matrix_shape": [int(hpa_matrix.shape[0]), int(hpa_matrix.shape[1])],
        "bgee_matrix_shape": [int(bgee_matrix.shape[0]), int(bgee_matrix.shape[1])],
        "hpa_log10_vmax": float(np.nanmax(np.log10(hpa_matrix.astype(float) + 1.0))),
        "bgee_log10_vmax": float(np.nanmax(np.log10(bgee_matrix.astype(float) + 1.0))),
        "hpa_raw_vmax": float(np.nanmax(hpa_matrix.astype(float))),
        "bgee_raw_vmax": float(np.nanmax(bgee_matrix.astype(float))),
        "zero_outline_style": {
            "enabled": True,
            "description": "Exact zero-valued cells are marked with a thin red outline on all compare heatmaps.",
            "edgecolor": ZERO_OUTLINE_COLOR,
            "linewidth": ZERO_OUTLINE_WIDTH,
            "hpa_exact_zero_count": int((hpa_matrix.astype(float) == 0.0).sum().sum()),
            "bgee_exact_zero_count": int((bgee_matrix.astype(float) == 0.0).sum().sum()),
        },
        "underlined_tissues": tissue_mapping.loc[
            tissue_mapping["is_underlined"], "compare_tissue_key"
        ].astype(str).tolist(),
        "alias_notes": alias_notes,
        "mapping_provenance": {
            "canonical_gene_names": (
                "Canonical gene identity comes from the Bgee human canonical map "
                "and the merged-v4 label refresh pipeline used by the existing heatmap builders."
            ),
            "tissue_mapping": (
                "Safe tissue synonym mapping is an explicit local compare-layer. "
                "The existing repo-level tissue intersection workflows use exact matches only."
            ),
            "display_labels": (
                "Presentation heatmaps use compact OX/OY codes. Source-native HPA and "
                "Bgee labels are preserved in the exported legends and map TSVs."
            ),
            "hpa_compare_value": (
                "HPA compare values use nTPM from the gene x tissue cells table. "
                "In this dataset each gene x tissue cell is unique, so the exported value "
                "matches the raw per-cell nTPM."
            ),
            "bgee_compare_value": "Bgee compare values use cell_mean_score from the processed expression table.",
        },
    }
    return metadata


def main() -> None:
    args = parse_args()
    hpa_tsv = Path(args.hpa_tsv)
    bgee_tsv = Path(args.bgee_tsv)
    map_tsv = Path(args.map_tsv)
    out_fig_dir = Path(args.out_fig_dir)
    out_data_dir = Path(args.out_data_dir)

    hpa_df = load_hpa_df(hpa_tsv)
    bgee_df = load_bgee_df(bgee_tsv)
    map_df = load_map_df(map_tsv)

    gene_mapping = build_gene_mapping(map_df, bgee_df)
    tissue_mapping = build_tissue_mapping(hpa_df, bgee_df)
    hpa_aligned = build_hpa_aligned_long(hpa_df, gene_mapping, tissue_mapping)
    bgee_aligned = build_bgee_aligned_long(bgee_df, gene_mapping, tissue_mapping)

    hpa_matrix = build_matrix(hpa_aligned, value_col="nTPM")
    bgee_matrix = build_matrix(bgee_aligned, value_col="cell_mean_score")
    validate_matrix_shape(hpa_matrix, name="HPA")
    validate_matrix_shape(bgee_matrix, name="Bgee")

    gene_legend_map = build_gene_legend_map(gene_mapping)
    hpa_tissue_legend_map = build_tissue_legend_map(
        tissue_mapping,
        source_col="hpa_tissue_label",
        label_col_name="hpa_tissue_label",
    )
    bgee_tissue_legend_map = build_tissue_legend_map(
        tissue_mapping,
        source_col="bgee_tissue_label",
        label_col_name="bgee_tissue_label",
    )

    compact_xticklabels = build_xticklabels(tissue_mapping)
    compact_yticklabels = build_yticklabels(gene_mapping)

    hpa_vmax = float(np.nanmax(np.log10(hpa_matrix.astype(float) + 1.0)))
    bgee_vmax = float(np.nanmax(np.log10(bgee_matrix.astype(float) + 1.0)))
    hpa_raw_vmax = float(np.nanmax(hpa_matrix.astype(float)))
    bgee_raw_vmax = float(np.nanmax(bgee_matrix.astype(float)))

    out_fig_dir.mkdir(parents=True, exist_ok=True)
    out_data_dir.mkdir(parents=True, exist_ok=True)

    save_table(
        gene_mapping[
            [
                "canonical_gene_name",
                "ensembl_gene_id",
                "hpa_gene_label",
                "bgee_gene_label",
                "label_match_type",
            ]
        ],
        out_data_dir / "hpa_vs_bgee_gene_mapping.tsv",
    )
    save_table(
        tissue_mapping[
            [
                "compare_tissue_key",
                "tissue_letter",
                "hpa_tissue_label",
                "bgee_tissue_label",
                "match_type",
                "is_underlined",
            ]
        ],
        out_data_dir / "hpa_vs_bgee_tissue_mapping.tsv",
    )
    save_table(hpa_aligned, out_data_dir / "hpa_vs_bgee_hpa_aligned_long.tsv")
    save_table(bgee_aligned, out_data_dir / "hpa_vs_bgee_bgee_aligned_long.tsv")
    save_table(gene_legend_map, out_fig_dir / "hpa_vs_bgee_gene_number_map.tsv")
    save_table(hpa_tissue_legend_map, out_fig_dir / "hpa_vs_bgee_hpa_tissue_letter_map.tsv")
    save_table(bgee_tissue_legend_map, out_fig_dir / "hpa_vs_bgee_bgee_tissue_letter_map.tsv")

    metadata = build_metadata(
        hpa_tsv=hpa_tsv,
        bgee_tsv=bgee_tsv,
        map_tsv=map_tsv,
        out_fig_dir=out_fig_dir,
        out_data_dir=out_data_dir,
        gene_mapping=gene_mapping,
        tissue_mapping=tissue_mapping,
        hpa_matrix=hpa_matrix,
        bgee_matrix=bgee_matrix,
    )
    (out_data_dir / "hpa_vs_bgee_metadata.json").write_text(
        json.dumps(metadata, indent=2),
        encoding="utf-8",
    )

    note = "Underlined OX codes = safe synonym mapping; red-outlined cells = exact zero"
    plot_heatmap(
        hpa_matrix,
        xticklabels=compact_xticklabels,
        yticklabels=compact_yticklabels,
        title="Human cH2A HPA (nTPM)",
        out_png=out_fig_dir / "hpa_vs_bgee_hpa_heatmap.png",
        out_svg=out_fig_dir / "hpa_vs_bgee_hpa_heatmap.svg",
        legend_note=note,
        vmax=hpa_vmax,
        transform_mode="log10p1",
    )
    plot_heatmap(
        bgee_matrix,
        xticklabels=compact_xticklabels,
        yticklabels=compact_yticklabels,
        title="Human cH2A Bgee (present_gold)",
        out_png=out_fig_dir / "hpa_vs_bgee_bgee_heatmap.png",
        out_svg=out_fig_dir / "hpa_vs_bgee_bgee_heatmap.svg",
        legend_note=note,
        vmax=bgee_vmax,
        transform_mode="log10p1",
    )
    plot_heatmap(
        hpa_matrix,
        xticklabels=compact_xticklabels,
        yticklabels=compact_yticklabels,
        title="Human cH2A HPA (nTPM, raw scale)",
        out_png=out_fig_dir / "hpa_vs_bgee_hpa_raw_heatmap.png",
        out_svg=out_fig_dir / "hpa_vs_bgee_hpa_raw_heatmap.svg",
        legend_note=note,
        vmax=hpa_raw_vmax,
        transform_mode="none",
    )
    plot_heatmap(
        bgee_matrix,
        xticklabels=compact_xticklabels,
        yticklabels=compact_yticklabels,
        title="Human cH2A Bgee (present_gold, raw scale)",
        out_png=out_fig_dir / "hpa_vs_bgee_bgee_raw_heatmap.png",
        out_svg=out_fig_dir / "hpa_vs_bgee_bgee_raw_heatmap.svg",
        legend_note=note,
        vmax=bgee_raw_vmax,
        transform_mode="none",
    )
    plot_gene_legend(
        gene_legend_map,
        out_png=out_fig_dir / "hpa_vs_bgee_gene_number_legend.png",
        out_svg=out_fig_dir / "hpa_vs_bgee_gene_number_legend.svg",
    )
    plot_tissue_legend(
        hpa_tissue_legend_map,
        title="OX: HPA tissue index",
        out_png=out_fig_dir / "hpa_vs_bgee_hpa_tissue_letter_legend.png",
        out_svg=out_fig_dir / "hpa_vs_bgee_hpa_tissue_letter_legend.svg",
    )
    plot_tissue_legend(
        bgee_tissue_legend_map,
        title="OX: Bgee tissue index",
        out_png=out_fig_dir / "hpa_vs_bgee_bgee_tissue_letter_legend.png",
        out_svg=out_fig_dir / "hpa_vs_bgee_bgee_tissue_letter_legend.svg",
    )
    plot_colorbar(
        title="HPA color scale: log10(nTPM + 1)",
        vmax=hpa_vmax,
        out_png=out_fig_dir / "hpa_vs_bgee_hpa_colorbar.png",
        out_svg=out_fig_dir / "hpa_vs_bgee_hpa_colorbar.svg",
    )
    plot_colorbar(
        title="Bgee color scale: log10(Expression score + 1)",
        vmax=bgee_vmax,
        out_png=out_fig_dir / "hpa_vs_bgee_bgee_colorbar.png",
        out_svg=out_fig_dir / "hpa_vs_bgee_bgee_colorbar.svg",
    )
    plot_colorbar(
        title="HPA color scale: nTPM",
        vmax=hpa_raw_vmax,
        out_png=out_fig_dir / "hpa_vs_bgee_hpa_raw_colorbar.png",
        out_svg=out_fig_dir / "hpa_vs_bgee_hpa_raw_colorbar.svg",
    )
    plot_colorbar(
        title="Bgee color scale: Expression score",
        vmax=bgee_raw_vmax,
        out_png=out_fig_dir / "hpa_vs_bgee_bgee_raw_colorbar.png",
        out_svg=out_fig_dir / "hpa_vs_bgee_bgee_raw_colorbar.svg",
    )

    print(f"HPA_MATRIX_SHAPE={hpa_matrix.shape[0]}x{hpa_matrix.shape[1]}")
    print(f"BGEE_MATRIX_SHAPE={bgee_matrix.shape[0]}x{bgee_matrix.shape[1]}")
    print(f"SAFE_UNDERLINED_TISSUES={int(tissue_mapping['is_underlined'].sum())}")
    print(f"GENE_ALIAS_COUNT={int(gene_mapping['label_match_type'].eq('alias').sum())}")
    print(f"PRESENTATION_TISSUE_CODES={len(tissue_mapping)}")


if __name__ == "__main__":
    main()
