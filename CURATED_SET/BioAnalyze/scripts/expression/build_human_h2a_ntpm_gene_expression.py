#!/usr/bin/env python
"""Build a gene-based human cH2A raw-nTPM workflow from HPA test.tsv."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[4]
GITHUB_ROOT = SCRIPT_PATH.parents[5]
DEFAULT_INPUT_TSV = (
    GITHUB_ROOT / "histonedb_external_storage" / "BioAnalyze" / "raw" / "expression_nTPM" / "human" / "test.tsv"
)
DEFAULT_OUT_DATA_DIR = REPO_ROOT / "CURATED_SET" / "BioAnalyze" / "data" / "expression_nTPM" / "human"
DEFAULT_OUT_FIG_DIR = REPO_ROOT / "CURATED_SET" / "BioAnalyze" / "figures" / "expression_nTPM" / "human"

REQUIRED_COLUMNS = ["Gene name", "Tissue", "nTPM", "variant"]
TRACEABILITY_INPUT_COLUMNS = [
    "Unnamed: 0",
    "Gene",
    "accession",
    "protein_accession",
]
HIGH_THRESHOLD = 10.0
MODERATE_THRESHOLD = 1.0
SUMMARY_CLASS_COLORS = {
    "high": "#bf4e30",
    "moderate": "#c9a227",
    "low": "#5a9e6f",
}
H2A5_SPLIT_GENES = ["H2AC11", "H2AC13", "H2AC15", "H2AC16", "H2AC17"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Rebuild the human cH2A HPA raw-nTPM workflow from test.tsv using "
            "Gene name (H2AC*) instead of legacy cH2A variant labels."
        )
    )
    parser.add_argument(
        "--input-tsv",
        default=str(DEFAULT_INPUT_TSV),
        help="Raw human HPA cH2A test.tsv source file.",
    )
    parser.add_argument(
        "--out-data-dir",
        default=str(DEFAULT_OUT_DATA_DIR),
        help="Output directory for derived TSV/JSON tables.",
    )
    parser.add_argument(
        "--out-fig-dir",
        default=str(DEFAULT_OUT_FIG_DIR),
        help="Output directory for figures.",
    )
    return parser.parse_args()


def require_columns(df: pd.DataFrame, columns: Sequence[str]) -> None:
    missing = [column for column in columns if column not in df.columns]
    if missing:
        raise RuntimeError(f"Missing required columns: {', '.join(missing)}")


def safe_slug(text: str) -> str:
    slug = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in text.strip())
    slug = "_".join(part for part in slug.split("_") if part)
    return slug or "item"


def clean_text_columns(df: pd.DataFrame, columns: Sequence[str]) -> pd.DataFrame:
    for column in columns:
        if column not in df.columns:
            df[column] = ""
        df[column] = df[column].fillna("").astype(str).str.strip()
    return df


def join_unique(values: Iterable[object]) -> str:
    items = sorted({str(value).strip() for value in values if str(value).strip()})
    return ";".join(items)


def first_nonempty(values: Iterable[object]) -> str:
    for value in values:
        text = str(value).strip()
        if text:
            return text
    return ""


def classify_gene(mean_ntpm: float) -> str:
    if mean_ntpm > HIGH_THRESHOLD:
        return "high"
    if mean_ntpm > MODERATE_THRESHOLD:
        return "moderate"
    return "low"


def load_source_df(input_tsv: Path) -> pd.DataFrame:
    if not input_tsv.exists():
        raise FileNotFoundError(input_tsv)

    source_df = pd.read_csv(input_tsv, sep="\t", low_memory=True)
    require_columns(source_df, REQUIRED_COLUMNS)
    source_df = clean_text_columns(
        source_df,
        ["Gene name", "Tissue", "variant", *TRACEABILITY_INPUT_COLUMNS],
    )
    source_df["nTPM"] = pd.to_numeric(source_df["nTPM"], errors="coerce").fillna(0.0).astype(float)
    source_df = source_df[
        source_df["Gene name"].ne("")
        & source_df["Tissue"].ne("")
        & source_df["variant"].ne("")
    ].copy()
    if source_df.empty:
        raise RuntimeError("No valid rows remain after cleaning the input TSV.")

    source_df["source_row_id"] = (
        source_df["Unnamed: 0"]
        .where(source_df["Unnamed: 0"].ne(""), source_df.index.astype(str))
        .astype(str)
    )
    return source_df


def build_gene_variant_map(source_df: pd.DataFrame) -> pd.DataFrame:
    mapping_df = (
        source_df.groupby("Gene name", as_index=False)
        .agg(
            legacy_variant=("variant", join_unique),
            source_gene_ids=("Gene", join_unique),
            source_accessions=("accession", join_unique),
            source_protein_accessions=("protein_accession", join_unique),
        )
        .rename(columns={"Gene name": "gene_name"})
        .sort_values("gene_name", ascending=True, kind="stable")
        .reset_index(drop=True)
    )
    return mapping_df


def build_cells_df(source_df: pd.DataFrame) -> pd.DataFrame:
    grouped = (
        source_df.groupby(["Gene name", "Tissue"], as_index=False)
        .agg(
            legacy_variant=("variant", join_unique),
            mean_nTPM=("nTPM", "mean"),
            median_nTPM=("nTPM", "median"),
            max_nTPM=("nTPM", "max"),
            min_nTPM=("nTPM", "min"),
            n_obs=("nTPM", "size"),
            source_row_ids=("source_row_id", join_unique),
            source_gene_ids=("Gene", join_unique),
            source_accessions=("accession", join_unique),
            source_protein_accessions=("protein_accession", join_unique),
        )
        .rename(columns={"Gene name": "gene_name", "Tissue": "tissue"})
    )
    numeric_cols = ["mean_nTPM", "median_nTPM", "max_nTPM", "min_nTPM"]
    for column in numeric_cols:
        grouped[column] = pd.to_numeric(grouped[column], errors="coerce").fillna(0.0).astype(float)
    grouped["n_obs"] = pd.to_numeric(grouped["n_obs"], errors="coerce").fillna(0).astype(int)
    grouped = grouped.sort_values(
        by=["gene_name", "tissue"],
        ascending=[True, True],
        kind="stable",
    ).reset_index(drop=True)
    return grouped


def build_nonzero_cells_df(source_df: pd.DataFrame) -> pd.DataFrame:
    nonzero_source_df = source_df[source_df["nTPM"].gt(0)].copy()
    if nonzero_source_df.empty:
        raise RuntimeError("No nTPM > 0 rows remain after applying the notebook filter.")
    return build_cells_df(nonzero_source_df)


def build_summary_df(
    all_cells_df: pd.DataFrame,
    nonzero_cells_df: pd.DataFrame,
    gene_variant_df: pd.DataFrame,
) -> pd.DataFrame:
    all_genes_df = gene_variant_df[["gene_name", "legacy_variant"]].copy()
    summary_nonzero_df = (
        nonzero_cells_df.groupby("gene_name", as_index=False)
        .agg(
            mean_nTPM=("median_nTPM", "mean"),
            median_nTPM=("median_nTPM", "median"),
            max_nTPM=("median_nTPM", "max"),
            expressed_tissue_count=("tissue", "nunique"),
        )
        .sort_values(["gene_name"], ascending=[True], kind="stable")
    )
    top_tissues_df = (
        nonzero_cells_df.sort_values(
            by=["gene_name", "median_nTPM", "mean_nTPM", "tissue"],
            ascending=[True, False, False, True],
            kind="stable",
        )
        .drop_duplicates(subset=["gene_name"], keep="first")
        .loc[:, ["gene_name", "tissue", "median_nTPM", "mean_nTPM"]]
        .rename(
            columns={
                "tissue": "top_tissue",
                "median_nTPM": "top_tissue_median_nTPM",
                "mean_nTPM": "top_tissue_mean_nTPM",
            }
        )
    )
    summary_df = all_genes_df.merge(summary_nonzero_df, on="gene_name", how="left")
    summary_df = summary_df.merge(top_tissues_df, on="gene_name", how="left")

    all_gene_totals_df = (
        all_cells_df.groupby("gene_name", as_index=False)
        .agg(
            source_tissue_count=("tissue", "nunique"),
            source_mean_nTPM=("median_nTPM", "mean"),
            source_max_nTPM=("median_nTPM", "max"),
        )
        .sort_values("gene_name", ascending=True, kind="stable")
    )
    summary_df = summary_df.merge(all_gene_totals_df, on="gene_name", how="left")

    for column in [
        "mean_nTPM",
        "median_nTPM",
        "max_nTPM",
        "top_tissue_median_nTPM",
        "top_tissue_mean_nTPM",
        "source_mean_nTPM",
        "source_max_nTPM",
    ]:
        summary_df[column] = pd.to_numeric(summary_df[column], errors="coerce").fillna(0.0).astype(float)
    for column in ["expressed_tissue_count", "source_tissue_count"]:
        summary_df[column] = pd.to_numeric(summary_df[column], errors="coerce").fillna(0).astype(int)
    summary_df["top_tissue"] = summary_df["top_tissue"].fillna("").astype(str).str.strip()
    summary_df["expression_class"] = summary_df["mean_nTPM"].apply(classify_gene)
    summary_df = summary_df.sort_values(
        by=["mean_nTPM", "median_nTPM", "gene_name"],
        ascending=[False, False, True],
        kind="stable",
    ).reset_index(drop=True)
    return summary_df


def build_tissue_counts_df(nonzero_cells_df: pd.DataFrame) -> pd.DataFrame:
    tissue_counts_df = (
        nonzero_cells_df.groupby("tissue", as_index=False)
        .agg(
            expressed_gene_count=("gene_name", "nunique"),
            genes=("gene_name", join_unique),
        )
        .sort_values(by=["expressed_gene_count", "tissue"], ascending=[False, True], kind="stable")
        .reset_index(drop=True)
    )
    return tissue_counts_df


def build_legacy_variant_audit_df(source_df: pd.DataFrame, nonzero_cells_df: pd.DataFrame) -> pd.DataFrame:
    all_variant_df = (
        source_df.groupby("variant", as_index=False)
        .agg(
            all_genes=("Gene name", join_unique),
            source_gene_count=("Gene name", "nunique"),
            source_row_count=("Gene name", "size"),
            source_tissue_count=("Tissue", "nunique"),
        )
        .rename(columns={"variant": "legacy_variant"})
    )
    nonzero_variant_df = (
        nonzero_cells_df.groupby("legacy_variant", as_index=False)
        .agg(
            expressed_genes=("gene_name", join_unique),
            expressed_gene_count=("gene_name", "nunique"),
            nonzero_tissue_count=("tissue", "nunique"),
        )
    )
    audit_df = all_variant_df.merge(nonzero_variant_df, on="legacy_variant", how="left")
    audit_df["expressed_genes"] = audit_df["expressed_genes"].fillna("").astype(str).str.strip()
    for column in ["expressed_gene_count", "nonzero_tissue_count"]:
        audit_df[column] = pd.to_numeric(audit_df[column], errors="coerce").fillna(0).astype(int)
    audit_df = audit_df.sort_values("legacy_variant", ascending=True, kind="stable").reset_index(drop=True)
    return audit_df


def ensure_parent_dirs(paths: Sequence[Path]) -> None:
    for path in paths:
        path.mkdir(parents=True, exist_ok=True)


def save_table(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def style_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", alpha=0.2, linewidth=0.8)


def save_figure(fig: plt.Figure, png_path: Path, svg_path: Path) -> None:
    png_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)


def plot_per_gene_barplots(nonzero_cells_df: pd.DataFrame, out_dir: Path) -> List[str]:
    generated: List[str] = []
    for gene_name in sorted(nonzero_cells_df["gene_name"].unique().tolist()):
        gene_df = nonzero_cells_df[nonzero_cells_df["gene_name"].eq(gene_name)].copy()
        gene_df = gene_df.sort_values(
            by=["median_nTPM", "mean_nTPM", "tissue"],
            ascending=[False, False, True],
            kind="stable",
        ).reset_index(drop=True)
        fig_w = max(10, gene_df.shape[0] * 0.35)
        fig, ax = plt.subplots(figsize=(fig_w, 5.5))
        ax.bar(gene_df["tissue"], gene_df["median_nTPM"], color="#3f7f93")
        ax.set_title(f"{gene_name} nTPM across tissues")
        ax.set_xlabel("Tissue")
        ax.set_ylabel("median nTPM")
        ax.tick_params(axis="x", rotation=75, labelsize=8)
        style_axis(ax)
        slug = safe_slug(gene_name)
        png_path = out_dir / f"{slug}_tissue_ntpm_barplot.png"
        svg_path = out_dir / f"{slug}_tissue_ntpm_barplot.svg"
        save_figure(fig, png_path, svg_path)
        generated.append(gene_name)
    return generated


def plot_per_tissue_barplots(nonzero_cells_df: pd.DataFrame, out_dir: Path) -> List[str]:
    generated: List[str] = []
    for tissue in sorted(nonzero_cells_df["tissue"].unique().tolist()):
        tissue_df = nonzero_cells_df[nonzero_cells_df["tissue"].eq(tissue)].copy()
        tissue_df = tissue_df.sort_values(
            by=["median_nTPM", "mean_nTPM", "gene_name"],
            ascending=[False, False, True],
            kind="stable",
        ).reset_index(drop=True)
        reference_median = float(tissue_df["median_nTPM"].median()) if not tissue_df.empty else 0.0
        fig_w = max(8, tissue_df.shape[0] * 0.55)
        fig, ax = plt.subplots(figsize=(fig_w, 5.2))
        ax.bar(tissue_df["gene_name"], tissue_df["median_nTPM"], color="#7a8f4c")
        ax.axhline(reference_median, linestyle="--", color="black", linewidth=1.1)
        ax.set_title(f"{tissue} cH2A gene ranking")
        ax.set_xlabel("Gene name")
        ax.set_ylabel("median nTPM")
        ax.tick_params(axis="x", rotation=65, labelsize=8)
        style_axis(ax)
        slug = safe_slug(tissue)
        png_path = out_dir / f"{slug}_gene_ntpm_ranked_barplot.png"
        svg_path = out_dir / f"{slug}_gene_ntpm_ranked_barplot.svg"
        save_figure(fig, png_path, svg_path)
        generated.append(tissue)
    return generated


def plot_mean_ranking(summary_df: pd.DataFrame, out_dir: Path) -> None:
    plot_df = summary_df[summary_df["expressed_tissue_count"].gt(0)].copy()
    plot_df = plot_df.sort_values(
        by=["mean_nTPM", "median_nTPM", "gene_name"],
        ascending=[True, True, False],
        kind="stable",
    ).reset_index(drop=True)
    colors = [SUMMARY_CLASS_COLORS.get(value, "#888888") for value in plot_df["expression_class"]]
    fig_h = max(5, plot_df.shape[0] * 0.48)
    fig, ax = plt.subplots(figsize=(8.5, fig_h))
    bars = ax.barh(plot_df["gene_name"], plot_df["mean_nTPM"], color=colors)
    ax.bar_label(bars, fmt="%.2f", fontsize=8, padding=3)
    ax.set_title("Human cH2A genes ranked by mean nonzero nTPM")
    ax.set_xlabel("mean nTPM")
    ax.set_ylabel("Gene name")
    style_axis(ax)
    png_path = out_dir / "h2a_human_gene_ntpm_mean_ranking.png"
    svg_path = out_dir / "h2a_human_gene_ntpm_mean_ranking.svg"
    save_figure(fig, png_path, svg_path)


def plot_tissue_counts(tissue_counts_df: pd.DataFrame, out_dir: Path) -> None:
    plot_df = tissue_counts_df.sort_values(
        by=["expressed_gene_count", "tissue"],
        ascending=[False, True],
        kind="stable",
    ).reset_index(drop=True)
    fig_w = max(11, plot_df.shape[0] * 0.35)
    fig, ax = plt.subplots(figsize=(fig_w, 5.5))
    bars = ax.bar(plot_df["tissue"], plot_df["expressed_gene_count"], color="#8363a8")
    ax.bar_label(bars, fontsize=7, padding=2)
    ax.set_title("Number of expressed cH2A genes by tissue")
    ax.set_xlabel("Tissue")
    ax.set_ylabel("Expressed gene count")
    ax.tick_params(axis="x", rotation=75, labelsize=8)
    style_axis(ax)
    png_path = out_dir / "h2a_human_gene_count_by_tissue.png"
    svg_path = out_dir / "h2a_human_gene_count_by_tissue.svg"
    save_figure(fig, png_path, svg_path)


def plot_h2a5_split_boxplot(nonzero_cells_df: pd.DataFrame, out_dir: Path) -> None:
    plot_df = nonzero_cells_df[nonzero_cells_df["gene_name"].isin(H2A5_SPLIT_GENES)].copy()
    ordered_genes = [gene for gene in H2A5_SPLIT_GENES if gene in set(plot_df["gene_name"])]
    if not ordered_genes:
        raise RuntimeError("No nonzero rows were found for the cH2A.5 split genes.")

    values = [
        plot_df.loc[plot_df["gene_name"].eq(gene_name), "median_nTPM"].astype(float).to_numpy()
        for gene_name in ordered_genes
    ]
    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    box = ax.boxplot(
        values,
        vert=False,
        patch_artist=True,
        tick_labels=ordered_genes,
    )
    palette = ["#6a8caf", "#8ab17d", "#d4a373", "#b56576", "#7b6d8d"]
    for patch, color in zip(box["boxes"], palette, strict=False):
        patch.set_facecolor(color)
        patch.set_alpha(0.9)
    ax.set_title("cH2A.5 split genes across tissues")
    ax.set_xlabel("nonzero median nTPM")
    ax.set_ylabel("Gene name")
    style_axis(ax)
    png_path = out_dir / "h2a5_split_gene_boxplot.png"
    svg_path = out_dir / "h2a5_split_gene_boxplot.svg"
    save_figure(fig, png_path, svg_path)


def build_metadata(
    *,
    input_tsv: Path,
    out_data_dir: Path,
    out_fig_dir: Path,
    source_df: pd.DataFrame,
    all_cells_df: pd.DataFrame,
    nonzero_cells_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    tissue_counts_df: pd.DataFrame,
    audit_df: pd.DataFrame,
) -> Dict[str, object]:
    expressed_genes = (
        summary_df.loc[summary_df["expressed_tissue_count"].gt(0), "gene_name"].astype(str).tolist()
    )
    zero_only_genes = (
        summary_df.loc[summary_df["expressed_tissue_count"].eq(0), "gene_name"].astype(str).tolist()
    )
    class_counts = {
        key: int(value)
        for key, value in summary_df["expression_class"].value_counts().sort_index().items()
    }
    top_tissue_row = tissue_counts_df.iloc[0].to_dict() if not tissue_counts_df.empty else {}
    metadata: Dict[str, object] = {
        "input_tsv": input_tsv.as_posix(),
        "out_data_dir": out_data_dir.as_posix(),
        "out_fig_dir": out_fig_dir.as_posix(),
        "required_columns": REQUIRED_COLUMNS,
        "filter_rule": "nTPM > 0",
        "thresholds": {
            "high_gt": HIGH_THRESHOLD,
            "moderate_gt": MODERATE_THRESHOLD,
            "low_lte": MODERATE_THRESHOLD,
        },
        "source_row_count": int(source_df.shape[0]),
        "source_gene_count": int(source_df["Gene name"].nunique()),
        "source_tissue_count": int(source_df["Tissue"].nunique()),
        "aggregated_cell_count": int(all_cells_df.shape[0]),
        "nonzero_cell_count": int(nonzero_cells_df.shape[0]),
        "expressed_gene_count": int(len(expressed_genes)),
        "nonzero_tissue_count": int(nonzero_cells_df["tissue"].nunique()),
        "expressed_genes": expressed_genes,
        "zero_only_genes": zero_only_genes,
        "class_counts": class_counts,
        "legacy_variant_gene_map": {
            row["legacy_variant"]: row["all_genes"].split(";") if row["all_genes"] else []
            for row in audit_df.to_dict(orient="records")
        },
        "top_tissue_by_expressed_gene_count": {
            "tissue": str(top_tissue_row.get("tissue", "")),
            "expressed_gene_count": int(top_tissue_row.get("expressed_gene_count", 0) or 0),
        },
        "notes": [
            "This workflow is intentionally separate from normalized Expression score outputs.",
            "Gene name is the primary analysis key; legacy variant labels are kept for traceability.",
            "Zero-only genes are preserved in raw summary/audit tables and excluded from expressed-gene plots.",
        ],
    }
    return metadata


def run_workflow(input_tsv: Path, out_data_dir: Path, out_fig_dir: Path) -> Dict[str, object]:
    source_df = load_source_df(input_tsv)
    gene_variant_df = build_gene_variant_map(source_df)
    all_cells_df = build_cells_df(source_df)
    nonzero_cells_df = build_nonzero_cells_df(source_df)
    summary_df = build_summary_df(all_cells_df, nonzero_cells_df, gene_variant_df)
    tissue_counts_df = build_tissue_counts_df(nonzero_cells_df)
    audit_df = build_legacy_variant_audit_df(source_df, nonzero_cells_df)

    ensure_parent_dirs(
        [
            out_data_dir,
            out_fig_dir / "per_gene",
            out_fig_dir / "per_tissue",
            out_fig_dir / "summary",
            out_fig_dir / "special_cases",
        ]
    )

    save_table(all_cells_df, out_data_dir / "h2a_human_gene_ntpm_cells.tsv")
    save_table(nonzero_cells_df, out_data_dir / "h2a_human_gene_ntpm_nonzero.tsv")
    save_table(summary_df, out_data_dir / "h2a_human_gene_ntpm_summary.tsv")
    save_table(tissue_counts_df, out_data_dir / "h2a_human_gene_ntpm_tissue_counts.tsv")
    save_table(audit_df, out_data_dir / "h2a_human_gene_ntpm_legacy_variant_audit.tsv")

    plot_per_gene_barplots(nonzero_cells_df, out_fig_dir / "per_gene")
    plot_per_tissue_barplots(nonzero_cells_df, out_fig_dir / "per_tissue")
    plot_mean_ranking(summary_df, out_fig_dir / "summary")
    plot_tissue_counts(tissue_counts_df, out_fig_dir / "summary")
    plot_h2a5_split_boxplot(nonzero_cells_df, out_fig_dir / "special_cases")

    metadata = build_metadata(
        input_tsv=input_tsv,
        out_data_dir=out_data_dir,
        out_fig_dir=out_fig_dir,
        source_df=source_df,
        all_cells_df=all_cells_df,
        nonzero_cells_df=nonzero_cells_df,
        summary_df=summary_df,
        tissue_counts_df=tissue_counts_df,
        audit_df=audit_df,
    )
    metadata_path = out_data_dir / "h2a_human_gene_ntpm_metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return metadata


def main() -> None:
    args = parse_args()
    input_tsv = Path(args.input_tsv).expanduser().resolve()
    out_data_dir = Path(args.out_data_dir).expanduser().resolve()
    out_fig_dir = Path(args.out_fig_dir).expanduser().resolve()
    metadata = run_workflow(input_tsv, out_data_dir, out_fig_dir)
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
