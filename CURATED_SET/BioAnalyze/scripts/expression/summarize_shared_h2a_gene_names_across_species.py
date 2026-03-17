#!/usr/bin/env python
"""Summarize shared H2A gene names across species and build a reusable index."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Iterable, List

import matplotlib.pyplot as plt
import pandas as pd


DEFAULT_HEATMAP_DIR = Path(r"CURATED_SET/BioAnalyze/figures/heatmaps")
DEFAULT_PROCESSED_DIR = Path(r"CURATED_SET/BioAnalyze/data/processed")
DEFAULT_STATS_DIR = Path(r"CURATED_SET/BioAnalyze/stats")
DEFAULT_INDEX_DIR = Path(r"CURATED_SET/BioAnalyze/data/gene_compare/index")
DETAIL_INDEX_FILENAME = "shared_h2a_gene_names_across_species_detail.csv"
HEATMAP_TO_PROCESSED_DIR = {"human": "homo_sapiens"}
HEATMAP_TO_SPECIES_NAME = {"human": "Homo sapiens"}
SKIP_HEATMAP_DIRS = {"alligned_human_pan", "gene_compare"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Count in how many species each H2A Gene name appears in the "
            "species-level *_present_gold.tsv files used for heatmaps and build "
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
        help="Directory with per-species processed *_present_gold.tsv files.",
    )
    parser.add_argument(
        "--out-dir",
        default=str(DEFAULT_STATS_DIR),
        help="Directory for summary table and plot.",
    )
    parser.add_argument(
        "--out-index-dir",
        default=str(DEFAULT_INDEX_DIR),
        help="Directory for the detailed reusable gene index.",
    )
    return parser.parse_args()


def path_str(path: Path) -> str:
    return path.as_posix()


def iter_species_dirs(heatmap_dir: Path) -> Iterable[str]:
    for child in sorted(heatmap_dir.iterdir()):
        if not child.is_dir():
            continue
        if child.name in SKIP_HEATMAP_DIRS:
            continue
        yield child.name


def resolve_processed_dir_name(heatmap_species_dir: str) -> str:
    return HEATMAP_TO_PROCESSED_DIR.get(heatmap_species_dir, heatmap_species_dir)


def infer_species_name(heatmap_species_dir: str, map_df: pd.DataFrame) -> str:
    if "species_name" in map_df.columns:
        species_names = (
            map_df["species_name"].fillna("").astype(str).str.strip().loc[lambda s: s.ne("")]
        )
        if not species_names.empty:
            return species_names.iloc[0]
    return HEATMAP_TO_SPECIES_NAME.get(
        heatmap_species_dir,
        heatmap_species_dir.replace("_", " ").title(),
    )


def find_present_gold_file(processed_species_dir: Path) -> Path:
    matches = sorted(processed_species_dir.glob("*_expr_advanced_H2A_present_gold.tsv"))
    if not matches:
        matches = sorted(processed_species_dir.glob("*present_gold.tsv"))
    if not matches:
        raise FileNotFoundError(f"No *_present_gold.tsv found in {processed_species_dir}")
    if len(matches) > 1:
        raise RuntimeError(f"Multiple *_present_gold.tsv files found in {processed_species_dir}")
    return matches[0]


def find_canonical_map_file(processed_species_dir: Path) -> Path:
    matches = sorted(processed_species_dir.glob("*canonical_variant_map.tsv"))
    if not matches:
        raise FileNotFoundError(f"No *canonical_variant_map.tsv found in {processed_species_dir}")
    if len(matches) > 1:
        raise RuntimeError(
            f"Multiple *canonical_variant_map.tsv files found in {processed_species_dir}"
        )
    return matches[0]


def normalize_map_df(map_tsv: Path) -> pd.DataFrame:
    map_df = pd.read_csv(map_tsv, sep="\t", dtype=str)
    for col in [
        "species_name",
        "ensembl_gene_id",
        "gene_name",
        "class",
        "label",
    ]:
        if col not in map_df.columns:
            map_df[col] = ""
        map_df[col] = map_df[col].fillna("").astype(str).str.strip()
    map_df = map_df[map_df["ensembl_gene_id"].ne("")].drop_duplicates(
        subset=["ensembl_gene_id"], keep="first"
    )
    return map_df


def load_present_gold(expr_tsv: Path) -> pd.DataFrame:
    expr_df = pd.read_csv(
        expr_tsv,
        sep="\t",
        dtype=str,
        usecols=[
            "Gene ID",
            "Gene name",
            "Anatomical entity name",
            "Expression",
            "Call quality",
            "Expression score",
        ],
        low_memory=True,
    )
    for col in [
        "Gene ID",
        "Gene name",
        "Anatomical entity name",
        "Expression",
        "Call quality",
    ]:
        expr_df[col] = expr_df[col].fillna("").astype(str).str.strip()
    expr_df["Expression score"] = pd.to_numeric(expr_df["Expression score"], errors="coerce")
    expr_df = expr_df[
        expr_df["Gene ID"].ne("")
        & expr_df["Anatomical entity name"].ne("")
        & expr_df["Expression"].eq("present")
        & expr_df["Call quality"].eq("gold quality")
        & expr_df["Expression score"].notna()
    ].copy()
    return expr_df


def build_detail_rows(
    heatmap_species_dir: str,
    expr_tsv: Path,
    map_tsv: Path,
) -> List[dict]:
    expr_df = load_present_gold(expr_tsv)
    map_df = normalize_map_df(map_tsv)
    species_name = infer_species_name(heatmap_species_dir, map_df)

    expr_gene_stats = (
        expr_df.groupby("Gene ID", as_index=False)
        .agg(
            row_count=("Gene ID", "size"),
            tissue_count=("Anatomical entity name", "nunique"),
        )
        .rename(columns={"Gene ID": "ensembl_gene_id"})
    )

    raw_name_map = (
        expr_df[["Gene ID", "Gene name"]]
        .assign(_len=lambda x: x["Gene name"].str.len())
        .sort_values(["Gene ID", "_len"], ascending=[True, False])
        .drop_duplicates(subset=["Gene ID"], keep="first")
        .rename(columns={"Gene ID": "ensembl_gene_id", "Gene name": "raw_gene_name"})
        .drop(columns=["_len"])
    )

    detail_df = expr_gene_stats.merge(raw_name_map, on="ensembl_gene_id", how="left")
    detail_df = detail_df.merge(
        map_df[["ensembl_gene_id", "gene_name", "class", "label"]],
        on="ensembl_gene_id",
        how="left",
    )

    detail_df["gene_name"] = (
        detail_df["gene_name"].fillna("").astype(str).str.strip()
    )
    detail_df["canonical_gene_name"] = detail_df["gene_name"]
    detail_df.loc[detail_df["canonical_gene_name"].eq(""), "canonical_gene_name"] = (
        detail_df.loc[detail_df["canonical_gene_name"].eq(""), "raw_gene_name"]
        .fillna("")
        .astype(str)
        .str.strip()
    )
    detail_df["gene_name"] = detail_df["canonical_gene_name"]
    detail_df["class"] = detail_df["class"].fillna("").astype(str).str.strip()
    detail_df["label"] = detail_df["label"].fillna("").astype(str).str.strip()
    detail_df = detail_df[detail_df["gene_name"].ne("")].copy()

    rows: List[dict] = []
    for row in detail_df.to_dict(orient="records"):
        rows.append(
            {
                "gene_name": row["gene_name"],
                "canonical_gene_name": row["canonical_gene_name"],
                "species_dir": heatmap_species_dir,
                "species_name": species_name,
                "ensembl_gene_id": row["ensembl_gene_id"],
                "map_label": row["label"],
                "class": row["class"],
                "present_gold_path": path_str(expr_tsv),
                "map_path": path_str(map_tsv),
                "tissue_count": int(row["tissue_count"]),
                "row_count": int(row["row_count"]),
            }
        )
    return rows


def build_detail_index(heatmap_dir: Path, processed_dir: Path) -> pd.DataFrame:
    detail_rows: List[dict] = []

    for heatmap_species_dir in iter_species_dirs(heatmap_dir):
        processed_species_dir = processed_dir / resolve_processed_dir_name(heatmap_species_dir)
        if not processed_species_dir.exists():
            raise FileNotFoundError(
                f"Processed directory is missing for heatmap species '{heatmap_species_dir}': "
                f"{processed_species_dir}"
            )

        expr_tsv = find_present_gold_file(processed_species_dir)
        map_tsv = find_canonical_map_file(processed_species_dir)
        detail_rows.extend(build_detail_rows(heatmap_species_dir, expr_tsv, map_tsv))

    detail_df = pd.DataFrame(detail_rows)
    if detail_df.empty:
        return pd.DataFrame(
            columns=[
                "gene_name",
                "canonical_gene_name",
                "species_dir",
                "species_name",
                "ensembl_gene_id",
                "map_label",
                "class",
                "present_gold_path",
                "map_path",
                "tissue_count",
                "row_count",
            ]
        )

    detail_df = detail_df.sort_values(
        by=["gene_name", "species_dir", "ensembl_gene_id"],
        ascending=[True, True, True],
    ).reset_index(drop=True)
    return detail_df


def summarize_species_ensembl_ids(group: pd.DataFrame) -> str:
    chunks: List[str] = []
    for species_dir, species_df in group.groupby("species_dir", sort=True):
        ids = sorted(species_df["ensembl_gene_id"].dropna().astype(str).unique().tolist())
        chunks.append(f"{species_dir}:{';'.join(ids)}")
    return "|".join(chunks)


def build_summary(detail_df: pd.DataFrame, detail_index_path: Path) -> pd.DataFrame:
    if detail_df.empty:
        return pd.DataFrame(
            columns=[
                "gene_name",
                "canonical_gene_name",
                "species_count",
                "species",
                "species_with_gene",
                "species_ensembl_ids",
                "detail_index_path",
            ]
        )

    rows: List[dict] = []
    for gene_name, group in detail_df.groupby("gene_name", sort=False):
        species_dirs = sorted(group["species_dir"].dropna().astype(str).unique().tolist())
        species_csv = ",".join(species_dirs)
        rows.append(
            {
                "gene_name": gene_name,
                "canonical_gene_name": gene_name,
                "species_count": len(species_dirs),
                "species": species_csv,
                "species_with_gene": species_csv,
                "species_ensembl_ids": summarize_species_ensembl_ids(group),
                "detail_index_path": path_str(detail_index_path),
            }
        )

    summary_df = pd.DataFrame(rows)
    summary_df = summary_df[summary_df["species_count"] > 1].copy()
    if summary_df.empty:
        return summary_df

    summary_df = summary_df.sort_values(
        by=["species_count", "gene_name"], ascending=[False, True]
    ).reset_index(drop=True)
    return summary_df


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
    ax.set_xlabel("Number of species with gene in *_present_gold.tsv")
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
    detail_df = build_detail_index(heatmap_dir, processed_dir)
    detail_df.to_csv(detail_index_path, index=False, encoding="utf-8")

    summary_df = build_summary(detail_df, detail_index_path)

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
