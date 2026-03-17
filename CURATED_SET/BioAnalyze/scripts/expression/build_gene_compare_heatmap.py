#!/usr/bin/env python
"""Build a cross-species heatmap for one canonical gene name."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Dict, Iterable, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


DEFAULT_HEATMAP_DIR = Path(r"CURATED_SET/BioAnalyze/figures/heatmaps")
DEFAULT_PROCESSED_DIR = Path(r"CURATED_SET/BioAnalyze/data/processed")
DEFAULT_OUT_FIG_ROOT = Path(r"CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare")
DEFAULT_OUT_DATA_ROOT = Path(r"CURATED_SET/BioAnalyze/data/gene_compare")
DEFAULT_DETAIL_INDEX = (
    DEFAULT_OUT_DATA_ROOT / "index" / "shared_h2a_gene_names_across_species_detail.csv"
)
HEATMAP_TO_PROCESSED_DIR = {"human": "homo_sapiens"}
HEATMAP_TO_SPECIES_NAME = {"human": "Homo sapiens"}
SKIP_HEATMAP_DIRS = {"alligned_human_pan", "gene_compare"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a cross-species heatmap for one canonical gene name using "
            "the present+gold H2A files and the canonical variant maps."
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
        default=str(DEFAULT_OUT_FIG_ROOT),
        help="Root directory for gene_compare heatmap figures.",
    )
    parser.add_argument(
        "--out-data-root",
        default=str(DEFAULT_OUT_DATA_ROOT),
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
        action="store_true",
        help="Include all known species as empty columns when the gene is absent.",
    )
    parser.add_argument(
        "--aggregate",
        default="mean",
        choices=["mean"],
        help="Aggregation for multiple rows per tissue/species.",
    )
    return parser.parse_args()


def path_str(path: Path) -> str:
    return path.as_posix()


def safe_slug(text: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9._-]+", "_", text.strip())
    slug = slug.strip("._")
    return slug or "gene"


def iter_species_dirs(heatmap_dir: Path) -> Iterable[str]:
    for child in sorted(heatmap_dir.iterdir()):
        if not child.is_dir():
            continue
        if child.name in SKIP_HEATMAP_DIRS:
            continue
        yield child.name


def resolve_processed_dir_name(heatmap_species_dir: str) -> str:
    return HEATMAP_TO_PROCESSED_DIR.get(heatmap_species_dir, heatmap_species_dir)


def infer_species_name(heatmap_species_dir: str) -> str:
    return HEATMAP_TO_SPECIES_NAME.get(
        heatmap_species_dir,
        heatmap_species_dir.replace("_", " ").title(),
    )


def load_present_gold(expr_tsv: Path) -> pd.DataFrame:
    expr_df = pd.read_csv(
        expr_tsv,
        sep="\t",
        dtype=str,
        usecols=[
            "Gene ID",
            "Anatomical entity name",
            "Expression",
            "Call quality",
            "Expression score",
        ],
        low_memory=True,
    )
    for col in ["Gene ID", "Anatomical entity name", "Expression", "Call quality"]:
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


def build_detail_index_fallback(heatmap_dir: Path, processed_dir: Path) -> pd.DataFrame:
    rows: List[dict] = []
    for heatmap_species_dir in iter_species_dirs(heatmap_dir):
        processed_species_dir = processed_dir / resolve_processed_dir_name(heatmap_species_dir)
        expr_matches = sorted(processed_species_dir.glob("*_expr_advanced_H2A_present_gold.tsv"))
        map_matches = sorted(processed_species_dir.glob("*canonical_variant_map.tsv"))
        if len(expr_matches) != 1 or len(map_matches) != 1:
            continue

        expr_tsv = expr_matches[0]
        map_tsv = map_matches[0]
        expr_df = load_present_gold(expr_tsv)
        map_df = pd.read_csv(map_tsv, sep="\t", dtype=str)
        for col in ["ensembl_gene_id", "gene_name", "class", "label"]:
            if col not in map_df.columns:
                map_df[col] = ""
            map_df[col] = map_df[col].fillna("").astype(str).str.strip()
        map_df = map_df[map_df["ensembl_gene_id"].ne("")].drop_duplicates(
            subset=["ensembl_gene_id"], keep="first"
        )

        expr_gene_stats = (
            expr_df.groupby("Gene ID", as_index=False)
            .agg(
                row_count=("Gene ID", "size"),
                tissue_count=("Anatomical entity name", "nunique"),
            )
            .rename(columns={"Gene ID": "ensembl_gene_id"})
        )
        joined = expr_gene_stats.merge(map_df, on="ensembl_gene_id", how="left")
        joined["gene_name"] = joined["gene_name"].fillna("").astype(str).str.strip()
        joined = joined[joined["gene_name"].ne("")].copy()

        species_name = infer_species_name(heatmap_species_dir)
        for row in joined.to_dict(orient="records"):
            rows.append(
                {
                    "gene_name": row["gene_name"],
                    "canonical_gene_name": row["gene_name"],
                    "species_dir": heatmap_species_dir,
                    "species_name": species_name,
                    "ensembl_gene_id": row["ensembl_gene_id"],
                    "map_label": row.get("label", ""),
                    "class": row.get("class", ""),
                    "present_gold_path": path_str(expr_tsv),
                    "map_path": path_str(map_tsv),
                    "tissue_count": int(row["tissue_count"]),
                    "row_count": int(row["row_count"]),
                }
            )

    return pd.DataFrame(rows)


def load_detail_index(detail_index_path: Path, heatmap_dir: Path, processed_dir: Path) -> pd.DataFrame:
    if detail_index_path.exists():
        detail_df = pd.read_csv(detail_index_path, dtype=str)
        if detail_df.empty:
            return detail_df
        for col in [
            "gene_name",
            "canonical_gene_name",
            "species_dir",
            "species_name",
            "ensembl_gene_id",
            "map_label",
            "class",
            "present_gold_path",
            "map_path",
        ]:
            if col not in detail_df.columns:
                detail_df[col] = ""
            detail_df[col] = detail_df[col].fillna("").astype(str).str.strip()
        for col in ["tissue_count", "row_count"]:
            if col not in detail_df.columns:
                detail_df[col] = 0
            detail_df[col] = pd.to_numeric(detail_df[col], errors="coerce").fillna(0).astype(int)
        return detail_df
    return build_detail_index_fallback(heatmap_dir, processed_dir)


def collect_gene_rows(
    gene_name: str,
    detail_df: pd.DataFrame,
    include_absent_species: bool,
    heatmap_dir: Path,
) -> tuple[pd.DataFrame, List[str]]:
    target_norm = gene_name.strip().casefold()
    gene_rows = detail_df[
        detail_df["canonical_gene_name"].fillna("").astype(str).str.casefold().eq(target_norm)
    ].copy()

    if gene_rows.empty:
        gene_rows = detail_df[
            detail_df["gene_name"].fillna("").astype(str).str.casefold().eq(target_norm)
        ].copy()

    if gene_rows.empty:
        raise RuntimeError(f"Gene name '{gene_name}' was not found in the detail index.")

    species_order = sorted(gene_rows["species_dir"].dropna().astype(str).unique().tolist())
    if include_absent_species:
        species_order = list(iter_species_dirs(heatmap_dir))
    return gene_rows, species_order


def build_long_dataframe(gene_name: str, gene_rows: pd.DataFrame) -> pd.DataFrame:
    long_frames: List[pd.DataFrame] = []
    for species_dir, species_df in gene_rows.groupby("species_dir", sort=True):
        expr_path = Path(species_df["present_gold_path"].iloc[0])
        expr_df = load_present_gold(expr_path)
        target_ids = sorted(species_df["ensembl_gene_id"].dropna().astype(str).unique().tolist())
        filtered = expr_df[expr_df["Gene ID"].isin(target_ids)].copy()
        if filtered.empty:
            continue

        filtered = filtered.rename(
            columns={
                "Gene ID": "ensembl_gene_id",
                "Anatomical entity name": "tissue",
                "Expression score": "expression_score",
            }
        )
        filtered["gene_name"] = gene_name
        filtered["species_dir"] = species_dir
        filtered["species_name"] = species_df["species_name"].iloc[0]

        meta_df = species_df[
            ["ensembl_gene_id", "map_label", "class"]
        ].drop_duplicates(subset=["ensembl_gene_id"], keep="first")
        filtered = filtered.merge(meta_df, on="ensembl_gene_id", how="left")
        long_frames.append(
            filtered[
                [
                    "gene_name",
                    "species_dir",
                    "species_name",
                    "ensembl_gene_id",
                    "tissue",
                    "expression_score",
                    "map_label",
                    "class",
                ]
            ]
        )

    if not long_frames:
        raise RuntimeError(
            f"Gene name '{gene_name}' exists in the index but has no present+gold rows."
        )

    long_df = pd.concat(long_frames, ignore_index=True)
    long_df["agg_expression_score"] = (
        long_df.groupby(["gene_name", "species_dir", "tissue"])["expression_score"]
        .transform("mean")
        .astype(float)
    )
    return long_df


def build_matrix(long_df: pd.DataFrame, species_order: List[str]) -> pd.DataFrame:
    matrix_df = (
        long_df.groupby(["tissue", "species_dir"], as_index=False)["expression_score"]
        .mean()
        .pivot(index="tissue", columns="species_dir", values="expression_score")
    )
    matrix_df = matrix_df.reindex(columns=species_order)

    tissue_presence = matrix_df.notna().sum(axis=1)
    matrix_df = matrix_df.assign(_species_presence=tissue_presence)
    matrix_df = matrix_df.sort_values(by=["_species_presence"], ascending=[False], kind="stable")
    ordered_tissues = sorted(
        matrix_df.index.tolist(),
        key=lambda tissue: (-int(matrix_df.loc[tissue, "_species_presence"]), tissue),
    )
    matrix_df = matrix_df.reindex(ordered_tissues).drop(columns="_species_presence")
    return matrix_df


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
    ax.set_title(f"{gene_name} Expression Across Species")
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


def main() -> None:
    args = parse_args()

    gene_name = args.gene_name.strip()
    heatmap_dir = Path(args.heatmap_dir)
    processed_dir = Path(args.processed_dir)
    out_fig_root = Path(args.out_fig_root)
    out_data_root = Path(args.out_data_root)
    detail_index_path = Path(args.detail_index)

    if not gene_name:
        raise RuntimeError("Gene name must not be empty.")
    if not heatmap_dir.exists():
        raise FileNotFoundError(heatmap_dir)
    if not processed_dir.exists():
        raise FileNotFoundError(processed_dir)

    detail_df = load_detail_index(detail_index_path, heatmap_dir, processed_dir)
    gene_rows, species_order = collect_gene_rows(
        gene_name,
        detail_df,
        args.include_absent_species,
        heatmap_dir,
    )
    canonical_gene_name = gene_rows["canonical_gene_name"].iloc[0]
    gene_slug = safe_slug(canonical_gene_name)

    long_df = build_long_dataframe(canonical_gene_name, gene_rows)
    matrix_df = build_matrix(long_df, species_order)

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
        tissue_mode=args.tissue_mode,
        aggregate=args.aggregate,
        include_absent_species=args.include_absent_species,
        detail_index_path=detail_index_path,
        metadata_path=out_metadata_json,
    )

    print(f"Saved long table to {out_long_csv}")
    print(f"Saved matrix table to {out_matrix_csv}")
    print(f"Saved heatmap to {out_png}")
    print(f"Saved heatmap to {out_svg}")
    print(f"Saved metadata to {out_metadata_json}")


if __name__ == "__main__":
    main()
