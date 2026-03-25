#!/usr/bin/env python
"""Shared helpers for cross-species H2A gene comparison workflows."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import pandas as pd

from normalized_expression_common import (
    CELL_STATUS_OBSERVED_ZERO,
    CELL_STATUS_PRESENT_GOLD,
    load_processed_expression_cells,
)


DEFAULT_HEATMAP_ROOT = Path(r"CURATED_SET/BioAnalyze/figures/heatmaps")
DEFAULT_HEATMAP_DIR = DEFAULT_HEATMAP_ROOT / "species"
DEFAULT_PROCESSED_DIR = Path(r"CURATED_SET/BioAnalyze/data/processed")
DEFAULT_STATS_DIR = Path(r"CURATED_SET/BioAnalyze/stats")
DEFAULT_SHARED_STATS_DIR = DEFAULT_STATS_DIR / "shared_genes"
DEFAULT_ACCESSION_STATS_DIR = DEFAULT_STATS_DIR / "accession_stats"
DEFAULT_RANKING_STATS_DIR = DEFAULT_STATS_DIR / "ranking"
DEFAULT_RANKING_TABLES_DIR = DEFAULT_RANKING_STATS_DIR / "tables"
DEFAULT_RANKING_REPORTS_DIR = DEFAULT_RANKING_STATS_DIR / "reports"
DEFAULT_RANKING_PLOTS_DIR = DEFAULT_RANKING_STATS_DIR / "plots"
DEFAULT_GENE_COMPARE_FIG_ROOT = DEFAULT_HEATMAP_ROOT / "gene_compare"
DEFAULT_GENE_COMPARE_DATA_ROOT = Path(r"CURATED_SET/BioAnalyze/data/gene_compare")
DEFAULT_RANKING_OUT_DIR = DEFAULT_RANKING_STATS_DIR
DEFAULT_INDEX_DIR = DEFAULT_GENE_COMPARE_DATA_ROOT / "index"
DEFAULT_SHARED_INDEX = DEFAULT_SHARED_STATS_DIR / "shared_h2a_gene_names_across_species.csv"
DETAIL_INDEX_FILENAME = "shared_h2a_gene_names_across_species_detail.csv"
DEFAULT_DETAIL_INDEX = DEFAULT_INDEX_DIR / DETAIL_INDEX_FILENAME
HEATMAP_TO_PROCESSED_DIR = {"human": "homo_sapiens"}
HEATMAP_TO_SPECIES_NAME = {"human": "Homo sapiens"}
SKIP_HEATMAP_DIRS = {"alligned_human_pan", "gene_compare"}
GENERIC_TISSUES = {
    "material anatomical entity",
    "anatomical system",
    "multicellular organism",
}


def path_str(path: Path) -> str:
    return path.as_posix()


def safe_slug(text: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9._-]+", "_", text.strip())
    slug = slug.strip("._")
    return slug or "gene"


def ensure_mean_aggregate(aggregate: str) -> str:
    if aggregate != "mean":
        raise ValueError(f"Unsupported aggregate: {aggregate}")
    return aggregate


def iter_species_dirs(heatmap_dir: Path) -> Iterable[str]:
    for child in sorted(heatmap_dir.iterdir()):
        if not child.is_dir():
            continue
        if child.name.startswith("_"):
            continue
        if child.name in SKIP_HEATMAP_DIRS:
            continue
        yield child.name


def resolve_processed_dir_name(heatmap_species_dir: str) -> str:
    return HEATMAP_TO_PROCESSED_DIR.get(heatmap_species_dir, heatmap_species_dir)


def infer_species_name(
    heatmap_species_dir: str,
    map_df: Optional[pd.DataFrame] = None,
) -> str:
    if map_df is not None and "species_name" in map_df.columns:
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
    for col in ["species_name", "ensembl_gene_id", "gene_name", "class", "label"]:
        if col not in map_df.columns:
            map_df[col] = ""
        map_df[col] = map_df[col].fillna("").astype(str).str.strip()
    map_df = map_df[map_df["ensembl_gene_id"].ne("")].drop_duplicates(
        subset=["ensembl_gene_id"], keep="first"
    )
    return map_df


def load_processed_expression(
    expr_tsv: Path,
    expr_cache: Optional[Dict[str, pd.DataFrame]] = None,
) -> pd.DataFrame:
    return load_processed_expression_cells(expr_tsv, expr_cache=expr_cache)


def build_detail_rows(
    heatmap_species_dir: str,
    expr_tsv: Path,
    map_tsv: Path,
    expr_cache: Optional[Dict[str, pd.DataFrame]] = None,
) -> List[dict]:
    expr_df = load_processed_expression(expr_tsv, expr_cache=expr_cache)
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

    detail_df["gene_name"] = detail_df["gene_name"].fillna("").astype(str).str.strip()
    detail_df["canonical_gene_name"] = detail_df["gene_name"]
    missing_mask = detail_df["canonical_gene_name"].eq("")
    detail_df.loc[missing_mask, "canonical_gene_name"] = (
        detail_df.loc[missing_mask, "raw_gene_name"].fillna("").astype(str).str.strip()
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


def normalize_detail_index(detail_df: pd.DataFrame) -> pd.DataFrame:
    if detail_df.empty:
        return detail_df

    detail_df = detail_df.copy()
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
    detail_df = detail_df.sort_values(
        by=["gene_name", "species_dir", "ensembl_gene_id"],
        ascending=[True, True, True],
    ).reset_index(drop=True)
    return detail_df


def build_detail_index(
    heatmap_dir: Path,
    processed_dir: Path,
    expr_cache: Optional[Dict[str, pd.DataFrame]] = None,
) -> pd.DataFrame:
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
        detail_rows.extend(
            build_detail_rows(
                heatmap_species_dir,
                expr_tsv,
                map_tsv,
                expr_cache=expr_cache,
            )
        )

    if not detail_rows:
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

    return normalize_detail_index(pd.DataFrame(detail_rows))


def load_detail_index(
    detail_index_path: Path,
    heatmap_dir: Optional[Path] = None,
    processed_dir: Optional[Path] = None,
    expr_cache: Optional[Dict[str, pd.DataFrame]] = None,
) -> pd.DataFrame:
    if detail_index_path.exists():
        return normalize_detail_index(pd.read_csv(detail_index_path, dtype=str))
    if heatmap_dir is None or processed_dir is None:
        raise FileNotFoundError(detail_index_path)
    return build_detail_index(heatmap_dir, processed_dir, expr_cache=expr_cache)


def summarize_species_ensembl_ids(group: pd.DataFrame) -> str:
    chunks: List[str] = []
    for species_dir, species_df in group.groupby("species_dir", sort=True):
        ids = sorted(species_df["ensembl_gene_id"].dropna().astype(str).unique().tolist())
        chunks.append(f"{species_dir}:{';'.join(ids)}")
    return "|".join(chunks)


def build_shared_gene_summary(detail_df: pd.DataFrame, detail_index_path: Path) -> pd.DataFrame:
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

    return summary_df.sort_values(
        by=["species_count", "gene_name"], ascending=[False, True]
    ).reset_index(drop=True)


def collect_gene_rows(
    gene_name: str,
    detail_df: pd.DataFrame,
    include_absent_species: bool = False,
    heatmap_dir: Optional[Path] = None,
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
        if heatmap_dir is None:
            raise ValueError("heatmap_dir is required when include_absent_species=True")
        species_order = list(iter_species_dirs(heatmap_dir))
    return gene_rows, species_order


def build_gene_long_dataframe(
    gene_name: str,
    gene_rows: pd.DataFrame,
    expr_cache: Optional[Dict[str, pd.DataFrame]] = None,
    aggregate: str = "mean",
) -> pd.DataFrame:
    ensure_mean_aggregate(aggregate)
    long_frames: List[pd.DataFrame] = []

    for species_dir, species_df in gene_rows.groupby("species_dir", sort=True):
        expr_path = Path(species_df["present_gold_path"].iloc[0])
        expr_df = load_processed_expression(expr_path, expr_cache=expr_cache)
        target_ids = sorted(species_df["ensembl_gene_id"].dropna().astype(str).unique().tolist())
        filtered = expr_df[expr_df["Gene ID"].isin(target_ids)].copy()
        if filtered.empty:
            continue

        filtered = filtered.rename(
            columns={
                "Gene ID": "ensembl_gene_id",
                "Anatomical entity name": "tissue",
                "cell_mean_score": "cell_mean_score",
                "cell_std_score": "cell_std_score",
                "cell_n": "cell_n",
                "cell_status": "cell_status",
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
                    "cell_mean_score",
                    "cell_std_score",
                    "cell_n",
                    "cell_status",
                    "map_label",
                    "class",
                ]
            ]
        )

    if not long_frames:
        raise RuntimeError(
            f"Gene name '{gene_name}' exists in the index but has no normalized rows."
        )

    long_df = pd.concat(long_frames, ignore_index=True)
    long_df["cell_mean_score"] = pd.to_numeric(long_df["cell_mean_score"], errors="coerce")
    long_df["cell_std_score"] = pd.to_numeric(long_df["cell_std_score"], errors="coerce").fillna(0.0)
    long_df["cell_n"] = pd.to_numeric(long_df["cell_n"], errors="coerce").fillna(0).astype(int)
    long_df["cell_status"] = long_df["cell_status"].fillna("").astype(str).str.strip()
    long_df["expression_score"] = long_df["cell_mean_score"]
    long_df["agg_expression_score"] = long_df["cell_mean_score"]
    return long_df


def summarize_species_gene_tissue(long_df: pd.DataFrame) -> pd.DataFrame:
    if long_df.empty:
        return pd.DataFrame(
            columns=[
                "gene_name",
                "species_dir",
                "species_name",
                "tissue",
                "cell_mean_score",
                "cell_std_score",
                "cell_n",
                "cell_status",
                "expression_score",
                "agg_expression_score",
            ]
        )

    summary_rows: List[dict] = []
    grouped = long_df.groupby(["gene_name", "species_dir", "species_name", "tissue"], sort=False)
    for (gene_name, species_dir, species_name, tissue), group in grouped:
        mean_score = float(group["cell_mean_score"].mean())
        if len(group) == 1:
            std_score = float(group["cell_std_score"].iloc[0])
            status = str(group["cell_status"].iloc[0])
        else:
            std_score = float(group["cell_mean_score"].std(ddof=1))
            if pd.isna(std_score):
                std_score = 0.0
            status = (
                CELL_STATUS_PRESENT_GOLD
                if group["cell_status"].eq(CELL_STATUS_PRESENT_GOLD).any()
                else CELL_STATUS_OBSERVED_ZERO
            )
        summary_rows.append(
            {
                "gene_name": gene_name,
                "species_dir": species_dir,
                "species_name": species_name,
                "tissue": tissue,
                "cell_mean_score": mean_score,
                "cell_std_score": std_score,
                "cell_n": int(group["cell_n"].sum()),
                "cell_status": status,
                "expression_score": mean_score,
                "agg_expression_score": mean_score,
            }
        )

    return pd.DataFrame(summary_rows)


def build_matrix(
    long_df: pd.DataFrame,
    species_order: List[str],
    aggregate: str = "mean",
) -> pd.DataFrame:
    ensure_mean_aggregate(aggregate)
    species_level_df = summarize_species_gene_tissue(long_df)
    matrix_df = (
        species_level_df.groupby(["tissue", "species_dir"], as_index=False)["cell_mean_score"]
        .mean()
        .pivot(index="tissue", columns="species_dir", values="cell_mean_score")
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
