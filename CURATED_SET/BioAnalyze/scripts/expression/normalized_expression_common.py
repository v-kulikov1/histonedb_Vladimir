"""Shared helpers for normalized H2A expression cell tables."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd


RAW_EXPR_USECOLS = [
    "Gene ID",
    "Gene name",
    "Anatomical entity name",
    "Expression",
    "Call quality",
    "Expression score",
]

CELL_STATUS_PRESENT_GOLD = "present_gold"
CELL_STATUS_OBSERVED_ZERO = "observed_zero"
CELL_STATUS_MISSING = "missing"

PROCESSED_EXPR_TEXT_COLS = [
    "Gene ID",
    "Gene name",
    "Anatomical entity name",
    "cell_status",
]

PROCESSED_EXPR_NUMERIC_COLS = [
    "cell_mean_score",
    "cell_std_score",
    "cell_n",
    "Expression score",
    "Expression std",
    "Expression n",
    "raw_row_count",
    "qualifying_row_count",
]


def _clean_text_columns(df: pd.DataFrame, columns: Sequence[str]) -> pd.DataFrame:
    for col in columns:
        if col not in df.columns:
            df[col] = ""
        df[col] = df[col].fillna("").astype(str).str.strip()
    return df


def _std_or_zero(values: pd.Series) -> float:
    std_value = values.astype(float).std(ddof=1)
    if pd.isna(std_value):
        return 0.0
    return float(std_value)


def _longest_name_per_gene(df: pd.DataFrame) -> pd.DataFrame:
    return (
        df[["Gene ID", "Gene name"]]
        .assign(_len=lambda x: x["Gene name"].str.len())
        .sort_values(["Gene ID", "_len"], ascending=[True, False])
        .drop_duplicates(subset=["Gene ID"], keep="first")
        .drop(columns=["_len"])
    )


def normalize_h2a_expression_cells(
    expr_path: Path,
    ensg_set: set[str],
    chunksize: int,
) -> Tuple[pd.DataFrame, int, int]:
    """Build one normalized row per Gene ID x tissue cell for H2A genes."""
    rows_total = 0
    rows_h2a = 0
    chunks: List[pd.DataFrame] = []

    for chunk in pd.read_csv(
        expr_path,
        sep="\t",
        dtype=str,
        usecols=RAW_EXPR_USECOLS,
        chunksize=chunksize,
        low_memory=True,
    ):
        rows_total += len(chunk)
        chunk = _clean_text_columns(
            chunk,
            ["Gene ID", "Gene name", "Anatomical entity name", "Expression", "Call quality"],
        )
        chunk = chunk[
            chunk["Gene ID"].isin(ensg_set)
            & chunk["Gene ID"].ne("")
            & chunk["Anatomical entity name"].ne("")
        ].copy()
        if chunk.empty:
            continue

        rows_h2a += len(chunk)
        chunk["Expression score"] = pd.to_numeric(chunk["Expression score"], errors="coerce")
        chunks.append(chunk)

    if not chunks:
        return (
            pd.DataFrame(
                columns=[
                    "Gene ID",
                    "Gene name",
                    "Anatomical entity name",
                    "cell_mean_score",
                    "cell_std_score",
                    "cell_n",
                    "cell_status",
                    "Expression score",
                    "Expression std",
                    "Expression n",
                    "raw_row_count",
                    "qualifying_row_count",
                ]
            ),
            rows_total,
            rows_h2a,
        )

    raw_df = pd.concat(chunks, ignore_index=True)
    raw_df = _clean_text_columns(
        raw_df,
        ["Gene ID", "Gene name", "Anatomical entity name", "Expression", "Call quality"],
    )

    observed_cells = (
        raw_df.groupby(["Gene ID", "Anatomical entity name"], as_index=False)
        .agg(raw_row_count=("Gene ID", "size"))
        .reset_index(drop=True)
    )

    qualifying_df = raw_df[
        raw_df["Expression"].eq("present")
        & raw_df["Call quality"].eq("gold quality")
        & raw_df["Expression score"].notna()
    ].copy()

    if qualifying_df.empty:
        qualifying_stats = pd.DataFrame(
            columns=[
                "Gene ID",
                "Anatomical entity name",
                "cell_mean_score",
                "cell_std_score",
                "cell_n",
                "qualifying_row_count",
            ]
        )
    else:
        qualifying_stats = (
            qualifying_df.groupby(["Gene ID", "Anatomical entity name"], as_index=False)
            .agg(
                cell_mean_score=("Expression score", "mean"),
                cell_std_score=("Expression score", _std_or_zero),
                cell_n=("Expression score", "size"),
                qualifying_row_count=("Expression score", "size"),
            )
            .reset_index(drop=True)
        )

    normalized_df = observed_cells.merge(
        qualifying_stats,
        on=["Gene ID", "Anatomical entity name"],
        how="left",
    )
    normalized_df = normalized_df.merge(
        _longest_name_per_gene(raw_df),
        on="Gene ID",
        how="left",
    )
    normalized_df = _clean_text_columns(normalized_df, ["Gene ID", "Gene name", "Anatomical entity name"])

    normalized_df["cell_n"] = pd.to_numeric(normalized_df["cell_n"], errors="coerce").fillna(0).astype(int)
    normalized_df["qualifying_row_count"] = (
        pd.to_numeric(normalized_df["qualifying_row_count"], errors="coerce").fillna(0).astype(int)
    )
    normalized_df["raw_row_count"] = (
        pd.to_numeric(normalized_df["raw_row_count"], errors="coerce").fillna(0).astype(int)
    )
    normalized_df["cell_mean_score"] = pd.to_numeric(
        normalized_df["cell_mean_score"], errors="coerce"
    )
    normalized_df["cell_std_score"] = pd.to_numeric(
        normalized_df["cell_std_score"], errors="coerce"
    )

    zero_mask = normalized_df["cell_n"].eq(0)
    normalized_df.loc[zero_mask, "cell_mean_score"] = 0.0
    normalized_df.loc[zero_mask, "cell_std_score"] = 0.0
    normalized_df.loc[~zero_mask & normalized_df["cell_std_score"].isna(), "cell_std_score"] = 0.0
    normalized_df["cell_status"] = CELL_STATUS_PRESENT_GOLD
    normalized_df.loc[zero_mask, "cell_status"] = CELL_STATUS_OBSERVED_ZERO

    normalized_df["Expression score"] = normalized_df["cell_mean_score"]
    normalized_df["Expression std"] = normalized_df["cell_std_score"]
    normalized_df["Expression n"] = normalized_df["cell_n"]

    normalized_df = normalized_df[
        [
            "Gene ID",
            "Gene name",
            "Anatomical entity name",
            "cell_mean_score",
            "cell_std_score",
            "cell_n",
            "cell_status",
            "Expression score",
            "Expression std",
            "Expression n",
            "raw_row_count",
            "qualifying_row_count",
        ]
    ].sort_values(
        by=["Gene ID", "Anatomical entity name"],
        ascending=[True, True],
    ).reset_index(drop=True)

    return normalized_df, rows_total, rows_h2a


def load_processed_expression_cells(
    expr_tsv: Path,
    expr_cache: Optional[Dict[str, pd.DataFrame]] = None,
) -> pd.DataFrame:
    cache_key = expr_tsv.as_posix()
    if expr_cache is not None and cache_key in expr_cache:
        return expr_cache[cache_key]

    expr_df = pd.read_csv(expr_tsv, sep="\t", low_memory=True)
    expr_df = _clean_text_columns(expr_df, PROCESSED_EXPR_TEXT_COLS)

    if "cell_mean_score" not in expr_df.columns and "Expression score" in expr_df.columns:
        expr_df["cell_mean_score"] = expr_df["Expression score"]
    if "cell_std_score" not in expr_df.columns and "Expression std" in expr_df.columns:
        expr_df["cell_std_score"] = expr_df["Expression std"]
    if "cell_n" not in expr_df.columns and "Expression n" in expr_df.columns:
        expr_df["cell_n"] = expr_df["Expression n"]

    for col in PROCESSED_EXPR_NUMERIC_COLS:
        if col not in expr_df.columns:
            expr_df[col] = 0
        expr_df[col] = pd.to_numeric(expr_df[col], errors="coerce")

    expr_df["cell_n"] = expr_df["cell_n"].fillna(0).astype(int)
    expr_df["Expression n"] = expr_df["Expression n"].fillna(expr_df["cell_n"]).astype(int)
    expr_df["raw_row_count"] = expr_df["raw_row_count"].fillna(0).astype(int)
    expr_df["qualifying_row_count"] = expr_df["qualifying_row_count"].fillna(0).astype(int)
    expr_df["cell_mean_score"] = expr_df["cell_mean_score"].astype(float)
    expr_df["cell_std_score"] = expr_df["cell_std_score"].fillna(0.0).astype(float)
    expr_df["Expression score"] = expr_df["Expression score"].fillna(expr_df["cell_mean_score"]).astype(float)
    expr_df["Expression std"] = expr_df["Expression std"].fillna(expr_df["cell_std_score"]).astype(float)

    if "cell_status" not in expr_df.columns:
        expr_df["cell_status"] = CELL_STATUS_PRESENT_GOLD
        expr_df.loc[expr_df["cell_mean_score"].eq(0) & expr_df["cell_n"].eq(0), "cell_status"] = (
            CELL_STATUS_OBSERVED_ZERO
        )

    expr_df = expr_df[
        expr_df["Gene ID"].ne("") & expr_df["Anatomical entity name"].ne("")
    ].copy()

    if expr_cache is not None:
        expr_cache[cache_key] = expr_df
    return expr_df


def sort_gene_labels(label_map: Dict[str, str], gene_ids: Iterable[str]) -> List[str]:
    labels = [label_map[gid] for gid in gene_ids if gid in label_map]
    return sorted(labels, key=lambda s: (s.split(":", 1)[0], s))


def build_tissue_coverage_table(
    expr_df: pd.DataFrame,
    gene_ids: Sequence[str],
    *,
    threshold: float,
    panel: str,
) -> pd.DataFrame:
    gene_ids = [str(gid).strip() for gid in gene_ids if str(gid).strip()]
    gene_ids = list(dict.fromkeys(gene_ids))

    coverage_columns = [
        "anatomical_entity_name",
        "observed_gene_count",
        "display_gene_count",
        "fill_rate",
        "threshold",
        "kept",
        "panel",
    ]
    if not gene_ids:
        return pd.DataFrame(columns=coverage_columns)

    expr_df = expr_df.copy()
    for col in ["Gene ID", "Anatomical entity name"]:
        if col not in expr_df.columns:
            expr_df[col] = ""
        expr_df[col] = expr_df[col].fillna("").astype(str).str.strip()

    expr_df = expr_df[
        expr_df["Gene ID"].isin(gene_ids) & expr_df["Anatomical entity name"].ne("")
    ].copy()

    display_gene_count = len(gene_ids)
    if expr_df.empty:
        return pd.DataFrame(columns=coverage_columns)

    coverage_df = (
        expr_df.groupby("Anatomical entity name", as_index=False)
        .agg(observed_gene_count=("Gene ID", "nunique"))
        .rename(columns={"Anatomical entity name": "anatomical_entity_name"})
    )
    coverage_df["display_gene_count"] = display_gene_count
    coverage_df["fill_rate"] = (
        coverage_df["observed_gene_count"] / coverage_df["display_gene_count"]
    ).astype(float)
    coverage_df["threshold"] = float(threshold)
    coverage_df["kept"] = coverage_df["fill_rate"].ge(float(threshold))
    coverage_df["panel"] = panel
    coverage_df = coverage_df.sort_values(
        by=["kept", "fill_rate", "anatomical_entity_name"],
        ascending=[False, False, True],
        kind="stable",
    ).reset_index(drop=True)
    return coverage_df[coverage_columns]


def build_species_heatmap_display_index(
    map_df: pd.DataFrame,
    expr_df: pd.DataFrame,
    *,
    preferred_id_col: str,
) -> pd.DataFrame:
    """Choose one representative ENSG row per displayed species heatmap entity."""
    if map_df.empty:
        return pd.DataFrame(
            columns=[
                "ensembl_gene_id",
                "gene_name",
                "class",
                "label",
                "display_key",
                "display_id",
                "has_any_rows",
                "has_present_gold",
                "non_missing_tissues",
                "present_gold_tissues",
                "total_cell_n",
            ]
        )

    map_df = map_df.copy()
    for col in ["ensembl_gene_id", "gene_name", "class", "label", preferred_id_col]:
        if col not in map_df.columns:
            map_df[col] = ""
        map_df[col] = map_df[col].fillna("").astype(str).str.strip()

    map_df = map_df[map_df["ensembl_gene_id"].ne("")].drop_duplicates(
        subset=["ensembl_gene_id"], keep="first"
    )

    expr_df = expr_df.copy()
    for col in ["Gene ID", "Anatomical entity name", "cell_status"]:
        if col not in expr_df.columns:
            expr_df[col] = ""
        expr_df[col] = expr_df[col].fillna("").astype(str).str.strip()
    if "cell_n" not in expr_df.columns:
        expr_df["cell_n"] = 0
    expr_df["cell_n"] = pd.to_numeric(expr_df["cell_n"], errors="coerce").fillna(0).astype(int)

    expr_stats = (
        expr_df.groupby("Gene ID", as_index=False)
        .agg(
            non_missing_tissues=("Anatomical entity name", "nunique"),
            total_cell_n=("cell_n", "sum"),
        )
        .rename(columns={"Gene ID": "ensembl_gene_id"})
    )
    present_gold_stats = (
        expr_df[expr_df["cell_status"].eq(CELL_STATUS_PRESENT_GOLD)]
        .groupby("Gene ID", as_index=False)
        .agg(present_gold_tissues=("Anatomical entity name", "nunique"))
        .rename(columns={"Gene ID": "ensembl_gene_id"})
    )

    display_df = map_df.merge(expr_stats, on="ensembl_gene_id", how="left")
    display_df = display_df.merge(present_gold_stats, on="ensembl_gene_id", how="left")
    for col in ["non_missing_tissues", "present_gold_tissues", "total_cell_n"]:
        display_df[col] = pd.to_numeric(display_df[col], errors="coerce").fillna(0).astype(int)

    display_df["has_any_rows"] = display_df["non_missing_tissues"].gt(0)
    display_df["has_present_gold"] = display_df["present_gold_tissues"].gt(0)

    display_df["display_id"] = display_df[preferred_id_col].where(display_df[preferred_id_col].ne(""), "")
    missing_display_id = display_df["display_id"].eq("")
    display_df.loc[missing_display_id, "display_id"] = display_df.loc[missing_display_id, "gene_name"]
    missing_display_id = display_df["display_id"].eq("")
    display_df.loc[missing_display_id, "display_id"] = display_df.loc[missing_display_id, "ensembl_gene_id"]
    display_df["display_key"] = display_df["display_id"]

    display_df = display_df.sort_values(
        by=[
            "display_key",
            "has_any_rows",
            "has_present_gold",
            "non_missing_tissues",
            "present_gold_tissues",
            "total_cell_n",
            "ensembl_gene_id",
        ],
        ascending=[True, False, False, False, False, False, True],
        kind="stable",
    )
    display_df = display_df[display_df["has_any_rows"]].drop_duplicates(
        subset=["display_key"], keep="first"
    )

    display_df["label"] = display_df.apply(
        lambda row: (
            f"{row['gene_name']}:{row[preferred_id_col]}"
            if row[preferred_id_col]
            else f"{row['gene_name']}:{row['ensembl_gene_id']}"
            if row["gene_name"]
            else row["ensembl_gene_id"]
        ),
        axis=1,
    )

    return display_df.reset_index(drop=True)
