#!/usr/bin/env python
"""Rank cross-species H2A expression differences across shared genes and tissues."""

from __future__ import annotations

import argparse
import itertools
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

import pandas as pd

from gene_compare_common import (
    DEFAULT_DETAIL_INDEX,
    DEFAULT_GENE_COMPARE_FIG_ROOT,
    DEFAULT_HEATMAP_DIR,
    DEFAULT_PROCESSED_DIR,
    DEFAULT_RANKING_OUT_DIR,
    DEFAULT_RANKING_PLOTS_DIR,
    DEFAULT_SHARED_INDEX,
    GENERIC_TISSUES,
    build_gene_long_dataframe,
    build_shared_gene_summary,
    collect_gene_rows,
    load_detail_index,
    safe_slug,
    summarize_species_gene_tissue,
)


DEFAULT_MIN_GENE_SPECIES_COUNT = 4
MIN_CONSERVED_TISSUES = 10
ALL_SHARED_SHORTLIST_LIMIT = 10
VARIANT_SHORTLIST_LIMIT = 10
CONSERVED_SHORTLIST_LIMIT = 3
DEFAULT_LOW_QUANTILES = [0.05, 0.10]
DEFAULT_REPORT_LANGUAGE = "ru"
DEFAULT_REPORT_STEM = "cross_species_expression_report_ru"
TOP_VARIABLE_GENE_LIMIT = 10
TOP_CONSERVED_GENE_LIMIT = 10
TOP_CASE_LIMIT = 10
TOP_SPECIES_LIMIT = 8
TOP_PAIR_LIMIT = 10


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Rank strong cross-species H2A expression differences across shared "
            "canonical gene names, produce summary tables, and write manuscript-ready reports."
        )
    )
    parser.add_argument("--shared-index", default=str(DEFAULT_SHARED_INDEX))
    parser.add_argument("--detail-index", default=str(DEFAULT_DETAIL_INDEX))
    parser.add_argument("--heatmap-dir", default=str(DEFAULT_HEATMAP_DIR))
    parser.add_argument("--processed-dir", default=str(DEFAULT_PROCESSED_DIR))
    parser.add_argument("--out-dir", default=str(DEFAULT_RANKING_OUT_DIR))
    parser.add_argument("--plots-dir", default=str(DEFAULT_RANKING_PLOTS_DIR))
    parser.add_argument("--min-species-per-tissue", default=4, type=int)
    parser.add_argument("--candidate-quantile", default=0.90, type=float)
    parser.add_argument("--high-confidence-quantile", default=0.95, type=float)
    parser.add_argument("--low-quantiles", nargs="*", type=float, default=DEFAULT_LOW_QUANTILES)
    parser.add_argument(
        "--class-specific-low-quantiles",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument(
        "--scope",
        default="two-tier",
        choices=["two-tier", "all-shared", "variants-only"],
    )
    parser.add_argument("--report-language", default=DEFAULT_REPORT_LANGUAGE)
    parser.add_argument(
        "--write-detailed-report",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument("--report-stem", default=DEFAULT_REPORT_STEM)
    return parser.parse_args()


def quantile_label(value: float) -> str:
    return f"p{int(round(float(value) * 100)):02d}"


def normalize_low_quantiles(values: Sequence[float]) -> List[float]:
    cleaned: List[float] = []
    for value in values:
        value = float(value)
        if not 0 < value < 0.5:
            raise ValueError("--low-quantiles values must be between 0 and 0.5")
        cleaned.append(value)
    return sorted(set(cleaned))


def load_shared_index(shared_index_path: Path, detail_df: pd.DataFrame, detail_index_path: Path) -> pd.DataFrame:
    if shared_index_path.exists():
        shared_df = pd.read_csv(shared_index_path, dtype=str)
    else:
        shared_df = build_shared_gene_summary(detail_df, detail_index_path)
    if shared_df.empty:
        return shared_df

    for col in ["gene_name", "canonical_gene_name"]:
        if col not in shared_df.columns:
            shared_df[col] = ""
        shared_df[col] = shared_df[col].fillna("").astype(str).str.strip()
    if "species_count" not in shared_df.columns:
        shared_df["species_count"] = 0
    shared_df["species_count"] = pd.to_numeric(shared_df["species_count"], errors="coerce").fillna(0)
    return shared_df


def choose_target_genes(shared_df: pd.DataFrame, min_gene_species_count: int) -> List[str]:
    filtered = shared_df[shared_df["species_count"] >= int(min_gene_species_count)].copy()
    if filtered.empty:
        return []
    filtered = filtered.sort_values(
        by=["species_count", "canonical_gene_name"],
        ascending=[False, True],
    )
    return filtered["canonical_gene_name"].tolist()


def summarize_gene_class(gene_rows: pd.DataFrame) -> str:
    classes = sorted(
        {
            value
            for value in gene_rows["class"].fillna("").astype(str).str.strip().tolist()
            if value
        }
    )
    return ",".join(classes) if classes else ""


def build_species_level_long(detail_df: pd.DataFrame, target_genes: Iterable[str]) -> pd.DataFrame:
    expr_cache: Dict[str, pd.DataFrame] = {}
    frames: List[pd.DataFrame] = []
    for gene_name in target_genes:
        gene_rows, _ = collect_gene_rows(gene_name, detail_df)
        gene_class = summarize_gene_class(gene_rows)
        gene_long_df = build_gene_long_dataframe(
            gene_name,
            gene_rows,
            expr_cache=expr_cache,
            aggregate="mean",
        )
        species_level = summarize_species_gene_tissue(gene_long_df)
        species_level = species_level[
            [
                "gene_name",
                "species_dir",
                "species_name",
                "tissue",
                "cell_mean_score",
                "cell_std_score",
                "cell_n",
                "cell_status",
                "expression_score",
            ]
        ].copy()
        species_level["gene_class"] = gene_class
        species_level["gene_species_count"] = int(gene_rows["species_dir"].nunique())
        frames.append(species_level)

    if not frames:
        return pd.DataFrame(
            columns=[
                "gene_name",
                "species_dir",
                "species_name",
                "tissue",
                "expression_score",
                "cell_mean_score",
                "cell_std_score",
                "cell_n",
                "cell_status",
                "gene_class",
                "gene_species_count",
            ]
        )

    long_df = pd.concat(frames, ignore_index=True)
    long_df["expression_score"] = pd.to_numeric(long_df["expression_score"], errors="coerce")
    long_df = long_df[long_df["expression_score"].notna()].copy()
    long_df["gene_species_count"] = pd.to_numeric(
        long_df["gene_species_count"], errors="coerce"
    ).fillna(0).astype(int)
    long_df["cell_mean_score"] = pd.to_numeric(long_df["cell_mean_score"], errors="coerce")
    long_df["cell_std_score"] = pd.to_numeric(long_df["cell_std_score"], errors="coerce").fillna(0.0)
    long_df["cell_n"] = pd.to_numeric(long_df["cell_n"], errors="coerce").fillna(0).astype(int)
    long_df = long_df[~long_df["tissue"].isin(GENERIC_TISSUES)].copy()
    return long_df


def build_gene_tissue_summary(species_level_long_df: pd.DataFrame) -> pd.DataFrame:
    comparable_df = species_level_long_df.copy()
    if comparable_df.empty:
        return pd.DataFrame()

    summary_df = (
        comparable_df.groupby(["gene_name", "tissue"], as_index=False)
        .agg(
            gene_class=("gene_class", "first"),
            gene_species_count=("gene_species_count", "first"),
            species_n=("species_dir", "nunique"),
            mean_score=("expression_score", "mean"),
            median_score=("expression_score", "median"),
            std_score=("expression_score", "std"),
            min_score=("expression_score", "min"),
            max_score=("expression_score", "max"),
        )
    )
    summary_df = summary_df[summary_df["species_n"] >= 2].copy()
    summary_df["range"] = summary_df["max_score"] - summary_df["min_score"]

    idx_max = comparable_df.groupby(["gene_name", "tissue"])["expression_score"].idxmax()
    idx_min = comparable_df.groupby(["gene_name", "tissue"])["expression_score"].idxmin()

    max_df = comparable_df.loc[
        idx_max,
        ["gene_name", "tissue", "species_dir", "species_name", "expression_score"],
    ].rename(
        columns={
            "species_dir": "max_species",
            "species_name": "max_species_name",
            "expression_score": "_max_score_check",
        }
    )
    min_df = comparable_df.loc[
        idx_min,
        ["gene_name", "tissue", "species_dir", "species_name", "expression_score"],
    ].rename(
        columns={
            "species_dir": "min_species",
            "species_name": "min_species_name",
            "expression_score": "_min_score_check",
        }
    )

    summary_df = summary_df.merge(max_df, on=["gene_name", "tissue"], how="left")
    summary_df = summary_df.merge(min_df, on=["gene_name", "tissue"], how="left")
    summary_df["std_score"] = summary_df["std_score"].fillna(0.0)
    return summary_df.sort_values(
        by=["range", "species_n", "gene_name", "tissue"],
        ascending=[False, False, True, True],
    ).reset_index(drop=True)


def assign_candidate_flags(
    summary_df: pd.DataFrame,
    min_species_per_tissue: int,
    candidate_quantile: float,
    high_confidence_quantile: float,
    low_quantiles: Sequence[float],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, Dict[str, float], Dict[str, pd.DataFrame]]:
    conservative_base = summary_df[summary_df["species_n"] >= int(min_species_per_tissue)].copy()
    if conservative_base.empty:
        raise RuntimeError("No gene/tissue combinations met the conservative species filter.")

    thresholds: Dict[str, float] = {
        quantile_label(candidate_quantile): float(conservative_base["range"].quantile(candidate_quantile)),
        quantile_label(high_confidence_quantile): float(
            conservative_base["range"].quantile(high_confidence_quantile)
        ),
    }
    for low_quantile in low_quantiles:
        thresholds[quantile_label(low_quantile)] = float(conservative_base["range"].quantile(low_quantile))

    summary_df = summary_df.copy()
    summary_df["passes_conservative_species_filter"] = summary_df["species_n"] >= int(
        min_species_per_tissue
    )
    summary_df["is_conservative_candidate"] = summary_df["passes_conservative_species_filter"] & (
        summary_df["range"] >= thresholds[quantile_label(candidate_quantile)]
    )
    summary_df["is_high_confidence_candidate"] = summary_df["passes_conservative_species_filter"] & (
        summary_df["range"] >= thresholds[quantile_label(high_confidence_quantile)]
    )
    summary_df["candidate_level"] = "background"
    summary_df.loc[summary_df["is_conservative_candidate"], "candidate_level"] = "conservative"
    summary_df.loc[
        summary_df["is_high_confidence_candidate"], "candidate_level"
    ] = "high_confidence"

    low_tables: Dict[str, pd.DataFrame] = {}
    summary_df["low_variability_level"] = "background"
    for low_quantile in low_quantiles:
        label = quantile_label(low_quantile)
        col_name = f"is_low_variability_{label}"
        summary_df[col_name] = summary_df["passes_conservative_species_filter"] & (
            summary_df["range"] <= thresholds[label]
        )
    for low_quantile in low_quantiles:
        label = quantile_label(low_quantile)
        low_mask = summary_df[f"is_low_variability_{label}"]
        summary_df.loc[low_mask, "low_variability_level"] = label
        low_df = summary_df[low_mask].copy()
        low_df["low_variability_level"] = label
        low_df = low_df.sort_values(
            by=["range", "species_n", "gene_name", "tissue"],
            ascending=[True, False, True, True],
        ).reset_index(drop=True)
        low_tables[label] = low_df

    conservative_df = summary_df[summary_df["is_conservative_candidate"]].copy()
    high_conf_df = summary_df[summary_df["is_high_confidence_candidate"]].copy()
    conservative_df = conservative_df.sort_values(
        by=["range", "species_n", "gene_name", "tissue"],
        ascending=[False, False, True, True],
    ).reset_index(drop=True)
    high_conf_df = high_conf_df.sort_values(
        by=["range", "species_n", "gene_name", "tissue"],
        ascending=[False, False, True, True],
    ).reset_index(drop=True)
    return summary_df, conservative_df, high_conf_df, thresholds, low_tables


def build_pairwise_contrasts(
    species_level_long_df: pd.DataFrame,
    candidate_df: pd.DataFrame,
) -> pd.DataFrame:
    empty_cols = [
        "gene_name",
        "gene_class",
        "tissue",
        "candidate_level",
        "species_n",
        "combo_range",
        "max_species",
        "max_species_name",
        "min_species",
        "min_species_name",
        "species_high",
        "species_high_name",
        "high_score",
        "species_low",
        "species_low_name",
        "low_score",
        "abs_diff",
    ]
    if candidate_df.empty:
        return pd.DataFrame(columns=empty_cols)

    candidate_lookup = {
        (row["gene_name"], row["tissue"]): row for row in candidate_df.to_dict(orient="records")
    }
    rows: List[dict] = []
    for (gene_name, tissue), combo_df in species_level_long_df.groupby(["gene_name", "tissue"]):
        candidate_row = candidate_lookup.get((gene_name, tissue))
        if candidate_row is None:
            continue

        records = combo_df[["species_dir", "species_name", "expression_score"]].to_dict(orient="records")
        for left, right in itertools.combinations(records, 2):
            if left["expression_score"] >= right["expression_score"]:
                high = left
                low = right
            else:
                high = right
                low = left
            rows.append(
                {
                    "gene_name": gene_name,
                    "gene_class": candidate_row["gene_class"],
                    "tissue": tissue,
                    "candidate_level": candidate_row["candidate_level"],
                    "species_n": int(candidate_row["species_n"]),
                    "combo_range": float(candidate_row["range"]),
                    "max_species": candidate_row["max_species"],
                    "max_species_name": candidate_row["max_species_name"],
                    "min_species": candidate_row["min_species"],
                    "min_species_name": candidate_row["min_species_name"],
                    "species_high": high["species_dir"],
                    "species_high_name": high["species_name"],
                    "high_score": float(high["expression_score"]),
                    "species_low": low["species_dir"],
                    "species_low_name": low["species_name"],
                    "low_score": float(low["expression_score"]),
                    "abs_diff": float(high["expression_score"] - low["expression_score"]),
                }
            )

    if not rows:
        return pd.DataFrame(columns=empty_cols)

    return pd.DataFrame(rows).sort_values(
        by=["abs_diff", "combo_range", "gene_name", "tissue"],
        ascending=[False, False, True, True],
    ).reset_index(drop=True)


def build_gene_variability_summary(
    summary_df: pd.DataFrame,
    min_species_per_tissue: int,
    low_quantiles: Sequence[float],
) -> pd.DataFrame:
    conservative_summary = summary_df[summary_df["species_n"] >= int(min_species_per_tissue)].copy()
    if conservative_summary.empty:
        return pd.DataFrame()

    agg_map = {
        "gene_class": ("gene_class", "first"),
        "gene_species_count": ("gene_species_count", "first"),
        "tissues_compared": ("tissue", "count"),
        "max_range": ("range", "max"),
        "median_range": ("range", "median"),
        "mean_range": ("range", "mean"),
        "p90_range": ("range", lambda s: float(s.quantile(0.90))),
        "p10_range": ("range", lambda s: float(s.quantile(0.10))),
        "global_p90_hits": ("is_conservative_candidate", "sum"),
        "global_p95_hits": ("is_high_confidence_candidate", "sum"),
        "gene_mean_score": ("mean_score", "mean"),
        "gene_median_of_tissue_medians": ("median_score", "median"),
    }
    for low_quantile in low_quantiles:
        label = quantile_label(low_quantile)
        agg_map[f"global_{label}_hits"] = (f"is_low_variability_{label}", "sum")

    gene_df = conservative_summary.groupby("gene_name", as_index=False).agg(**agg_map)
    count_cols = [
        "gene_species_count",
        "tissues_compared",
        "global_p90_hits",
        "global_p95_hits",
    ] + [f"global_{quantile_label(q)}_hits" for q in low_quantiles]
    for col in count_cols:
        gene_df[col] = pd.to_numeric(gene_df[col], errors="coerce").fillna(0).astype(int)

    return gene_df.sort_values(
        by=["median_range", "max_range", "gene_name"],
        ascending=[True, True, True],
    ).reset_index(drop=True)


def build_gene_expression_overall_summary(species_level_long_df: pd.DataFrame) -> pd.DataFrame:
    if species_level_long_df.empty:
        return pd.DataFrame()

    gene_df = (
        species_level_long_df.groupby("gene_name", as_index=False)
        .agg(
            gene_class=("gene_class", "first"),
            gene_species_count=("gene_species_count", "first"),
            observations=("expression_score", "size"),
            species_observed=("species_dir", "nunique"),
            tissues_observed=("tissue", "nunique"),
            overall_mean_score=("expression_score", "mean"),
            overall_median_score=("expression_score", "median"),
            overall_std_score=("expression_score", "std"),
        )
    )
    gene_df["gene_species_count"] = pd.to_numeric(gene_df["gene_species_count"], errors="coerce").fillna(0).astype(int)
    gene_df["observations"] = pd.to_numeric(gene_df["observations"], errors="coerce").fillna(0).astype(int)
    gene_df["species_observed"] = pd.to_numeric(gene_df["species_observed"], errors="coerce").fillna(0).astype(int)
    gene_df["tissues_observed"] = pd.to_numeric(gene_df["tissues_observed"], errors="coerce").fillna(0).astype(int)
    gene_df["overall_std_score"] = gene_df["overall_std_score"].fillna(0.0)
    return gene_df.sort_values(by=["gene_name"], ascending=[True]).reset_index(drop=True)


def build_class_variability_summary(
    summary_df: pd.DataFrame,
    min_species_per_tissue: int,
    low_quantiles: Sequence[float],
) -> pd.DataFrame:
    class_df = summary_df[summary_df["species_n"] >= int(min_species_per_tissue)].copy()
    if class_df.empty:
        return pd.DataFrame()

    agg_map = {
        "genes": ("gene_name", "nunique"),
        "rows": ("range", "size"),
        "mean_range": ("range", "mean"),
        "median_range": ("range", "median"),
        "p05_range": ("range", lambda s: float(s.quantile(0.05))),
        "p10_range": ("range", lambda s: float(s.quantile(0.10))),
        "p90_range": ("range", lambda s: float(s.quantile(0.90))),
        "p95_range": ("range", lambda s: float(s.quantile(0.95))),
        "min_range": ("range", "min"),
        "max_range": ("range", "max"),
        "global_p90_hits": ("is_conservative_candidate", "sum"),
        "global_p95_hits": ("is_high_confidence_candidate", "sum"),
    }
    for low_quantile in low_quantiles:
        label = quantile_label(low_quantile)
        agg_map[f"global_{label}_hits"] = (f"is_low_variability_{label}", "sum")

    class_df = class_df.groupby("gene_class", as_index=False).agg(**agg_map)
    int_cols = ["genes", "rows", "global_p90_hits", "global_p95_hits"] + [
        f"global_{quantile_label(q)}_hits" for q in low_quantiles
    ]
    for col in int_cols:
        class_df[col] = pd.to_numeric(class_df[col], errors="coerce").fillna(0).astype(int)
    return class_df.sort_values(by=["gene_class"], ascending=[True]).reset_index(drop=True)


def build_species_extrema_summary(
    summary_df: pd.DataFrame,
    min_species_per_tissue: int,
    low_quantiles: Sequence[float],
) -> pd.DataFrame:
    base_df = summary_df[summary_df["species_n"] >= int(min_species_per_tissue)].copy()
    if base_df.empty:
        return pd.DataFrame()

    species = sorted(
        set(base_df["min_species"].dropna().astype(str)) | set(base_df["max_species"].dropna().astype(str))
    )
    subset_defs: List[tuple[str, pd.Series]] = [
        ("overall", pd.Series(True, index=base_df.index)),
        ("high_p90", base_df["is_conservative_candidate"]),
        ("high_p95", base_df["is_high_confidence_candidate"]),
    ]
    for low_quantile in low_quantiles:
        label = quantile_label(low_quantile)
        subset_defs.append((f"low_{label}", base_df[f"is_low_variability_{label}"]))

    rows: List[dict] = []
    for species_dir in species:
        row = {"species_dir": species_dir}
        for subset_name, mask in subset_defs:
            subset = base_df[mask].copy()
            total = int(len(subset))
            min_hits = int(subset["min_species"].eq(species_dir).sum())
            max_hits = int(subset["max_species"].eq(species_dir).sum())
            row[f"{subset_name}_rows"] = total
            row[f"{subset_name}_min_hits"] = min_hits
            row[f"{subset_name}_max_hits"] = max_hits
            row[f"{subset_name}_min_share"] = float(min_hits / total) if total else 0.0
            row[f"{subset_name}_max_share"] = float(max_hits / total) if total else 0.0
        rows.append(row)

    return pd.DataFrame(rows).sort_values(
        by=["overall_min_hits", "overall_max_hits", "species_dir"],
        ascending=[False, False, True],
    ).reset_index(drop=True)


def build_species_pair_summary(
    summary_df: pd.DataFrame,
    min_species_per_tissue: int,
    low_quantiles: Sequence[float],
) -> pd.DataFrame:
    base_df = summary_df[summary_df["species_n"] >= int(min_species_per_tissue)].copy()
    if base_df.empty:
        return pd.DataFrame()

    subset_defs: List[tuple[str, pd.Series]] = [
        ("overall", pd.Series(True, index=base_df.index)),
        ("high_p90", base_df["is_conservative_candidate"]),
        ("high_p95", base_df["is_high_confidence_candidate"]),
    ]
    for low_quantile in low_quantiles:
        label = quantile_label(low_quantile)
        subset_defs.append((f"low_{label}", base_df[f"is_low_variability_{label}"]))

    frames: List[pd.DataFrame] = []
    for subset_name, mask in subset_defs:
        subset = base_df[mask].copy()
        if subset.empty:
            continue
        pair_df = (
            subset.groupby(["min_species", "max_species"], as_index=False)
            .agg(count=("gene_name", "size"))
            .sort_values(by=["count", "min_species", "max_species"], ascending=[False, True, True])
            .reset_index(drop=True)
        )
        pair_df["subset_label"] = subset_name
        pair_df["subset_rows"] = int(len(subset))
        pair_df["share"] = pair_df["count"] / float(len(subset))
        frames.append(pair_df)

    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def build_class_low_variability_candidates(
    summary_df: pd.DataFrame,
    min_species_per_tissue: int,
    low_quantiles: Sequence[float],
) -> pd.DataFrame:
    base_df = summary_df[summary_df["species_n"] >= int(min_species_per_tissue)].copy()
    if base_df.empty:
        return pd.DataFrame()

    rows: List[pd.DataFrame] = []
    for gene_class in sorted(base_df["gene_class"].dropna().astype(str).unique().tolist()):
        class_df = base_df[base_df["gene_class"].eq(gene_class)].copy()
        if class_df.empty:
            continue
        for low_quantile in low_quantiles:
            label = quantile_label(low_quantile)
            threshold = float(class_df["range"].quantile(low_quantile))
            selected = class_df[class_df["range"] <= threshold].copy()
            selected["class_specific_low_level"] = label
            selected["class_specific_low_threshold"] = threshold
            rows.append(selected)

    if not rows:
        return pd.DataFrame()

    low_df = pd.concat(rows, ignore_index=True)
    low_df = low_df.sort_values(
        by=["gene_class", "class_specific_low_level", "range", "species_n", "gene_name", "tissue"],
        ascending=[True, True, True, False, True, True],
    ).reset_index(drop=True)
    return low_df


def take_shortlist_rows(
    high_conf_df: pd.DataFrame,
    conservative_df: pd.DataFrame,
    limit: int,
    gene_class: Optional[str] = None,
) -> pd.DataFrame:
    frames: List[pd.DataFrame] = []
    if gene_class is None:
        high_subset = high_conf_df
        cons_subset = conservative_df
    else:
        high_subset = high_conf_df[high_conf_df["gene_class"].eq(gene_class)].copy()
        cons_subset = conservative_df[conservative_df["gene_class"].eq(gene_class)].copy()
    frames.extend([high_subset, cons_subset])

    picked_keys = set()
    rows: List[dict] = []
    for frame in frames:
        for record in frame.to_dict(orient="records"):
            key = (record["gene_name"], record["tissue"])
            if key in picked_keys:
                continue
            picked_keys.add(key)
            rows.append(record)
            if len(rows) >= limit:
                return pd.DataFrame(rows)
    return pd.DataFrame(rows)


def take_ranked_summary_rows(
    summary_df: pd.DataFrame,
    limit: int,
    min_species_per_tissue: int,
    gene_class: Optional[str] = None,
) -> pd.DataFrame:
    ranked_df = summary_df[summary_df["species_n"] >= int(min_species_per_tissue)].copy()
    if gene_class is not None:
        ranked_df = ranked_df[ranked_df["gene_class"].eq(gene_class)].copy()
    ranked_df = ranked_df.sort_values(
        by=["range", "species_n", "gene_name", "tissue"],
        ascending=[False, False, True, True],
    )
    return ranked_df.head(limit).reset_index(drop=True)


def format_candidate_sentence(row: pd.Series) -> str:
    return (
        f'In {row["tissue"]}, {row["gene_name"]} expression was higher in '
        f'{row["max_species_name"]} ({row["max_score"]:.2f}) than in '
        f'{row["min_species_name"]} ({row["min_score"]:.2f}), yielding a '
        f'cross-species range of {row["range"]:.2f} across {int(row["species_n"])} species.'
    )


def write_candidate_section(lines: List[str], title: str, section_df: pd.DataFrame) -> None:
    lines.append(f"## {title}")
    if section_df.empty:
        lines.append("")
        lines.append("- No candidates met the current filtering rules.")
        lines.append("")
        return

    lines.append("")
    for row in section_df.to_dict(orient="records"):
        label = (
            f'- `{row["gene_name"]}` in `{row["tissue"]}`: '
            f'{row["min_species"]} {row["min_score"]:.2f} vs '
            f'{row["max_species"]} {row["max_score"]:.2f} '
            f'(range {row["range"]:.2f}, {int(row["species_n"])} species, '
            f'{row["candidate_level"]}).'
        )
        lines.append(label)
        lines.append(f'  Template: "{format_candidate_sentence(pd.Series(row))}"')
    lines.append("")


def write_conserved_section(lines: List[str], gene_variability_df: pd.DataFrame) -> None:
    lines.append("## Conserved genes")
    lines.append("")
    conserved_df = gene_variability_df[gene_variability_df["tissues_compared"] >= MIN_CONSERVED_TISSUES]
    conserved_df = conserved_df.head(CONSERVED_SHORTLIST_LIMIT)
    if conserved_df.empty:
        lines.append("- No genes met the minimum tissue-count rule for the conserved block.")
        lines.append("")
        return

    for row in conserved_df.to_dict(orient="records"):
        sentence = (
            f'{row["gene_name"]} showed comparatively stable expression across species '
            f'(median range {row["median_range"]:.2f}, max range {row["max_range"]:.2f}, '
            f'{int(row["tissues_compared"])} tissues).'
        )
        lines.append(
            f'- `{row["gene_name"]}` ({row["gene_class"] or "unclassified"}): '
            f'median range {row["median_range"]:.2f}, max range {row["max_range"]:.2f}, '
            f'{int(row["tissues_compared"])} tissues.'
        )
        lines.append(f'  Template: "{sentence}"')
    lines.append("")


def build_manuscript_shortlist(
    *,
    out_path: Path,
    scope: str,
    summary_df: pd.DataFrame,
    conservative_df: pd.DataFrame,
    high_conf_df: pd.DataFrame,
    gene_variability_df: pd.DataFrame,
    candidate_threshold: float,
    high_conf_threshold: float,
    min_species_per_tissue: int,
    target_gene_count: int,
) -> None:
    lines: List[str] = []
    lines.append("# Cross-species H2A manuscript shortlist")
    lines.append("")
    lines.append(
        f"- Genes analyzed: {target_gene_count} shared genes with species_count >= {DEFAULT_MIN_GENE_SPECIES_COUNT}."
    )
    lines.append(
        f"- Tissue filter for candidate selection: species_n >= {min_species_per_tissue}."
    )
    lines.append(
        f"- Generic tissues removed: {', '.join(sorted(GENERIC_TISSUES))}."
    )
    lines.append(
        f"- Conservative threshold: range >= global p90 = {candidate_threshold:.2f}."
    )
    lines.append(
        f"- High-confidence threshold: range >= global p95 = {high_conf_threshold:.2f}."
    )
    lines.append(
        f"- Candidate counts: {len(conservative_df)} conservative, {len(high_conf_df)} high-confidence."
    )
    lines.append("")

    if scope in {"two-tier", "all-shared"}:
        all_shared_df = take_shortlist_rows(
            high_conf_df,
            conservative_df,
            limit=ALL_SHARED_SHORTLIST_LIMIT,
        )
        write_candidate_section(lines, "All shared genes", all_shared_df)

    if scope in {"two-tier", "variants-only"}:
        variant_df = take_ranked_summary_rows(
            summary_df,
            limit=VARIANT_SHORTLIST_LIMIT,
            min_species_per_tissue=min_species_per_tissue,
            gene_class="variant",
        )
        write_candidate_section(lines, "Variant genes", variant_df)

    write_conserved_section(lines, gene_variability_df)
    out_path.write_text("\n".join(lines), encoding="utf-8")


def resolve_md_link(path: Path, label: str) -> str:
    return f"[{label}]({path.resolve().as_posix()})"


def maybe_resolve_existing(*paths: Path) -> Optional[Path]:
    for path in paths:
        if path.exists():
            return path
    return None


def load_panel_index(plots_dir: Path) -> pd.DataFrame:
    panel_index_path = plots_dir / "panel_membership.csv"
    if not panel_index_path.exists():
        return pd.DataFrame(
            columns=[
                "gene_name",
                "tissue",
                "panel_family",
                "panel_label",
                "panel_scope",
                "panel_direction",
                "panel_quantile",
                "page",
                "file_png",
                "file_svg",
            ]
        )
    return pd.read_csv(panel_index_path)


def case_panel_links(gene_name: str, tissue: str, panel_index_df: pd.DataFrame) -> str:
    if panel_index_df.empty:
        return "not available"
    case_df = panel_index_df[
        panel_index_df["gene_name"].eq(gene_name) & panel_index_df["tissue"].eq(tissue)
    ].copy()
    if case_df.empty:
        return "not available"

    case_df = case_df.sort_values(
        by=["panel_direction", "panel_quantile", "panel_scope", "page"],
        ascending=[True, True, True, True],
    )
    links: List[str] = []
    for row in case_df.drop_duplicates(subset=["panel_family", "page"], keep="first").to_dict(orient="records"):
        file_path = Path(str(row["file_png"]))
        if not file_path.exists():
            file_path = Path(str(row["file_svg"]))
        if not file_path.exists():
            continue
        links.append(resolve_md_link(file_path, str(row["panel_family"])))
    return ", ".join(links) if links else "not available"


def gene_heatmap_link(gene_name: str) -> str:
    gene_slug = safe_slug(gene_name)
    png_path = DEFAULT_GENE_COMPARE_FIG_ROOT / gene_slug / f"{gene_slug}_gene_compare_heatmap.png"
    svg_path = DEFAULT_GENE_COMPARE_FIG_ROOT / gene_slug / f"{gene_slug}_gene_compare_heatmap.svg"
    existing = maybe_resolve_existing(png_path, svg_path)
    if existing is None:
        return "not available"
    return resolve_md_link(existing, f"{gene_name} heatmap")


def format_float(value: object, digits: int = 2) -> str:
    if pd.isna(value):
        return "NA"
    return f"{float(value):.{digits}f}"


def format_int(value: object) -> str:
    if pd.isna(value):
        return "0"
    return str(int(value))


def markdown_table(df: pd.DataFrame, columns: Sequence[str], headers: Sequence[str]) -> List[str]:
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
    for row in df.to_dict(orient="records"):
        values = [str(row[col]) for col in columns]
        lines.append("| " + " | ".join(values) + " |")
    return lines


def build_top_gene_table(
    gene_variability_df: pd.DataFrame,
    gene_expression_df: pd.DataFrame,
    *,
    sort_by: List[str],
    ascending: List[bool],
    limit: int,
) -> pd.DataFrame:
    merged = gene_variability_df.merge(
        gene_expression_df,
        on=["gene_name", "gene_class", "gene_species_count"],
        how="left",
    )
    merged = merged.sort_values(by=sort_by, ascending=ascending).head(limit).copy()
    merged = merged.assign(
        gene=lambda df: df["gene_name"],
        gene_class=lambda df: df["gene_class"],
        tissues=lambda df: df["tissues_compared"].map(format_int),
        max_range=lambda df: df["max_range"].map(format_float),
        median_range=lambda df: df["median_range"].map(format_float),
        mean_range=lambda df: df["mean_range"].map(format_float),
        p90_hits=lambda df: df["global_p90_hits"].map(format_int),
        p95_hits=lambda df: df["global_p95_hits"].map(format_int),
        p10_hits=lambda df: df["global_p10_hits"].map(format_int),
        p05_hits=lambda df: df["global_p05_hits"].map(format_int),
        mean_expr=lambda df: df["overall_mean_score"].map(format_float),
        median_expr=lambda df: df["overall_median_score"].map(format_float),
        heatmap=lambda df: df["gene_name"].map(gene_heatmap_link),
    )
    return merged[
        [
            "gene",
            "gene_class",
            "tissues",
            "max_range",
            "median_range",
            "mean_range",
            "p90_hits",
            "p95_hits",
            "p10_hits",
            "p05_hits",
            "mean_expr",
            "median_expr",
            "heatmap",
        ]
    ]


def build_top_case_table(case_df: pd.DataFrame, panel_index_df: pd.DataFrame, limit: int) -> pd.DataFrame:
    case_df = case_df.head(limit).copy()
    case_df = case_df.assign(
        gene=lambda df: df["gene_name"],
        gene_class=lambda df: df["gene_class"],
        tissue=lambda df: df["tissue"],
        species_n=lambda df: df["species_n"].map(format_int),
        range=lambda df: df["range"].map(format_float),
        mean_score=lambda df: df["mean_score"].map(format_float),
        median_score=lambda df: df["median_score"].map(format_float),
        min_case=lambda df: df.apply(
            lambda row: f'{row["min_species"]} {format_float(row["min_score"])}', axis=1
        ),
        max_case=lambda df: df.apply(
            lambda row: f'{row["max_species"]} {format_float(row["max_score"])}', axis=1
        ),
        heatmap=lambda df: df["gene_name"].map(gene_heatmap_link),
        panels=lambda df: df.apply(
            lambda row: case_panel_links(row["gene_name"], row["tissue"], panel_index_df),
            axis=1,
        ),
    )
    return case_df[
        [
            "gene",
            "gene_class",
            "tissue",
            "species_n",
            "range",
            "mean_score",
            "median_score",
            "min_case",
            "max_case",
            "heatmap",
            "panels",
        ]
    ]


def top_species_table(
    species_extrema_df: pd.DataFrame,
    *,
    subset_prefix: str,
    mode: str,
    limit: int,
) -> pd.DataFrame:
    count_col = f"{subset_prefix}_{mode}_hits"
    share_col = f"{subset_prefix}_{mode}_share"
    rows_col = f"{subset_prefix}_rows"
    df = species_extrema_df.sort_values(
        by=[count_col, share_col, "species_dir"],
        ascending=[False, False, True],
    ).head(limit).copy()
    df = df.assign(
        species=lambda d: d["species_dir"],
        hits=lambda d: d[count_col].map(format_int),
        share=lambda d: d[share_col].map(lambda x: format_float(100 * x)),
        rows=lambda d: d[rows_col].map(format_int),
    )
    return df[["species", "hits", "share", "rows"]]


def top_pair_table(species_pair_df: pd.DataFrame, subset_label: str, limit: int) -> pd.DataFrame:
    subset_df = species_pair_df[species_pair_df["subset_label"].eq(subset_label)].copy()
    subset_df = subset_df.sort_values(
        by=["count", "share", "min_species", "max_species"],
        ascending=[False, False, True, True],
    ).head(limit)
    subset_df = subset_df.assign(
        pair=lambda d: d["min_species"] + " -> " + d["max_species"],
        count=lambda d: d["count"].map(format_int),
        share=lambda d: d["share"].map(lambda x: format_float(100 * x)),
        rows=lambda d: d["subset_rows"].map(format_int),
    )
    return subset_df[["pair", "count", "share", "rows"]]


def overall_stats_bullet(class_variability_df: pd.DataFrame, gene_class: str) -> Optional[str]:
    row_df = class_variability_df[class_variability_df["gene_class"].eq(gene_class)].copy()
    if row_df.empty:
        return None
    row = row_df.iloc[0]
    return (
        f"- `{gene_class}`: {int(row['genes'])} genes, {int(row['rows'])} comparable gene*tissue rows, "
        f"median range {row['median_range']:.2f}, mean range {row['mean_range']:.2f}, "
        f"p90 {row['p90_range']:.2f}, p95 {row['p95_range']:.2f}, max range {row['max_range']:.2f}."
    )


def interpretation_bullets(
    class_variability_df: pd.DataFrame,
    gene_variability_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    class_low_df: pd.DataFrame,
    species_extrema_df: pd.DataFrame,
    panel_index_df: pd.DataFrame,
) -> List[str]:
    bullets: List[str] = []
    clustered_df = class_variability_df[class_variability_df["gene_class"].eq("clustered")]
    variant_df = class_variability_df[class_variability_df["gene_class"].eq("variant")]
    if not clustered_df.empty and not variant_df.empty:
        c_row = clustered_df.iloc[0]
        v_row = variant_df.iloc[0]
        bullets.append(
            f"- В целом clustered-гены варьируют сильнее, чем variant-гены: медиана range {c_row['median_range']:.2f} против {v_row['median_range']:.2f}, а p95 {c_row['p95_range']:.2f} против {v_row['p95_range']:.2f}."
        )

    top_max_gene = gene_variability_df.sort_values(
        by=["max_range", "median_range", "gene_name"], ascending=[False, False, True]
    ).head(1)
    if not top_max_gene.empty:
        row = top_max_gene.iloc[0]
        bullets.append(
            f"- `{row['gene_name']}` имеет самый большой наблюдаемый cross-species range на уровне gene*tissue: максимальный range {row['max_range']:.2f}, медианный range по тканям {row['median_range']:.2f}."
        )

    top_var_variant = gene_variability_df[gene_variability_df["gene_class"].eq("variant")].sort_values(
        by=["max_range", "median_range", "gene_name"], ascending=[False, False, True]
    ).head(1)
    if not top_var_variant.empty:
        row = top_var_variant.iloc[0]
        bullets.append(
            f"- Среди variant-генов наиболее сильную вариабельность показывает `{row['gene_name']}`: максимальный range {row['max_range']:.2f}, медианный range {row['median_range']:.2f}."
        )

    conserved_df = gene_variability_df.sort_values(
        by=["median_range", "max_range", "gene_name"], ascending=[True, True, True]
    ).head(3)
    if not conserved_df.empty:
        genes = ", ".join(
            f"`{row.gene_name}` (median range {row.median_range:.2f})" for row in conserved_df.itertuples()
        )
        bullets.append(
            f"- Наиболее консервативный блок формируют {genes}; эти гены системно попадают в low-tail по межвидовому разбросу."
        )

    top_variant_cases = summary_df[summary_df["gene_class"].eq("variant")].sort_values(
        by=["range", "species_n", "gene_name", "tissue"], ascending=[False, False, True, True]
    ).head(3)
    for row in top_variant_cases.to_dict(orient="records"):
        panels = case_panel_links(row["gene_name"], row["tissue"], panel_index_df)
        bullets.append(
            f"- Вариантный кейс `{row['gene_name']}` / `{row['tissue']}`: range {row['range']:.2f}, "
            f"среднее {row['mean_score']:.2f}, медиана {row['median_score']:.2f}; минимум у {row['min_species']} ({row['min_score']:.2f}), максимум у {row['max_species']} ({row['max_score']:.2f}). Panels: {panels}."
        )

    high_min_df = top_species_table(species_extrema_df, subset_prefix="high_p90", mode="min", limit=2)
    high_max_df = top_species_table(species_extrema_df, subset_prefix="high_p90", mode="max", limit=2)
    if not high_min_df.empty and not high_max_df.empty:
        low_species = ", ".join(f"{row.species} ({row.hits} hits)" for row in high_min_df.itertuples())
        high_species = ", ".join(f"{row.species} ({row.hits} hits)" for row in high_max_df.itertuples())
        bullets.append(
            f"- В high-tail (`p90`) нижние значения чаще всего дают {low_species}, тогда как верхние значения чаще наблюдаются у {high_species}."
        )

    if not class_low_df.empty:
        clustered_low = class_low_df[
            class_low_df["gene_class"].eq("clustered") & class_low_df["class_specific_low_level"].eq("p10")
        ].head(1)
        if not clustered_low.empty:
            row = clustered_low.iloc[0]
            bullets.append(
                f"- Даже внутри clustered-класса есть сравнительно стабильные сочетания, например `{row['gene_name']}` / `{row['tissue']}` (class-specific {row['class_specific_low_level']}, range {row['range']:.2f})."
            )

    return bullets


def build_detailed_report(
    *,
    out_path: Path,
    report_language: str,
    min_species_per_tissue: int,
    target_gene_count: int,
    summary_df: pd.DataFrame,
    high_conf_df: pd.DataFrame,
    low_tables: Dict[str, pd.DataFrame],
    class_low_df: pd.DataFrame,
    gene_variability_df: pd.DataFrame,
    gene_expression_df: pd.DataFrame,
    class_variability_df: pd.DataFrame,
    species_extrema_df: pd.DataFrame,
    species_pair_df: pd.DataFrame,
    thresholds: Dict[str, float],
    panel_index_df: pd.DataFrame,
) -> None:
    if report_language != "ru":
        raise ValueError("Currently only --report-language ru is supported.")

    summary_base = summary_df[summary_df["species_n"] >= int(min_species_per_tissue)].copy()
    lines: List[str] = []
    lines.append("# Подробный отчёт по cross-species сравнению экспрессии H2A на уровне gene*tissue")
    lines.append("")
    lines.append("## Overview")
    lines.append("")
    lines.append(
        f"- В анализ включены {target_gene_count} shared H2A genes с `species_count >= {DEFAULT_MIN_GENE_SPECIES_COUNT}`."
    )
    lines.append(
        f"- Основной фильтр для интерпретации gene*tissue: `species_n >= {min_species_per_tissue}`."
    )
    lines.append(
        f"- После удаления generic tissues осталось {len(summary_base)} сравнимых gene*tissue сочетаний."
    )
    lines.append(f"- Исключённые generic tissues: {', '.join(sorted(GENERIC_TISSUES))}.")
    lines.append(
        f"- Глобальные пороги: p05 = {thresholds['p05']:.2f}, p10 = {thresholds['p10']:.2f}, p90 = {thresholds['p90']:.2f}, p95 = {thresholds['p95']:.2f}."
    )
    lines.append("")

    lines.append("## Class-level patterns")
    lines.append("")
    class_bullets = [overall_stats_bullet(class_variability_df, gene_class) for gene_class in ["clustered", "variant"]]
    lines.extend([bullet for bullet in class_bullets if bullet])
    lines.append("")
    class_table = class_variability_df.copy()
    class_table = class_table.assign(
        genes=class_table["genes"].map(format_int),
        rows=class_table["rows"].map(format_int),
        mean_range=class_table["mean_range"].map(format_float),
        median_range=class_table["median_range"].map(format_float),
        p90_range=class_table["p90_range"].map(format_float),
        p95_range=class_table["p95_range"].map(format_float),
        max_range=class_table["max_range"].map(format_float),
        global_p90_hits=class_table["global_p90_hits"].map(format_int),
        global_p95_hits=class_table["global_p95_hits"].map(format_int),
        global_p10_hits=class_table["global_p10_hits"].map(format_int),
        global_p05_hits=class_table["global_p05_hits"].map(format_int),
    )
    lines.extend(
        markdown_table(
            class_table,
            [
                "gene_class",
                "genes",
                "rows",
                "median_range",
                "mean_range",
                "p90_range",
                "p95_range",
                "max_range",
                "global_p90_hits",
                "global_p95_hits",
                "global_p10_hits",
                "global_p05_hits",
            ],
            [
                "Class",
                "Genes",
                "Rows",
                "Median range",
                "Mean range",
                "p90",
                "p95",
                "Max",
                "Global p90 hits",
                "Global p95 hits",
                "Global p10 hits",
                "Global p05 hits",
            ],
        )
    )
    lines.append("")

    lines.append("## Most variable genes")
    lines.append("")
    lines.append("- Рейтинг по максимальному наблюдаемому range:")
    lines.append("")
    lines.extend(
        markdown_table(
            build_top_gene_table(
                gene_variability_df,
                gene_expression_df,
                sort_by=["max_range", "median_range", "gene_name"],
                ascending=[False, False, True],
                limit=TOP_VARIABLE_GENE_LIMIT,
            ),
            [
                "gene",
                "gene_class",
                "tissues",
                "max_range",
                "median_range",
                "mean_range",
                "p90_hits",
                "p95_hits",
                "mean_expr",
                "median_expr",
                "heatmap",
            ],
            [
                "Gene",
                "Class",
                "Tissues",
                "Max range",
                "Median range",
                "Mean range",
                "p90 hits",
                "p95 hits",
                "Overall mean expr",
                "Overall median expr",
                "Heatmap",
            ],
        )
    )
    lines.append("")
    lines.append("- Рейтинг по медианному range across tissues:")
    lines.append("")
    lines.extend(
        markdown_table(
            build_top_gene_table(
                gene_variability_df,
                gene_expression_df,
                sort_by=["median_range", "max_range", "gene_name"],
                ascending=[False, False, True],
                limit=TOP_VARIABLE_GENE_LIMIT,
            ),
            [
                "gene",
                "gene_class",
                "tissues",
                "median_range",
                "max_range",
                "mean_range",
                "p90_hits",
                "p95_hits",
                "mean_expr",
                "median_expr",
                "heatmap",
            ],
            [
                "Gene",
                "Class",
                "Tissues",
                "Median range",
                "Max range",
                "Mean range",
                "p90 hits",
                "p95 hits",
                "Overall mean expr",
                "Overall median expr",
                "Heatmap",
            ],
        )
    )
    lines.append("")

    lines.append("## Most conserved genes")
    lines.append("")
    lines.extend(
        markdown_table(
            build_top_gene_table(
                gene_variability_df,
                gene_expression_df,
                sort_by=["median_range", "max_range", "gene_name"],
                ascending=[True, True, True],
                limit=TOP_CONSERVED_GENE_LIMIT,
            ),
            [
                "gene",
                "gene_class",
                "tissues",
                "median_range",
                "max_range",
                "mean_range",
                "p10_hits",
                "p05_hits",
                "mean_expr",
                "median_expr",
                "heatmap",
            ],
            [
                "Gene",
                "Class",
                "Tissues",
                "Median range",
                "Max range",
                "Mean range",
                "p10 hits",
                "p05 hits",
                "Overall mean expr",
                "Overall median expr",
                "Heatmap",
            ],
        )
    )
    lines.append("")

    lines.append("## Most variable gene*tissue combinations")
    lines.append("")
    lines.append("- Global p95 cases:")
    lines.append("")
    lines.extend(
        markdown_table(
            build_top_case_table(high_conf_df, panel_index_df, TOP_CASE_LIMIT),
            [
                "gene",
                "gene_class",
                "tissue",
                "species_n",
                "range",
                "mean_score",
                "median_score",
                "min_case",
                "max_case",
                "heatmap",
                "panels",
            ],
            [
                "Gene",
                "Class",
                "Tissue",
                "Species n",
                "Range",
                "Mean score",
                "Median score",
                "Min species",
                "Max species",
                "Heatmap",
                "Panels",
            ],
        )
    )
    lines.append("")
    lines.append("- Variant-biased high-tail cases:")
    lines.append("")
    variant_high_df = summary_base[summary_base["gene_class"].eq("variant")].sort_values(
        by=["range", "species_n", "gene_name", "tissue"], ascending=[False, False, True, True]
    )
    lines.extend(
        markdown_table(
            build_top_case_table(variant_high_df, panel_index_df, TOP_CASE_LIMIT),
            [
                "gene",
                "tissue",
                "species_n",
                "range",
                "mean_score",
                "median_score",
                "min_case",
                "max_case",
                "heatmap",
                "panels",
            ],
            [
                "Gene",
                "Tissue",
                "Species n",
                "Range",
                "Mean score",
                "Median score",
                "Min species",
                "Max species",
                "Heatmap",
                "Panels",
            ],
        )
    )
    lines.append("")

    lines.append("## Most conserved gene*tissue combinations")
    lines.append("")
    for label in ["p05", "p10"]:
        low_df = low_tables.get(label, pd.DataFrame())
        if low_df.empty:
            continue
        lines.append(f"- Global low-tail `{label}` cases:")
        lines.append("")
        lines.extend(
            markdown_table(
                build_top_case_table(low_df, panel_index_df, TOP_CASE_LIMIT),
                [
                    "gene",
                    "gene_class",
                    "tissue",
                    "species_n",
                    "range",
                    "mean_score",
                    "median_score",
                    "min_case",
                    "max_case",
                    "heatmap",
                    "panels",
                ],
                [
                    "Gene",
                    "Class",
                    "Tissue",
                    "Species n",
                    "Range",
                    "Mean score",
                    "Median score",
                    "Min species",
                    "Max species",
                    "Heatmap",
                    "Panels",
                ],
            )
        )
        lines.append("")

    if not class_low_df.empty:
        lines.append("- Class-specific low-tail cases:")
        lines.append("")
        class_low_table = class_low_df.head(TOP_CASE_LIMIT).copy()
        class_low_table = class_low_table.assign(
            gene=class_low_table["gene_name"],
            class_low=class_low_table["class_specific_low_level"],
            species_n=class_low_table["species_n"].map(format_int),
            range=class_low_table["range"].map(format_float),
            mean_score=class_low_table["mean_score"].map(format_float),
            median_score=class_low_table["median_score"].map(format_float),
            min_case=class_low_table.apply(
                lambda row: f'{row["min_species"]} {format_float(row["min_score"])}', axis=1
            ),
            max_case=class_low_table.apply(
                lambda row: f'{row["max_species"]} {format_float(row["max_score"])}', axis=1
            ),
            heatmap=class_low_table["gene_name"].map(gene_heatmap_link),
            panels=class_low_table.apply(
                lambda row: case_panel_links(row["gene_name"], row["tissue"], panel_index_df),
                axis=1,
            ),
        )
        lines.extend(
            markdown_table(
                class_low_table,
                [
                    "gene",
                    "gene_class",
                    "class_low",
                    "tissue",
                    "species_n",
                    "range",
                    "mean_score",
                    "median_score",
                    "min_case",
                    "max_case",
                    "heatmap",
                    "panels",
                ],
                [
                    "Gene",
                    "Class",
                    "Class low",
                    "Tissue",
                    "Species n",
                    "Range",
                    "Mean score",
                    "Median score",
                    "Min species",
                    "Max species",
                    "Heatmap",
                    "Panels",
                ],
            )
        )
        lines.append("")

    lines.append("## Species tendencies")
    lines.append("")
    species_sections = [
        ("overall", "min", "Наиболее частые `min_species` в полном comparable set"),
        ("overall", "max", "Наиболее частые `max_species` в полном comparable set"),
        ("high_p90", "min", "Наиболее частые `min_species` в high-tail p90"),
        ("high_p90", "max", "Наиболее частые `max_species` в high-tail p90"),
        ("low_p10", "min", "Наиболее частые `min_species` в low-tail p10"),
        ("low_p10", "max", "Наиболее частые `max_species` в low-tail p10"),
    ]
    for subset_prefix, mode, title in species_sections:
        table_df = top_species_table(
            species_extrema_df,
            subset_prefix=subset_prefix,
            mode=mode,
            limit=TOP_SPECIES_LIMIT,
        )
        if table_df.empty:
            continue
        lines.append(f"- {title}:")
        lines.append("")
        lines.extend(
            markdown_table(
                table_df,
                ["species", "hits", "share", "rows"],
                ["Species", "Hits", "Share %", "Subset rows"],
            )
        )
        lines.append("")

    if not species_pair_df.empty:
        lines.append("- Наиболее частые пары `min_species -> max_species`:")
        lines.append("")
        for subset_label in ["overall", "high_p90", "low_p10"]:
            pair_table = top_pair_table(species_pair_df, subset_label, TOP_PAIR_LIMIT)
            if pair_table.empty:
                continue
            lines.append(f"  - `{subset_label}`:")
            lines.append("")
            lines.extend(
                markdown_table(
                    pair_table,
                    ["pair", "count", "share", "rows"],
                    ["Pair", "Count", "Share %", "Subset rows"],
                )
            )
            lines.append("")

    lines.append("## Interpretation-ready findings")
    lines.append("")
    for bullet in interpretation_bullets(
        class_variability_df,
        gene_variability_df,
        summary_base,
        class_low_df,
        species_extrema_df,
        panel_index_df,
    ):
        lines.append(bullet)
    lines.append("")

    lines.append("## Caveats")
    lines.append("")
    lines.append("- Все сопоставления тканей основаны на exact string matching по `Anatomical entity name`; синонимы и онтологическое слияние не применялись.")
    lines.append("- Основной narrative строится только на сочетаниях с `species_n >= 4`, поэтому менее покрытые ткани в отчёт не включены как ключевые доказательства.")
    lines.append("- Species-level выводы следует трактовать как статистические tendencies внутри текущего набора сравнимых тканей, а не как абсолютные свойства вида.")
    lines.append("- Общие средние и медианы по гену рассчитаны по всем доступным `species × tissue` наблюдениям данного гена после удаления generic tissues.")
    lines.append("")

    out_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    if not 0 < args.candidate_quantile < 1:
        raise ValueError("--candidate-quantile must be between 0 and 1")
    if not 0 < args.high_confidence_quantile < 1:
        raise ValueError("--high-confidence-quantile must be between 0 and 1")
    if args.candidate_quantile > args.high_confidence_quantile:
        raise ValueError("--candidate-quantile must be <= --high-confidence-quantile")

    low_quantiles = normalize_low_quantiles(args.low_quantiles)
    shared_index_path = Path(args.shared_index)
    detail_index_path = Path(args.detail_index)
    heatmap_dir = Path(args.heatmap_dir)
    processed_dir = Path(args.processed_dir)
    out_dir = Path(args.out_dir)
    plots_dir = Path(args.plots_dir)
    tables_dir = out_dir / "tables"
    reports_dir = out_dir / "reports"
    tables_dir.mkdir(parents=True, exist_ok=True)
    reports_dir.mkdir(parents=True, exist_ok=True)

    detail_df = load_detail_index(detail_index_path, heatmap_dir=heatmap_dir, processed_dir=processed_dir)
    shared_df = load_shared_index(shared_index_path, detail_df, detail_index_path)
    target_genes = choose_target_genes(shared_df, DEFAULT_MIN_GENE_SPECIES_COUNT)
    if not target_genes:
        raise RuntimeError("No shared genes met the conservative gene-count filter.")

    species_level_long_df = build_species_level_long(detail_df, target_genes)
    summary_df = build_gene_tissue_summary(species_level_long_df)
    summary_df, conservative_df, high_conf_df, thresholds, low_tables = assign_candidate_flags(
        summary_df,
        min_species_per_tissue=args.min_species_per_tissue,
        candidate_quantile=args.candidate_quantile,
        high_confidence_quantile=args.high_confidence_quantile,
        low_quantiles=low_quantiles,
    )
    candidate_union_df = pd.concat([conservative_df, high_conf_df], ignore_index=True)
    candidate_union_df = candidate_union_df.drop_duplicates(subset=["gene_name", "tissue"], keep="last").reset_index(drop=True)
    pairwise_df = build_pairwise_contrasts(species_level_long_df, candidate_union_df)
    gene_variability_df = build_gene_variability_summary(summary_df, args.min_species_per_tissue, low_quantiles)
    gene_expression_df = build_gene_expression_overall_summary(species_level_long_df)
    class_variability_df = build_class_variability_summary(summary_df, args.min_species_per_tissue, low_quantiles)
    species_extrema_df = build_species_extrema_summary(summary_df, args.min_species_per_tissue, low_quantiles)
    species_pair_df = build_species_pair_summary(summary_df, args.min_species_per_tissue, low_quantiles)
    class_low_df = build_class_low_variability_candidates(summary_df, args.min_species_per_tissue, low_quantiles) if args.class_specific_low_quantiles else pd.DataFrame()

    summary_df.to_csv(tables_dir / "gene_tissue_summary.csv", index=False, encoding="utf-8")
    conservative_df.to_csv(tables_dir / "conservative_candidates.csv", index=False, encoding="utf-8")
    high_conf_df.to_csv(tables_dir / "high_confidence_candidates.csv", index=False, encoding="utf-8")
    pairwise_df.to_csv(tables_dir / "pairwise_contrasts.csv", index=False, encoding="utf-8")
    gene_variability_df.to_csv(tables_dir / "gene_variability_summary.csv", index=False, encoding="utf-8")
    gene_expression_df.to_csv(tables_dir / "gene_expression_overall_summary.csv", index=False, encoding="utf-8")
    class_variability_df.to_csv(tables_dir / "class_variability_summary.csv", index=False, encoding="utf-8")
    species_extrema_df.to_csv(tables_dir / "species_extrema_summary.csv", index=False, encoding="utf-8")
    species_pair_df.to_csv(tables_dir / "species_pair_summary.csv", index=False, encoding="utf-8")
    if args.class_specific_low_quantiles:
        class_low_df.to_csv(tables_dir / "class_low_variability_candidates.csv", index=False, encoding="utf-8")
    for label, low_df in low_tables.items():
        low_df.to_csv(tables_dir / f"low_variability_candidates_{label}.csv", index=False, encoding="utf-8")

    out_shortlist_md = reports_dir / "manuscript_shortlist.md"
    build_manuscript_shortlist(
        out_path=out_shortlist_md,
        scope=args.scope,
        summary_df=summary_df,
        conservative_df=conservative_df,
        high_conf_df=high_conf_df,
        gene_variability_df=gene_variability_df,
        candidate_threshold=thresholds[quantile_label(args.candidate_quantile)],
        high_conf_threshold=thresholds[quantile_label(args.high_confidence_quantile)],
        min_species_per_tissue=args.min_species_per_tissue,
        target_gene_count=len(target_genes),
    )

    if args.write_detailed_report:
        build_detailed_report(
            out_path=reports_dir / f"{args.report_stem}.md",
            report_language=args.report_language,
            min_species_per_tissue=args.min_species_per_tissue,
            target_gene_count=len(target_genes),
            summary_df=summary_df,
            high_conf_df=high_conf_df,
            low_tables=low_tables,
            class_low_df=class_low_df,
            gene_variability_df=gene_variability_df,
            gene_expression_df=gene_expression_df,
            class_variability_df=class_variability_df,
            species_extrema_df=species_extrema_df,
            species_pair_df=species_pair_df,
            thresholds=thresholds,
            panel_index_df=load_panel_index(plots_dir),
        )


if __name__ == "__main__":
    main()
