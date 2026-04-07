#!/usr/bin/env python
"""Build statistics and a difference report for the human HPA vs Bgee compare layer."""

from __future__ import annotations

import argparse
import json
from datetime import datetime
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

import numpy as np
import pandas as pd

BIOANALYZE_SCRIPTS_ROOT = Path(__file__).resolve().parents[1]
if str(BIOANALYZE_SCRIPTS_ROOT) not in sys.path:
    sys.path.insert(0, str(BIOANALYZE_SCRIPTS_ROOT))

from bioanalyze_paths import get_bioanalyze_data_root, get_bioanalyze_figures_root, get_bioanalyze_stats_root


DEFAULT_HPA_TSV = get_bioanalyze_data_root() / "processed" / "intersections" / "human_hpa_bgee" / "hpa_vs_bgee_hpa_aligned_long.tsv"
DEFAULT_BGEE_TSV = get_bioanalyze_data_root() / "processed" / "intersections" / "human_hpa_bgee" / "hpa_vs_bgee_bgee_aligned_long.tsv"
DEFAULT_COMPARE_METADATA = get_bioanalyze_data_root() / "processed" / "intersections" / "human_hpa_bgee" / "hpa_vs_bgee_metadata.json"
DEFAULT_OUT_DIR = get_bioanalyze_stats_root() / "compare_nTPM_bgee"
DEFAULT_TOP_N = 20

MERGE_KEYS = [
    "gene_order",
    "gene_number",
    "tissue_order",
    "tissue_letter",
    "canonical_gene_name",
    "compare_tissue_key",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build statistics and a markdown report describing agreement and "
            "differences between the human HPA-vs-Bgee compare heatmaps."
        )
    )
    parser.add_argument("--hpa-tsv", default=str(DEFAULT_HPA_TSV))
    parser.add_argument("--bgee-tsv", default=str(DEFAULT_BGEE_TSV))
    parser.add_argument("--compare-metadata", default=str(DEFAULT_COMPARE_METADATA))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument("--top-n", default=DEFAULT_TOP_N, type=int)
    return parser.parse_args()


def require_columns(df: pd.DataFrame, columns: Sequence[str], *, label: str) -> None:
    missing = [column for column in columns if column not in df.columns]
    if missing:
        raise RuntimeError(f"{label} is missing required columns: {', '.join(missing)}")


def load_hpa_df(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path, sep="\t", low_memory=True)
    require_columns(
        df,
        MERGE_KEYS
        + [
            "ensembl_gene_id",
            "hpa_gene_label",
            "hpa_tissue_label",
            "bgee_tissue_label",
            "match_type",
            "is_underlined",
            "nTPM",
        ],
        label="HPA aligned TSV",
    )
    df["nTPM"] = pd.to_numeric(df["nTPM"], errors="coerce")
    df["is_underlined"] = df["is_underlined"].fillna(False).astype(bool)
    return df


def load_bgee_df(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path, sep="\t", low_memory=True)
    require_columns(
        df,
        MERGE_KEYS
        + [
            "ensembl_gene_id",
            "bgee_gene_label",
            "hpa_tissue_label",
            "bgee_tissue_label",
            "match_type",
            "is_underlined",
            "cell_mean_score",
        ],
        label="Bgee aligned TSV",
    )
    df["cell_mean_score"] = pd.to_numeric(df["cell_mean_score"], errors="coerce")
    df["is_underlined"] = df["is_underlined"].fillna(False).astype(bool)
    return df


def load_optional_metadata(path: Path) -> Dict[str, object]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def _validate_equal_columns(merged: pd.DataFrame, base_column: str) -> pd.DataFrame:
    left_col = f"{base_column}_hpa"
    right_col = f"{base_column}_bgee"
    if left_col not in merged.columns or right_col not in merged.columns:
        return merged

    left = merged[left_col]
    right = merged[right_col]
    if left.fillna("").astype(str).ne(right.fillna("").astype(str)).any():
        raise RuntimeError(f"Merged compare tables disagree on column {base_column}.")
    merged[base_column] = left
    return merged.drop(columns=[left_col, right_col])


def build_paired_df(hpa_df: pd.DataFrame, bgee_df: pd.DataFrame) -> pd.DataFrame:
    merged = hpa_df.merge(
        bgee_df,
        on=MERGE_KEYS,
        how="inner",
        suffixes=("_hpa", "_bgee"),
    )
    if merged.empty:
        raise RuntimeError("Merged compare table is empty.")

    for column in ["ensembl_gene_id", "hpa_tissue_label", "bgee_tissue_label", "match_type", "is_underlined"]:
        merged = _validate_equal_columns(merged, column)

    merged["nTPM"] = pd.to_numeric(merged["nTPM"], errors="coerce")
    merged["cell_mean_score"] = pd.to_numeric(merged["cell_mean_score"], errors="coerce")
    merged["paired_is_missing"] = merged["nTPM"].isna() | merged["cell_mean_score"].isna()

    paired = merged.loc[~merged["paired_is_missing"]].copy()
    paired["log_hpa"] = np.log10(paired["nTPM"] + 1.0)
    paired["log_bgee"] = np.log10(paired["cell_mean_score"] + 1.0)
    paired["signed_diff_log"] = paired["log_hpa"] - paired["log_bgee"]
    paired["abs_diff_log"] = paired["signed_diff_log"].abs()
    paired["signed_diff_raw"] = paired["nTPM"] - paired["cell_mean_score"]
    paired["abs_diff_raw"] = paired["signed_diff_raw"].abs()
    paired["dominance_label"] = np.select(
        [
            np.isclose(paired["signed_diff_raw"], 0.0),
            paired["signed_diff_raw"] > 0.0,
            paired["signed_diff_raw"] < 0.0,
        ],
        ["equal", "hpa_higher", "bgee_higher"],
        default="equal",
    )
    paired["pair_status"] = np.where(
        (paired["nTPM"] > 0.0) | (paired["cell_mean_score"] > 0.0),
        "nonzero_any",
        "both_zero",
    )

    ordered_columns = [
        "gene_order",
        "gene_number",
        "canonical_gene_name",
        "ensembl_gene_id",
        "hpa_gene_label",
        "bgee_gene_label",
        "tissue_order",
        "tissue_letter",
        "compare_tissue_key",
        "hpa_tissue_label",
        "bgee_tissue_label",
        "match_type",
        "is_underlined",
        "nTPM",
        "cell_mean_score",
        "log_hpa",
        "log_bgee",
        "signed_diff_log",
        "abs_diff_log",
        "signed_diff_raw",
        "abs_diff_raw",
        "dominance_label",
        "pair_status",
    ]
    return paired[ordered_columns].copy()


def summarize_correlations(df: pd.DataFrame, *, subset_label: str, subset_description: str) -> Dict[str, object]:
    result: Dict[str, object] = {
        "subset_label": subset_label,
        "subset_description": subset_description,
        "pair_count": int(len(df)),
        "hpa_zero_count": int((df["nTPM"] <= 0.0).sum()),
        "bgee_zero_count": int((df["cell_mean_score"] <= 0.0).sum()),
    }
    if len(df) < 2:
        result.update(
            {
                "pearson_log": np.nan,
                "spearman_log": np.nan,
                "pearson_raw": np.nan,
                "spearman_raw": np.nan,
            }
        )
        return result

    result.update(
        {
            "pearson_log": float(df["log_hpa"].corr(df["log_bgee"], method="pearson")),
            "spearman_log": float(df["log_hpa"].corr(df["log_bgee"], method="spearman")),
            "pearson_raw": float(df["nTPM"].corr(df["cell_mean_score"], method="pearson")),
            "spearman_raw": float(df["nTPM"].corr(df["cell_mean_score"], method="spearman")),
        }
    )
    return result


def build_correlation_summary(paired_df: pd.DataFrame) -> pd.DataFrame:
    rows = [
        summarize_correlations(
            paired_df,
            subset_label="all_non_missing",
            subset_description="All paired cells with non-missing values in both sources.",
        ),
        summarize_correlations(
            paired_df.loc[(paired_df["nTPM"] > 0.0) | (paired_df["cell_mean_score"] > 0.0)].copy(),
            subset_label="nonzero_any",
            subset_description="Paired cells where either HPA or Bgee is > 0.",
        ),
    ]
    return pd.DataFrame(rows)


def _top_rows(df: pd.DataFrame, *, top_n: int, subset_label: str) -> pd.DataFrame:
    ranked = (
        df.sort_values(
            by=["abs_diff_log", "abs_diff_raw", "canonical_gene_name", "tissue_order"],
            ascending=[False, False, True, True],
            kind="stable",
        )
        .head(top_n)
        .copy()
    )
    ranked.insert(0, "subset_label", subset_label)
    ranked.insert(1, "rank_within_subset", np.arange(1, len(ranked) + 1, dtype=int))
    return ranked


def build_top_differing_cells(paired_df: pd.DataFrame, *, top_n: int) -> pd.DataFrame:
    frames = [
        _top_rows(paired_df, top_n=top_n, subset_label="overall"),
        _top_rows(
            paired_df.loc[paired_df["dominance_label"].eq("hpa_higher")].copy(),
            top_n=top_n,
            subset_label="hpa_higher",
        ),
        _top_rows(
            paired_df.loc[paired_df["dominance_label"].eq("bgee_higher")].copy(),
            top_n=top_n,
            subset_label="bgee_higher",
        ),
    ]
    return pd.concat(frames, ignore_index=True)


def build_gene_difference_summary(paired_df: pd.DataFrame) -> pd.DataFrame:
    grouped = (
        paired_df.groupby(
            ["gene_order", "gene_number", "canonical_gene_name", "hpa_gene_label", "bgee_gene_label"],
            as_index=False,
        )
        .agg(
            comparable_cell_count=("compare_tissue_key", "size"),
            mean_abs_log_diff=("abs_diff_log", "mean"),
            median_abs_log_diff=("abs_diff_log", "median"),
            max_abs_log_diff=("abs_diff_log", "max"),
            mean_abs_raw_diff=("abs_diff_raw", "mean"),
            median_abs_raw_diff=("abs_diff_raw", "median"),
            max_abs_raw_diff=("abs_diff_raw", "max"),
            safe_synonym_cell_count=("is_underlined", "sum"),
        )
    )

    idx = paired_df.groupby("canonical_gene_name")["abs_diff_log"].idxmax()
    top_rows = paired_df.loc[idx, [
        "canonical_gene_name",
        "compare_tissue_key",
        "hpa_tissue_label",
        "bgee_tissue_label",
        "match_type",
        "is_underlined",
        "abs_diff_log",
        "abs_diff_raw",
        "dominance_label",
    ]].rename(
        columns={
            "compare_tissue_key": "top_compare_tissue_key",
            "hpa_tissue_label": "top_hpa_tissue_label",
            "bgee_tissue_label": "top_bgee_tissue_label",
            "match_type": "top_match_type",
            "is_underlined": "top_is_underlined",
            "abs_diff_log": "top_abs_diff_log",
            "abs_diff_raw": "top_abs_diff_raw",
            "dominance_label": "top_dominance_label",
        }
    )
    summary = grouped.merge(top_rows, on="canonical_gene_name", how="left")
    return summary.sort_values(
        by=["mean_abs_log_diff", "max_abs_log_diff", "canonical_gene_name"],
        ascending=[False, False, True],
        kind="stable",
    ).reset_index(drop=True)


def build_tissue_difference_summary(paired_df: pd.DataFrame) -> pd.DataFrame:
    grouped = (
        paired_df.groupby(
            [
                "tissue_order",
                "tissue_letter",
                "compare_tissue_key",
                "hpa_tissue_label",
                "bgee_tissue_label",
                "match_type",
                "is_underlined",
            ],
            as_index=False,
        )
        .agg(
            comparable_cell_count=("canonical_gene_name", "size"),
            mean_abs_log_diff=("abs_diff_log", "mean"),
            median_abs_log_diff=("abs_diff_log", "median"),
            max_abs_log_diff=("abs_diff_log", "max"),
            mean_abs_raw_diff=("abs_diff_raw", "mean"),
            median_abs_raw_diff=("abs_diff_raw", "median"),
            max_abs_raw_diff=("abs_diff_raw", "max"),
        )
    )

    idx = paired_df.groupby("compare_tissue_key")["abs_diff_log"].idxmax()
    top_rows = paired_df.loc[idx, [
        "compare_tissue_key",
        "canonical_gene_name",
        "hpa_gene_label",
        "bgee_gene_label",
        "abs_diff_log",
        "abs_diff_raw",
        "dominance_label",
    ]].rename(
        columns={
            "canonical_gene_name": "top_canonical_gene_name",
            "hpa_gene_label": "top_hpa_gene_label",
            "bgee_gene_label": "top_bgee_gene_label",
            "abs_diff_log": "top_abs_diff_log",
            "abs_diff_raw": "top_abs_diff_raw",
            "dominance_label": "top_dominance_label",
        }
    )
    summary = grouped.merge(top_rows, on="compare_tissue_key", how="left")
    return summary.sort_values(
        by=["mean_abs_log_diff", "max_abs_log_diff", "tissue_order"],
        ascending=[False, False, True],
        kind="stable",
    ).reset_index(drop=True)


def save_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def fmt_value(value: object, *, digits: int = 4) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return "-"
    if isinstance(value, (np.floating, float)):
        return f"{float(value):.{digits}f}"
    if isinstance(value, (np.integer, int)):
        return str(int(value))
    if isinstance(value, (bool, np.bool_)):
        return "True" if bool(value) else "False"
    return str(value)


def markdown_table(
    df: pd.DataFrame,
    *,
    columns: Sequence[str],
    headers: Sequence[str],
    digits_map: Dict[str, int] | None = None,
) -> List[str]:
    digits_map = digits_map or {}
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for _, row in df.iterrows():
        values = [fmt_value(row[column], digits=digits_map.get(column, 4)) for column in columns]
        lines.append("| " + " | ".join(values) + " |")
    return lines


def write_report(
    out_md: Path,
    *,
    paired_df: pd.DataFrame,
    correlation_df: pd.DataFrame,
    top_diff_df: pd.DataFrame,
    gene_summary_df: pd.DataFrame,
    tissue_summary_df: pd.DataFrame,
    compare_metadata: Dict[str, object],
    stats_metadata: Dict[str, object],
) -> None:
    generated_at = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    figure_dir = str(
        compare_metadata.get("out_fig_dir", str(get_bioanalyze_figures_root() / "heatmaps" / "compare_nTPM_bgee"))
    )
    total_slots = int(stats_metadata["total_compare_slots"])
    non_missing_slots = int(stats_metadata["non_missing_pair_count"])
    nonzero_any_slots = int(stats_metadata["nonzero_any_pair_count"])

    overall_top = top_diff_df.loc[top_diff_df["subset_label"].eq("overall")].head(10).copy()
    hpa_top = top_diff_df.loc[top_diff_df["subset_label"].eq("hpa_higher")].head(8).copy()
    bgee_top = top_diff_df.loc[top_diff_df["subset_label"].eq("bgee_higher")].head(8).copy()
    top_genes = gene_summary_df.head(10).copy()
    top_tissues = tissue_summary_df.head(10).copy()
    safe_hits = overall_top.loc[overall_top["match_type"].eq("safe_synonym")].copy()

    lines: List[str] = [
        "# HPA vs Bgee Compare Statistics",
        "",
        f"Generated: {generated_at}",
        "",
        f"- Compare figures: `{figure_dir}`",
        f"- HPA compare table: `{stats_metadata['hpa_tsv']}`",
        f"- Bgee compare table: `{stats_metadata['bgee_tsv']}`",
        f"- Paired stats table: `{out_md.parent / 'paired_cells.tsv'}`",
        "",
        "## Overview",
        "",
        f"- Total compare slots: {total_slots}",
        f"- Non-missing paired cells: {non_missing_slots}",
        f"- Nonzero-any paired cells: {nonzero_any_slots}",
        "- Primary comparison space: `log10(value + 1)`",
        "- HPA compare value: `nTPM`",
        "- Bgee compare value: `cell_mean_score`",
        "",
        "## Correlation Summary",
        "",
    ]
    lines.extend(
        markdown_table(
            correlation_df,
            columns=[
                "subset_label",
                "pair_count",
                "pearson_log",
                "spearman_log",
                "pearson_raw",
                "spearman_raw",
            ],
            headers=[
                "Subset",
                "Pairs",
                "Pearson log",
                "Spearman log",
                "Pearson raw",
                "Spearman raw",
            ],
            digits_map={
                "pearson_log": 6,
                "spearman_log": 6,
                "pearson_raw": 6,
                "spearman_raw": 6,
            },
        )
    )

    lines.extend(["", "## Most Different Cells", ""])
    lines.extend(
        markdown_table(
            overall_top,
            columns=[
                "canonical_gene_name",
                "hpa_tissue_label",
                "bgee_tissue_label",
                "match_type",
                "nTPM",
                "cell_mean_score",
                "abs_diff_log",
                "dominance_label",
            ],
            headers=[
                "Gene",
                "HPA tissue",
                "Bgee tissue",
                "Match",
                "HPA nTPM",
                "Bgee score",
                "Abs log diff",
                "Dominance",
            ],
            digits_map={
                "nTPM": 4,
                "cell_mean_score": 4,
                "abs_diff_log": 6,
            },
        )
    )

    if not safe_hits.empty:
        lines.extend(["", "Safe-synonym cases among the strongest differences:", ""])
        for _, row in safe_hits.iterrows():
            lines.append(
                f"- `{row['canonical_gene_name']}`: `{row['hpa_tissue_label']}` vs `{row['bgee_tissue_label']}` "
                f"(`abs_diff_log={row['abs_diff_log']:.6f}`)"
            )

    lines.extend(["", "## Most Different HPA-Higher Cells", ""])
    lines.extend(
        markdown_table(
            hpa_top,
            columns=[
                "canonical_gene_name",
                "hpa_tissue_label",
                "bgee_tissue_label",
                "nTPM",
                "cell_mean_score",
                "abs_diff_log",
            ],
            headers=["Gene", "HPA tissue", "Bgee tissue", "HPA nTPM", "Bgee score", "Abs log diff"],
            digits_map={"nTPM": 4, "cell_mean_score": 4, "abs_diff_log": 6},
        )
    )

    lines.extend(["", "## Most Different Bgee-Higher Cells", ""])
    lines.extend(
        markdown_table(
            bgee_top,
            columns=[
                "canonical_gene_name",
                "hpa_tissue_label",
                "bgee_tissue_label",
                "nTPM",
                "cell_mean_score",
                "abs_diff_log",
            ],
            headers=["Gene", "HPA tissue", "Bgee tissue", "HPA nTPM", "Bgee score", "Abs log diff"],
            digits_map={"nTPM": 4, "cell_mean_score": 4, "abs_diff_log": 6},
        )
    )

    lines.extend(["", "## Genes With The Strongest Overall Disagreement", ""])
    lines.extend(
        markdown_table(
            top_genes,
            columns=[
                "canonical_gene_name",
                "comparable_cell_count",
                "mean_abs_log_diff",
                "median_abs_log_diff",
                "max_abs_log_diff",
                "top_hpa_tissue_label",
                "top_bgee_tissue_label",
            ],
            headers=[
                "Gene",
                "Cells",
                "Mean abs log diff",
                "Median abs log diff",
                "Max abs log diff",
                "Top HPA tissue",
                "Top Bgee tissue",
            ],
            digits_map={
                "mean_abs_log_diff": 6,
                "median_abs_log_diff": 6,
                "max_abs_log_diff": 6,
            },
        )
    )

    lines.extend(["", "## Tissues With The Strongest Overall Disagreement", ""])
    lines.extend(
        markdown_table(
            top_tissues,
            columns=[
                "compare_tissue_key",
                "hpa_tissue_label",
                "bgee_tissue_label",
                "match_type",
                "mean_abs_log_diff",
                "median_abs_log_diff",
                "max_abs_log_diff",
                "top_canonical_gene_name",
            ],
            headers=[
                "Compare tissue",
                "HPA tissue",
                "Bgee tissue",
                "Match",
                "Mean abs log diff",
                "Median abs log diff",
                "Max abs log diff",
                "Top gene",
            ],
            digits_map={
                "mean_abs_log_diff": 6,
                "median_abs_log_diff": 6,
                "max_abs_log_diff": 6,
            },
        )
    )

    lines.extend(
        [
            "",
            "## Interpretation Note",
            "",
            "- Correlation and difference ranking are computed in the same space shown on the heatmaps: `log10(value + 1)`.",
            "- Raw HPA `nTPM` and raw Bgee `cell_mean_score` are reported next to the log-scale differences so large disagreements remain biologically interpretable.",
        ]
    )

    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_stats_metadata(
    *,
    hpa_tsv: Path,
    bgee_tsv: Path,
    compare_metadata_path: Path,
    out_dir: Path,
    top_n: int,
    merged_slot_count: int,
    paired_df: pd.DataFrame,
    top_diff_df: pd.DataFrame,
) -> Dict[str, object]:
    return {
        "hpa_tsv": hpa_tsv.as_posix(),
        "bgee_tsv": bgee_tsv.as_posix(),
        "compare_metadata_path": compare_metadata_path.as_posix(),
        "out_dir": out_dir.as_posix(),
        "top_n": int(top_n),
        "total_compare_slots": int(merged_slot_count),
        "non_missing_pair_count": int(len(paired_df)),
        "nonzero_any_pair_count": int(((paired_df["nTPM"] > 0.0) | (paired_df["cell_mean_score"] > 0.0)).sum()),
        "both_zero_pair_count": int(paired_df["pair_status"].eq("both_zero").sum()),
        "safe_synonym_pair_count": int(paired_df["is_underlined"].sum()),
        "max_abs_log_diff": float(paired_df["abs_diff_log"].max()) if not paired_df.empty else np.nan,
        "max_abs_raw_diff": float(paired_df["abs_diff_raw"].max()) if not paired_df.empty else np.nan,
        "top_difference_examples": top_diff_df.head(5)[
            [
                "subset_label",
                "canonical_gene_name",
                "hpa_tissue_label",
                "bgee_tissue_label",
                "abs_diff_log",
            ]
        ].to_dict(orient="records"),
    }


def main() -> None:
    args = parse_args()
    hpa_tsv = Path(args.hpa_tsv)
    bgee_tsv = Path(args.bgee_tsv)
    compare_metadata_path = Path(args.compare_metadata)
    out_dir = Path(args.out_dir)
    top_n = max(int(args.top_n), 1)

    hpa_df = load_hpa_df(hpa_tsv)
    bgee_df = load_bgee_df(bgee_tsv)
    compare_metadata = load_optional_metadata(compare_metadata_path)

    merged_slot_count = int(
        hpa_df.merge(
            bgee_df,
            on=MERGE_KEYS,
            how="inner",
        ).shape[0]
    )
    paired_df = build_paired_df(hpa_df, bgee_df)
    correlation_df = build_correlation_summary(paired_df)
    top_diff_df = build_top_differing_cells(paired_df, top_n=top_n)
    gene_summary_df = build_gene_difference_summary(paired_df)
    tissue_summary_df = build_tissue_difference_summary(paired_df)

    out_dir.mkdir(parents=True, exist_ok=True)
    save_tsv(paired_df, out_dir / "paired_cells.tsv")
    save_tsv(correlation_df, out_dir / "correlation_summary.tsv")
    save_tsv(top_diff_df, out_dir / "top_differing_cells.tsv")
    save_tsv(gene_summary_df, out_dir / "gene_difference_summary.tsv")
    save_tsv(tissue_summary_df, out_dir / "tissue_difference_summary.tsv")

    stats_metadata = build_stats_metadata(
        hpa_tsv=hpa_tsv,
        bgee_tsv=bgee_tsv,
        compare_metadata_path=compare_metadata_path,
        out_dir=out_dir,
        top_n=top_n,
        merged_slot_count=merged_slot_count,
        paired_df=paired_df,
        top_diff_df=top_diff_df,
    )
    (out_dir / "metadata.json").write_text(json.dumps(stats_metadata, indent=2), encoding="utf-8")
    write_report(
        out_md=out_dir / "report_md.md",
        paired_df=paired_df,
        correlation_df=correlation_df,
        top_diff_df=top_diff_df,
        gene_summary_df=gene_summary_df,
        tissue_summary_df=tissue_summary_df,
        compare_metadata=compare_metadata,
        stats_metadata=stats_metadata,
    )

    print(f"TOTAL_COMPARE_SLOTS={merged_slot_count}")
    print(f"NON_MISSING_PAIRED_CELLS={len(paired_df)}")
    print(
        "LOG_PEARSON_ALL="
        + f"{correlation_df.loc[correlation_df['subset_label'].eq('all_non_missing'), 'pearson_log'].iloc[0]:.6f}"
    )
    print(
        "LOG_SPEARMAN_ALL="
        + f"{correlation_df.loc[correlation_df['subset_label'].eq('all_non_missing'), 'spearman_log'].iloc[0]:.6f}"
    )
    print(f"OUT_DIR={out_dir}")


if __name__ == "__main__":
    main()
