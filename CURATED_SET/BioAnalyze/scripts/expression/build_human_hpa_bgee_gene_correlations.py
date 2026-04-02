#!/usr/bin/env python
"""Build per-gene HPA vs Bgee raw correlation barplots for the human compare layer."""

from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Sequence

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DEFAULT_PAIRED_TSV = Path(
    r"CURATED_SET/BioAnalyze/stats/compare_nTPM_bgee/paired_cells.tsv"
)
DEFAULT_COMPARE_METADATA = Path(
    r"CURATED_SET/BioAnalyze/data/processed/intersections/human_hpa_bgee/hpa_vs_bgee_metadata.json"
)
DEFAULT_OUT_DIR = Path(
    r"CURATED_SET/BioAnalyze/stats/compare_nTPM_bgee/per_gene_correlations"
)

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
POSITIVE_COLOR = "#6ea6c8"
NEGATIVE_COLOR = "#d7886a"
MISSING_COLOR = "#b3b3b3"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute per-gene raw Pearson/Spearman correlations between the "
            "human HPA and Bgee compare tables and build barplots."
        )
    )
    parser.add_argument("--paired-tsv", default=str(DEFAULT_PAIRED_TSV))
    parser.add_argument("--compare-metadata", default=str(DEFAULT_COMPARE_METADATA))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    return parser.parse_args()


def require_columns(df: pd.DataFrame, columns: Sequence[str], *, label: str) -> None:
    missing = [column for column in columns if column not in df.columns]
    if missing:
        raise RuntimeError(f"{label} is missing required columns: {', '.join(missing)}")


def load_paired_df(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path, sep="\t", low_memory=True)
    require_columns(
        df,
        [
            "gene_order",
            "gene_number",
            "canonical_gene_name",
            "hpa_gene_label",
            "bgee_gene_label",
            "compare_tissue_key",
            "hpa_tissue_label",
            "bgee_tissue_label",
            "match_type",
            "is_underlined",
            "nTPM",
            "cell_mean_score",
        ],
        label="paired_cells.tsv",
    )
    df["gene_order"] = pd.to_numeric(df["gene_order"], errors="coerce").fillna(0).astype(int)
    df["gene_number"] = pd.to_numeric(df["gene_number"], errors="coerce").fillna(0).astype(int)
    df["nTPM"] = pd.to_numeric(df["nTPM"], errors="coerce")
    df["cell_mean_score"] = pd.to_numeric(df["cell_mean_score"], errors="coerce")
    df["is_underlined"] = df["is_underlined"].fillna(False).astype(bool)
    return df


def load_optional_metadata(path: Path) -> Dict[str, object]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def filter_mode_rows(gene_df: pd.DataFrame, mode: str) -> pd.DataFrame:
    if mode == "keep_zeros":
        filtered = gene_df.copy()
    elif mode == "drop_zeros":
        filtered = gene_df[
            gene_df["nTPM"].gt(0.0) & gene_df["cell_mean_score"].gt(0.0)
        ].copy()
    else:
        raise RuntimeError(f"Unsupported mode: {mode}")
    return filtered


def correlation_status(filtered_df: pd.DataFrame) -> str:
    if len(filtered_df) < 2:
        return "insufficient_pairs"
    if (
        filtered_df["nTPM"].nunique(dropna=True) < 2
        or filtered_df["cell_mean_score"].nunique(dropna=True) < 2
    ):
        return "constant_input"
    return "ok"


def compute_gene_correlations(paired_df: pd.DataFrame, *, mode: str) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    grouped = paired_df.groupby(
        ["gene_order", "gene_number", "canonical_gene_name", "hpa_gene_label", "bgee_gene_label"],
        sort=False,
    )
    for (gene_order, gene_number, canonical_gene_name, hpa_gene_label, bgee_gene_label), gene_df in grouped:
        filtered_df = filter_mode_rows(gene_df, mode)
        status = correlation_status(filtered_df)
        pearson_raw = np.nan
        spearman_raw = np.nan
        if status == "ok":
            pearson_raw = float(filtered_df["nTPM"].corr(filtered_df["cell_mean_score"], method="pearson"))
            spearman_raw = float(filtered_df["nTPM"].corr(filtered_df["cell_mean_score"], method="spearman"))

        rows.append(
            {
                "mode": mode,
                "gene_order": int(gene_order),
                "gene_number": int(gene_number),
                "canonical_gene_name": canonical_gene_name,
                "hpa_gene_label": hpa_gene_label,
                "bgee_gene_label": bgee_gene_label,
                "paired_tissue_count": int(len(filtered_df)),
                "safe_synonym_tissue_count": int(filtered_df["is_underlined"].sum()),
                "pearson_raw": pearson_raw,
                "spearman_raw": spearman_raw,
                "status": status,
                "total_tissue_count_before_filter": int(len(gene_df)),
            }
        )

    result = pd.DataFrame(rows)
    result["sort_key"] = result["canonical_gene_name"].map(
        {gene_name: idx for idx, gene_name in enumerate(GENE_ORDER, start=1)}
    )
    result = result.sort_values(
        by=["sort_key", "gene_order", "canonical_gene_name"],
        ascending=[True, True, True],
        kind="stable",
    ).drop(columns=["sort_key"])
    return result.reset_index(drop=True)


def save_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def style_axis(ax: plt.Axes) -> None:
    for spine in ax.spines.values():
        spine.set_visible(False)


def correlation_bar_colors(values: pd.Series) -> List[str]:
    colors: List[str] = []
    for value in values:
        if pd.isna(value):
            colors.append(MISSING_COLOR)
        elif float(value) >= 0.0:
            colors.append(POSITIVE_COLOR)
        else:
            colors.append(NEGATIVE_COLOR)
    return colors


def plot_correlation_bars(
    df: pd.DataFrame,
    *,
    metric_col: str,
    title: str,
    ylabel: str,
    out_png: Path,
    out_svg: Path,
) -> None:
    fig_w = max(10.5, len(df) * 0.58)
    fig, ax = plt.subplots(figsize=(fig_w, 5.8))
    plot_df = df.copy()
    plot_df["plot_value"] = plot_df[metric_col].fillna(0.0)
    colors = correlation_bar_colors(plot_df[metric_col])
    bars = ax.bar(plot_df["canonical_gene_name"], plot_df["plot_value"], color=colors)
    ax.axhline(0.0, color="#333333", linewidth=1.0)
    ax.set_ylim(-1.0, 1.0)
    ax.set_title(title)
    ax.set_xlabel("Gene")
    ax.set_ylabel(ylabel)
    ax.tick_params(axis="x", rotation=65, labelsize=9)
    ax.tick_params(axis="y", labelsize=10)
    style_axis(ax)

    for bar, raw_value in zip(bars, plot_df[metric_col].tolist()):
        x = bar.get_x() + bar.get_width() / 2.0
        if pd.isna(raw_value):
            ax.text(x, 0.03, "NA", ha="center", va="bottom", fontsize=8, color="#555555", rotation=90)
        else:
            y = float(raw_value)
            offset = 0.03 if y >= 0 else -0.03
            va = "bottom" if y >= 0 else "top"
            ax.text(x, y + offset, f"{y:.2f}", ha="center", va=va, fontsize=8, color="#333333", rotation=90)

    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    plt.close(fig)


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


def strongest_and_weakest(df: pd.DataFrame, metric_col: str, limit: int = 5) -> tuple[pd.DataFrame, pd.DataFrame]:
    valid = df[df[metric_col].notna()].copy()
    strongest = valid.sort_values(
        by=[metric_col, "canonical_gene_name"],
        ascending=[False, True],
        kind="stable",
    ).head(limit)
    weakest = valid.sort_values(
        by=[metric_col, "canonical_gene_name"],
        ascending=[True, True],
        kind="stable",
    ).head(limit)
    return strongest, weakest


def write_report(
    out_md: Path,
    *,
    keep_df: pd.DataFrame,
    drop_df: pd.DataFrame,
    compare_metadata: Dict[str, object],
) -> None:
    generated_at = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    figure_dir = str(compare_metadata.get("out_fig_dir", "CURATED_SET/BioAnalyze/figures/heatmaps/compare_nTPM_bgee"))

    lines: List[str] = [
        "# Per-Gene HPA vs Bgee Correlations",
        "",
        f"Generated: {generated_at}",
        "",
        f"- Source paired table: `CURATED_SET/BioAnalyze/stats/compare_nTPM_bgee/paired_cells.tsv`",
        f"- Compare figures: `{figure_dir}`",
        "",
        "## Modes",
        "",
        "- `keep_zeros`: all non-missing paired tissues are kept; zeroes are preserved.",
        "- `drop_zeros`: only tissues with both `nTPM > 0` and `cell_mean_score > 0` are used.",
        "- Correlations are raw-only here: `pearson_raw` and `spearman_raw`.",
        "",
    ]

    mode_configs = [
        ("keep_zeros", keep_df),
        ("drop_zeros", drop_df),
    ]
    for mode_name, mode_df in mode_configs:
        lines.extend([f"## {mode_name}", ""])
        for metric_col, metric_label in [("pearson_raw", "Pearson raw"), ("spearman_raw", "Spearman raw")]:
            strongest, weakest = strongest_and_weakest(mode_df, metric_col)
            lines.extend([f"### {metric_label}: strongest", ""])
            lines.extend(
                markdown_table(
                    strongest,
                    columns=[
                        "canonical_gene_name",
                        "bgee_gene_label",
                        "paired_tissue_count",
                        metric_col,
                        "status",
                    ],
                    headers=["Gene", "Bgee label", "Tissues", metric_label, "Status"],
                    digits_map={metric_col: 6},
                )
            )
            lines.extend(["", f"### {metric_label}: weakest", ""])
            lines.extend(
                markdown_table(
                    weakest,
                    columns=[
                        "canonical_gene_name",
                        "bgee_gene_label",
                        "paired_tissue_count",
                        metric_col,
                        "status",
                    ],
                    headers=["Gene", "Bgee label", "Tissues", metric_label, "Status"],
                    digits_map={metric_col: 6},
                )
            )
            lines.append("")

        insufficient = mode_df[mode_df["status"].ne("ok")].copy()
        lines.append("### Genes with insufficient or constant input")
        lines.append("")
        if insufficient.empty:
            lines.append("- None.")
        else:
            for _, row in insufficient.iterrows():
                lines.append(
                    f"- `{row['canonical_gene_name']}` (`Bgee: {row['bgee_gene_label']}`): "
                    f"status=`{row['status']}`, usable tissues={int(row['paired_tissue_count'])}"
                )
        lines.append("")

    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    paired_tsv = Path(args.paired_tsv)
    compare_metadata_path = Path(args.compare_metadata)
    out_dir = Path(args.out_dir)

    paired_df = load_paired_df(paired_tsv)
    compare_metadata = load_optional_metadata(compare_metadata_path)

    keep_df = compute_gene_correlations(paired_df, mode="keep_zeros")
    drop_df = compute_gene_correlations(paired_df, mode="drop_zeros")

    out_dir.mkdir(parents=True, exist_ok=True)
    save_tsv(keep_df, out_dir / "gene_correlation_keep_zeros.tsv")
    save_tsv(drop_df, out_dir / "gene_correlation_drop_zeros.tsv")

    plot_correlation_bars(
        keep_df,
        metric_col="pearson_raw",
        title="Per-gene HPA vs Bgee Pearson correlation (keep zeros)",
        ylabel="Pearson raw",
        out_png=out_dir / "pearson_raw_keep_zeros.png",
        out_svg=out_dir / "pearson_raw_keep_zeros.svg",
    )
    plot_correlation_bars(
        keep_df,
        metric_col="spearman_raw",
        title="Per-gene HPA vs Bgee Spearman correlation (keep zeros)",
        ylabel="Spearman raw",
        out_png=out_dir / "spearman_raw_keep_zeros.png",
        out_svg=out_dir / "spearman_raw_keep_zeros.svg",
    )
    plot_correlation_bars(
        drop_df,
        metric_col="pearson_raw",
        title="Per-gene HPA vs Bgee Pearson correlation (drop zeros)",
        ylabel="Pearson raw",
        out_png=out_dir / "pearson_raw_drop_zeros.png",
        out_svg=out_dir / "pearson_raw_drop_zeros.svg",
    )
    plot_correlation_bars(
        drop_df,
        metric_col="spearman_raw",
        title="Per-gene HPA vs Bgee Spearman correlation (drop zeros)",
        ylabel="Spearman raw",
        out_png=out_dir / "spearman_raw_drop_zeros.png",
        out_svg=out_dir / "spearman_raw_drop_zeros.svg",
    )

    write_report(
        out_md=out_dir / "report_md.md",
        keep_df=keep_df,
        drop_df=drop_df,
        compare_metadata=compare_metadata,
    )

    print(f"KEEP_ZERO_GENES={len(keep_df)}")
    print(f"DROP_ZERO_GENES={len(drop_df)}")
    print(f"OUT_DIR={out_dir}")


if __name__ == "__main__":
    main()
