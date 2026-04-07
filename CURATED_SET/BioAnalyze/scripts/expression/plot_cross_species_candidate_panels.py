#!/usr/bin/env python
"""Plot overview and paginated panels for cross-species H2A candidates."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Sequence, Tuple

import matplotlib.pyplot as plt
import pandas as pd

from gene_compare_common import (
    CROSS_SPECIES_CLASS_COLORS,
    DEFAULT_GENE_COMPARE_DATA_ROOT,
    DEFAULT_RANKING_PLOTS_DIR,
    DEFAULT_RANKING_TABLES_DIR,
    safe_slug,
    summarize_species_gene_tissue,
)


DEFAULT_OUT_DIR = DEFAULT_RANKING_PLOTS_DIR
DEFAULT_TOP_OVERALL = 2
DEFAULT_TOP_VARIANT = 4
DEFAULT_MIN_SPECIES = 4
DEFAULT_QUANTILE = 0.95
DEFAULT_HIGH_QUANTILES = [0.90, 0.95]
DEFAULT_LOW_QUANTILES = [0.05, 0.10]
DEFAULT_PANELS_PER_PAGE = 6
DEFAULT_PREFERRED_VARIANTS = [
    ("H2AX", "adipose tissue"),
    ("MACROH2A2", "zone of skin"),
    ("H2AJ", "testis"),
    ("H2AJ", "pituitary gland"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create overview plots and paginated species-level panels for "
            "cross-species H2A candidates."
        )
    )
    parser.add_argument(
        "--ranking-dir",
        default=str(DEFAULT_RANKING_TABLES_DIR),
        help="Directory with gene_tissue_summary.csv and candidate tables.",
    )
    parser.add_argument(
        "--gene-compare-data-root",
        default=str(DEFAULT_GENE_COMPARE_DATA_ROOT),
        help="Root directory with per-gene gene_compare long CSV files.",
    )
    parser.add_argument(
        "--out-dir",
        default=str(DEFAULT_OUT_DIR),
        help="Directory for output figures.",
    )
    parser.add_argument(
        "--top-overall",
        default=DEFAULT_TOP_OVERALL,
        type=int,
        help="Number of top overall candidates to include when --candidates is not used.",
    )
    parser.add_argument(
        "--top-variant",
        default=DEFAULT_TOP_VARIANT,
        type=int,
        help="Number of top variant candidates to include when --candidates is not used.",
    )
    parser.add_argument(
        "--candidates",
        nargs="*",
        default=[],
        help="Optional explicit candidates in GENE::tissue format.",
    )
    parser.add_argument(
        "--include-class-high-panels",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Build class-specific high-tail panel pages for each high quantile.",
    )
    parser.add_argument(
        "--include-class-p95",
        action="store_true",
        help="Backward-compatible alias for class-specific high-tail panel pages.",
    )
    parser.add_argument(
        "--include-global-low-panels",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Build global low-tail panel pages for each low quantile.",
    )
    parser.add_argument(
        "--include-class-low-panels",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Build class-specific low-tail panel pages for each low quantile.",
    )
    parser.add_argument(
        "--min-species",
        default=DEFAULT_MIN_SPECIES,
        type=int,
        help="Minimum species_n for candidate selection in overview and paginated panels.",
    )
    parser.add_argument(
        "--quantile",
        default=DEFAULT_QUANTILE,
        type=float,
        help="Backward-compatible high quantile for --include-class-p95.",
    )
    parser.add_argument(
        "--high-quantiles",
        nargs="*",
        type=float,
        default=DEFAULT_HIGH_QUANTILES,
        help="High-tail quantiles for class-specific candidate pages.",
    )
    parser.add_argument(
        "--low-quantiles",
        nargs="*",
        type=float,
        default=DEFAULT_LOW_QUANTILES,
        help="Low-tail quantiles for global/class-specific candidate pages.",
    )
    parser.add_argument(
        "--panels-per-page",
        default=DEFAULT_PANELS_PER_PAGE,
        type=int,
        help="Maximum number of candidate panels per page for paginated outputs.",
    )
    return parser.parse_args()


def parse_candidate_values(values: List[str]) -> List[Tuple[str, str]]:
    pairs: List[Tuple[str, str]] = []
    for value in values:
        if "::" not in value:
            raise ValueError(f"Invalid candidate format: {value}. Use GENE::tissue")
        gene_name, tissue = value.split("::", 1)
        gene_name = gene_name.strip()
        tissue = tissue.strip()
        if not gene_name or not tissue:
            raise ValueError(f"Invalid candidate format: {value}. Use GENE::tissue")
        pairs.append((gene_name, tissue))
    return pairs


def normalize_quantiles(values: Sequence[float], *, low_tail: bool) -> List[float]:
    cleaned: List[float] = []
    for value in values:
        if not 0 < float(value) < 1:
            raise ValueError("All quantiles must be between 0 and 1")
        if low_tail and float(value) >= 0.5:
            raise ValueError("Low-tail quantiles must be below 0.5")
        cleaned.append(float(value))
    order = cleaned if low_tail else sorted(cleaned, reverse=False)
    return sorted(set(order))


def quantile_label(value: float) -> str:
    return f"p{int(round(value * 100)):02d}"


def choose_default_candidates(
    summary_df: pd.DataFrame,
    high_conf_df: pd.DataFrame,
    top_overall: int,
    top_variant: int,
) -> pd.DataFrame:
    overall_df = high_conf_df.sort_values(
        by=["range", "species_n", "gene_name", "tissue"],
        ascending=[False, False, True, True],
    ).head(top_overall)
    variant_df = (
        summary_df[(summary_df["gene_class"] == "variant") & (summary_df["species_n"] >= 4)]
        .sort_values(
            by=["range", "species_n", "gene_name", "tissue"],
            ascending=[False, False, True, True],
        )
    )

    preferred_frames: List[pd.DataFrame] = []
    for gene_name, tissue in DEFAULT_PREFERRED_VARIANTS:
        preferred = variant_df[
            variant_df["gene_name"].eq(gene_name) & variant_df["tissue"].eq(tissue)
        ].copy()
        if not preferred.empty:
            preferred_frames.append(preferred.head(1))
    preferred_df = pd.concat(preferred_frames, ignore_index=True) if preferred_frames else pd.DataFrame()

    chosen_variant_df = pd.concat([preferred_df, variant_df], ignore_index=True)
    chosen_variant_df = chosen_variant_df.drop_duplicates(
        subset=["gene_name", "tissue"], keep="first"
    ).head(top_variant)

    chosen_df = pd.concat([overall_df, chosen_variant_df], ignore_index=True)
    chosen_df = chosen_df.drop_duplicates(subset=["gene_name", "tissue"], keep="first")
    return chosen_df.reset_index(drop=True)


def select_quantile_candidates(
    summary_df: pd.DataFrame,
    *,
    min_species: int,
    quantile: float,
    direction: str,
    gene_class: str | None = None,
) -> tuple[pd.DataFrame, float]:
    candidate_df = summary_df[summary_df["species_n"] >= int(min_species)].copy()
    if gene_class is not None:
        candidate_df = candidate_df[candidate_df["gene_class"].eq(gene_class)].copy()
    if candidate_df.empty:
        return candidate_df, float("nan")

    threshold = float(candidate_df["range"].quantile(quantile))
    if direction == "high":
        candidate_df = candidate_df[candidate_df["range"] >= threshold].copy()
        candidate_df = candidate_df.sort_values(
            by=["range", "species_n", "gene_name", "tissue"],
            ascending=[False, False, True, True],
        )
    elif direction == "low":
        candidate_df = candidate_df[candidate_df["range"] <= threshold].copy()
        candidate_df = candidate_df.sort_values(
            by=["range", "species_n", "gene_name", "tissue"],
            ascending=[True, False, True, True],
        )
    else:
        raise ValueError(f"Unsupported direction: {direction}")

    candidate_df = candidate_df.reset_index(drop=True)
    candidate_df["panel_direction"] = direction
    candidate_df["panel_quantile"] = float(quantile)
    candidate_df["panel_threshold"] = threshold
    candidate_df["panel_scope"] = gene_class or "global"
    return candidate_df, threshold


def load_candidate_scores(
    gene_compare_data_root: Path,
    gene_name: str,
    tissue: str,
) -> pd.DataFrame:
    gene_slug = safe_slug(gene_name)
    long_path = gene_compare_data_root / gene_slug / f"{gene_slug}_gene_compare_long.csv"
    if not long_path.exists():
        raise FileNotFoundError(long_path)

    long_df = pd.read_csv(long_path)
    if "cell_mean_score" not in long_df.columns and "expression_score" in long_df.columns:
        long_df["cell_mean_score"] = long_df["expression_score"]
    if "cell_std_score" not in long_df.columns:
        long_df["cell_std_score"] = 0.0
    if "cell_n" not in long_df.columns:
        long_df["cell_n"] = 0
    if "cell_status" not in long_df.columns:
        long_df["cell_status"] = ""
    filtered = long_df[long_df["tissue"].eq(tissue)].copy()
    if filtered.empty:
        raise RuntimeError(f"No rows found for {gene_name} / {tissue} in {long_path}")

    species_level_df = summarize_species_gene_tissue(filtered)
    plot_df = species_level_df[
        ["species_dir", "species_name", "cell_mean_score", "cell_std_score", "cell_n", "cell_status"]
    ].copy()
    plot_df = plot_df.rename(columns={"cell_mean_score": "expression_score"})
    plot_df = plot_df.sort_values(by=["expression_score", "species_name"], ascending=[True, True]).reset_index(
        drop=True
    )
    return plot_df


def plot_overview(
    summary_df: pd.DataFrame,
    highlight_df: pd.DataFrame,
    out_png: Path,
    out_svg: Path,
    min_species: int = DEFAULT_MIN_SPECIES,
) -> None:
    plot_df = summary_df[summary_df["species_n"] >= int(min_species)].copy()
    color_map = CROSS_SPECIES_CLASS_COLORS

    fig, ax = plt.subplots(figsize=(11, 7))
    for gene_class, class_df in plot_df.groupby("gene_class", dropna=False):
        ax.scatter(
            class_df["species_n"],
            class_df["range"],
            s=class_df["species_n"] * 28,
            alpha=0.65,
            color=color_map.get(gene_class, "#888888"),
            label=gene_class or "unclassified",
            edgecolor="white",
            linewidth=0.6,
        )

    for row in highlight_df.to_dict(orient="records"):
        ax.scatter(
            row["species_n"],
            row["range"],
            s=row["species_n"] * 42,
            color="#B22222",
            edgecolor="black",
            linewidth=0.8,
            zorder=5,
        )
        ax.annotate(
            f'{row["gene_name"]}\n{row["tissue"]}',
            (row["species_n"], row["range"]),
            xytext=(6, 4),
            textcoords="offset points",
            fontsize=8,
        )

    ax.set_title("Cross-species H2A candidates by range and species coverage")
    ax.set_xlabel("Number of species with tissue-level score")
    ax.set_ylabel("Expression score range across species")
    ax.legend(title="Gene class")
    ax.grid(alpha=0.2)
    plt.tight_layout()

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    plt.close(fig)


def panel_title(row: dict) -> str:
    subtitle = f'range {row["range"]:.2f}, {int(row["species_n"])} species'
    if "panel_direction" in row and "panel_quantile" in row:
        direction = str(row.get("panel_direction", ""))
        quantile = quantile_label(float(row.get("panel_quantile", 0.0)))
        if direction == "low":
            subtitle = f"{quantile} low-tail | {subtitle}"
        elif direction == "high":
            subtitle = f"{quantile} high-tail | {subtitle}"
    return f'{row["gene_name"]} | {row["tissue"]}\n{subtitle}'


def plot_candidate_panels(
    candidate_df: pd.DataFrame,
    gene_compare_data_root: Path,
    out_png: Path,
    out_svg: Path,
) -> None:
    if candidate_df.empty:
        raise RuntimeError("No candidates available for focused panels.")

    n_panels = len(candidate_df)
    n_cols = 2
    n_rows = (n_panels + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(14, max(4 * n_rows, 6)))
    axes_flat = axes.flatten() if hasattr(axes, "flatten") else [axes]
    color_map = CROSS_SPECIES_CLASS_COLORS

    for ax, row in zip(axes_flat, candidate_df.to_dict(orient="records")):
        plot_df = load_candidate_scores(
            gene_compare_data_root,
            row["gene_name"],
            row["tissue"],
        )
        ax.barh(
            plot_df["species_name"],
            plot_df["expression_score"],
            xerr=plot_df["cell_std_score"],
            color=color_map.get(row["gene_class"], "#888888"),
            alpha=0.9,
            ecolor="#222222",
            capsize=3,
        )
        ax.set_title(panel_title(row))
        ax.set_xlabel("Expression score")
        ax.set_ylabel("")
        ax.set_xlim(
            0,
            max(
                100,
                float((plot_df["expression_score"] + plot_df["cell_std_score"]).max()) + 5,
            ),
        )
        ax.grid(axis="x", alpha=0.2)

    for ax in axes_flat[n_panels:]:
        ax.axis("off")

    plt.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    plt.close(fig)


def build_panel_rows(
    page_df: pd.DataFrame,
    out_png: Path,
    out_svg: Path,
    stem_prefix: str,
    panel_scope: str,
    panel_direction: str,
    quantile: float,
    page_idx: int,
) -> List[dict]:
    rows: List[dict] = []
    for row in page_df.to_dict(orient="records"):
        rows.append(
            {
                "gene_name": row["gene_name"],
                "tissue": row["tissue"],
                "gene_class": row.get("gene_class", ""),
                "panel_family": stem_prefix,
                "panel_label": stem_prefix.replace("_", " "),
                "panel_scope": panel_scope,
                "panel_direction": panel_direction,
                "panel_quantile": float(quantile),
                "page": page_idx,
                "file_png": out_png.as_posix(),
                "file_svg": out_svg.as_posix(),
            }
        )
    return rows


def write_paginated_panels(
    candidate_df: pd.DataFrame,
    gene_compare_data_root: Path,
    out_dir: Path,
    stem_prefix: str,
    panels_per_page: int,
    *,
    panel_scope: str,
    panel_direction: str,
    quantile: float,
) -> tuple[List[Path], List[dict]]:
    if candidate_df.empty:
        return [], []

    written_paths: List[Path] = []
    panel_rows: List[dict] = []
    total = len(candidate_df)
    for page_start in range(0, total, panels_per_page):
        page_idx = page_start // panels_per_page + 1
        page_df = candidate_df.iloc[page_start : page_start + panels_per_page].reset_index(drop=True)
        out_png = out_dir / f"{stem_prefix}_page{page_idx}.png"
        out_svg = out_dir / f"{stem_prefix}_page{page_idx}.svg"
        plot_candidate_panels(page_df, gene_compare_data_root, out_png, out_svg)
        written_paths.extend([out_png, out_svg])
        panel_rows.extend(
            build_panel_rows(
                page_df,
                out_png,
                out_svg,
                stem_prefix=stem_prefix,
                panel_scope=panel_scope,
                panel_direction=panel_direction,
                quantile=quantile,
                page_idx=page_idx,
            )
        )
    return written_paths, panel_rows


def append_focus_panel_rows(
    panel_rows: List[dict],
    candidate_df: pd.DataFrame,
    out_png: Path,
    out_svg: Path,
) -> None:
    for row in candidate_df.to_dict(orient="records"):
        panel_rows.append(
            {
                "gene_name": row["gene_name"],
                "tissue": row["tissue"],
                "gene_class": row.get("gene_class", ""),
                "panel_family": "candidate_focus_panels",
                "panel_label": "candidate focus panels",
                "panel_scope": "focus",
                "panel_direction": "mixed",
                "panel_quantile": float("nan"),
                "page": 1,
                "file_png": out_png.as_posix(),
                "file_svg": out_svg.as_posix(),
            }
        )


def write_panel_index(panel_rows: List[dict], out_csv: Path) -> None:
    if not panel_rows:
        pd.DataFrame(
            columns=[
                "gene_name",
                "tissue",
                "gene_class",
                "panel_family",
                "panel_label",
                "panel_scope",
                "panel_direction",
                "panel_quantile",
                "page",
                "file_png",
                "file_svg",
            ]
        ).to_csv(out_csv, index=False, encoding="utf-8")
        return

    panel_df = pd.DataFrame(panel_rows).drop_duplicates(
        subset=["gene_name", "tissue", "panel_family", "page"], keep="first"
    )
    panel_df = panel_df.sort_values(
        by=["panel_direction", "panel_quantile", "panel_scope", "gene_name", "tissue", "page"],
        ascending=[True, True, True, True, True, True],
    ).reset_index(drop=True)
    panel_df.to_csv(out_csv, index=False, encoding="utf-8")


def main() -> None:
    args = parse_args()

    if args.panels_per_page < 1:
        raise ValueError("--panels-per-page must be >= 1")

    high_quantiles = normalize_quantiles(args.high_quantiles, low_tail=False)
    if args.include_class_p95 and args.quantile not in high_quantiles:
        high_quantiles = sorted(set(high_quantiles + [float(args.quantile)]))
    low_quantiles = normalize_quantiles(args.low_quantiles, low_tail=True)

    ranking_dir = Path(args.ranking_dir)
    gene_compare_data_root = Path(args.gene_compare_data_root)
    out_dir = Path(args.out_dir)

    summary_path = ranking_dir / "gene_tissue_summary.csv"
    high_conf_path = ranking_dir / "high_confidence_candidates.csv"
    if not summary_path.exists():
        raise FileNotFoundError(summary_path)
    if not high_conf_path.exists():
        raise FileNotFoundError(high_conf_path)

    summary_df = pd.read_csv(summary_path)
    high_conf_df = pd.read_csv(high_conf_path)

    explicit_pairs = parse_candidate_values(args.candidates)
    if explicit_pairs:
        key_df = pd.DataFrame(explicit_pairs, columns=["gene_name", "tissue"])
        candidate_df = key_df.merge(summary_df, on=["gene_name", "tissue"], how="left")
        if candidate_df.isna().any(axis=None):
            missing = candidate_df[candidate_df["range"].isna()][["gene_name", "tissue"]]
            raise RuntimeError(
                "Missing candidate(s) in summary: "
                + ", ".join(f"{r.gene_name}::{r.tissue}" for r in missing.itertuples())
            )
    else:
        candidate_df = choose_default_candidates(
            summary_df,
            high_conf_df,
            top_overall=args.top_overall,
            top_variant=args.top_variant,
        )

    out_dir.mkdir(parents=True, exist_ok=True)
    out_overview_png = out_dir / "candidate_range_overview.png"
    out_overview_svg = out_dir / "candidate_range_overview.svg"
    out_panels_png = out_dir / "candidate_focus_panels.png"
    out_panels_svg = out_dir / "candidate_focus_panels.svg"
    panel_index_csv = out_dir / "panel_membership.csv"

    plot_overview(
        summary_df,
        candidate_df,
        out_overview_png,
        out_overview_svg,
        min_species=args.min_species,
    )
    plot_candidate_panels(candidate_df, gene_compare_data_root, out_panels_png, out_panels_svg)

    print(f"Saved overview plot to {out_overview_png}")
    print(f"Saved overview plot to {out_overview_svg}")
    print(f"Saved candidate panels to {out_panels_png}")
    print(f"Saved candidate panels to {out_panels_svg}")

    panel_rows: List[dict] = []
    append_focus_panel_rows(panel_rows, candidate_df, out_panels_png, out_panels_svg)

    if args.include_class_high_panels or args.include_class_p95:
        for quantile in high_quantiles:
            for gene_class in ["clustered", "variant"]:
                class_df, threshold = select_quantile_candidates(
                    summary_df,
                    min_species=args.min_species,
                    quantile=quantile,
                    direction="high",
                    gene_class=gene_class,
                )
                stem_prefix = f"{gene_class}_{quantile_label(quantile)}_panels"
                written, rows = write_paginated_panels(
                    class_df,
                    gene_compare_data_root=gene_compare_data_root,
                    out_dir=out_dir,
                    stem_prefix=stem_prefix,
                    panels_per_page=args.panels_per_page,
                    panel_scope=gene_class,
                    panel_direction="high",
                    quantile=quantile,
                )
                panel_rows.extend(rows)
                if written:
                    print(
                        f"Saved {len(class_df)} {gene_class} candidate panel(s) across "
                        f"{len(written) // 2} page(s) for {quantile_label(quantile)} "
                        f"(threshold {threshold:.2f})."
                    )
                else:
                    print(
                        f"No {gene_class} candidates met the class-specific "
                        f"{quantile_label(quantile)} threshold."
                    )

    if args.include_global_low_panels:
        for quantile in low_quantiles:
            low_df, threshold = select_quantile_candidates(
                summary_df,
                min_species=args.min_species,
                quantile=quantile,
                direction="low",
                gene_class=None,
            )
            stem_prefix = f"global_{quantile_label(quantile)}_low_panels"
            written, rows = write_paginated_panels(
                low_df,
                gene_compare_data_root=gene_compare_data_root,
                out_dir=out_dir,
                stem_prefix=stem_prefix,
                panels_per_page=args.panels_per_page,
                panel_scope="global",
                panel_direction="low",
                quantile=quantile,
            )
            panel_rows.extend(rows)
            if written:
                print(
                    f"Saved {len(low_df)} global low-tail panel(s) across "
                    f"{len(written) // 2} page(s) for {quantile_label(quantile)} "
                    f"(threshold {threshold:.2f})."
                )
            else:
                print(f"No global candidates met the {quantile_label(quantile)} low-tail threshold.")

    if args.include_class_low_panels:
        for quantile in low_quantiles:
            for gene_class in ["clustered", "variant"]:
                class_df, threshold = select_quantile_candidates(
                    summary_df,
                    min_species=args.min_species,
                    quantile=quantile,
                    direction="low",
                    gene_class=gene_class,
                )
                stem_prefix = f"{gene_class}_{quantile_label(quantile)}_low_panels"
                written, rows = write_paginated_panels(
                    class_df,
                    gene_compare_data_root=gene_compare_data_root,
                    out_dir=out_dir,
                    stem_prefix=stem_prefix,
                    panels_per_page=args.panels_per_page,
                    panel_scope=gene_class,
                    panel_direction="low",
                    quantile=quantile,
                )
                panel_rows.extend(rows)
                if written:
                    print(
                        f"Saved {len(class_df)} {gene_class} low-tail panel(s) across "
                        f"{len(written) // 2} page(s) for {quantile_label(quantile)} "
                        f"(threshold {threshold:.2f})."
                    )
                else:
                    print(
                        f"No {gene_class} candidates met the class-specific "
                        f"{quantile_label(quantile)} low-tail threshold."
                    )

    write_panel_index(panel_rows, panel_index_csv)
    print(f"Saved panel membership index to {panel_index_csv}")


if __name__ == "__main__":
    main()
