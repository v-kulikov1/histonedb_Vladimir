#!/usr/bin/env python
# Build Bgee H2A expression heatmaps split into clustered (canonical) vs variants.
# Rows: unique Ensembl Gene ID (labels: GeneName:HGNC). Columns: anatomical entities.

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

BIOANALYZE_SCRIPTS_ROOT = Path(__file__).resolve().parents[1]
if str(BIOANALYZE_SCRIPTS_ROOT) not in sys.path:
    sys.path.insert(0, str(BIOANALYZE_SCRIPTS_ROOT))

from bioanalyze_paths import get_bioanalyze_data_root, get_bioanalyze_figures_root
from normalized_expression_common import (
    build_species_heatmap_display_index,
    build_tissue_coverage_table,
    load_processed_expression_cells,
    sort_gene_labels,
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Build human H2A heatmaps from normalized H2A cell tables: "
            "all, clustered (canonical), and variants. "
            "Rows are unique Ensembl IDs; labels are merged gene_name:HGNC."
        )
    )
    p.add_argument(
        "--expr",
        default=str(
            get_bioanalyze_data_root()
            / "processed"
            / "homo_sapiens"
            / "Homo_sapiens_expr_advanced_H2A_present_gold.tsv"
        ),
        help="Normalized H2A cell TSV for Homo sapiens.",
    )
    p.add_argument(
        "--h2a-merged",
        default=str(get_bioanalyze_data_root() / "merged" / "mammalia_H2A_merged_with_taxonomy_v4.csv"),
        help="Merged v4 H2A dataset with taxonomy, gene_name, and HGNC IDs.",
    )
    p.add_argument(
        "--out-dir",
        default=str(get_bioanalyze_figures_root() / "heatmaps" / "species" / "human"),
        help="Output directory for heatmap images.",
    )
    p.add_argument(
        "--out-map",
        default=str(
            get_bioanalyze_data_root() / "processed" / "homo_sapiens" / "h2a_hs_canonical_variant_map.tsv"
        ),
        help="Output TSV for ENSG->label/class mapping.",
    )
    p.add_argument(
        "--show-all-xticks",
        action="store_true",
        help="Force all anatomical entity labels on X axis (default: True).",
    )
    p.add_argument(
        "--min-tissue-fill-rate",
        type=float,
        default=0.0,
        help=(
            "Keep only anatomical entities with non-missing coverage >= this fraction "
            "within each panel. observed_zero counts as filled."
        ),
    )
    return p.parse_args()


def build_label_maps(
    expr_df: pd.DataFrame, h2a_df: pd.DataFrame
) -> Tuple[Dict[str, str], Dict[str, str]]:
    # Canonical gene_name from merged v4 (stable per ENSG: longest non-empty name)
    name_map = (
        h2a_df[["ensembl_gene_id", "gene_name"]]
        .assign(_len=lambda x: x["gene_name"].str.len())
        .sort_values(["ensembl_gene_id", "_len"], ascending=[True, False])
        .drop_duplicates("ensembl_gene_id", keep="first")
        .set_index("ensembl_gene_id")["gene_name"]
        .to_dict()
    )

    # HGNC map from v3
    hgnc_map = (
        h2a_df[["ensembl_gene_id", "hgnc_id"]]
        .sort_values(["ensembl_gene_id", "hgnc_id"])
        .drop_duplicates("ensembl_gene_id", keep="first")
        .set_index("ensembl_gene_id")["hgnc_id"]
        .to_dict()
    )

    label_map: Dict[str, str] = {}
    for gid in sorted(
        {
            gid.strip()
            for gid in h2a_df["ensembl_gene_id"].dropna().astype(str).tolist()
            if gid.strip()
        }
    ):
        gname = (name_map.get(gid, "") or "").strip() or gid
        hgnc = (hgnc_map.get(gid, "") or "").strip()
        label_map[gid] = f"{gname}:{hgnc}" if hgnc else f"{gname}:{gid}"

    # Ensure uniqueness
    seen: Dict[str, int] = {}
    for gid, lab in list(label_map.items()):
        if lab in seen:
            seen[lab] += 1
            label_map[gid] = f"{lab}#{seen[lab]}"
        else:
            seen[lab] = 1

    return label_map, name_map


def classify_ensg(h2a_hs: pd.DataFrame) -> Dict[str, str]:
    # canonical = variant == "clustered H2A"
    cls_map: Dict[str, str] = {}
    for gid, grp in h2a_hs.groupby("ensembl_gene_id", dropna=False):
        if not isinstance(gid, str) or not gid.strip():
            continue
        variants = {v.strip() for v in grp["variant"].dropna().astype(str).tolist()}
        if "clustered H2A" in variants:
            cls = "clustered"
        else:
            cls = "variant"
        cls_map[gid] = cls
    return cls_map


def heatmap(
    df: pd.DataFrame,
    label_map: Dict[str, str],
    gene_ids: list[str],
    tissue_names: list[str],
    title: str,
    out_png: Path,
    out_svg: Path,
    show_all_xticks: bool,
) -> Tuple[int, int]:
    df = df.copy()
    df["Gene row label"] = df["Gene ID"].map(label_map)
    mat = (
        df.pivot_table(
            index="Gene row label",
            columns="Anatomical entity name",
            values="cell_mean_score",
            aggfunc="mean",
        )
    )
    mat = mat.reindex(index=sort_gene_labels(label_map, gene_ids), columns=tissue_names)
    if mat.shape[0] == 0 or mat.shape[1] == 0 or mat.notna().sum().sum() == 0:
        raise RuntimeError(f"Empty matrix for {title}: {mat.shape}")

    mat_log = np.log10(mat + 1)

    sns.set(style="whitegrid")
    fig_w = max(48, mat_log.shape[1] * 0.25)
    fig_h = max(10, mat_log.shape[0] * 0.5)

    plt.figure(figsize=(fig_w, fig_h))
    ax = sns.heatmap(
        mat_log,
        cmap="viridis",
        mask=mat_log.isna(),
        linewidths=0.1,
        cbar_kws={"label": "log10(Expression score + 1)"},
        xticklabels=1 if show_all_xticks else "auto",
        yticklabels=1,
    )
    ax.set_title(title)
    ax.set_xlabel("Anatomical entity name")
    ax.set_ylabel("GeneName:HGNC (one row per displayed gene)")
    plt.xticks(rotation=90, fontsize=5 if show_all_xticks else 6)
    plt.yticks(rotation=0, fontsize=8)
    plt.tight_layout()

    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_svg, bbox_inches="tight")
    plt.close()

    return mat.shape[0], mat.shape[1]


def main() -> None:
    args = parse_args()

    expr_path = Path(args.expr)
    h2a_path = Path(args.h2a_merged)
    out_dir = Path(args.out_dir)
    out_map = Path(args.out_map)
    show_all_xticks = True if args.show_all_xticks else True
    min_tissue_fill_rate = float(args.min_tissue_fill_rate)

    if not 0.0 <= min_tissue_fill_rate <= 1.0:
        raise ValueError("min_tissue_fill_rate must be between 0.0 and 1.0")

    if not expr_path.exists():
        raise FileNotFoundError(expr_path)
    if not h2a_path.exists():
        raise FileNotFoundError(h2a_path)

    expr = load_processed_expression_cells(expr_path)

    h2a = pd.read_csv(h2a_path, dtype=str)
    required_h2a_cols = {"species_name", "ensembl_gene_id", "hgnc_id", "variant"}
    required_h2a_cols.add("gene_name")
    missing = required_h2a_cols - set(h2a.columns)
    if missing:
        raise RuntimeError(f"Missing columns in h2a merged: {missing}")

    h2a_hs = h2a[h2a["species_name"].eq("Homo sapiens")].copy()
    h2a_hs["ensembl_gene_id"] = h2a_hs["ensembl_gene_id"].fillna("").astype(str).str.strip()
    h2a_hs["gene_name"] = h2a_hs["gene_name"].fillna("").astype(str).str.strip()
    h2a_hs["hgnc_id"] = h2a_hs["hgnc_id"].fillna("").astype(str).str.strip()
    h2a_hs["variant"] = h2a_hs["variant"].fillna("").astype(str).str.strip()
    h2a_hs = h2a_hs[h2a_hs["ensembl_gene_id"] != ""].copy()

    cls_map = classify_ensg(h2a_hs)

    label_map, name_map = build_label_maps(expr, h2a_hs)

    # Build mapping table (for reproducibility)
    map_rows = []
    for gid, cls in sorted(cls_map.items()):
        map_rows.append(
            {
                "ensembl_gene_id": gid,
                "gene_name": name_map.get(gid, ""),
                "hgnc_id": h2a_hs.loc[
                    h2a_hs["ensembl_gene_id"] == gid, "hgnc_id"
                ].dropna().astype(str).head(1).squeeze()
                if gid in h2a_hs["ensembl_gene_id"].values
                else "",
                "class": cls,
                "label": label_map.get(gid, ""),
            }
        )
    out_map.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(map_rows).to_csv(out_map, sep="\t", index=False)

    display_map_rows = pd.DataFrame(map_rows)
    display_df = build_species_heatmap_display_index(
        display_map_rows,
        expr,
        preferred_id_col="hgnc_id",
    )
    display_gene_ids = display_df["ensembl_gene_id"].dropna().astype(str).tolist()
    if not display_gene_ids:
        raise RuntimeError("No displayable human heatmap rows remained after deduplication.")

    display_label_map = dict(
        zip(
            display_df["ensembl_gene_id"].astype(str),
            display_df["label"].astype(str),
        )
    )
    expr = expr[expr["Gene ID"].isin(display_gene_ids)].copy()
    expr["class"] = expr["Gene ID"].map(display_df.set_index("ensembl_gene_id")["class"].to_dict())
    canonical = expr[expr["class"].eq("clustered")].copy()
    variants = expr[expr["class"].eq("variant")].copy()
    tissue_names = sorted(expr["Anatomical entity name"].dropna().astype(str).unique().tolist())
    canonical_gene_ids = display_df.loc[
        display_df["class"].eq("clustered"), "ensembl_gene_id"
    ].astype(str).tolist()
    variant_gene_ids = display_df.loc[
        display_df["class"].eq("variant"), "ensembl_gene_id"
    ].astype(str).tolist()

    if canonical.empty:
        raise RuntimeError("No clustered (canonical) rows after filtering.")
    if variants.empty:
        raise RuntimeError("No variant rows after filtering.")

    def panel_tissues_and_report(
        panel_df: pd.DataFrame,
        panel_gene_ids: list[str],
        panel_key: str,
    ) -> list[str]:
        coverage_df = build_tissue_coverage_table(
            panel_df,
            panel_gene_ids,
            threshold=min_tissue_fill_rate,
            panel=panel_key,
        )
        if min_tissue_fill_rate > 0:
            coverage_tsv = out_dir / f"h2a_{panel_key}_tissue_coverage.tsv"
            out_dir.mkdir(parents=True, exist_ok=True)
            coverage_df.to_csv(coverage_tsv, sep="\t", index=False)
        kept_tissues = coverage_df.loc[coverage_df["kept"], "anatomical_entity_name"].astype(str).tolist()
        if min_tissue_fill_rate > 0 and not kept_tissues:
            raise RuntimeError(
                f"No tissues passed min_tissue_fill_rate={min_tissue_fill_rate:.2f} for panel '{panel_key}'."
            )
        if min_tissue_fill_rate > 0:
            return kept_tissues
        return coverage_df["anatomical_entity_name"].astype(str).tolist()

    all_gene_ids = display_df["ensembl_gene_id"].astype(str).tolist()
    all_tissue_names = panel_tissues_and_report(expr, all_gene_ids, "all")
    canonical_tissue_names = panel_tissues_and_report(canonical, canonical_gene_ids, "clustered")
    variant_tissue_names = panel_tissues_and_report(variants, variant_gene_ids, "variants")

    all_png = out_dir / "h2a_all.png"
    all_svg = out_dir / "h2a_all.svg"
    can_png = out_dir / "h2a_clustered.png"
    can_svg = out_dir / "h2a_clustered.svg"
    var_png = out_dir / "h2a_variants.png"
    var_svg = out_dir / "h2a_variants.svg"

    all_rows, all_cols = heatmap(
        expr,
        display_label_map,
        all_gene_ids,
        all_tissue_names,
        "H2A Human Expression (normalized cells) - all",
        all_png,
        all_svg,
        show_all_xticks,
    )
    can_rows, can_cols = heatmap(
        canonical,
        display_label_map,
        canonical_gene_ids,
        canonical_tissue_names,
        "H2A Human Expression (normalized cells) - clustered H2A",
        can_png,
        can_svg,
        show_all_xticks,
    )
    var_rows, var_cols = heatmap(
        variants,
        display_label_map,
        variant_gene_ids,
        variant_tissue_names,
        "H2A Human Expression (normalized cells) - variants",
        var_png,
        var_svg,
        show_all_xticks,
    )

    print(f"ALL_ROWS={all_rows} ALL_COLS={all_cols}")
    print(f"CANONICAL_ROWS={can_rows} CANONICAL_COLS={can_cols}")
    print(f"VARIANT_ROWS={var_rows} VARIANT_COLS={var_cols}")
    print(f"MAP_TSV={out_map}")
    print(f"ALL_SVG={all_svg}")
    print(f"CANONICAL_SVG={can_svg}")
    print(f"VARIANT_SVG={var_svg}")


if __name__ == "__main__":
    main()
