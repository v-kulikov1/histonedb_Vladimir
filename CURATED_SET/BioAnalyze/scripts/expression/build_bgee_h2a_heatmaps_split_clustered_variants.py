#!/usr/bin/env python
# Build Bgee H2A expression heatmaps split into clustered (canonical) vs variants.
# Rows: unique Ensembl Gene ID (labels: GeneName:HGNC). Columns: anatomical entities.

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Build two heatmaps from Bgee advanced H2A data: "
            "clustered (canonical) vs variants. "
            "Rows are unique Ensembl IDs; labels are merged gene_name:HGNC."
        )
    )
    p.add_argument(
        "--expr",
        default=r"CURATED_SET/BioAnalyze/data/processed/homo_sapiens/Homo_sapiens_expr_advanced_H2A_present_gold.tsv",
        help="H2A-only Bgee advanced TSV (preferably present+gold filtered).",
    )
    p.add_argument(
        "--h2a-merged",
        default=r"CURATED_SET/BioAnalyze/data/merged/mammalia_H2A_merged_with_taxonomy_v4.csv",
        help="Merged v4 H2A dataset with taxonomy, gene_name, and HGNC IDs.",
    )
    p.add_argument(
        "--out-dir",
        default=r"CURATED_SET/BioAnalyze/figures/heatmaps/human",
        help="Output directory for heatmap images.",
    )
    p.add_argument(
        "--out-map",
        default=r"CURATED_SET/BioAnalyze/data/processed/homo_sapiens/h2a_hs_canonical_variant_map.tsv",
        help="Output TSV for ENSG->label/class mapping.",
    )
    p.add_argument(
        "--show-all-xticks",
        action="store_true",
        help="Force all anatomical entity labels on X axis (default: True).",
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
            values="Expression score",
            aggfunc="mean",
        )
        .dropna(axis=0, how="all")
        .dropna(axis=1, how="all")
    )
    if mat.shape[0] == 0 or mat.shape[1] == 0:
        raise RuntimeError(f"Empty matrix for {title}: {mat.shape}")

    # Sort rows by Gene name prefix (before ':')
    labels_sorted = sorted(mat.index.tolist(), key=lambda s: (s.split(":", 1)[0], s))
    mat = mat.reindex(labels_sorted)

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
    ax.set_ylabel("GeneName:HGNC (one row per Ensembl ID)")
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

    if not expr_path.exists():
        raise FileNotFoundError(expr_path)
    if not h2a_path.exists():
        raise FileNotFoundError(h2a_path)

    expr = pd.read_csv(expr_path, sep="\t", dtype=str, low_memory=True)
    for c in ["Gene ID", "Gene name", "Anatomical entity name"]:
        if c not in expr.columns:
            raise RuntimeError(f"Missing column in expr: {c}")
        expr[c] = expr[c].fillna("").astype(str).str.strip()
    if "Expression score" not in expr.columns:
        raise RuntimeError("Missing column in expr: Expression score")
    expr["Expression score"] = pd.to_numeric(expr["Expression score"], errors="coerce")
    expr = expr[expr["Expression score"].notna()].copy()

    # If present/gold columns exist, enforce filter (safe even if already filtered)
    if "Expression" in expr.columns and "Call quality" in expr.columns:
        expr["Expression"] = expr["Expression"].fillna("").astype(str).str.strip()
        expr["Call quality"] = expr["Call quality"].fillna("").astype(str).str.strip()
        expr = expr[
            expr["Expression"].eq("present") & expr["Call quality"].eq("gold quality")
        ].copy()

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

    expr["class"] = expr["Gene ID"].map(cls_map)
    canonical = expr[expr["class"].eq("clustered")].copy()
    variants = expr[expr["class"].eq("variant")].copy()

    if canonical.empty:
        raise RuntimeError("No clustered (canonical) rows after filtering.")
    if variants.empty:
        raise RuntimeError("No variant rows after filtering.")

    can_png = out_dir / "h2a_clustered.png"
    can_svg = out_dir / "h2a_clustered.svg"
    var_png = out_dir / "h2a_variants.png"
    var_svg = out_dir / "h2a_variants.svg"

    can_rows, can_cols = heatmap(
        canonical,
        label_map,
        "H2A Human Expression (present + gold) - clustered H2A",
        can_png,
        can_svg,
        show_all_xticks,
    )
    var_rows, var_cols = heatmap(
        variants,
        label_map,
        "H2A Human Expression (present + gold) - variants",
        var_png,
        var_svg,
        show_all_xticks,
    )

    print(f"CANONICAL_ROWS={can_rows} CANONICAL_COLS={can_cols}")
    print(f"VARIANT_ROWS={var_rows} VARIANT_COLS={var_cols}")
    print(f"MAP_TSV={out_map}")
    print(f"CANONICAL_SVG={can_svg}")
    print(f"VARIANT_SVG={var_svg}")


if __name__ == "__main__":
    main()
