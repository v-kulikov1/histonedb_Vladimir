#!/usr/bin/env python
"""
Build heatmaps for the dataset with more rows, using only intersections across:
- tissues (Anatomical entity name)
- histones (Gene name)

No normalization is applied. Intersections are exact string matches.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from normalized_expression_common import load_processed_expression_cells


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Build heatmaps for the larger of two normalized H2A datasets, "
            "keeping only tissue and merged gene_name intersections."
        )
    )
    p.add_argument(
        "--a-present-gold",
        default=r"CURATED_SET/BioAnalyze/data/processed/homo_sapiens/Homo_sapiens_expr_advanced_H2A_present_gold.tsv",
        help="Normalized H2A TSV for dataset A.",
    )
    p.add_argument(
        "--b-present-gold",
        default=r"CURATED_SET/BioAnalyze/data/processed/pan_troglodytes/pan_troglodytes_expr_advanced_H2A_present_gold.tsv",
        help="Normalized H2A TSV for dataset B.",
    )
    p.add_argument(
        "--a-species",
        default="Homo sapiens",
        help="Species name for dataset A (as in merged v3).",
    )
    p.add_argument(
        "--b-species",
        default="Pan troglodytes",
        help="Species name for dataset B (as in merged v3).",
    )
    p.add_argument(
        "--merged",
        default=r"CURATED_SET/BioAnalyze/data/merged/mammalia_H2A_merged_with_taxonomy_v4.csv",
        help="Merged v4 dataset with gene_name, variant, and IDs.",
    )
    p.add_argument(
        "--out-dir",
        default=r"CURATED_SET/BioAnalyze/figures/heatmaps/alligned_human_pan",
        help="Output directory for heatmaps.",
    )
    p.add_argument(
        "--out-processed-dir",
        default=r"CURATED_SET/BioAnalyze/data/processed/intersections",
        help="Output directory for processed TSVs.",
    )
    p.add_argument(
        "--out-prefix",
        default="",
        help="Output prefix for heatmap files (default: intersection_<primary_slug>).",
    )
    p.add_argument(
        "--id-col",
        default="auto",
        choices=["auto", "hgnc_id", "vgnc_id"],
        help="ID column for row labels. auto = HGNC for human, otherwise VGNC.",
    )
    return p.parse_args()


def slugify_species(species: str) -> str:
    return species.strip().lower().replace(" ", "_")


def classify_ensg(h2a_sp: pd.DataFrame) -> Dict[str, str]:
    cls_map: Dict[str, str] = {}
    for gid, grp in h2a_sp.groupby("ensembl_gene_id", dropna=False):
        if not isinstance(gid, str) or not gid.strip():
            continue
        variants = {v.strip() for v in grp["variant"].dropna().astype(str).tolist()}
        cls_map[gid] = "clustered" if "clustered H2A" in variants else "variant"
    return cls_map


def build_label_maps(
    expr_df: pd.DataFrame, h2a_df: pd.DataFrame, id_col: str
) -> Tuple[Dict[str, str], Dict[str, str]]:
    name_map = (
        h2a_df[["ensembl_gene_id", "gene_name"]]
        .assign(_len=lambda x: x["gene_name"].str.len())
        .sort_values(["ensembl_gene_id", "_len"], ascending=[True, False])
        .drop_duplicates("ensembl_gene_id", keep="first")
        .set_index("ensembl_gene_id")["gene_name"]
        .to_dict()
    )
    id_map = (
        h2a_df[["ensembl_gene_id", id_col]]
        .sort_values(["ensembl_gene_id", id_col])
        .drop_duplicates("ensembl_gene_id", keep="first")
        .set_index("ensembl_gene_id")[id_col]
        .to_dict()
    )

    label_map: Dict[str, str] = {}
    for gid in sorted(expr_df["Gene ID"].unique()):
        gname = (name_map.get(gid, "") or "").strip() or gid
        id_val = (id_map.get(gid, "") or "").strip()
        label_map[gid] = f"{gname}:{id_val}" if id_val else f"{gname}:{gid}"

    seen: Dict[str, int] = {}
    for gid, lab in list(label_map.items()):
        if lab in seen:
            seen[lab] += 1
            label_map[gid] = f"{lab}#{seen[lab]}"
        else:
            seen[lab] = 1

    return label_map, name_map


def heatmap(
    df: pd.DataFrame,
    label_map: Dict[str, str],
    target_names: List[str],
    title: str,
    out_png: Path,
    out_svg: Path,
    id_label: str,
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
        .reindex(columns=target_names)
    )
    if mat.shape[0] == 0 or mat.shape[1] == 0 or mat.notna().sum().sum() == 0:
        raise RuntimeError(f"Empty matrix for {title}: {mat.shape}")

    labels_sorted = sorted(mat.index.tolist(), key=lambda s: (s.split(":", 1)[0], s))
    mat = mat.reindex(labels_sorted)
    mat_log = np.log10(mat + 1)

    sns.set(style="whitegrid")
    fig_w = max(24, len(target_names) * 0.6)
    fig_h = max(10, mat_log.shape[0] * 0.5)

    plt.figure(figsize=(fig_w, fig_h))
    ax = sns.heatmap(
        mat_log,
        cmap="viridis",
        mask=mat_log.isna(),
        linewidths=0.1,
        cbar_kws={"label": "log10(Expression score + 1)"},
        xticklabels=1,
        yticklabels=1,
    )
    ax.set_title(title)
    ax.set_xlabel("Anatomical entity name (intersection)")
    ax.set_ylabel(f"GeneName:{id_label} (one row per Ensembl ID)")
    plt.xticks(rotation=90, fontsize=6)
    plt.yticks(rotation=0, fontsize=8)
    plt.tight_layout()

    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_svg, bbox_inches="tight")
    plt.close()

    return mat.shape[0], mat.shape[1]


def main() -> None:
    args = parse_args()
    a_path = Path(args.a_present_gold)
    b_path = Path(args.b_present_gold)
    merged_path = Path(args.merged)
    out_dir = Path(args.out_dir)
    out_processed_dir = Path(args.out_processed_dir)

    for p in [a_path, b_path, merged_path]:
        if not p.exists():
            raise FileNotFoundError(p)

    a_df = load_processed_expression_cells(a_path)
    b_df = load_processed_expression_cells(b_path)

    a_rows = len(a_df)
    b_rows = len(b_df)
    if a_rows >= b_rows:
        primary_df, secondary_df = a_df, b_df
        primary_species, secondary_species = args.a_species, args.b_species
    else:
        primary_df, secondary_df = b_df, a_df
        primary_species, secondary_species = args.b_species, args.a_species

    h2a = pd.read_csv(merged_path, dtype=str)
    required_h2a_cols = {
        "species_name",
        "ensembl_gene_id",
        "gene_name",
        "variant",
        "hgnc_id",
        "vgnc_id",
    }
    missing = required_h2a_cols - set(h2a.columns)
    if missing:
        raise RuntimeError(f"Missing columns in merged v4: {missing}")

    for col in ["species_name", "ensembl_gene_id", "gene_name", "variant", "hgnc_id", "vgnc_id"]:
        h2a[col] = h2a[col].fillna("").astype(str).str.strip()

    def attach_merged_gene_name(df: pd.DataFrame, species: str) -> pd.DataFrame:
        name_map_df = (
            h2a[h2a["species_name"].eq(species)][["ensembl_gene_id", "gene_name"]]
            .assign(_len=lambda x: x["gene_name"].str.len())
            .sort_values(["ensembl_gene_id", "_len"], ascending=[True, False])
            .drop_duplicates("ensembl_gene_id", keep="first")
            .drop(columns=["_len"])
        )
        out = df.merge(
            name_map_df,
            left_on="Gene ID",
            right_on="ensembl_gene_id",
            how="left",
        )
        out = out.rename(columns={"gene_name": "merged_gene_name"})
        out["merged_gene_name"] = (
            out["merged_gene_name"].fillna("").astype(str).str.strip()
        )
        return out.drop(columns=["ensembl_gene_id"])

    primary_df = attach_merged_gene_name(primary_df, primary_species)
    secondary_df = attach_merged_gene_name(secondary_df, secondary_species)

    tissues = sorted(
        set(primary_df["Anatomical entity name"]) & set(secondary_df["Anatomical entity name"])
    )
    genes = set(primary_df["merged_gene_name"]) & set(secondary_df["merged_gene_name"])
    genes.discard("")

    primary_df = primary_df[
        primary_df["Anatomical entity name"].isin(tissues)
        & primary_df["merged_gene_name"].isin(genes)
    ].copy()

    primary_df["Expression score"] = pd.to_numeric(
        primary_df["cell_mean_score"], errors="coerce"
    )
    primary_df["cell_mean_score"] = primary_df["Expression score"]
    primary_df = primary_df[primary_df["cell_mean_score"].notna()].copy()
    if primary_df.empty or not tissues or not genes:
        raise RuntimeError("No overlapping tissues/merged gene_name values after intersection.")

    h2a_sp = h2a[h2a["species_name"].eq(primary_species)].copy()
    h2a_sp = h2a_sp[h2a_sp["ensembl_gene_id"] != ""].copy()

    id_col = args.id_col
    if id_col == "auto":
        id_col = "hgnc_id" if primary_species == "Homo sapiens" else "vgnc_id"

    label_map, _ = build_label_maps(primary_df, h2a_sp, id_col)
    cls_map = classify_ensg(h2a_sp)
    primary_df["class"] = primary_df["Gene ID"].map(cls_map)

    canonical = primary_df[primary_df["class"].eq("clustered")].copy()
    variants = primary_df[primary_df["class"].eq("variant")].copy()
    if canonical.empty:
        raise RuntimeError("No clustered (canonical) rows after filtering.")
    if variants.empty:
        raise RuntimeError("No variant rows after filtering.")

    primary_slug = slugify_species(primary_species)
    out_prefix = args.out_prefix.strip() or f"intersection_{primary_slug}"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_processed_dir.mkdir(parents=True, exist_ok=True)

    processed_tsv = (
        out_processed_dir / f"{primary_slug}_expr_advanced_H2A_present_gold_intersection.tsv"
    )
    primary_df.to_csv(processed_tsv, sep="\t", index=False)

    all_png = out_dir / f"{out_prefix}_all.png"
    all_svg = out_dir / f"{out_prefix}_all.svg"
    can_png = out_dir / f"{out_prefix}_clustered.png"
    can_svg = out_dir / f"{out_prefix}_clustered.svg"
    var_png = out_dir / f"{out_prefix}_variants.png"
    var_svg = out_dir / f"{out_prefix}_variants.svg"

    all_rows, all_cols = heatmap(
        primary_df,
        label_map,
        tissues,
        f"H2A Expression (normalized cells) - {primary_species} (intersection with {secondary_species})",
        all_png,
        all_svg,
        "HGNC" if id_col == "hgnc_id" else "VGNC",
    )
    can_rows, can_cols = heatmap(
        canonical,
        label_map,
        tissues,
        f"H2A Expression (normalized cells) - {primary_species} (intersection, clustered H2A)",
        can_png,
        can_svg,
        "HGNC" if id_col == "hgnc_id" else "VGNC",
    )
    var_rows, var_cols = heatmap(
        variants,
        label_map,
        tissues,
        f"H2A Expression (normalized cells) - {primary_species} (intersection, variants)",
        var_png,
        var_svg,
        "HGNC" if id_col == "hgnc_id" else "VGNC",
    )

    print(f"PRIMARY_SPECIES={primary_species}")
    print(f"SECONDARY_SPECIES={secondary_species}")
    print(f"A_ROWS={a_rows} B_ROWS={b_rows}")
    print(f"TISSUE_INTERSECTION={len(tissues)}")
    print(f"GENE_INTERSECTION={len(genes)}")
    print(f"ALL_ROWS={all_rows} ALL_COLS={all_cols}")
    print(f"CANONICAL_ROWS={can_rows} CANONICAL_COLS={can_cols}")
    print(f"VARIANT_ROWS={var_rows} VARIANT_COLS={var_cols}")
    print(f"PROCESSED_TSV={processed_tsv}")
    print(f"ALL_SVG={all_svg}")
    print(f"CANONICAL_SVG={can_svg}")
    print(f"VARIANT_SVG={var_svg}")


if __name__ == "__main__":
    main()
