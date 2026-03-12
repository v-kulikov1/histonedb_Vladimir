#!/usr/bin/env python
# Build H2A expression heatmaps for any species from Bgee advanced TSV.
# Outputs: all H2A, clustered-only, variants-only heatmaps and mapping table.

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
            "Build H2A expression heatmaps for a given species from Bgee advanced TSV. "
            "Filters to present+gold, rows=Ensembl ID, labels=GeneName:<ID> (fallback ENSG)."
        )
    )
    p.add_argument(
        "--species",
        required=True,
        help='Species name as in v3 file, e.g. "Homo sapiens".',
    )
    p.add_argument(
        "--expr",
        required=True,
        help="Bgee advanced TSV for the species (all conditions).",
    )
    p.add_argument(
        "--merged",
        default=r"CURATED_SET/BioAnalyze/data/merged/mammalia_H2A_merged_with_taxonomy_v3.csv",
        help="Merged v3 H2A dataset with taxonomy, HGNC IDs, and variant.",
    )
    p.add_argument(
        "--out-dir",
        default=r"CURATED_SET/BioAnalyze/figures/heatmaps",
        help="Base output directory for heatmaps (species subfolder is appended).",
    )
    p.add_argument(
        "--out-processed-dir",
        default=r"CURATED_SET/BioAnalyze/data/processed",
        help="Output directory for processed TSVs.",
    )
    p.add_argument(
        "--id-col",
        default="auto",
        choices=["auto", "hgnc_id", "vgnc_id"],
        help=(
            "ID column to use in row labels. "
            "auto = hgnc_id for Homo sapiens, otherwise vgnc_id."
        ),
    )
    p.add_argument(
        "--chunksize",
        type=int,
        default=200000,
        help="Chunk size for streaming Bgee TSV.",
    )
    p.add_argument(
        "--square-cells",
        action="store_true",
        help="Use square-like cells by sizing figure from rows/cols.",
    )
    p.add_argument(
        "--cell-size",
        type=float,
        default=0.7,
        help="Cell size (inches) when --square-cells is set.",
    )
    p.add_argument(
        "--min-width",
        type=float,
        default=12.0,
        help="Minimum figure width (inches) when --square-cells is set.",
    )
    p.add_argument(
        "--min-height",
        type=float,
        default=8.0,
        help="Minimum figure height (inches) when --square-cells is set.",
    )
    return p.parse_args()


def slugify_species(species: str) -> str:
    return species.strip().lower().replace(" ", "_")


def classify_ensg(h2a_hs: pd.DataFrame) -> Dict[str, str]:
    cls_map: Dict[str, str] = {}
    for gid, grp in h2a_hs.groupby("ensembl_gene_id", dropna=False):
        if not isinstance(gid, str) or not gid.strip():
            continue
        variants = {v.strip() for v in grp["variant"].dropna().astype(str).tolist()}
        cls_map[gid] = "clustered" if "clustered H2A" in variants else "variant"
    return cls_map


def build_label_maps(
    expr_df: pd.DataFrame, h2a_df: pd.DataFrame, id_col: str
) -> Tuple[Dict[str, str], Dict[str, str]]:
    name_map = (
        expr_df[["Gene ID", "Gene name"]]
        .assign(_len=lambda x: x["Gene name"].str.len())
        .sort_values(["Gene ID", "_len"], ascending=[True, False])
        .drop_duplicates("Gene ID", keep="first")
        .set_index("Gene ID")["Gene name"]
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

    # Ensure uniqueness
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
    title: str,
    out_png: Path,
    out_svg: Path,
    id_label: str,
    square_cells: bool,
    cell_size: float,
    min_width: float,
    min_height: float,
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

    labels_sorted = sorted(mat.index.tolist(), key=lambda s: (s.split(":", 1)[0], s))
    mat = mat.reindex(labels_sorted)
    mat_log = np.log10(mat + 1)

    sns.set(style="whitegrid")
    if square_cells:
        fig_w = max(min_width, mat_log.shape[1] * cell_size)
        fig_h = max(min_height, mat_log.shape[0] * cell_size)
    else:
        fig_w = max(48, mat_log.shape[1] * 0.25)
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
    ax.set_xlabel("Anatomical entity name")
    ax.set_ylabel(f"GeneName:{id_label} (one row per Ensembl ID)")
    plt.xticks(rotation=90, fontsize=5)
    plt.yticks(rotation=0, fontsize=8)
    plt.tight_layout()

    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_svg, bbox_inches="tight")
    plt.close()

    return mat.shape[0], mat.shape[1]


def main() -> None:
    args = parse_args()
    species = args.species.strip()
    expr_path = Path(args.expr)
    merged_path = Path(args.merged)
    out_dir = Path(args.out_dir)
    out_processed_dir = Path(args.out_processed_dir)
    id_col = args.id_col.strip().lower()
    if id_col == "auto":
        id_col = "hgnc_id" if species == "Homo sapiens" else "vgnc_id"

    if not expr_path.exists():
        raise FileNotFoundError(expr_path)
    if not merged_path.exists():
        raise FileNotFoundError(merged_path)

    h2a = pd.read_csv(merged_path, dtype=str)
    required_cols = {"species_name", "ensembl_gene_id", "hgnc_id", "vgnc_id", "variant"}
    missing = required_cols - set(h2a.columns)
    if missing:
        raise RuntimeError(f"Missing columns in merged v3: {missing}")

    h2a_sp = h2a[h2a["species_name"].eq(species)].copy()
    h2a_sp["ensembl_gene_id"] = h2a_sp["ensembl_gene_id"].fillna("").astype(str).str.strip()
    h2a_sp["hgnc_id"] = h2a_sp["hgnc_id"].fillna("").astype(str).str.strip()
    h2a_sp["vgnc_id"] = h2a_sp["vgnc_id"].fillna("").astype(str).str.strip()
    h2a_sp["variant"] = h2a_sp["variant"].fillna("").astype(str).str.strip()
    h2a_sp = h2a_sp[h2a_sp["ensembl_gene_id"] != ""].copy()

    ensg_set = set(h2a_sp["ensembl_gene_id"].unique())
    if not ensg_set:
        raise RuntimeError(f"No ENSG for species: {species}")

    usecols = [
        "Gene ID",
        "Gene name",
        "Anatomical entity name",
        "Expression",
        "Call quality",
        "Expression score",
    ]
    chunks = []
    rows_total = 0
    rows_kept = 0
    for ch in pd.read_csv(expr_path, sep="\t", dtype=str, usecols=usecols, chunksize=args.chunksize, low_memory=True):
        rows_total += len(ch)
        ch["Gene ID"] = ch["Gene ID"].fillna("").astype(str).str.strip()
        ch["Gene name"] = ch["Gene name"].fillna("").astype(str).str.strip()
        ch["Expression"] = ch["Expression"].fillna("").astype(str).str.strip()
        ch["Call quality"] = ch["Call quality"].fillna("").astype(str).str.strip()

        keep = ch[
            ch["Gene ID"].isin(ensg_set)
            & ch["Expression"].eq("present")
            & ch["Call quality"].eq("gold quality")
        ].copy()
        if not keep.empty:
            rows_kept += len(keep)
            chunks.append(keep)

    if not chunks:
        raise RuntimeError("No rows after filtering to H2A + present + gold quality.")

    expr = pd.concat(chunks, ignore_index=True)
    expr["Expression score"] = pd.to_numeric(expr["Expression score"], errors="coerce")
    expr = expr[expr["Expression score"].notna()].copy()
    if expr.empty:
        raise RuntimeError("All Expression score are NaN after conversion.")

    slug = slugify_species(species)
    out_processed_dir.mkdir(parents=True, exist_ok=True)
    processed_tsv = out_processed_dir / f"{slug}_expr_advanced_H2A_present_gold.tsv"
    expr.to_csv(processed_tsv, sep="\t", index=False)

    cls_map = classify_ensg(h2a_sp)
    label_map, name_map = build_label_maps(expr, h2a_sp, id_col)

    # Mapping table for reproducibility
    map_rows = []
    for gid, cls in sorted(cls_map.items()):
        id_val = (
            h2a_sp.loc[h2a_sp["ensembl_gene_id"] == gid, id_col]
            .dropna()
            .astype(str)
            .head(1)
            .squeeze()
            if gid in h2a_sp["ensembl_gene_id"].values
            else ""
        )
        map_rows.append(
            {
                "species_name": species,
                "ensembl_gene_id": gid,
                "gene_name": name_map.get(gid, ""),
                "hgnc_id": (
                    h2a_sp.loc[h2a_sp["ensembl_gene_id"] == gid, "hgnc_id"]
                    .dropna()
                    .astype(str)
                    .head(1)
                    .squeeze()
                    if gid in h2a_sp["ensembl_gene_id"].values
                    else ""
                ),
                "vgnc_id": (
                    h2a_sp.loc[h2a_sp["ensembl_gene_id"] == gid, "vgnc_id"]
                    .dropna()
                    .astype(str)
                    .head(1)
                    .squeeze()
                    if gid in h2a_sp["ensembl_gene_id"].values
                    else ""
                ),
                "id_type": id_col,
                "id_value": id_val,
                "class": cls,
                "label": label_map.get(gid, ""),
            }
        )
    map_tsv = out_processed_dir / f"{slug}_h2a_canonical_variant_map.tsv"
    pd.DataFrame(map_rows).to_csv(map_tsv, sep="\t", index=False)

    expr["class"] = expr["Gene ID"].map(cls_map)
    canonical = expr[expr["class"].eq("clustered")].copy()
    variants = expr[expr["class"].eq("variant")].copy()

    if canonical.empty:
        raise RuntimeError("No clustered (canonical) rows after filtering.")
    if variants.empty:
        raise RuntimeError("No variant rows after filtering.")

    out_dir = out_dir if out_dir.name == slug else out_dir / slug
    out_dir.mkdir(parents=True, exist_ok=True)
    all_png = out_dir / f"h2a_{slug}_all.png"
    all_svg = out_dir / f"h2a_{slug}_all.svg"
    can_png = out_dir / f"h2a_{slug}_clustered.png"
    can_svg = out_dir / f"h2a_{slug}_clustered.svg"
    var_png = out_dir / f"h2a_{slug}_variants.png"
    var_svg = out_dir / f"h2a_{slug}_variants.svg"

    all_rows, all_cols = heatmap(
        expr,
        label_map,
        f"H2A Expression (present + gold) - {species} (all)",
        all_png,
        all_svg,
        "HGNC" if id_col == "hgnc_id" else "VGNC",
        args.square_cells,
        args.cell_size,
        args.min_width,
        args.min_height,
    )
    can_rows, can_cols = heatmap(
        canonical,
        label_map,
        f"H2A Expression (present + gold) - {species} (clustered H2A)",
        can_png,
        can_svg,
        "HGNC" if id_col == "hgnc_id" else "VGNC",
        args.square_cells,
        args.cell_size,
        args.min_width,
        args.min_height,
    )
    var_rows, var_cols = heatmap(
        variants,
        label_map,
        f"H2A Expression (present + gold) - {species} (variants)",
        var_png,
        var_svg,
        "HGNC" if id_col == "hgnc_id" else "VGNC",
        args.square_cells,
        args.cell_size,
        args.min_width,
        args.min_height,
    )

    print(f"SPECIES={species}")
    print(f"ID_COL={id_col}")
    print(f"ROWS_TOTAL_SCANNED={rows_total}")
    print(f"ROWS_AFTER_FILTER={rows_kept}")
    print(f"ALL_ROWS={all_rows} ALL_COLS={all_cols}")
    print(f"CANONICAL_ROWS={can_rows} CANONICAL_COLS={can_cols}")
    print(f"VARIANT_ROWS={var_rows} VARIANT_COLS={var_cols}")
    print(f"PROCESSED_TSV={processed_tsv}")
    print(f"MAP_TSV={map_tsv}")
    print(f"ALL_SVG={all_svg}")
    print(f"CANONICAL_SVG={can_svg}")
    print(f"VARIANT_SVG={var_svg}")


if __name__ == "__main__":
    main()
