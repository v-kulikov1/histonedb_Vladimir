#!/usr/bin/env python
# Build H2A expression heatmaps for any species from Bgee advanced TSV.
# Outputs: all H2A, clustered-only, variants-only heatmaps and mapping table.

from __future__ import annotations

import argparse
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

BIOANALYZE_SCRIPTS_ROOT = Path(__file__).resolve().parents[1]
if str(BIOANALYZE_SCRIPTS_ROOT) not in sys.path:
    sys.path.insert(0, str(BIOANALYZE_SCRIPTS_ROOT))

from bioanalyze_paths import get_bioanalyze_data_root, get_bioanalyze_figures_root
from normalized_expression_common import (
    build_species_heatmap_display_index,
    build_tissue_coverage_table,
    load_processed_expression_cells,
    normalize_h2a_expression_cells,
    sort_gene_labels,
)


DEFAULT_MERGED = get_bioanalyze_data_root() / "merged" / "mammalia_H2A_merged_with_taxonomy_v4.csv"
DEFAULT_OUT_DIR = get_bioanalyze_figures_root() / "heatmaps" / "species"
DEFAULT_OUT_PROCESSED_DIR = get_bioanalyze_data_root() / "processed"
CANONICAL_RULES = ("legacy", "canonical_like")


@dataclass
class BuildResult:
    species: str
    slug: str
    status: str
    reason: str
    id_col: str
    canonical_rule: str
    allow_partial_splits: bool
    rows_total_scanned: int = 0
    rows_after_filter: int = 0
    matched_genes_total: int = 0
    matched_genes_clustered: int = 0
    matched_genes_variant: int = 0
    all_rows: int = 0
    all_cols: int = 0
    canonical_rows: int = 0
    canonical_cols: int = 0
    variant_rows: int = 0
    variant_cols: int = 0
    processed_tsv: str = ""
    map_tsv: str = ""
    all_png: str = ""
    all_svg: str = ""
    canonical_png: str = ""
    canonical_svg: str = ""
    variant_png: str = ""
    variant_svg: str = ""

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Build H2A expression heatmaps for a given species from Bgee advanced TSV. "
            "Normalizes H2A Gene ID x tissue cells so present+gold keeps the real score, "
            "observed-but-nonqualifying cells become 0, and truly missing cells stay empty."
        )
    )
    p.add_argument(
        "--species",
        required=True,
        help='Species name as in merged v4 file, e.g. "Homo sapiens".',
    )
    p.add_argument(
        "--expr",
        required=True,
        help="Bgee advanced TSV for the species (all conditions).",
    )
    p.add_argument(
        "--merged",
        default=str(DEFAULT_MERGED),
        help="Merged v4 H2A dataset with taxonomy, gene_name, HGNC/VGNC IDs, and variant.",
    )
    p.add_argument(
        "--out-dir",
        default=str(DEFAULT_OUT_DIR),
        help="Base output directory for heatmaps (species subfolder is appended).",
    )
    p.add_argument(
        "--out-processed-dir",
        default=str(DEFAULT_OUT_PROCESSED_DIR),
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
        "--canonical-rule",
        default="legacy",
        choices=list(CANONICAL_RULES),
        help=(
            "legacy = canonical only when variant == clustered H2A. "
            "canonical_like = canonical when variant == clustered H2A or startswith cH2A."
        ),
    )
    p.add_argument(
        "--allow-partial-splits",
        action="store_true",
        help=(
            "Write all heatmap plus whichever split is available instead of failing "
            "when one of clustered/variants is empty."
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


def slugify_species(species: str) -> str:
    return species.strip().lower().replace(" ", "_")


def resolve_id_col(species: str, id_col: str) -> str:
    id_col = id_col.strip().lower()
    if id_col == "auto":
        return "hgnc_id" if species == "Homo sapiens" else "vgnc_id"
    return id_col


def classify_ensg(h2a_df: pd.DataFrame, canonical_rule: str) -> Dict[str, str]:
    cls_map: Dict[str, str] = {}
    for gid, grp in h2a_df.groupby("ensembl_gene_id", dropna=False):
        if not isinstance(gid, str) or not gid.strip():
            continue
        variants = {v.strip() for v in grp["variant"].dropna().astype(str).tolist()}
        is_clustered = "clustered H2A" in variants
        if canonical_rule == "canonical_like":
            is_clustered = is_clustered or any(v.startswith("cH2A") for v in variants)
        cls_map[gid] = "clustered" if is_clustered else "variant"
    return cls_map


def build_label_maps(
    h2a_df: pd.DataFrame, id_col: str
) -> Tuple[Dict[str, str], Dict[str, str]]:
    merged_name_map = (
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
    for gid in sorted(
        {
            gid.strip()
            for gid in h2a_df["ensembl_gene_id"].dropna().astype(str).tolist()
            if gid.strip()
        }
    ):
        gname = (merged_name_map.get(gid, "") or "").strip() or gid
        id_val = (id_map.get(gid, "") or "").strip()
        label_map[gid] = f"{gname}:{id_val}" if id_val else f"{gname}:{gid}"

    seen: Dict[str, int] = {}
    for gid, lab in list(label_map.items()):
        if lab in seen:
            seen[lab] += 1
            label_map[gid] = f"{lab}#{seen[lab]}"
        else:
            seen[lab] = 1

    return label_map, merged_name_map


def heatmap(
    df: pd.DataFrame,
    label_map: Dict[str, str],
    gene_ids: List[str],
    tissue_names: List[str],
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
            values="cell_mean_score",
            aggfunc="mean",
        )
    )
    ordered_labels = sort_gene_labels(label_map, gene_ids)
    mat = mat.reindex(index=ordered_labels, columns=tissue_names)
    if mat.shape[0] == 0 or mat.shape[1] == 0 or mat.notna().sum().sum() == 0:
        raise RuntimeError(f"Empty matrix for {title}: {mat.shape}")

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
    ax.set_ylabel(f"GeneName:{id_label} (one row per displayed gene)")
    plt.xticks(rotation=90, fontsize=5)
    plt.yticks(rotation=0, fontsize=8)
    plt.tight_layout()

    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_svg, bbox_inches="tight")
    plt.close()

    return mat.shape[0], mat.shape[1]


def ensure_species_dir(base_dir: Path, slug: str) -> Path:
    if base_dir.name == slug:
        return base_dir
    if base_dir.parent.name == slug:
        return base_dir
    return base_dir / slug


def load_species_h2a(merged_path: Path, species: str) -> pd.DataFrame:
    h2a = pd.read_csv(merged_path, dtype=str)
    required_cols = {
        "species_name",
        "ensembl_gene_id",
        "gene_name",
        "hgnc_id",
        "vgnc_id",
        "variant",
    }
    missing = required_cols - set(h2a.columns)
    if missing:
        raise RuntimeError(f"Missing columns in merged v4: {missing}")

    h2a_sp = h2a[h2a["species_name"].eq(species)].copy()
    for col in ["ensembl_gene_id", "gene_name", "hgnc_id", "vgnc_id", "variant"]:
        h2a_sp[col] = h2a_sp[col].fillna("").astype(str).str.strip()
    h2a_sp = h2a_sp[h2a_sp["ensembl_gene_id"] != ""].copy()
    return h2a_sp


def build_map_rows(
    species: str,
    h2a_sp: pd.DataFrame,
    cls_map: Dict[str, str],
    label_map: Dict[str, str],
    merged_name_map: Dict[str, str],
    id_col: str,
) -> pd.DataFrame:
    map_rows = []
    for gid, cls in sorted(cls_map.items()):
        subset = h2a_sp.loc[h2a_sp["ensembl_gene_id"] == gid]
        id_val = subset[id_col].dropna().astype(str).head(1).squeeze() if not subset.empty else ""
        map_rows.append(
            {
                "species_name": species,
                "ensembl_gene_id": gid,
                "gene_name": merged_name_map.get(gid, ""),
                "hgnc_id": subset["hgnc_id"].dropna().astype(str).head(1).squeeze() if not subset.empty else "",
                "vgnc_id": subset["vgnc_id"].dropna().astype(str).head(1).squeeze() if not subset.empty else "",
                "id_type": id_col,
                "id_value": id_val,
                "class": cls,
                "label": label_map.get(gid, ""),
            }
        )
    return pd.DataFrame(map_rows)


def build_species_heatmaps(
    *,
    species: str,
    expr_path: Path,
    merged_path: Path,
    out_dir: Path,
    out_processed_dir: Path,
    id_col: str = "auto",
    canonical_rule: str = "legacy",
    allow_partial_splits: bool = False,
    chunksize: int = 200000,
    square_cells: bool = False,
    cell_size: float = 0.7,
    min_width: float = 12.0,
    min_height: float = 8.0,
    min_tissue_fill_rate: float = 0.0,
) -> BuildResult:
    species = species.strip()
    slug = slugify_species(species)
    resolved_id_col = resolve_id_col(species, id_col)
    result = BuildResult(
        species=species,
        slug=slug,
        status="success",
        reason="",
        id_col=resolved_id_col,
        canonical_rule=canonical_rule,
        allow_partial_splits=allow_partial_splits,
    )

    if canonical_rule not in CANONICAL_RULES:
        raise ValueError(f"Unsupported canonical rule: {canonical_rule}")
    if not 0.0 <= min_tissue_fill_rate <= 1.0:
        raise ValueError("min_tissue_fill_rate must be between 0.0 and 1.0")
    if not expr_path.exists():
        raise FileNotFoundError(expr_path)
    if not merged_path.exists():
        raise FileNotFoundError(merged_path)

    h2a_sp = load_species_h2a(merged_path, species)
    ensg_set = set(h2a_sp["ensembl_gene_id"].unique())
    if not ensg_set:
        result.status = "skipped"
        result.reason = "no_ensembl_ids_in_merged_v4"
        if allow_partial_splits:
            return result
        raise RuntimeError(f"No ENSG for species: {species}")

    expr, rows_total, rows_kept = normalize_h2a_expression_cells(expr_path, ensg_set, chunksize)
    result.rows_total_scanned = rows_total
    result.rows_after_filter = rows_kept

    if expr.empty:
        result.status = "skipped"
        result.reason = "no_h2a_rows_after_join"
        if allow_partial_splits:
            return result
        raise RuntimeError("No raw H2A rows were found after joining the species Ensembl IDs.")

    species_processed_dir = ensure_species_dir(out_processed_dir, slug)
    species_processed_dir.mkdir(parents=True, exist_ok=True)
    processed_tsv = species_processed_dir / f"{slug}_expr_advanced_H2A_present_gold.tsv"
    expr.to_csv(processed_tsv, sep="\t", index=False)
    result.processed_tsv = str(processed_tsv)

    cls_map = classify_ensg(h2a_sp, canonical_rule)
    label_map, merged_name_map = build_label_maps(h2a_sp, resolved_id_col)
    map_tsv = species_processed_dir / f"{slug}_h2a_canonical_variant_map.tsv"
    map_df = build_map_rows(
        species,
        h2a_sp,
        cls_map,
        label_map,
        merged_name_map,
        resolved_id_col,
    )
    map_df.to_csv(map_tsv, sep="\t", index=False)
    result.map_tsv = str(map_tsv)

    expr = load_processed_expression_cells(processed_tsv)
    expr["class"] = expr["Gene ID"].map(cls_map)
    display_df = build_species_heatmap_display_index(
        map_df,
        expr,
        preferred_id_col=resolved_id_col,
    )
    display_gene_ids = display_df["ensembl_gene_id"].dropna().astype(str).tolist()
    if not display_gene_ids:
        result.status = "skipped"
        result.reason = "no_displayable_genes_after_dedup"
        if allow_partial_splits:
            return result
        raise RuntimeError("No displayable species heatmap rows remained after deduplication.")

    display_label_map = dict(
        zip(
            display_df["ensembl_gene_id"].astype(str),
            display_df["label"].astype(str),
        )
    )
    expr = expr[expr["Gene ID"].isin(display_gene_ids)].copy()
    expr["class"] = expr["Gene ID"].map(
        display_df.set_index("ensembl_gene_id")["class"].to_dict()
    )
    canonical = expr[expr["class"].eq("clustered")].copy()
    variants = expr[expr["class"].eq("variant")].copy()
    result.matched_genes_total = int(len(display_df))
    result.matched_genes_clustered = int(display_df["class"].eq("clustered").sum())
    result.matched_genes_variant = int(display_df["class"].eq("variant").sum())
    tissue_names = sorted(expr["Anatomical entity name"].dropna().astype(str).unique().tolist())
    all_gene_ids = display_df["ensembl_gene_id"].astype(str).tolist()
    canonical_gene_ids = display_df.loc[
        display_df["class"].eq("clustered"), "ensembl_gene_id"
    ].astype(str).tolist()
    variant_gene_ids = display_df.loc[
        display_df["class"].eq("variant"), "ensembl_gene_id"
    ].astype(str).tolist()

    if canonical.empty and variants.empty:
        result.status = "skipped"
        result.reason = "no_rows_classified_for_heatmap"
        if allow_partial_splits:
            return result
        raise RuntimeError("No classified rows available for heatmap generation.")

    if (canonical.empty or variants.empty) and not allow_partial_splits:
        if canonical.empty:
            raise RuntimeError("No clustered (canonical) rows after filtering.")
        raise RuntimeError("No variant rows after filtering.")

    species_out_dir = ensure_species_dir(out_dir, slug)
    species_out_dir.mkdir(parents=True, exist_ok=True)
    id_label = "HGNC" if resolved_id_col == "hgnc_id" else "VGNC"

    def panel_tissues_and_report(
        panel_df: pd.DataFrame,
        panel_gene_ids: List[str],
        panel_key: str,
    ) -> List[str]:
        coverage_df = build_tissue_coverage_table(
            panel_df,
            panel_gene_ids,
            threshold=min_tissue_fill_rate,
            panel=panel_key,
        )
        if min_tissue_fill_rate > 0:
            coverage_tsv = species_out_dir / f"h2a_{slug}_{panel_key}_tissue_coverage.tsv"
            coverage_df.to_csv(coverage_tsv, sep="\t", index=False)
        kept_tissues = coverage_df.loc[coverage_df["kept"], "anatomical_entity_name"].astype(str).tolist()
        if min_tissue_fill_rate > 0 and not kept_tissues:
            raise RuntimeError(
                f"No tissues passed min_tissue_fill_rate={min_tissue_fill_rate:.2f} for panel '{panel_key}'."
            )
        if min_tissue_fill_rate > 0:
            return kept_tissues
        return coverage_df["anatomical_entity_name"].astype(str).tolist()

    all_png = species_out_dir / f"h2a_{slug}_all.png"
    all_svg = species_out_dir / f"h2a_{slug}_all.svg"
    all_tissue_names = panel_tissues_and_report(expr, all_gene_ids, "all")
    result.all_rows, result.all_cols = heatmap(
        expr,
        display_label_map,
        all_gene_ids,
        all_tissue_names,
        f"H2A Expression (normalized cells) - {species} (all)",
        all_png,
        all_svg,
        id_label,
        square_cells,
        cell_size,
        min_width,
        min_height,
    )
    result.all_png = str(all_png)
    result.all_svg = str(all_svg)

    if not canonical.empty:
        can_png = species_out_dir / f"h2a_{slug}_clustered.png"
        can_svg = species_out_dir / f"h2a_{slug}_clustered.svg"
        canonical_tissue_names = panel_tissues_and_report(canonical, canonical_gene_ids, "clustered")
        result.canonical_rows, result.canonical_cols = heatmap(
            canonical,
            display_label_map,
            canonical_gene_ids,
            canonical_tissue_names,
            f"H2A Expression (normalized cells) - {species} (clustered H2A)",
            can_png,
            can_svg,
            id_label,
            square_cells,
            cell_size,
            min_width,
            min_height,
        )
        result.canonical_png = str(can_png)
        result.canonical_svg = str(can_svg)

    if not variants.empty:
        var_png = species_out_dir / f"h2a_{slug}_variants.png"
        var_svg = species_out_dir / f"h2a_{slug}_variants.svg"
        variant_tissue_names = panel_tissues_and_report(variants, variant_gene_ids, "variants")
        result.variant_rows, result.variant_cols = heatmap(
            variants,
            display_label_map,
            variant_gene_ids,
            variant_tissue_names,
            f"H2A Expression (normalized cells) - {species} (variants)",
            var_png,
            var_svg,
            id_label,
            square_cells,
            cell_size,
            min_width,
            min_height,
        )
        result.variant_png = str(var_png)
        result.variant_svg = str(var_svg)

    if canonical.empty or variants.empty:
        missing = "clustered" if canonical.empty else "variants"
        result.status = "partial"
        result.reason = f"missing_{missing}_split"

    return result


def print_result(result: BuildResult) -> None:
    print(f"STATUS={result.status}")
    print(f"REASON={result.reason}")
    print(f"SPECIES={result.species}")
    print(f"SLUG={result.slug}")
    print(f"ID_COL={result.id_col}")
    print(f"CANONICAL_RULE={result.canonical_rule}")
    print(f"ALLOW_PARTIAL_SPLITS={result.allow_partial_splits}")
    print(f"ROWS_TOTAL_SCANNED={result.rows_total_scanned}")
    print(f"ROWS_AFTER_FILTER={result.rows_after_filter}")
    print(f"MATCHED_GENES_TOTAL={result.matched_genes_total}")
    print(f"MATCHED_GENES_CLUSTERED={result.matched_genes_clustered}")
    print(f"MATCHED_GENES_VARIANT={result.matched_genes_variant}")
    print(f"ALL_ROWS={result.all_rows} ALL_COLS={result.all_cols}")
    print(f"CANONICAL_ROWS={result.canonical_rows} CANONICAL_COLS={result.canonical_cols}")
    print(f"VARIANT_ROWS={result.variant_rows} VARIANT_COLS={result.variant_cols}")
    print(f"PROCESSED_TSV={result.processed_tsv}")
    print(f"MAP_TSV={result.map_tsv}")
    print(f"ALL_PNG={result.all_png}")
    print(f"ALL_SVG={result.all_svg}")
    print(f"CANONICAL_PNG={result.canonical_png}")
    print(f"CANONICAL_SVG={result.canonical_svg}")
    print(f"VARIANT_PNG={result.variant_png}")
    print(f"VARIANT_SVG={result.variant_svg}")


def main() -> None:
    args = parse_args()
    result = build_species_heatmaps(
        species=args.species,
        expr_path=Path(args.expr),
        merged_path=Path(args.merged),
        out_dir=Path(args.out_dir),
        out_processed_dir=Path(args.out_processed_dir),
        id_col=args.id_col,
        canonical_rule=args.canonical_rule,
        allow_partial_splits=args.allow_partial_splits,
        chunksize=args.chunksize,
        square_cells=args.square_cells,
        cell_size=args.cell_size,
        min_width=args.min_width,
        min_height=args.min_height,
        min_tissue_fill_rate=args.min_tissue_fill_rate,
    )
    print_result(result)


if __name__ == "__main__":
    main()
