
#!/usr/bin/env python
"""Build single-figure cross-species gene:anatomical_entity barplots with shared legends."""

from __future__ import annotations

import argparse
import itertools
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
import numpy as np
import pandas as pd

from gene_compare_common import (
    CROSS_SPECIES_CLASS_COLORS,
    DEFAULT_DETAIL_INDEX,
    DEFAULT_GENE_COMPARE_DATA_ROOT,
    DEFAULT_GENE_TISSUE_BARPLOT_DIR,
    DEFAULT_GENE_TISSUE_BARPLOT_TABLES_DIR,
    DEFAULT_HEATMAP_DIR,
    DEFAULT_PROCESSED_DIR,
    GENERIC_TISSUES,
    build_gene_long_dataframe,
    collect_gene_rows,
    load_detail_index,
    summarize_species_gene_tissue,
)
from normalized_expression_common import sort_gene_labels


DEFAULT_SPECIES_COUNTS = [4, 5, 6, 7]
DEFAULT_TABLE_STEM = "cross_species_gene_tissue_barplot"
SUMMARY_FILENAME = f"{DEFAULT_TABLE_STEM}_coverage_summary.tsv"
CLASS_ORDER = ["clustered", "variant"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build cross-species gene:anatomical_entity barplots from ranking-compatible "
            "species-level mean scores using fixed comparable species subsets."
        )
    )
    parser.add_argument("--detail-index", default=str(DEFAULT_DETAIL_INDEX))
    parser.add_argument("--heatmap-dir", default=str(DEFAULT_HEATMAP_DIR))
    parser.add_argument("--processed-dir", default=str(DEFAULT_PROCESSED_DIR))
    parser.add_argument("--gene-compare-data-root", default=str(DEFAULT_GENE_COMPARE_DATA_ROOT))
    parser.add_argument("--out-dir", default=str(DEFAULT_GENE_TISSUE_BARPLOT_DIR))
    parser.add_argument("--tables-dir", default=str(DEFAULT_GENE_TISSUE_BARPLOT_TABLES_DIR))
    parser.add_argument(
        "--species-counts",
        nargs="*",
        type=int,
        default=DEFAULT_SPECIES_COUNTS,
        help="Species-set sizes to optimize and plot.",
    )
    return parser.parse_args()


def normalize_species_counts(values: Sequence[int]) -> List[int]:
    cleaned: List[int] = []
    for value in values:
        value = int(value)
        if value < 2:
            raise ValueError("--species-counts values must be >= 2")
        cleaned.append(value)
    return sorted(set(cleaned))


def prefer_species_order(species_dirs: Iterable[str]) -> List[str]:
    unique = sorted({str(value).strip() for value in species_dirs if str(value).strip()})
    if "human" in unique:
        return ["human"] + [value for value in unique if value != "human"]
    return unique


def first_nonempty(values: Iterable[object]) -> str:
    for value in values:
        text = str(value).strip()
        if text:
            return text
    return ""


def axis_code(index: int) -> str:
    if index < 0:
        raise ValueError("Axis index must be non-negative.")
    label = ""
    value = index
    while True:
        value, remainder = divmod(value, 26)
        label = chr(ord("A") + remainder) + label
        if value == 0:
            return label
        value -= 1


def clean_text_columns(df: pd.DataFrame, columns: Sequence[str]) -> pd.DataFrame:
    for col in columns:
        if col not in df.columns:
            df[col] = ""
        df[col] = df[col].fillna("").astype(str).str.strip()
    return df


def build_presence_summary(
    gene_compare_data_root: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, Dict[str, str]]:
    frames: List[pd.DataFrame] = []
    for long_path in sorted(gene_compare_data_root.glob("*/*_gene_compare_long.csv")):
        long_df = pd.read_csv(
            long_path,
            usecols=["gene_name", "species_dir", "species_name", "tissue", "class"],
        )
        long_df["tissue"] = long_df["tissue"].fillna("").astype(str).str.strip()
        long_df = long_df[long_df["tissue"].ne("") & ~long_df["tissue"].isin(GENERIC_TISSUES)].copy()
        if long_df.empty:
            continue
        frames.append(
            long_df[["gene_name", "species_dir", "species_name", "tissue", "class"]].drop_duplicates()
        )

    if not frames:
        raise RuntimeError(f"No gene_compare long tables found under {gene_compare_data_root}")

    all_df = pd.concat(frames, ignore_index=True)
    all_df["gene_name"] = all_df["gene_name"].fillna("").astype(str).str.strip()
    all_df["species_dir"] = all_df["species_dir"].fillna("").astype(str).str.strip()
    all_df["species_name"] = all_df["species_name"].fillna("").astype(str).str.strip()
    all_df["class"] = all_df["class"].fillna("").astype(str).str.strip()
    all_df = all_df[
        all_df["gene_name"].ne("")
        & all_df["species_dir"].ne("")
        & all_df["tissue"].ne("")
    ].copy()

    species_name_map = (
        all_df[["species_dir", "species_name"]]
        .sort_values(["species_dir", "species_name"], ascending=[True, True], kind="stable")
        .drop_duplicates(subset=["species_dir"], keep="first")
        .set_index("species_dir")["species_name"]
        .to_dict()
    )

    combo_rows: List[dict] = []
    for (gene_name, tissue), group in all_df.groupby(["gene_name", "tissue"], sort=False):
        species_dirs = sorted(group["species_dir"].dropna().astype(str).unique().tolist())
        combo_rows.append(
            {
                "gene_name": gene_name,
                "tissue": tissue,
                "gene_tissue_label": f"{gene_name}:{tissue}",
                "gene_class": first_nonempty(group["class"].tolist()),
                "species_dirs": tuple(species_dirs),
                "species_n": int(len(species_dirs)),
            }
        )

    combo_df = pd.DataFrame(combo_rows)
    if combo_df.empty:
        raise RuntimeError("No comparable gene:tissue combinations were found.")

    species_rows = []
    for species_dir in sorted(species_name_map):
        combo_count = int(combo_df["species_dirs"].apply(lambda values: species_dir in values).sum())
        species_rows.append(
            {
                "species_dir": species_dir,
                "species_name": species_name_map.get(species_dir, species_dir),
                "combo_count": combo_count,
            }
        )
    species_df = pd.DataFrame(species_rows).sort_values(
        by=["combo_count", "species_dir"],
        ascending=[False, True],
        kind="stable",
    ).reset_index(drop=True)
    return combo_df, species_df, species_name_map


def choose_optimal_species_subsets(
    combo_df: pd.DataFrame,
    species_df: pd.DataFrame,
    species_name_map: Dict[str, str],
    species_counts: Sequence[int],
) -> pd.DataFrame:
    all_species = sorted(species_df["species_dir"].astype(str).tolist())
    combo_species_sets = [set(values) for values in combo_df["species_dirs"].tolist()]
    freq_map = dict(zip(species_df["species_dir"].astype(str), species_df["combo_count"].astype(int)))

    rows: List[dict] = []
    for species_count in species_counts:
        if species_count > len(all_species):
            raise RuntimeError(
                f"Requested species_count={species_count}, but only {len(all_species)} comparable species exist."
            )

        best_rows: List[tuple[tuple[str, ...], int]] = []
        best_count = -1
        for subset in itertools.combinations(all_species, species_count):
            subset_set = set(subset)
            retained_count = sum(1 for combo_set in combo_species_sets if subset_set.issubset(combo_set))
            if retained_count > best_count:
                best_count = retained_count
                best_rows = [(subset, retained_count)]
            elif retained_count == best_count:
                best_rows.append((subset, retained_count))

        def sort_key(item: tuple[tuple[str, ...], int]) -> tuple[int, tuple[str, ...]]:
            subset = item[0]
            return (-int(sum(freq_map[species] for species in subset)), tuple(subset))

        best_rows = sorted(best_rows, key=sort_key)
        chosen_subset = best_rows[0][0]
        ordered_subset = prefer_species_order(chosen_subset)
        tied_subsets = [
            ",".join(prefer_species_order(candidate_subset))
            for candidate_subset, _ in best_rows
        ]

        rows.append(
            {
                "species_count": int(species_count),
                "selected_species_dirs": ",".join(ordered_subset),
                "selected_species_names": ",".join(
                    species_name_map.get(species, species) for species in ordered_subset
                ),
                "retained_combo_count": int(best_count),
                "tied_optimal_subset_count": int(len(best_rows)),
                "tied_optimal_subsets": " | ".join(tied_subsets),
            }
        )

    return pd.DataFrame(rows).sort_values(by=["species_count"], ascending=[True]).reset_index(drop=True)


def build_coverage_summary_table(
    combo_df: pd.DataFrame,
    species_df: pd.DataFrame,
    subset_df: pd.DataFrame,
    out_dir: Path,
) -> pd.DataFrame:
    rows: List[dict] = []
    total_combos = int(len(combo_df))
    total_species = int(len(species_df))
    for row in species_df.to_dict(orient="records"):
        rows.append(
            {
                "summary_type": "species_frequency",
                "total_combos": total_combos,
                "total_species": total_species,
                "species_count": "",
                "species_dir": row["species_dir"],
                "species_name": row["species_name"],
                "combo_count": int(row["combo_count"]),
                "selected_species_dirs": "",
                "selected_species_names": "",
                "retained_combo_count": "",
                "tied_optimal_subset_count": "",
                "tied_optimal_subsets": "",
                "output_dir": "",
            }
        )
    for row in subset_df.to_dict(orient="records"):
        species_count = int(row["species_count"])
        rows.append(
            {
                "summary_type": "optimal_subset",
                "total_combos": total_combos,
                "total_species": total_species,
                "species_count": species_count,
                "species_dir": "",
                "species_name": "",
                "combo_count": "",
                "selected_species_dirs": row["selected_species_dirs"],
                "selected_species_names": row["selected_species_names"],
                "retained_combo_count": int(row["retained_combo_count"]),
                "tied_optimal_subset_count": int(row["tied_optimal_subset_count"]),
                "tied_optimal_subsets": row["tied_optimal_subsets"],
                "output_dir": (out_dir / f"k{species_count}").as_posix(),
            }
        )
    return pd.DataFrame(rows)


def compute_barplot_stats(values: Sequence[float]) -> Dict[str, float]:
    arr = np.sort(np.asarray(values, dtype=float))
    if arr.size == 0:
        raise RuntimeError("Cannot compute barplot statistics for an empty value set.")

    mean_score = float(arr.mean())
    std_score = float(np.std(arr, ddof=1)) if arr.size > 1 else 0.0

    return {
        "n_species": int(arr.size),
        "mean_score": mean_score,
        "std_score": std_score,
        "min_score": float(arr.min()),
        "max_score": float(arr.max()),
        "range_score": float(arr.max() - arr.min()),
        "error_bar_low": max(0.0, mean_score - std_score),
        "error_bar_high": mean_score + std_score,
    }


def build_ranking_compatible_species_level_cache(
    detail_df: pd.DataFrame,
    gene_names: Sequence[str],
) -> Dict[str, pd.DataFrame]:
    expr_cache: Dict[str, pd.DataFrame] = {}
    cache: Dict[str, pd.DataFrame] = {}
    for gene_name in sorted(set(gene_names)):
        gene_rows, _ = collect_gene_rows(gene_name, detail_df)
        gene_long_df = build_gene_long_dataframe(
            gene_name,
            gene_rows,
            expr_cache=expr_cache,
            aggregate="mean",
        )
        species_level_df = summarize_species_gene_tissue(gene_long_df)
        if species_level_df.empty:
            continue
        species_level_df = species_level_df.copy()
        species_level_df["tissue"] = species_level_df["tissue"].fillna("").astype(str).str.strip()
        species_level_df = species_level_df[
            species_level_df["tissue"].ne("")
            & ~species_level_df["tissue"].isin(GENERIC_TISSUES)
        ].copy()
        cache[gene_name] = species_level_df
    return cache


def build_gene_meta(detail_df: pd.DataFrame, gene_names: Sequence[str]) -> pd.DataFrame:
    meta_df = detail_df[detail_df["canonical_gene_name"].isin(gene_names)].copy()
    meta_df = clean_text_columns(
        meta_df,
        ["canonical_gene_name", "class", "map_label", "species_dir", "species_name"],
    )
    rows: List[dict] = []
    for gene_name, group in meta_df.groupby("canonical_gene_name", sort=False):
        class_values = [value for value in group["class"].tolist() if value]
        if "clustered" in class_values:
            gene_class = "clustered"
        elif "variant" in class_values:
            gene_class = "variant"
        else:
            gene_class = class_values[0] if class_values else "variant"

        label_values = sorted({value for value in group["map_label"].tolist() if value})
        canonical_label = label_values[0] if label_values else gene_name
        rows.append(
            {
                "gene_name": gene_name,
                "gene_class": gene_class,
                "canonical_label": canonical_label,
            }
        )

    result = pd.DataFrame(rows)
    result = result.sort_values(
        by=["gene_class", "canonical_label", "gene_name"],
        ascending=[True, True, True],
        kind="stable",
    ).reset_index(drop=True)
    return result


def build_heatmap_gene_order(
    gene_meta_df: pd.DataFrame,
    present_gene_names: Sequence[str],
) -> List[str]:
    present_set = {str(value).strip() for value in present_gene_names if str(value).strip()}
    ordered: List[str] = []
    for gene_class in CLASS_ORDER:
        class_df = gene_meta_df[
            gene_meta_df["gene_class"].eq(gene_class) & gene_meta_df["gene_name"].isin(present_set)
        ].copy()
        if class_df.empty:
            continue

        label_map = dict(
            zip(class_df["gene_name"].astype(str), class_df["canonical_label"].astype(str))
        )
        ordered_labels = sort_gene_labels(label_map, class_df["gene_name"].astype(str).tolist())
        label_to_genes: Dict[str, List[str]] = {}
        for row in class_df.to_dict(orient="records"):
            label_to_genes.setdefault(str(row["canonical_label"]), []).append(str(row["gene_name"]))
        for label in ordered_labels:
            for gene_name in sorted(label_to_genes.get(label, [])):
                if gene_name not in ordered:
                    ordered.append(gene_name)

        leftovers = [
            gene_name
            for gene_name in class_df["gene_name"].astype(str).tolist()
            if gene_name not in ordered
        ]
        ordered.extend(sorted(leftovers))
    return ordered


def build_subset_tables(
    combo_df: pd.DataFrame,
    subset_row: pd.Series,
    species_level_cache: Dict[str, pd.DataFrame],
    species_name_map: Dict[str, str],
    gene_meta_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    selected_species_dirs = [
        value.strip()
        for value in str(subset_row["selected_species_dirs"]).split(",")
        if value.strip()
    ]
    selected_species_set = set(selected_species_dirs)
    species_order_map = {species_dir: idx + 1 for idx, species_dir in enumerate(selected_species_dirs)}

    retained_combo_df = combo_df[
        combo_df["species_dirs"].apply(lambda values: selected_species_set.issubset(set(values)))
    ].copy()
    if retained_combo_df.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    gene_order = build_heatmap_gene_order(
        gene_meta_df,
        retained_combo_df["gene_name"].drop_duplicates().astype(str).tolist(),
    )
    gene_number_map = {gene_name: idx + 1 for idx, gene_name in enumerate(gene_order)}

    anatomical_entity_names = sorted(retained_combo_df["tissue"].drop_duplicates().astype(str).tolist())
    anatomical_entity_letter_map = {
        anatomical_entity: axis_code(idx) for idx, anatomical_entity in enumerate(anatomical_entity_names)
    }

    gene_meta_subset_df = gene_meta_df[gene_meta_df["gene_name"].isin(gene_order)].copy()
    gene_meta_subset_df["gene_number"] = gene_meta_subset_df["gene_name"].map(gene_number_map)
    gene_meta_subset_df = gene_meta_subset_df.sort_values("gene_number", kind="stable").reset_index(drop=True)

    gene_map_df = gene_meta_subset_df[
        ["gene_number", "gene_name", "gene_class", "canonical_label"]
    ].copy()
    gene_map_df["axis"] = "OY"
    gene_map_df = gene_map_df[
        ["axis", "gene_number", "gene_name", "gene_class", "canonical_label"]
    ].copy()

    tissue_map_df = pd.DataFrame(
        {
            "axis": "OY",
            "anatomical_entity_letter": [
                anatomical_entity_letter_map[anatomical_entity]
                for anatomical_entity in anatomical_entity_names
            ],
            "anatomical_entity_name": anatomical_entity_names,
        }
    )

    gene_class_map = dict(zip(gene_meta_subset_df["gene_name"], gene_meta_subset_df["gene_class"]))
    canonical_label_map = dict(zip(gene_meta_subset_df["gene_name"], gene_meta_subset_df["canonical_label"]))

    source_rows: List[dict] = []
    stats_rows: List[dict] = []
    for combo in retained_combo_df.to_dict(orient="records"):
        gene_name = str(combo["gene_name"])
        tissue = str(combo["tissue"])
        species_level_df = species_level_cache.get(gene_name)
        if species_level_df is None or species_level_df.empty:
            continue

        subset_df = species_level_df[
            species_level_df["tissue"].eq(tissue)
            & species_level_df["species_dir"].isin(selected_species_dirs)
        ].copy()
        subset_df = subset_df.sort_values(
            by=["species_dir", "species_name"],
            ascending=[True, True],
            kind="stable",
        ).drop_duplicates(subset=["species_dir"], keep="first")
        if int(subset_df["species_dir"].nunique()) != len(selected_species_dirs):
            continue

        subset_df["gene_name"] = gene_name
        subset_df["anatomical_entity"] = tissue
        subset_df["gene_tissue_label"] = combo["gene_tissue_label"]
        subset_df["gene_class"] = gene_class_map.get(gene_name, str(combo["gene_class"]).strip() or "variant")
        subset_df["canonical_label"] = canonical_label_map.get(gene_name, gene_name)
        subset_df["species_count"] = int(subset_row["species_count"])
        subset_df["selected_species_dirs"] = ",".join(selected_species_dirs)
        subset_df["selected_species_names"] = ",".join(
            species_name_map.get(species_dir, species_dir) for species_dir in selected_species_dirs
        )
        subset_df["species_order"] = subset_df["species_dir"].map(species_order_map)
        subset_df["gene_number"] = int(gene_number_map[gene_name])
        subset_df["anatomical_entity_letter"] = anatomical_entity_letter_map[tissue]
        subset_df["plot_label"] = (
            f"{gene_number_map[gene_name]}:{anatomical_entity_letter_map[tissue]}"
        )
        source_rows.extend(
            subset_df[
                [
                    "species_count",
                    "selected_species_dirs",
                    "selected_species_names",
                    "gene_name",
                    "anatomical_entity",
                    "gene_tissue_label",
                    "gene_class",
                    "canonical_label",
                    "gene_number",
                    "anatomical_entity_letter",
                    "plot_label",
                    "species_order",
                    "species_dir",
                    "species_name",
                    "cell_mean_score",
                    "cell_std_score",
                    "cell_n",
                    "cell_status",
                    "expression_score",
                ]
            ].to_dict(orient="records")
        )

        values = subset_df.sort_values("species_order", kind="stable")["expression_score"].astype(float).tolist()
        stats_rows.append(
            {
                "species_count": int(subset_row["species_count"]),
                "selected_species_dirs": ",".join(selected_species_dirs),
                "selected_species_names": ",".join(
                    species_name_map.get(species_dir, species_dir) for species_dir in selected_species_dirs
                ),
                "gene_name": gene_name,
                "anatomical_entity": tissue,
                "gene_tissue_label": combo["gene_tissue_label"],
                "gene_class": subset_df["gene_class"].iloc[0],
                "canonical_label": subset_df["canonical_label"].iloc[0],
                "gene_number": int(gene_number_map[gene_name]),
                "anatomical_entity_letter": anatomical_entity_letter_map[tissue],
                "plot_label": (
                    f"{gene_number_map[gene_name]}:{anatomical_entity_letter_map[tissue]}"
                ),
                "selected_values_csv": ",".join(f"{value:.6f}" for value in values),
                **compute_barplot_stats(values),
            }
        )

    source_df = pd.DataFrame(source_rows)
    stats_df = pd.DataFrame(stats_rows)
    if source_df.empty or stats_df.empty:
        return source_df, stats_df, gene_map_df, tissue_map_df

    class_rank_map = {gene_class: idx for idx, gene_class in enumerate(CLASS_ORDER)}
    stats_df["class_rank"] = stats_df["gene_class"].map(class_rank_map).fillna(99).astype(int)
    stats_df["gene_rank"] = stats_df["gene_name"].map(gene_number_map).fillna(10**9).astype(int)
    stats_df = stats_df.sort_values(
        by=["class_rank", "gene_rank", "anatomical_entity", "gene_name"],
        ascending=[True, True, True, True],
        kind="stable",
    ).reset_index(drop=True)
    stats_df["plot_order"] = np.arange(1, len(stats_df) + 1, dtype=int)

    plot_order_map = dict(zip(stats_df["plot_label"], stats_df["plot_order"]))
    source_df["plot_order"] = source_df["plot_label"].map(plot_order_map).fillna(10**9).astype(int)
    source_df = source_df.sort_values(
        by=["plot_order", "species_order", "gene_name", "anatomical_entity", "species_name"],
        ascending=[True, True, True, True, True],
        kind="stable",
    ).reset_index(drop=True)

    stats_df = stats_df.drop(columns=["class_rank", "gene_rank"])
    return source_df, stats_df, gene_map_df, tissue_map_df


def plot_single_figure(
    source_df: pd.DataFrame,
    stats_df: pd.DataFrame,
    out_png: Path,
    out_svg: Path,
    *,
    species_count: int,
    class_name: str | None = None,
) -> List[Path]:
    if source_df.empty or stats_df.empty:
        return []

    labels = stats_df["plot_label"].astype(str).tolist()
    means = stats_df["mean_score"].astype(float).to_numpy()
    stds = stats_df["std_score"].astype(float).to_numpy()
    positions = np.arange(len(labels), dtype=float)
    max_score = float(
        max(
            100.0,
            np.nanmax(means + stds) + 5.0 if len(means) else 100.0,
        )
    )

    fig_h = max(9.0, len(labels) * 0.24 + 2.0)
    fig, ax = plt.subplots(figsize=(15.0, fig_h))
    bar_colors = [
        CROSS_SPECIES_CLASS_COLORS.get(str(row["gene_class"]), "#b7c0c8")
        for row in stats_df.to_dict(orient="records")
    ]
    ax.barh(
        positions,
        means,
        xerr=stds,
        color=bar_colors,
        edgecolor="#444444",
        linewidth=1.0,
        alpha=0.95,
        height=0.62,
        error_kw={"ecolor": "#444444", "elinewidth": 1.1, "capsize": 3.5, "capthick": 1.1},
    )

    selected_species_names = stats_df["selected_species_names"].iloc[0]
    ax.invert_yaxis()
    ax.set_yticks(positions)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Expression score", fontsize=13)
    ax.set_ylabel("gene_number:anatomical_entity_letter", fontsize=13)
    ax.set_xlim(0.0, max_score)
    class_title = f" | class={class_name}" if class_name else ""
    ax.set_title(
        f"Cross-species gene:anatomical_entity barplots | k={species_count}{class_title}\n"
        f"Species: {selected_species_names}",
        fontsize=15,
        fontweight="bold",
        pad=12,
    )
    ax.grid(axis="x", alpha=0.25)
    ax.grid(axis="y", alpha=0.08)
    ax.set_axisbelow(True)
    ax.tick_params(axis="x", labelsize=10)
    ax.tick_params(axis="y", labelsize=8)

    present_classes = [
        gene_class
        for gene_class in CLASS_ORDER
        if stats_df["gene_class"].astype(str).eq(gene_class).any()
    ]
    legend_handles = [
        Patch(
            facecolor=CROSS_SPECIES_CLASS_COLORS[gene_class],
            edgecolor="#444444",
            label=gene_class,
        )
        for gene_class in present_classes
    ]
    if legend_handles:
        ax.legend(handles=legend_handles, title="Gene class", loc="upper right", frameon=True)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=(0.03, 0.02, 0.99, 0.98))
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    plt.close(fig)
    return [out_png, out_svg]


def plot_gene_legend(gene_map_df: pd.DataFrame, out_png: Path, out_svg: Path) -> List[Path]:
    if gene_map_df.empty:
        return []

    fig_h = max(7.0, len(gene_map_df) * 0.50 + 1.2)
    fig, ax = plt.subplots(figsize=(9.4, fig_h))
    ax.axis("off")

    ax.text(
        0.05,
        0.985,
        "OY: Gene index",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=20,
        fontweight="bold",
        color="#111111",
    )

    top = 0.91
    bottom = 0.05
    step = (top - bottom) / max(len(gene_map_df) - 1, 1)
    x_number = 0.14
    x_swatch = 0.19
    swatch_size = 0.03
    x_gene = 0.255

    for idx, row in enumerate(gene_map_df.to_dict(orient="records")):
        y = top - idx * step
        color = CROSS_SPECIES_CLASS_COLORS.get(str(row["gene_class"]), "#333333")
        ax.text(
            x_number,
            y,
            str(int(row["gene_number"])),
            transform=ax.transAxes,
            va="center",
            ha="right",
            fontsize=18,
            fontweight="bold",
            color="#111111",
        )
        ax.add_patch(
            Rectangle(
                (x_swatch, y - swatch_size / 2.0),
                swatch_size,
                swatch_size,
                transform=ax.transAxes,
                facecolor=color,
                edgecolor="#444444",
                linewidth=0.8,
            )
        )
        ax.text(
            x_gene,
            y,
            str(row["gene_name"]),
            transform=ax.transAxes,
            va="center",
            ha="left",
            fontsize=18,
            color="#111111",
        )

    fig.tight_layout(rect=(0.02, 0.01, 0.995, 0.995))
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    plt.close(fig)
    return [out_png, out_svg]


def plot_anatomical_entity_legend(
    anatomical_entity_map_df: pd.DataFrame,
    out_png: Path,
    out_svg: Path,
) -> List[Path]:
    if anatomical_entity_map_df.empty:
        return []

    fig_h = max(8.0, len(anatomical_entity_map_df) * 0.50 + 1.2)
    fig, ax = plt.subplots(figsize=(10.8, fig_h))
    ax.axis("off")

    ax.text(
        0.05,
        0.985,
        "OY: Anatomical entity",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=24,
        fontweight="bold",
        color="#111111",
    )

    top = 0.92
    bottom = 0.04
    step = (top - bottom) / max(len(anatomical_entity_map_df) - 1, 1)
    x_letter = 0.13
    x_name = 0.21

    for idx, row in enumerate(anatomical_entity_map_df.to_dict(orient="records")):
        y = top - idx * step
        ax.text(
            x_letter,
            y,
            str(row["anatomical_entity_letter"]),
            transform=ax.transAxes,
            va="center",
            ha="right",
            fontsize=18,
            fontweight="bold",
            color="#111111",
        )
        ax.text(
            x_name,
            y,
            str(row["anatomical_entity_name"]),
            transform=ax.transAxes,
            va="center",
            ha="left",
            fontsize=18,
            color="#111111",
        )

    fig.tight_layout(rect=(0.02, 0.01, 0.995, 0.995))
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    plt.close(fig)
    return [out_png, out_svg]


def filter_class_outputs(
    source_df: pd.DataFrame,
    stats_df: pd.DataFrame,
    gene_map_df: pd.DataFrame,
    anatomical_entity_map_df: pd.DataFrame,
    class_name: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    class_source_df = source_df[source_df["gene_class"].eq(class_name)].copy()
    class_stats_df = stats_df[stats_df["gene_class"].eq(class_name)].copy()
    if class_source_df.empty or class_stats_df.empty:
        return class_source_df, class_stats_df, pd.DataFrame(), pd.DataFrame()

    class_gene_numbers = class_stats_df["gene_number"].drop_duplicates().astype(int).tolist()
    class_letters = (
        class_stats_df["anatomical_entity_letter"].drop_duplicates().astype(str).tolist()
    )
    class_gene_map_df = gene_map_df[gene_map_df["gene_number"].isin(class_gene_numbers)].copy()
    class_anatomical_entity_map_df = anatomical_entity_map_df[
        anatomical_entity_map_df["anatomical_entity_letter"].isin(class_letters)
    ].copy()
    return class_source_df, class_stats_df, class_gene_map_df, class_anatomical_entity_map_df


def write_output_set(
    output_dir: Path,
    *,
    plot_stem: str,
    table_stem: str,
    source_df: pd.DataFrame,
    stats_df: pd.DataFrame,
    gene_map_df: pd.DataFrame,
    anatomical_entity_map_df: pd.DataFrame,
    species_count: int,
    class_name: str | None = None,
) -> List[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    source_path = output_dir / f"{table_stem}_source.tsv"
    stats_path = output_dir / f"{table_stem}_stats.tsv"
    gene_map_path = output_dir / f"{table_stem}_gene_number_map.tsv"
    anatomical_entity_map_path = output_dir / f"{table_stem}_anatomical_entity_letter_map.tsv"
    source_df.to_csv(source_path, sep="\t", index=False)
    stats_df.to_csv(stats_path, sep="\t", index=False)
    gene_map_df.to_csv(gene_map_path, sep="\t", index=False)
    anatomical_entity_map_df.to_csv(anatomical_entity_map_path, sep="\t", index=False)

    written_paths: List[Path] = [source_path, stats_path, gene_map_path, anatomical_entity_map_path]
    written_paths.extend(
        plot_single_figure(
            source_df,
            stats_df,
            output_dir / f"{plot_stem}.png",
            output_dir / f"{plot_stem}.svg",
            species_count=species_count,
            class_name=class_name,
        )
    )
    written_paths.extend(
        plot_gene_legend(
            gene_map_df,
            output_dir / f"{plot_stem}_gene_legend.png",
            output_dir / f"{plot_stem}_gene_legend.svg",
        )
    )
    written_paths.extend(
        plot_anatomical_entity_legend(
            anatomical_entity_map_df,
            output_dir / f"{plot_stem}_anatomical_entity_legend.png",
            output_dir / f"{plot_stem}_anatomical_entity_legend.svg",
        )
    )
    return written_paths


def run(args: argparse.Namespace) -> Dict[str, object]:
    species_counts = normalize_species_counts(args.species_counts)
    out_dir = Path(args.out_dir)
    tables_dir = Path(args.tables_dir)
    tables_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    combo_df, species_df, species_name_map = build_presence_summary(Path(args.gene_compare_data_root))
    subset_df = choose_optimal_species_subsets(combo_df, species_df, species_name_map, species_counts)
    coverage_summary_df = build_coverage_summary_table(combo_df, species_df, subset_df, out_dir)
    coverage_summary_path = tables_dir / SUMMARY_FILENAME
    coverage_summary_df.to_csv(coverage_summary_path, sep="\t", index=False)

    retained_gene_names = sorted(
        {
            row["gene_name"]
            for row in combo_df[
                combo_df["species_dirs"].apply(
                    lambda values: any(
                        set(
                            species_dir.strip()
                            for species_dir in selected_species_dirs.split(",")
                            if species_dir.strip()
                        ).issubset(set(values))
                        for selected_species_dirs in subset_df["selected_species_dirs"].astype(str)
                    )
                )
            ].to_dict(orient="records")
        }
    )
    detail_df = load_detail_index(
        Path(args.detail_index),
        heatmap_dir=Path(args.heatmap_dir),
        processed_dir=Path(args.processed_dir),
    )
    gene_meta_df = build_gene_meta(detail_df, retained_gene_names)
    species_level_cache = build_ranking_compatible_species_level_cache(detail_df, retained_gene_names)

    outputs: Dict[str, object] = {
        "coverage_summary_path": coverage_summary_path,
        "subset_counts": {},
        "written_paths": [],
    }
    for subset_row in subset_df.to_dict(orient="records"):
        species_count = int(subset_row["species_count"])
        k_dir = out_dir / f"k{species_count}"
        k_dir.mkdir(parents=True, exist_ok=True)

        source_df, stats_df, gene_map_df, anatomical_entity_map_df = build_subset_tables(
            combo_df,
            pd.Series(subset_row),
            species_level_cache,
            species_name_map,
            gene_meta_df,
        )
        if source_df.empty or stats_df.empty:
            raise RuntimeError(f"No source/stat rows were produced for k={species_count}")

        written_paths = write_output_set(
            k_dir,
            plot_stem=f"{DEFAULT_TABLE_STEM}_k{species_count}",
            table_stem=f"{DEFAULT_TABLE_STEM}_k{species_count}",
            source_df=source_df,
            stats_df=stats_df,
            gene_map_df=gene_map_df,
            anatomical_entity_map_df=anatomical_entity_map_df,
            species_count=species_count,
        )

        class_output_dirs: Dict[str, Path] = {}
        class_counts: Dict[str, int] = {}
        for class_name in CLASS_ORDER:
            (
                class_source_df,
                class_stats_df,
                class_gene_map_df,
                class_anatomical_entity_map_df,
            ) = filter_class_outputs(
                source_df,
                stats_df,
                gene_map_df,
                anatomical_entity_map_df,
                class_name,
            )
            if class_source_df.empty or class_stats_df.empty:
                continue

            class_dir = k_dir / "by_class" / class_name
            class_output_dirs[class_name] = class_dir
            class_counts[class_name] = int(len(class_stats_df))
            written_paths.extend(
                write_output_set(
                    class_dir,
                    plot_stem=f"{DEFAULT_TABLE_STEM}_{class_name}_k{species_count}",
                    table_stem=f"{DEFAULT_TABLE_STEM}_{class_name}_k{species_count}",
                    source_df=class_source_df,
                    stats_df=class_stats_df,
                    gene_map_df=class_gene_map_df,
                    anatomical_entity_map_df=class_anatomical_entity_map_df,
                    species_count=species_count,
                    class_name=class_name,
                )
            )

        outputs["subset_counts"][species_count] = {
            "source_rows": int(len(source_df)),
            "stats_rows": int(len(stats_df)),
            "selected_species_dirs": str(subset_row["selected_species_dirs"]),
            "output_dir": k_dir,
            "class_counts": class_counts,
            "class_output_dirs": {name: path.as_posix() for name, path in class_output_dirs.items()},
        }
        outputs["written_paths"].extend(written_paths)

    return outputs


def main() -> None:
    args = parse_args()
    outputs = run(args)
    print(f"Saved coverage summary to {outputs['coverage_summary_path']}")
    for species_count, payload in outputs["subset_counts"].items():
        print(
            f"k={species_count}: {payload['stats_rows']} barplots | "
            f"{payload['selected_species_dirs']} | {payload['output_dir']}"
        )
    for path in outputs["written_paths"]:
        print(f"Saved asset to {path}")


if __name__ == "__main__":
    main()
