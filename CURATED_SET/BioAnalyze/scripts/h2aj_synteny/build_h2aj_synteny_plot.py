#!/usr/bin/env python
"""Build a presentation-friendly H2A.J synteny plot across mammals."""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List, Sequence

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd


BIOANALYZE_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT_DIR = Path(
    r"C:\Users\USER\Documents\GitHub\histonedb_external_storage\BioAnalyze\synteny"
)
DEFAULT_OUTPUT_DIR = BIOANALYZE_ROOT / "figures" / "h2aj_synteny"
DEFAULT_OUTPUT_STEM = "h2aj_synteny"
DEFAULT_NEIGHBORS = 1500

FILE_ORDER = [
    "Homo-sapiens.tsv",
    "Pan-troglodytes.tsv",
    "Mus-musculus.tsv",
    "Oryctolagus-cuniculus.tsv",
    "Erinaceus-europaeus.tsv",
    "Eptesicus-fuscus.tsv",
    "Canis-lupus.tsv",
    "Sus-scrofa.tsv",
    "Equus-caballus.tsv",
    "Manis-pentadactyla.tsv",
    "Loxodonta-africana.tsv",
    "Choloepus-didactylus.tsv",
    "Dasypus-novemcinctus.tsv",
]

ORDER_MAP = {
    "Homo sapiens": "Primates",
    "Pan troglodytes": "Primates",
    "Mus musculus": "Rodentia",
    "Oryctolagus cuniculus": "Lagomorpha",
    "Erinaceus europaeus": "Eulipotyphla",
    "Eptesicus fuscus": "Chiroptera",
    "Canis lupus": "Carnivora",
    "Sus scrofa": "Artiodactyla",
    "Equus caballus": "Perissodactyla",
    "Manis pentadactyla": "Pholidota",
    "Loxodonta africana": "Proboscidea",
    "Choloepus didactylus": "Pilosa",
    "Dasypus novemcinctus": "Cingulata",
}

H2AJ_NAME_PATTERN = re.compile(r"\bh2a[\s\-_.]*j\b|\bj[\s\-_.]*h2a\b", re.IGNORECASE)
SYMBOL_H2AJ_PATTERN = re.compile(r"^h2aj$", re.IGNORECASE)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a presentation-oriented H2A.J synteny plot from external TSV files."
        )
    )
    parser.add_argument(
        "--input-dir",
        default=str(DEFAULT_INPUT_DIR),
        help="Directory with mammalian synteny TSV files.",
    )
    parser.add_argument(
        "--output-dir",
        default=str(DEFAULT_OUTPUT_DIR),
        help="Directory for output figures.",
    )
    parser.add_argument(
        "--output-stem",
        default=DEFAULT_OUTPUT_STEM,
        help="Base filename for saved figures.",
    )
    parser.add_argument(
        "--neighbors",
        type=int,
        default=DEFAULT_NEIGHBORS,
        help="Number of flanking genes to inspect on each side before filtering.",
    )
    parser.add_argument(
        "--min-shared-species",
        type=int,
        default=len(FILE_ORDER),
        help="Minimum species count for a neighbor gene to remain in the plot.",
    )
    return parser.parse_args()


def normalize_name(name: object) -> str:
    words = re.sub(r"[^a-z0-9]", " ", str(name).lower()).split()
    return " ".join(sorted(words))


def invert_orientation(orientation: str) -> str:
    return "minus" if orientation == "plus" else "plus"


def species_name_from_path(path: Path) -> str:
    return path.stem.replace("-", " ")


def require_columns(df: pd.DataFrame, path: Path, columns: Sequence[str]) -> None:
    missing = [column for column in columns if column not in df.columns]
    if missing:
        raise ValueError(f"{path.name} is missing required columns: {', '.join(missing)}")


def resolve_input_paths(input_dir: Path) -> List[Path]:
    paths: List[Path] = []
    missing: List[str] = []
    for file_name in FILE_ORDER:
        path = input_dir / file_name
        if not path.exists():
            missing.append(file_name)
        else:
            paths.append(path)
    if missing:
        raise FileNotFoundError(
            "Missing expected synteny TSV files: " + ", ".join(sorted(missing))
        )
    return paths


def find_h2aj_row(df: pd.DataFrame) -> pd.Series:
    name_clean = df["Name"].fillna("").astype(str).str.lower()
    name_clean = name_clean.str.replace(r"[^a-z0-9\s\-_.]", " ", regex=True)
    h2aj_mask = name_clean.str.contains(H2AJ_NAME_PATTERN, regex=True)

    if not h2aj_mask.any() and "Symbol" in df.columns:
        symbol_clean = df["Symbol"].fillna("").astype(str).str.strip().str.lower()
        h2aj_mask = symbol_clean.str.match(SYMBOL_H2AJ_PATTERN)

    if not h2aj_mask.any():
        raise RuntimeError("H2A.J was not found in Name or Symbol columns.")
    return df.loc[h2aj_mask].iloc[0]


def build_species_row(path: Path, neighbors: int) -> Dict[str, object]:
    df = pd.read_csv(path, sep="\t")
    require_columns(df, path, ["Name", "Orientation", "Chromosome"])

    coord_col = "Begin" if "Begin" in df.columns else "Start"
    end_col = "End" if "End" in df.columns else "Stop"
    require_columns(df, path, [coord_col, end_col])

    h2aj_row = find_h2aj_row(df)
    chromosome = h2aj_row["Chromosome"]
    strand = str(h2aj_row["Orientation"]).strip().lower()
    if strand not in {"plus", "minus"}:
        raise RuntimeError(f"Unsupported H2A.J orientation '{strand}' in {path.name}.")

    chromosome_df = df.loc[df["Chromosome"].eq(chromosome)].copy()
    chromosome_df = chromosome_df.sort_values(by=coord_col).reset_index(drop=True)

    h2aj_idx = chromosome_df.index[chromosome_df["Name"].eq(h2aj_row["Name"])]
    if h2aj_idx.empty:
        raise RuntimeError(f"H2A.J index not found after chromosome sort in {path.name}.")

    idx = int(h2aj_idx[0])
    left_df = chromosome_df.iloc[max(0, idx - neighbors) : idx]
    right_df = chromosome_df.iloc[idx + 1 : idx + 1 + neighbors]

    left_genes = [
        {
            "name": str(row["Name"]),
            "norm": normalize_name(row["Name"]),
            "orient": str(row["Orientation"]).strip().lower(),
        }
        for _, row in left_df.iterrows()
    ]
    right_genes = [
        {
            "name": str(row["Name"]),
            "norm": normalize_name(row["Name"]),
            "orient": str(row["Orientation"]).strip().lower(),
        }
        for _, row in right_df.iterrows()
    ]

    if strand == "minus":
        left_genes, right_genes = right_genes[::-1], left_genes[::-1]
        for gene in left_genes + right_genes:
            gene["orient"] = invert_orientation(str(gene["orient"]))

    return {
        "species": species_name_from_path(path),
        "left": left_genes,
        "right": right_genes,
        "strand": "→",
    }


def build_filtered_rows(
    input_paths: Sequence[Path],
    *,
    neighbors: int,
    min_shared_species: int,
) -> List[Dict[str, object]]:
    all_rows: List[Dict[str, object]] = []
    gene_occurrences: Dict[str, set[str]] = {}

    for path in input_paths:
        row = build_species_row(path, neighbors)
        all_rows.append(row)
        for gene in row["left"] + row["right"]:
            gene_occurrences.setdefault(str(gene["norm"]), set()).add(str(row["species"]))

    frequent_genes = {
        gene_name
        for gene_name, species_set in gene_occurrences.items()
        if len(species_set) >= min_shared_species
    }
    if not frequent_genes:
        raise RuntimeError(
            f"No shared neighbor genes passed min_shared_species={min_shared_species}."
        )

    filtered_rows: List[Dict[str, object]] = []
    for row in all_rows:
        left_filtered = [gene for gene in row["left"] if gene["norm"] in frequent_genes]
        right_filtered = [gene for gene in row["right"] if gene["norm"] in frequent_genes]
        if not left_filtered and not right_filtered:
            continue
        filtered_rows.append(
            {
                "species": row["species"],
                "order": ORDER_MAP.get(str(row["species"]), "Unknown"),
                "left": left_filtered,
                "right": right_filtered,
                "strand": row["strand"],
            }
        )

    filtered_rows.sort(
        key=lambda row: FILE_ORDER.index(f"{str(row['species']).replace(' ', '-')}.tsv")
    )
    if not filtered_rows:
        raise RuntimeError("All rows were filtered out; nothing remains to plot.")
    return filtered_rows


def build_gene_color_map(rows: Sequence[Dict[str, object]]) -> Dict[str, object]:
    all_genes = sorted(
        {
            str(gene["norm"])
            for row in rows
            for gene in list(row["left"]) + list(row["right"])
        }
    )
    cmap = plt.get_cmap("tab20", max(1, len(all_genes)))
    return {gene_name: cmap(index) for index, gene_name in enumerate(all_genes)}


def build_gene_label_map(rows: Sequence[Dict[str, object]]) -> Dict[str, str]:
    labels: Dict[str, str] = {}
    for row in rows:
        for gene in list(row["left"]) + list(row["right"]):
            norm = str(gene["norm"])
            labels.setdefault(norm, str(gene["name"]))
    return labels


def draw_arrow(
    ax: plt.Axes,
    *,
    x_center: float,
    y: float,
    orientation: str,
    color: object,
    arrow_dx: float,
    arrow_width: float,
    arrow_head_width: float,
    arrow_head_length: float,
) -> None:
    direction = arrow_dx if orientation == "plus" else -arrow_dx
    ax.add_patch(
        mpatches.FancyArrow(
            x_center - direction / 2,
            y,
            direction,
            0,
            width=arrow_width,
            head_width=arrow_head_width,
            head_length=arrow_head_length,
            length_includes_head=True,
            color=color,
        )
    )


def plot_rows(rows: Sequence[Dict[str, object]], output_png: Path, output_svg: Path) -> None:
    gene_color_map = build_gene_color_map(rows)
    gene_label_map = build_gene_label_map(rows)
    max_left = max(len(row["left"]) for row in rows)
    max_right = max(len(row["right"]) for row in rows)
    n_rows = len(rows)

    x_spacing = 2.25
    arrow_dx = 1.18
    arrow_width = 0.22
    arrow_head_width = 0.56
    arrow_head_length = 0.42
    label_x = -(max_left + 1.7) * x_spacing
    species_x = label_x + 0.75
    fig_w = max(33, (max_left + max_right) * 1.45 + 24.0)
    fig_h = max(15.5, n_rows * 1.28 + 2.3)

    plt.close("all")
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    for idx, row in enumerate(rows):
        y = n_rows - idx
        species = str(row["species"])
        order = str(row["order"])

        ax.text(
            species_x,
            y + 0.18,
            species,
            ha="right",
            va="center",
            fontsize=21.1,
            fontstyle="italic",
            color="#222222",
        )
        ax.text(
            species_x,
            y - 0.2,
            order,
            ha="right",
            va="center",
            fontsize=21.6,
            fontweight="bold",
            color="black",
        )

        draw_arrow(
            ax,
            x_center=0,
            y=y,
            orientation="plus",
            color="black",
            arrow_dx=arrow_dx,
            arrow_width=0.28,
            arrow_head_width=0.68,
            arrow_head_length=0.50,
        )

        for gene_index, gene in enumerate(reversed(row["left"])):
            x = -(gene_index + 1) * x_spacing
            draw_arrow(
                ax,
                x_center=x,
                y=y,
                orientation=str(gene["orient"]),
                color=gene_color_map[str(gene["norm"])],
                arrow_dx=arrow_dx,
                arrow_width=arrow_width,
                arrow_head_width=arrow_head_width,
                arrow_head_length=arrow_head_length,
            )

        for gene_index, gene in enumerate(row["right"]):
            x = (gene_index + 1) * x_spacing
            draw_arrow(
                ax,
                x_center=x,
                y=y,
                orientation=str(gene["orient"]),
                color=gene_color_map[str(gene["norm"])],
                arrow_dx=arrow_dx,
                arrow_width=arrow_width,
                arrow_head_width=arrow_head_width,
                arrow_head_length=arrow_head_length,
            )

    ax.set_xlim(label_x - 0.8, (max_right + 1.8) * x_spacing)
    ax.set_ylim(0.15, n_rows + 1.05)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    legend_handles = [
        mpatches.Patch(color="black", label="H2A.J"),
    ]
    legend_handles.extend(
        [
        mpatches.Patch(color=color, label=gene_label_map[gene_name])
        for gene_name, color in gene_color_map.items()
        ]
    )
    ax.legend(
        handles=legend_handles,
        title="Genes",
        loc="center left",
        bbox_to_anchor=(1.08, 0.5),
        frameon=False,
        fontsize=18,
        title_fontsize=22,
        borderaxespad=0.0,
        handlelength=1.7,
        handleheight=1.3,
        handletextpad=0.9,
        labelspacing=0.95,
    )

    ax.text(
        0,
        0.32,
        "H2A.J",
        ha="center",
        va="top",
        fontsize=21,
        fontweight="bold",
        color="black",
    )

    ax.set_title("H2A.J Gene Synteny Across Mammals", fontsize=18, fontweight="bold", pad=16)
    plt.tight_layout(rect=(0, 0.03, 0.70, 1))

    output_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.savefig(output_svg, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)

    input_paths = resolve_input_paths(input_dir)
    rows = build_filtered_rows(
        input_paths,
        neighbors=args.neighbors,
        min_shared_species=args.min_shared_species,
    )

    output_png = output_dir / f"{args.output_stem}.png"
    output_svg = output_dir / f"{args.output_stem}.svg"
    plot_rows(rows, output_png, output_svg)

    print(f"Built H2A.J synteny plot for {len(rows)} species.")
    print(f"PNG: {output_png}")
    print(f"SVG: {output_svg}")


if __name__ == "__main__":
    main()
