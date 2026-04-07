#!/usr/bin/env python
"""Build annotated codon entropy heatmaps for SQK H2A datasets."""

from __future__ import annotations

import argparse
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import SeqIO
from Bio.Data import CodonTable
from matplotlib import cm, colors

BIOANALYZE_SCRIPTS_ROOT = Path(__file__).resolve().parents[1]
if str(BIOANALYZE_SCRIPTS_ROOT) not in sys.path:
    sys.path.insert(0, str(BIOANALYZE_SCRIPTS_ROOT))

from bioanalyze_paths import get_bioanalyze_raw_root


BIOANALYZE_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT_DIR = get_bioanalyze_raw_root() / "codons"
DEFAULT_OUTPUT_DIR = BIOANALYZE_ROOT / "figures" / "codons"

STANDARD_TABLE = CodonTable.unambiguous_dna_by_id[1]
AA_MAX_SYNONYMOUS_CODONS = {
    # Standard genetic code (NCBI translation table 1), excluding stop codons.
    "A": 4,
    "C": 2,
    "D": 2,
    "E": 2,
    "F": 2,
    "G": 4,
    "H": 2,
    "I": 3,
    "K": 2,
    "L": 6,
    "M": 1,
    "N": 2,
    "P": 4,
    "Q": 2,
    "R": 6,
    "S": 6,
    "T": 4,
    "V": 4,
    "W": 1,
    "Y": 2,
}

DATASETS = {
    "without-short": {
        "protein": "protein_from_SQK_nuc(without short).fasta",
        "cds": "SQK_nuc(without short).fasta",
        "stem": "sqk_nuc_without_short_codon_entropy_annotated",
        "presentation_mode": True,
        "cell_mode": "reference-blank-nonsyn",
    },
    "full": {
        "protein": "protein_from_SQK_nuc.fasta",
        "cds": "SQK_nuc.fasta",
        "stem": "sqk_nuc_full_codon_entropy_annotated",
        "presentation_mode": True,
        "cell_mode": "majority-star-nonsyn",
    },
}

ORDER_COLORS = {
    "Artiodactyla": "#1f77b4",
    "Carnivora": "#ff7f0e",
    "Rodentia": "#2ca02c",
    "Primates": "#d62728",
    "Chiroptera": "#9467bd",
    "Eulipotyphla": "#8c564b",
    "Perissodactyla": "#e377c2",
    "Lagomorpha": "#7f7f7f",
}

DEFAULT_REGIONS = [
    (17, 23, "α1ext"),
    (27, 38, "α1"),
    (40, 45, "loopL1"),
    (47, 74, "α2"),
    (76, 78, "loopL2"),
    (80, 89, "α3"),
    (100, 102, "β3"),
]

FULL_PRESENTATION_REGIONS = [
    (17, 23, "α1e"),
    (27, 38, "α1"),
    (40, 45, "L1"),
    (47, 74, "α2"),
    (76, 78, "L2"),
    (80, 89, "α3"),
    (100, 102, "β3"),
]

HIGHLIGHTED_POSITIONS = [43, 57, 62, 65, 78, 91, 92, 93, 123, 124, 125]
BLUE_POSITIONS = {43, 78}
GREEN_POSITIONS = {123, 124, 125}


def build_synonymous_codon_sets() -> Dict[str, frozenset[str]]:
    synonymous_codons: Dict[str, set[str]] = defaultdict(set)
    for codon, amino_acid in STANDARD_TABLE.forward_table.items():
        synonymous_codons[amino_acid].add(codon)
    return {
        amino_acid: frozenset(codons)
        for amino_acid, codons in sorted(synonymous_codons.items())
    }


AA_SYNONYMOUS_CODONS = build_synonymous_codon_sets()


def validate_synonymous_codon_counts() -> None:
    derived_counts = {
        amino_acid: len(codons)
        for amino_acid, codons in AA_SYNONYMOUS_CODONS.items()
    }
    if derived_counts != AA_MAX_SYNONYMOUS_CODONS:
        raise ValueError(
            "AA_MAX_SYNONYMOUS_CODONS does not match the standard genetic code: "
            f"{derived_counts}"
        )


validate_synonymous_codon_counts()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build annotated codon-entropy heatmaps from paired protein/CDS FASTA "
            "inputs for the SQK H2A datasets."
        )
    )
    parser.add_argument(
        "--input-dir",
        default=str(DEFAULT_INPUT_DIR),
        help="Directory with the four raw SQK FASTA inputs.",
    )
    parser.add_argument(
        "--output-dir",
        default=str(DEFAULT_OUTPUT_DIR),
        help="Directory for PNG/SVG outputs.",
    )
    parser.add_argument(
        "--dataset",
        default="all",
        choices=["without-short", "full", "all"],
        help="Dataset variant to build.",
    )
    parser.add_argument(
        "--with-legend",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Whether to write the standalone full-dataset legend assets.",
    )
    return parser.parse_args()


def resolve_dataset_names(selected: str) -> List[str]:
    if selected == "all":
        return ["without-short", "full"]
    return [selected]


def require_file(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Missing required file: {path}")


def read_reference_sequences(
    protein_path: Path,
    cds_path: Path,
) -> Tuple[str, List[str], str]:
    prot_iter = SeqIO.parse(str(protein_path), "fasta")
    cds_iter = SeqIO.parse(str(cds_path), "fasta")

    first_prot = next(prot_iter)
    first_cds = next(cds_iter)

    aa_ref = str(first_prot.seq)
    cds_ref = str(first_cds.seq)
    if cds_ref[-3:] in {"TAA", "TAG", "TGA"}:
        cds_ref = cds_ref[:-3]

    required_nt = len(aa_ref) * 3
    if len(cds_ref) < required_nt:
        raise ValueError(
            f"Reference CDS is shorter than the reference protein: {cds_path}"
        )

    cds_ref = cds_ref[:required_nt]
    codons_ref = [cds_ref[i : i + 3] for i in range(0, required_nt, 3)]
    reference_species = first_prot.description.split("|")[0].strip()
    return aa_ref, codons_ref, reference_species


def shannon_entropy(codons: Sequence[str]) -> float:
    counts = Counter(codons)
    freqs = np.array(list(counts.values()), dtype=float)
    probabilities = freqs / freqs.sum()
    return float(-np.sum(probabilities * np.log2(probabilities)))


def basis_amino_acid_for_cell(
    codons: Sequence[str], reference_aa: str
) -> str | None:
    amino_acid_counts: Counter[str] = Counter()
    for codon in codons:
        if codon == "---":
            continue
        amino_acid = STANDARD_TABLE.forward_table.get(codon)
        if amino_acid is not None:
            amino_acid_counts[amino_acid] += 1

    if not amino_acid_counts:
        return None

    top_count = max(amino_acid_counts.values())
    leaders = sorted(
        amino_acid
        for amino_acid, count in amino_acid_counts.items()
        if count == top_count
    )
    # Prefer the reference amino acid on ties to keep the basis biologically stable.
    if reference_aa in leaders:
        return reference_aa
    return leaders[0]


def synonymous_codons_for_amino_acid(
    codons: Sequence[str], basis_aa: str
) -> Tuple[List[str], bool]:
    synonymous_codons: List[str] = []
    has_alternative_amino_acid = False
    for codon in codons:
        if codon == "---":
            continue
        amino_acid = STANDARD_TABLE.forward_table.get(codon)
        if amino_acid == basis_aa:
            synonymous_codons.append(codon)
        elif amino_acid is not None:
            has_alternative_amino_acid = True
    return synonymous_codons, has_alternative_amino_acid


def normalized_entropy_majority_basis(codons: Sequence[str], reference_aa: str) -> float:
    """Normalize synonymous codon diversity relative to the most frequent amino acid in the cell."""
    basis_aa = basis_amino_acid_for_cell(codons, reference_aa)
    if basis_aa is None:
        return float("nan")

    synonymous_codons, _ = synonymous_codons_for_amino_acid(codons, basis_aa)
    max_synonymous_count = AA_MAX_SYNONYMOUS_CODONS.get(basis_aa)
    if max_synonymous_count is None:
        raise ValueError(f"Unsupported amino acid for entropy: {basis_aa}")
    if max_synonymous_count <= 1 or len(set(synonymous_codons)) <= 1:
        return 0.0
    return float(shannon_entropy(synonymous_codons) / np.log2(max_synonymous_count))


def normalized_entropy_reference_basis_or_nan(
    codons: Sequence[str], reference_aa: str
) -> float:
    """Return NaN when any translated codon is nonsynonymous to the reference AA."""
    synonymous_codons: List[str] = []
    for codon in codons:
        if codon == "---":
            continue
        amino_acid = STANDARD_TABLE.forward_table.get(codon)
        if amino_acid is None:
            continue
        if amino_acid != reference_aa:
            return float("nan")
        synonymous_codons.append(codon)

    if not synonymous_codons:
        return float("nan")

    max_synonymous_count = AA_MAX_SYNONYMOUS_CODONS.get(reference_aa)
    if max_synonymous_count is None:
        raise ValueError(f"Unsupported amino acid for entropy: {reference_aa}")
    if max_synonymous_count <= 1 or len(set(synonymous_codons)) <= 1:
        return 0.0
    return float(shannon_entropy(synonymous_codons) / np.log2(max_synonymous_count))


def collect_position_codons(
    protein_path: Path,
    cds_path: Path,
    *,
    cell_mode: str,
) -> Tuple[pd.DataFrame, pd.DataFrame | None, str]:
    aa_ref, codons_ref, reference_species = read_reference_sequences(protein_path, cds_path)
    cds_dict = {
        record.description.split("|")[0].strip(): str(record.seq)
        for record in SeqIO.parse(str(cds_path), "fasta")
    }

    order_counts: Counter[str] = Counter()
    for record in SeqIO.parse(str(protein_path), "fasta"):
        parts = record.description.split("|")
        order_counts[parts[-1].strip()] += 1

    valid_orders = sorted(order for order, count in order_counts.items() if count > 3)
    valid_order_set = set(valid_orders)

    position_codons: Dict[int, Dict[str, List[str]]] = defaultdict(lambda: defaultdict(list))
    for record in SeqIO.parse(str(protein_path), "fasta"):
        parts = record.description.split("|")
        species = parts[0].strip()
        order = parts[-1].strip()

        if species == reference_species:
            continue
        if order not in valid_order_set:
            continue
        if species not in cds_dict:
            continue

        aa_seq = str(record.seq)
        cds_seq = cds_dict[species]
        if cds_seq[-3:] in {"TAA", "TAG", "TGA"}:
            cds_seq = cds_seq[:-3]
        if len(cds_seq) < len(aa_seq) * 3:
            continue

        codon_list = [cds_seq[i : i + 3] for i in range(0, len(aa_seq) * 3, 3)]
        usable_length = min(len(codon_list), len(codons_ref), len(aa_seq))
        for position in range(usable_length):
            codon = codon_list[position]
            if codon == "---" or codon in STANDARD_TABLE.forward_table:
                position_codons[position][order].append(codon)

    positions = sorted(position_codons.keys())
    matrix: List[List[float]] = []
    non_syn_mask: pd.DataFrame | None = None
    if cell_mode == "majority-star-nonsyn":
        non_syn_mask = pd.DataFrame(
            False,
            index=[position + 1 for position in positions],
            columns=valid_orders,
        )

    for position in positions:
        row: List[float] = []
        human_aa = aa_ref[position]
        for order in valid_orders:
            codons = position_codons[position].get(order, [])
            if not codons:
                row.append(np.nan)
                continue

            if cell_mode == "reference-blank-nonsyn":
                row.append(normalized_entropy_reference_basis_or_nan(codons, human_aa))
                continue

            basis_aa = basis_amino_acid_for_cell(codons, human_aa)
            if basis_aa is None:
                row.append(np.nan)
                continue

            _, has_non_synonymous = synonymous_codons_for_amino_acid(codons, basis_aa)
            row.append(normalized_entropy_majority_basis(codons, human_aa))
            if has_non_synonymous and non_syn_mask is not None:
                non_syn_mask.at[position + 1, order] = True
        matrix.append(row)

    entropy_df = pd.DataFrame(
        matrix,
        index=[position + 1 for position in positions],
        columns=valid_orders,
    )
    return entropy_df, non_syn_mask, aa_ref


def draw_region_annotations(
    ax: plt.Axes,
    regions: Sequence[Tuple[int, int, str]],
    *,
    x: float,
    tick_width: float,
    font_size: int,
    line_width: float,
) -> None:
    for start, end, label in regions:
        y0 = start - 0.5
        y1 = end + 0.5
        ax.plot([x, x], [y0, y1], color="black", lw=line_width, clip_on=False)
        ax.plot([x, x + tick_width], [y0, y0], color="black", lw=line_width, clip_on=False)
        ax.plot([x, x + tick_width], [y1, y1], color="black", lw=line_width, clip_on=False)
        ax.text(
            x - 0.12,
            (y0 + y1) / 2,
            label,
            ha="right",
            va="center",
            fontsize=font_size,
            clip_on=False,
            fontweight="bold",
        )


def highlight_color_for_position(position: int) -> str:
    if position in BLUE_POSITIONS:
        return "blue"
    if position in GREEN_POSITIONS:
        return "green"
    return "red"


def draw_highlighted_positions(ax: plt.Axes, aa_ref: str) -> None:
    for position in HIGHLIGHTED_POSITIONS:
        y = position - 1
        aa = aa_ref[y]
        color = highlight_color_for_position(position)

        circle = patches.Circle(
            (-0.3, y + 0.5),
            radius=0.3,
            fill=False,
            edgecolor=color,
            linewidth=1.5,
            clip_on=False,
            zorder=10,
        )
        ax.add_patch(circle)
        ax.text(
            -0.4,
            y + 0.5,
            aa,
            color=color,
            fontsize=7,
            ha="center",
            va="center",
            zorder=11,
            clip_on=False,
            fontweight="bold",
        )


def draw_order_swatches(ax: plt.Axes, orders: Sequence[str]) -> None:
    swatch_height = 3
    gap_top = 0.6
    for index, order in enumerate(orders):
        color = ORDER_COLORS.get(order, "lightgray")
        rectangle = patches.Rectangle(
            (index, -(gap_top + swatch_height)),
            1.0,
            swatch_height,
            linewidth=0,
            facecolor=color,
            edgecolor="none",
            clip_on=False,
            zorder=15,
        )
        ax.add_patch(rectangle)

    ymin, _ = ax.get_ylim()
    ax.set_ylim(ymin, -gap_top)


def draw_non_syn_markers(ax: plt.Axes, non_syn_mask: pd.DataFrame) -> None:
    for row_index, position in enumerate(non_syn_mask.index):
        for column_index, order in enumerate(non_syn_mask.columns):
            if non_syn_mask.at[position, order]:
                ax.text(
                    column_index + 0.5,
                    row_index + 0.5,
                    "*",
                    ha="center",
                    va="center",
                    fontsize=20,
                    color="white",
                    weight="bold",
                    zorder=12,
                )


def build_heatmap(
    *,
    entropy_df: pd.DataFrame,
    non_syn_mask: pd.DataFrame | None,
    aa_ref: str,
    output_stem: Path,
    presentation_mode: bool,
) -> None:
    if entropy_df.empty:
        raise RuntimeError("Entropy matrix is empty; no heatmap can be written.")

    sns.set_theme(style="white")
    fig, ax = plt.subplots(figsize=(14, 14))
    sns.heatmap(
        entropy_df,
        cmap="viridis",
        cbar=not presentation_mode,
        cbar_kws={"label": "Норм. энтропия Шеннона (0–1)"},
        linewidths=0.5,
        linecolor="gray",
        yticklabels=1,
        ax=ax,
    )

    aa_labels = [f"{position} {aa_ref[position - 1]}" for position in entropy_df.index]
    fontsize = 6 if len(entropy_df) > 50 else 8
    ax.set_yticklabels(aa_labels, fontsize=fontsize, rotation=0)

    draw_order_swatches(ax, list(entropy_df.columns))
    if presentation_mode:
        draw_region_annotations(
            ax,
            FULL_PRESENTATION_REGIONS,
            x=-0.72,
            tick_width=0.14,
            font_size=16,
            line_width=2.2,
        )
    else:
        draw_region_annotations(
            ax,
            DEFAULT_REGIONS,
            x=-0.6,
            tick_width=0.1,
            font_size=8,
            line_width=1.5,
        )
    draw_highlighted_positions(ax, aa_ref)
    if non_syn_mask is not None:
        draw_non_syn_markers(ax, non_syn_mask)

    ax.set_xlabel("")
    ax.set_ylabel("")
    if presentation_mode:
        plt.subplots_adjust(left=0.155, right=0.98, top=0.98, bottom=0.04)
    else:
        plt.tight_layout()

    output_stem.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.savefig(output_stem.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def build_order_legend(output_dir: Path) -> None:
    handles = [
        patches.Patch(facecolor=ORDER_COLORS[order], edgecolor="none", label=order)
        for order in ORDER_COLORS
    ]

    fig, ax = plt.subplots(figsize=(6.5, 4.8))
    ax.axis("off")
    ax.legend(
        handles=handles,
        loc="center left",
        ncol=1,
        frameon=False,
        fontsize=16,
        handlelength=1.8,
        handleheight=1.2,
        labelspacing=1.1,
        borderpad=0.4,
    )
    plt.tight_layout()

    legend_stem = output_dir / "sqk_nuc_full_order_legend"
    plt.savefig(legend_stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.savefig(legend_stem.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def build_shannon_scale(output_dir: Path) -> None:
    fig = plt.figure(figsize=(3.2, 7.2))
    color_ax = fig.add_axes([0.42, 0.23, 0.18, 0.66])
    norm = colors.Normalize(vmin=0.0, vmax=1.0)
    scalar = cm.ScalarMappable(norm=norm, cmap="viridis")
    scalar.set_array([])
    colorbar = fig.colorbar(scalar, cax=color_ax, orientation="vertical")
    colorbar.set_ticks([0.0, 0.5, 1.0])
    colorbar.ax.tick_params(labelsize=15, width=1.2, length=6)
    fig.text(
        0.5,
        0.08,
        "Норм. энтропия Шеннона",
        ha="center",
        va="center",
        fontsize=18,
        fontweight="bold",
    )
    fig.text(
        0.5,
        0.035,
        "(0–1)",
        ha="center",
        va="center",
        fontsize=16,
    )

    scale_stem = output_dir / "sqk_nuc_full_shannon_scale"
    plt.savefig(scale_stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.savefig(scale_stem.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def build_dataset(dataset_name: str, input_dir: Path, output_dir: Path) -> pd.DataFrame:
    dataset_config = DATASETS[dataset_name]
    protein_path = input_dir / dataset_config["protein"]
    cds_path = input_dir / dataset_config["cds"]
    require_file(protein_path)
    require_file(cds_path)

    entropy_df, non_syn_mask, aa_ref = collect_position_codons(
        protein_path,
        cds_path,
        cell_mode=str(dataset_config["cell_mode"]),
    )
    build_heatmap(
        entropy_df=entropy_df,
        non_syn_mask=non_syn_mask,
        aa_ref=aa_ref,
        output_stem=output_dir / dataset_config["stem"],
        presentation_mode=bool(dataset_config["presentation_mode"]),
    )
    print(
        f"[{dataset_name}] wrote heatmap with shape "
        f"{entropy_df.shape[0]}x{entropy_df.shape[1]}"
    )
    return entropy_df


def main() -> None:
    args = parse_args()
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    dataset_names = resolve_dataset_names(args.dataset)
    for dataset_name in dataset_names:
        build_dataset(dataset_name, input_dir, output_dir)

    if args.with_legend and "full" in dataset_names:
        build_order_legend(output_dir)
        build_shannon_scale(output_dir)
        print("[legend] wrote full order legend")
        print("[legend] wrote full Shannon scale")


if __name__ == "__main__":
    main()
