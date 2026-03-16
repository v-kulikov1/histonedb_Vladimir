#!/usr/bin/env python
# Batch-build H2A processed tables and heatmaps for the remaining species.

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List

import pandas as pd

from build_bgee_h2a_heatmaps_any_species import BuildResult, build_species_heatmaps


DEFAULT_RAW_DIR = Path(r"C:\Users\USER\Documents\GitHub\histonedb_external_storage\BioAnalyze\raw")
DEFAULT_MERGED = Path(r"CURATED_SET/BioAnalyze/data/merged/mammalia_H2A_merged_with_taxonomy_v4.csv")
DEFAULT_OUT_DIR = Path(r"CURATED_SET/BioAnalyze/figures/heatmaps")
DEFAULT_OUT_PROCESSED_DIR = Path(r"CURATED_SET/BioAnalyze/data/processed")
DEFAULT_AUDIT_OUT = Path(r"CURATED_SET/BioAnalyze/audits/audit_h2a_remaining_species_batch_v4.tsv")

TARGET_SPECIES = [
    "Callithrix jacchus",
    "Canis lupus familiaris",
    "Cavia porcellus",
    "Equus caballus",
    "Felis catus",
    "Heterocephalus glaber",
    "Macaca mulatta",
    "Mus musculus",
    "Oryctolagus cuniculus",
    "Rattus norvegicus",
    "Sus scrofa",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Batch-build H2A present+gold tables, maps, and heatmaps for the "
            "remaining mammalian species using merged v4 labels."
        )
    )
    p.add_argument(
        "--raw-dir",
        default=str(DEFAULT_RAW_DIR),
        help="Directory with per-species Bgee advanced TSV files.",
    )
    p.add_argument(
        "--merged",
        default=str(DEFAULT_MERGED),
        help="Merged v4 H2A dataset with taxonomy, gene_name, IDs, and variant.",
    )
    p.add_argument(
        "--out-dir",
        default=str(DEFAULT_OUT_DIR),
        help="Base output directory for heatmaps.",
    )
    p.add_argument(
        "--out-processed-dir",
        default=str(DEFAULT_OUT_PROCESSED_DIR),
        help="Base output directory for processed TSVs.",
    )
    p.add_argument(
        "--audit-out",
        default=str(DEFAULT_AUDIT_OUT),
        help="Output TSV for batch status audit.",
    )
    p.add_argument(
        "--chunksize",
        type=int,
        default=200000,
        help="Chunk size for streaming Bgee TSV.",
    )
    p.add_argument(
        "--cell-size",
        type=float,
        default=0.7,
        help="Cell size (inches) for square-cell heatmaps.",
    )
    p.add_argument(
        "--min-width",
        type=float,
        default=12.0,
        help="Minimum figure width (inches) for square-cell heatmaps.",
    )
    p.add_argument(
        "--min-height",
        type=float,
        default=8.0,
        help="Minimum figure height (inches) for square-cell heatmaps.",
    )
    return p.parse_args()


def make_error_result(species: str, exc: Exception) -> BuildResult:
    slug = species.strip().lower().replace(" ", "_")
    return BuildResult(
        species=species,
        slug=slug,
        status="error",
        reason=str(exc),
        id_col="hgnc_id" if species == "Homo sapiens" else "vgnc_id",
        canonical_rule="canonical_like",
        allow_partial_splits=True,
    )


def run_batch(args: argparse.Namespace) -> List[BuildResult]:
    raw_dir = Path(args.raw_dir)
    merged_path = Path(args.merged)
    out_dir = Path(args.out_dir)
    out_processed_dir = Path(args.out_processed_dir)

    results: List[BuildResult] = []
    for species in TARGET_SPECIES:
        expr_path = raw_dir / f"{species.replace(' ', '_')}_expr_advanced_all_conditions.tsv"
        try:
            result = build_species_heatmaps(
                species=species,
                expr_path=expr_path,
                merged_path=merged_path,
                out_dir=out_dir,
                out_processed_dir=out_processed_dir,
                canonical_rule="canonical_like",
                allow_partial_splits=True,
                chunksize=args.chunksize,
                square_cells=True,
                cell_size=args.cell_size,
                min_width=args.min_width,
                min_height=args.min_height,
            )
        except Exception as exc:
            result = make_error_result(species, exc)

        results.append(result)
        print(f"{result.species}\t{result.status}\t{result.reason}")

    return results


def write_audit(results: List[BuildResult], audit_out: Path, raw_dir: Path) -> None:
    rows = []
    for result in results:
        row = result.to_dict()
        row["raw_expr_tsv"] = str(raw_dir / f"{result.species.replace(' ', '_')}_expr_advanced_all_conditions.tsv")
        rows.append(row)

    audit_out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(audit_out, sep="\t", index=False)


def main() -> None:
    args = parse_args()
    results = run_batch(args)
    audit_out = Path(args.audit_out)
    write_audit(results, audit_out, Path(args.raw_dir))
    print(f"AUDIT_TSV={audit_out}")


if __name__ == "__main__":
    main()
