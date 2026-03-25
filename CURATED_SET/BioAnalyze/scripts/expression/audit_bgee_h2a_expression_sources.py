#!/usr/bin/env python
"""Audit which Bgee source-specific tracks underlie H2A expression scores."""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Sequence

import pandas as pd


DEFAULT_RAW_DIR = Path(r"C:\Users\USER\Documents\GitHub\histonedb_external_storage\BioAnalyze\raw")
DEFAULT_MERGED = Path(r"CURATED_SET/BioAnalyze/data/merged/mammalia_H2A_merged_with_taxonomy_v4.csv")
DEFAULT_OUT_TSV = Path(r"CURATED_SET/BioAnalyze/audits/bgee_h2a_expression_source_audit_3species.tsv")
DEFAULT_OUT_MD = Path(
    r"CURATED_SET/BioAnalyze/stats/ranking/reports/bgee_h2a_expression_source_audit_3species.md"
)
DEFAULT_SPECIES = ["Homo sapiens", "Bos taurus", "Mus musculus"]

OVERALL_RAW_COLUMNS = {
    "gene_id": "Gene ID",
    "gene_name": "Gene name",
    "anatomical_entity_name": "Anatomical entity name",
    "expression": "Expression",
    "call_quality": "Call quality",
    "expression_score": "Expression score",
}

SOURCE_SCORE_COLUMNS = {
    "affymetrix_expression_score": "Affymetrix expression score",
    "est_expression_score": "EST expression score",
    "in_situ_hybridization_expression_score": "in situ hybridization expression score",
    "rna_seq_expression_score": "RNA-Seq expression score",
    "single_cell_rna_seq_expression_score": "single-cell RNA-Seq expression score",
}

SOURCE_LABELS = {
    "affymetrix_expression_score": "Affymetrix",
    "est_expression_score": "EST",
    "in_situ_hybridization_expression_score": "in situ hybridization",
    "rna_seq_expression_score": "RNA-Seq",
    "single_cell_rna_seq_expression_score": "single-cell RNA-Seq",
}

AUDIT_OUTPUT_COLUMNS = [
    "species",
    "species_slug",
    "gene_id",
    "gene_name",
    "anatomical_entity_name",
    "expression",
    "call_quality",
    "expression_score",
    "affymetrix_expression_score",
    "est_expression_score",
    "in_situ_hybridization_expression_score",
    "rna_seq_expression_score",
    "single_cell_rna_seq_expression_score",
    "active_sources",
    "active_source_count",
    "source_class",
    "heatmap_qualifying",
    "overall_equals_rnaseq_score",
]

@dataclass
class ExampleRow:
    gene_id: str
    gene_name: str
    anatomical_entity_name: str
    source_class: str
    active_sources: str
    expression_score: str
    rna_seq_expression_score: str


@dataclass
class SpeciesAuditSummary:
    species: str
    species_slug: str
    gene_count: int
    matched_raw_rows: int = 0
    qualifying_rows: int = 0
    source_class_counts_all: Counter = field(default_factory=Counter)
    source_class_counts_qualifying: Counter = field(default_factory=Counter)
    non_rna_examples: List[ExampleRow] = field(default_factory=list)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Audit Bgee H2A raw expression rows to determine whether heatmap-relevant "
            "Expression score values are supported by RNA-Seq only, multiple sources, "
            "non-RNA sources, or cannot be resolved."
        )
    )
    p.add_argument(
        "--species",
        action="append",
        default=None,
        help=(
            "Species name as in merged v4 file, e.g. 'Homo sapiens'. "
            "Repeat to audit multiple species. Defaults to Homo sapiens, Bos taurus, Mus musculus."
        ),
    )
    p.add_argument(
        "--raw-dir",
        default=str(DEFAULT_RAW_DIR),
        help="Directory containing per-species Bgee advanced all-conditions TSV files.",
    )
    p.add_argument(
        "--merged",
        default=str(DEFAULT_MERGED),
        help="Merged H2A v4 CSV with species_name and ensembl_gene_id columns.",
    )
    p.add_argument(
        "--out-tsv",
        default=str(DEFAULT_OUT_TSV),
        help="Detailed audit TSV output path.",
    )
    p.add_argument(
        "--out-md",
        default=str(DEFAULT_OUT_MD),
        help="Human-readable Markdown summary output path.",
    )
    return p.parse_args()


def slugify_species(species: str) -> str:
    return species.strip().lower().replace(" ", "_")


def clean_text(value: str) -> str:
    text = (value or "").strip()
    if len(text) >= 2 and text[0] == text[-1] == '"':
        text = text[1:-1].strip()
    if text.upper() == "NA":
        return ""
    return text


def parse_numeric(value: str) -> float | None:
    text = clean_text(value)
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def is_heatmap_qualifying(expression: str, call_quality: str, expression_score: str) -> bool:
    return (
        clean_text(expression) == "present"
        and clean_text(call_quality) == "gold quality"
        and parse_numeric(expression_score) is not None
    )


def active_source_labels(row: Dict[str, str]) -> List[str]:
    active: List[str] = []
    for out_col in SOURCE_SCORE_COLUMNS:
        if clean_text(row[out_col]):
            active.append(SOURCE_LABELS[out_col])
    return active


def classify_sources(active_sources: Sequence[str]) -> str:
    if not active_sources:
        return "unresolved"
    if active_sources == ["RNA-Seq"]:
        return "rna_seq_only"
    if "RNA-Seq" in active_sources:
        return "multi_source"
    return "non_rna_only"


def overall_equals_rnaseq_score(expression_score: str, rna_seq_expression_score: str) -> str:
    overall_num = parse_numeric(expression_score)
    rnaseq_num = parse_numeric(rna_seq_expression_score)
    if overall_num is None or rnaseq_num is None:
        return "not_applicable"
    return "yes" if math.isclose(overall_num, rnaseq_num, rel_tol=1e-9, abs_tol=1e-9) else "no"


def load_h2a_gene_sets(merged_path: Path, species_names: Sequence[str]) -> Dict[str, set[str]]:
    merged = pd.read_csv(merged_path, dtype=str).fillna("")
    required_cols = {"species_name", "ensembl_gene_id"}
    missing = sorted(required_cols.difference(merged.columns))
    if missing:
        raise RuntimeError(f"Missing required merged columns: {', '.join(missing)}")

    gene_sets: Dict[str, set[str]] = {}
    for species in species_names:
        subset = merged[
            merged["species_name"].eq(species) & merged["ensembl_gene_id"].str.strip().ne("")
        ]
        gene_set = {
            gid.strip()
            for gid in subset["ensembl_gene_id"].astype(str).tolist()
            if gid and gid.strip()
        }
        if not gene_set:
            raise RuntimeError(f"No non-empty ensembl_gene_id entries found for species '{species}'")
        gene_sets[species] = gene_set
    return gene_sets


def raw_path_for_species(raw_dir: Path, species: str) -> Path:
    return raw_dir / f"{species.replace(' ', '_')}_expr_advanced_all_conditions.tsv"


def resolve_header_indices(header: str) -> Dict[str, int]:
    header_fields = [clean_text(field) for field in header.rstrip("\r\n").split("\t")]
    raw_name_to_index = {name: idx for idx, name in enumerate(header_fields)}
    indices: Dict[str, int] = {}
    for out_col, raw_col in {**OVERALL_RAW_COLUMNS, **SOURCE_SCORE_COLUMNS}.items():
        if raw_col not in raw_name_to_index:
            raise RuntimeError(f"Missing raw column: {raw_col}")
        indices[out_col] = raw_name_to_index[raw_col]
    return indices


def read_selected_fields(parts: List[str], indices: Dict[str, int]) -> Dict[str, str]:
    row: Dict[str, str] = {}
    for out_col, idx in indices.items():
        row[out_col] = clean_text(parts[idx]) if idx < len(parts) else ""
    return row


def record_example(summary: SpeciesAuditSummary, row: Dict[str, str]) -> None:
    if len(summary.non_rna_examples) >= 5:
        return
    summary.non_rna_examples.append(
        ExampleRow(
            gene_id=row["gene_id"],
            gene_name=row["gene_name"],
            anatomical_entity_name=row["anatomical_entity_name"],
            source_class=row["source_class"],
            active_sources=row["active_sources"],
            expression_score=row["expression_score"],
            rna_seq_expression_score=row["rna_seq_expression_score"],
        )
    )


def audit_species(
    species: str,
    raw_dir: Path,
    gene_ids: set[str],
    writer,
) -> SpeciesAuditSummary:
    raw_path = raw_path_for_species(raw_dir, species)
    if not raw_path.exists():
        raise FileNotFoundError(f"Raw Bgee TSV not found: {raw_path}")

    summary = SpeciesAuditSummary(
        species=species,
        species_slug=slugify_species(species),
        gene_count=len(gene_ids),
    )

    with raw_path.open("r", encoding="utf-8", newline="") as handle:
        header = handle.readline()
        if not header:
            raise RuntimeError(f"Empty raw file: {raw_path}")
        indices = resolve_header_indices(header)

        for line_number, line in enumerate(handle, start=2):
            parts = line.rstrip("\n").split("\t")
            gene_id = clean_text(parts[indices["gene_id"]]) if indices["gene_id"] < len(parts) else ""
            if not gene_id or gene_id not in gene_ids:
                continue

            row = read_selected_fields(parts, indices)
            active_sources = active_source_labels(row)
            row["species"] = species
            row["species_slug"] = summary.species_slug
            row["active_sources"] = "|".join(active_sources) if active_sources else "none"
            row["active_source_count"] = str(len(active_sources))
            row["source_class"] = classify_sources(active_sources)
            row["heatmap_qualifying"] = (
                "yes"
                if is_heatmap_qualifying(
                    row["expression"],
                    row["call_quality"],
                    row["expression_score"],
                )
                else "no"
            )
            row["overall_equals_rnaseq_score"] = overall_equals_rnaseq_score(
                row["expression_score"], row["rna_seq_expression_score"]
            )

            summary.matched_raw_rows += 1
            summary.source_class_counts_all[row["source_class"]] += 1
            if row["heatmap_qualifying"] == "yes":
                summary.qualifying_rows += 1
                summary.source_class_counts_qualifying[row["source_class"]] += 1
                if row["source_class"] != "rna_seq_only":
                    record_example(summary, row)

            writer.writerow({col: row.get(col, "") for col in AUDIT_OUTPUT_COLUMNS})

            if line_number % 2_000_000 == 0:
                print(
                    f"{species}\tprocessed_line={line_number}\t"
                    f"matched_rows={summary.matched_raw_rows}\tqualifying_rows={summary.qualifying_rows}"
                )

    return summary


def qualifying_pct(summary: SpeciesAuditSummary, source_class: str) -> float:
    if summary.qualifying_rows == 0:
        return 0.0
    return (summary.source_class_counts_qualifying[source_class] / summary.qualifying_rows) * 100.0


def conclusion_for_species(summary: SpeciesAuditSummary) -> str:
    qualifying_total = summary.qualifying_rows
    rnaseq_only = summary.source_class_counts_qualifying["rna_seq_only"]
    if qualifying_total == 0:
        return f"For qualifying H2A rows in {summary.species}, no rows met the present+gold filter."
    if rnaseq_only == qualifying_total:
        return (
            f"For qualifying H2A rows in {summary.species}, expression score is exclusively RNA-Seq-based."
        )
    if rnaseq_only > (qualifying_total / 2.0):
        return (
            f"For qualifying H2A rows in {summary.species}, expression score is predominantly RNA-Seq-based."
        )
    return (
        f"For qualifying H2A rows in {summary.species}, expression score is not exclusively RNA-Seq-based."
    )


def build_summary_table(summaries: Sequence[SpeciesAuditSummary]) -> List[str]:
    lines = [
        "| Species | H2A genes scanned | Matched raw rows | Qualifying rows | RNA-Seq only | Multi-source | Non-RNA only | Unresolved |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for summary in summaries:
        lines.append(
            "| {species} | {gene_count} | {matched_raw_rows} | {qualifying_rows} | {rna_seq_only} | {multi_source} | {non_rna_only} | {unresolved} |".format(
                species=summary.species,
                gene_count=summary.gene_count,
                matched_raw_rows=summary.matched_raw_rows,
                qualifying_rows=summary.qualifying_rows,
                rna_seq_only=summary.source_class_counts_qualifying["rna_seq_only"],
                multi_source=summary.source_class_counts_qualifying["multi_source"],
                non_rna_only=summary.source_class_counts_qualifying["non_rna_only"],
                unresolved=summary.source_class_counts_qualifying["unresolved"],
            )
        )
    return lines


def example_table_lines(examples: Sequence[ExampleRow]) -> List[str]:
    lines = [
        "| Gene ID | Gene name | Anatomical entity | Source class | Active sources | Expression score | RNA-Seq expression score |",
        "| --- | --- | --- | --- | --- | ---: | ---: |",
    ]
    for example in examples:
        lines.append(
            "| {gene_id} | {gene_name} | {anatomical_entity_name} | {source_class} | {active_sources} | {expression_score} | {rna_seq_expression_score} |".format(
                gene_id=example.gene_id,
                gene_name=example.gene_name or "-",
                anatomical_entity_name=example.anatomical_entity_name or "-",
                source_class=example.source_class,
                active_sources=example.active_sources,
                expression_score=example.expression_score or "-",
                rna_seq_expression_score=example.rna_seq_expression_score or "-",
            )
        )
    return lines


def write_markdown_report(
    out_md: Path,
    out_tsv: Path,
    merged_path: Path,
    raw_dir: Path,
    summaries: Sequence[SpeciesAuditSummary],
) -> None:
    generated_at = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    lines: List[str] = [
        "# Bgee H2A Expression Source Audit",
        "",
        f"Generated: {generated_at}",
        "",
        f"- Merged H2A file: `{merged_path}`",
        f"- Raw Bgee dir: `{raw_dir}`",
        f"- Detailed audit TSV: `{out_tsv}`",
        "",
        "## Overview",
        "",
    ]
    lines.extend(build_summary_table(summaries))

    for summary in summaries:
        lines.extend(
            [
                "",
                f"## {summary.species}",
                "",
                f"- H2A genes scanned: {summary.gene_count}",
                f"- Matched raw rows: {summary.matched_raw_rows}",
                f"- Heatmap-qualifying rows (`present` + `gold quality` + non-empty overall `Expression score`): {summary.qualifying_rows}",
                (
                    "- Qualifying source classes: "
                    f"`rna_seq_only={summary.source_class_counts_qualifying['rna_seq_only']}` "
                    f"({qualifying_pct(summary, 'rna_seq_only'):.2f}%), "
                    f"`multi_source={summary.source_class_counts_qualifying['multi_source']}` "
                    f"({qualifying_pct(summary, 'multi_source'):.2f}%), "
                    f"`non_rna_only={summary.source_class_counts_qualifying['non_rna_only']}` "
                    f"({qualifying_pct(summary, 'non_rna_only'):.2f}%), "
                    f"`unresolved={summary.source_class_counts_qualifying['unresolved']}` "
                    f"({qualifying_pct(summary, 'unresolved'):.2f}%)"
                ),
                f"- Conclusion: {conclusion_for_species(summary)}",
            ]
        )
        if summary.non_rna_examples:
            lines.extend(["", "Non-`rna_seq_only` qualifying examples:", ""])
            lines.extend(example_table_lines(summary.non_rna_examples))
        else:
            lines.extend(["", "No non-`rna_seq_only` qualifying examples found."])

    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_audit(
    species_names: Sequence[str],
    raw_dir: Path,
    merged_path: Path,
    out_tsv: Path,
    out_md: Path,
) -> List[SpeciesAuditSummary]:
    gene_sets = load_h2a_gene_sets(merged_path, species_names)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    summaries: List[SpeciesAuditSummary] = []
    with out_tsv.open("w", encoding="utf-8", newline="") as handle:
        dict_writer = csv.DictWriter(handle, fieldnames=AUDIT_OUTPUT_COLUMNS, delimiter="\t")
        dict_writer.writeheader()
        for species in species_names:
            summary = audit_species(
                species=species,
                raw_dir=raw_dir,
                gene_ids=gene_sets[species],
                writer=dict_writer,
            )
            summaries.append(summary)
            print(
                f"{summary.species}\tgenes={summary.gene_count}\tmatched_rows={summary.matched_raw_rows}\t"
                f"qualifying_rows={summary.qualifying_rows}"
            )

    write_markdown_report(
        out_md=out_md,
        out_tsv=out_tsv,
        merged_path=merged_path,
        raw_dir=raw_dir,
        summaries=summaries,
    )
    return summaries


def main() -> None:
    args = parse_args()
    species_names = args.species or list(DEFAULT_SPECIES)
    summaries = run_audit(
        species_names=species_names,
        raw_dir=Path(args.raw_dir),
        merged_path=Path(args.merged),
        out_tsv=Path(args.out_tsv),
        out_md=Path(args.out_md),
    )
    print(f"AUDIT_TSV={Path(args.out_tsv)}")
    print(f"SUMMARY_MD={Path(args.out_md)}")
    print(f"SPECIES={','.join(summary.species for summary in summaries)}")


if __name__ == "__main__":
    main()
