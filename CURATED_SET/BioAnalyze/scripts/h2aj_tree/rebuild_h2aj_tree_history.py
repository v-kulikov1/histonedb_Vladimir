#!/usr/bin/env python
"""Rebuild the historical H2A.J tree-preparation workflow and evidence trail."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shutil
import sys
import zipfile
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Sequence
from xml.etree import ElementTree as ET

from Bio import Phylo, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

BIOANALYZE_SCRIPTS_ROOT = Path(__file__).resolve().parents[1]
if str(BIOANALYZE_SCRIPTS_ROOT) not in sys.path:
    sys.path.insert(0, str(BIOANALYZE_SCRIPTS_ROOT))

from bioanalyze_paths import (
    get_bioanalyze_raw_root,
    get_documents_root,
    get_downloads_root,
)


SCRIPT_DIR = Path(__file__).resolve().parent
BIOANALYZE_ROOT = SCRIPT_DIR.parents[1]
DEFAULT_OUTPUT_ROOT = BIOANALYZE_ROOT / "data" / "h2aj_tree"
DEFAULT_EXTERNAL_ROOT = get_bioanalyze_raw_root() / "tree"
DEFAULT_CODONS_ROOT = get_bioanalyze_raw_root() / "codons"

DEFAULT_HISTONES_ROOT = get_documents_root() / "Гистоны"
DEFAULT_MY_TEST_ROOT = get_documents_root() / "My_test"
DEFAULT_GRANT_ROOT = get_documents_root() / "Работа над грантом HistoneJ"
DEFAULT_CHAT_EXPORT = (
    get_downloads_root() / "Telegram Desktop" / "ChatExport_2026-03-29" / "result.json"
)
DEFAULT_DRAFT_DOCX = get_downloads_root() / "Telegram Desktop" / "Draft2.docx"
DEFAULT_DROP_CLEAN_IDS = [
    "Myotis-lucifugus|XM_006084274",
    "Homo-sapiens|AK303301",
    "Homo-sapiens|AL133626",
]

MANIFEST_COLUMNS = [
    "category",
    "logical_name",
    "source_path",
    "copied_to",
    "sha256",
    "mtime",
    "records",
    "sites",
    "note",
]

CHECKPOINT_COLUMNS = [
    "logical_name",
    "archived_path",
    "record_count",
    "site_count",
    "source_mtime",
    "note",
]

NOTEBOOK_EVIDENCE_COLUMNS = [
    "source_name",
    "source_path",
    "cell_index",
    "cell_type",
    "match_tag",
    "snippet",
]

CHAT_EVIDENCE_COLUMNS = [
    "date",
    "match_tag",
    "text",
]

DRAFT2_COLUMNS = [
    "slide",
    "text",
]

CYRILLIC_SMALL_ES = "с"
CYRILLIC_CAPITAL_ES = "С"


@dataclass(frozen=True)
class SourcePaths:
    histones_root: Path
    my_test_root: Path
    grant_root: Path
    chat_export_path: Path
    draft_docx_path: Path


@dataclass(frozen=True)
class SourceSpec:
    logical_name: str
    category: str
    source_path: Path
    destination_subdir: str
    note: str
    required: bool = True
    reuse_codons_archive: bool = False


@dataclass(frozen=True)
class ArchiveRecord:
    logical_name: str
    category: str
    source_path: Path
    archived_path: Path
    sha256: str
    mtime: str
    records: str
    sites: str
    note: str


NOTEBOOK_PATTERNS: Sequence[tuple[str, str]] = [
    ("blast_hsp_fragment", r"hsp\.sbjct\.replace\('-', ''\)"),
    ("partial_codon_warning", r"partial codon"),
    ("histones_csv_ch2a", r"histones\.csv|variant_group\s*==\s*['\"]cH2A['\"]"),
    ("phylip_conversion", r"All_AA_0206\.phy|All_NUC_1606\.phy|Converted\s+\d+\s+sequences"),
    ("rooting_monophyly_check", r"is_monophyletic|root_with_outgroup|common_ancestor"),
    ("coded_by_fetch", r"coded_by|nucleotide_sequences\.fasta"),
    ("translate_no_gaps", r"translate_no_gaps|protein_from_SQK_nuc|cleaned_nucleotide_output\.fasta"),
    ("motif_position_131", r"target_position\s*=\s*131"),
    ("sqk_shk_filter", r"SQK|SHK"),
]

CHAT_PATTERNS: Sequence[tuple[str, str]] = [
    ("phyml_hky85_memory", r"PhyML|HKY85|bootstrap"),
    ("mammal_outgroup_refine", r"млекоп|mammal|cH2A"),
    ("nucleotide_tree_wrong", r"неправильно построил дерево|белков.*вместо нуклеотид"),
    ("sqk_cleanup", r"нет мотива SQK|SQK|SHK"),
    ("short_sequences", r"коротк|128"),
    ("aligned_fragments", r"aligned sequences|aligned sequence|обрезал|C-кон"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Archive the historical H2A.J tree materials, rebuild the sequence "
            "inputs, and export evidence about why the nucleotide tree failed."
        )
    )
    parser.add_argument(
        "--phase",
        default="all",
        choices=["archive", "evidence", "aa", "nuc", "postprocess", "all"],
        help="Pipeline phase to execute.",
    )
    parser.add_argument(
        "--profile",
        default="both",
        choices=["historical", "clean", "both"],
        help="Output profile to write.",
    )
    parser.add_argument(
        "--external-root",
        default=str(DEFAULT_EXTERNAL_ROOT),
        help="External archive directory for copied historical evidence.",
    )
    parser.add_argument(
        "--output-root",
        default=str(DEFAULT_OUTPUT_ROOT),
        help="Data output root for reconstructed tree inputs and evidence.",
    )
    parser.add_argument(
        "--histones-root",
        default=str(DEFAULT_HISTONES_ROOT),
        help="Root directory that contains the historical tree files from the histone workspace.",
    )
    parser.add_argument(
        "--my-test-root",
        default=str(DEFAULT_MY_TEST_ROOT),
        help="Root directory with the local Jupyter notebooks used for extra preparation.",
    )
    parser.add_argument(
        "--grant-root",
        default=str(DEFAULT_GRANT_ROOT),
        help="Root directory with the grant-era notebooks (BLAST / SQK-SHK filtering).",
    )
    parser.add_argument(
        "--chat-export",
        default=str(DEFAULT_CHAT_EXPORT),
        help="Path to Telegram ChatExport result.json.",
    )
    parser.add_argument(
        "--draft-docx",
        default=str(DEFAULT_DRAFT_DOCX),
        help="Path to the Draft2.docx source.",
    )
    parser.add_argument(
        "--drop-clean-id",
        action="append",
        default=list(DEFAULT_DROP_CLEAN_IDS),
        help=(
            "Identifier prefix to remove from clean-profile outputs. "
            "May be repeated; applies only to clean outputs."
        ),
    )
    return parser.parse_args()


def resolve_histones_base(root: Path) -> Path:
    summer_dir = root / "Лето2025"
    return summer_dir if summer_dir.exists() else root


def resolve_my_test_tree_base(root: Path) -> Path:
    tree_dir = root / "Древо 21.05"
    return tree_dir if tree_dir.exists() else root


def build_source_specs(paths: SourcePaths) -> list[SourceSpec]:
    histones_base = resolve_histones_base(paths.histones_root)
    my_test_tree = resolve_my_test_tree_base(paths.my_test_root)

    histones_specs = [
        SourceSpec("H2AJ.fasta", "histones", histones_base / "H2AJ.fasta", "histones", "Early aligned H2A.J protein set."),
        SourceSpec("cH2A.fasta", "histones", histones_base / "cH2A.fasta", "histones", "Canonical H2A outgroup protein set."),
        SourceSpec("platypus.fasta", "histones", histones_base / "platypus.fasta", "histones", "Historical nonplacental protein set, labelled as platypus."),
        SourceSpec("All_AA_0206.fasta", "histones", histones_base / "All_AA_0206.fasta", "histones", "02 June 2025 broad outgroup amino-acid alignment."),
        SourceSpec("All_AA_0206.phy", "histones", histones_base / "All_AA_0206.phy", "histones", "PhyML-ready PHYLIP for the 02 June 2025 amino-acid alignment."),
        SourceSpec("all_aa_0206_1__phy_phyml.zip", "histones", histones_base / "all_aa_0206_1__phy_phyml.zip", "histones", "Archived PhyML web result bundle for the confirmed AA run."),
        SourceSpec("All_AA_1006.fasta", "histones", histones_base / "All_AA_1006.fasta", "histones", "10 June 2025 mammal cH2A + nonplacental amino-acid alignment."),
        SourceSpec("All_AA_1006.phy", "histones", histones_base / "All_AA_1006.phy", "histones", "PHYLIP for the 10 June 2025 mammal/nonplacental AA set."),
        SourceSpec("Nuc_merged_sequences.fasta", "histones", histones_base / "Nuc_merged_sequences.fasta", "histones", "Merged nucleotide FASTA before the 16 June alignment."),
        SourceSpec("All_NUC_1606.fasta", "histones", histones_base / "All_NUC_1606.fasta", "histones", "16 June 2025 nucleotide alignment."),
        SourceSpec("All_NUC_1606.phy", "histones", histones_base / "All_NUC_1606.phy", "histones", "PHYLIP for the 16 June 2025 nucleotide alignment."),
        SourceSpec("nuc_tree_bootstrap.nwk", "histones", histones_base / "nuc_tree_bootstrap.nwk", "histones", "Historical unrooted/bootstrapped nucleotide tree."),
        SourceSpec("rooted_tree.nwk", "histones", histones_base / "rooted_tree.nwk", "histones", "Historical rooted nucleotide tree written on 17 June 2025."),
        SourceSpec("SQK_nuc.fasta", "histones", histones_base / "SQK_nuc.fasta", "histones", "08 July 2025 SQK nucleotide set.", reuse_codons_archive=True),
        SourceSpec("SQK_nuc(without short).fasta", "histones", histones_base / "SQK_nuc(without short).fasta", "histones", "08 July 2025 SQK nucleotide set without short sequences.", reuse_codons_archive=True),
        SourceSpec("protein_from_SQK_nuc.fasta", "histones", histones_base / "protein_from_SQK_nuc.fasta", "histones", "Translated proteins from SQK_nuc.", reuse_codons_archive=True),
        SourceSpec("protein_from_SQK_nuc(without short).fasta", "histones", histones_base / "protein_from_SQK_nuc(without short).fasta", "histones", "Translated proteins from SQK_nuc(without short).", reuse_codons_archive=True),
        SourceSpec("cH2A_nuc + platypus_nuc + SQK_nuc.fasta", "histones", histones_base / "cH2A_nuc + platypus_nuc + SQK_nuc.fasta", "histones", "09 July 2025 merged nucleotide set."),
        SourceSpec("All_NUC_0907.phy", "histones", histones_base / "All_NUC_0907.phy", "histones", "PHYLIP for the July merged nucleotide set."),
    ]

    my_test_specs = [
        SourceSpec("Phylogenetic tree.ipynb", "my_test", paths.my_test_root / "Phylogenetic tree.ipynb", "my_test", "Early local FastTree notebook."),
        SourceSpec("02_06 H2A.J + H2A + platypus .ipynb", "my_test", my_test_tree / "02_06 H2A.J + H2A + platypus .ipynb", "my_test", "Later notebook with MUSCLE alignment and protein tree preparation."),
        SourceSpec("All_NUC_0907.fasta", "my_test", my_test_tree / "All_NUC_0907.fasta", "my_test", "July 2025 aligned nucleotide FASTA present only in My_test."),
        SourceSpec("Protein_merged_sequences2907.fasta", "my_test", my_test_tree / "Protein_merged_sequences2907.fasta", "my_test", "Late-July merged protein FASTA."),
        SourceSpec("Protein_merged.fasta", "my_test", my_test_tree / "Protein_merged.fasta", "my_test", "Late merged protein source file."),
        SourceSpec("Nuc_merged_sequences (2).fasta", "my_test", my_test_tree / "Nuc_merged_sequences (2).fasta", "my_test", "My_test copy of the July merged nucleotide set."),
    ]

    grant_specs = [
        SourceSpec("02_06 File_for_alignment_outgroup.ipynb", "grant", paths.grant_root / "02_06 File_for_alignment_outgroup.ipynb", "grant", "Notebook with cH2A extraction, PHYLIP conversion, and rooting logic."),
        SourceSpec("15_05 BLAST N .ipynb", "grant", paths.grant_root / "15_05 BLAST N .ipynb", "grant", "Key notebook showing BLAST-HSP fragment handling for nucleotide hits."),
        SourceSpec("07_07 Фильтруем нуклеотиды из BLAST.ipynb", "grant", paths.grant_root / "07_07 Фильтруем нуклеотиды из BLAST.ipynb", "grant", "Notebook for July cleanup of BLAST-derived nucleotide hits."),
        SourceSpec("Копия Фильтрация SHK_SQK.ipynb", "grant", paths.grant_root / "Копия Фильтрация SHK_SQK.ipynb", "grant", "Notebook that classifies SQK/SHK by aligned position."),
        SourceSpec("Важные скрипты.ipynb", "grant", paths.grant_root / "Важные скрипты.ipynb", "grant", "Notebook with nucleotide cleanup, translation, coded_by fetch, and FASTA reordering."),
        SourceSpec("Optimize tree 11_06.ipynb", "grant", paths.grant_root / "Optimize tree 11_06.ipynb", "grant", "Tree optimization notebook from the mammal/nonplacental stage.", required=False),
        SourceSpec("BLAST .ipynb", "grant", paths.grant_root / "BLAST .ipynb", "grant", "Earlier BLAST screening notebook.", required=False),
        SourceSpec("Копия H2AJ_ филогенетическое древо.ipynb", "grant", paths.grant_root / "Копия H2AJ_ филогенетическое древо.ipynb", "grant", "Additional phylogenetic notebook snapshot.", required=False),
    ]

    telegram_specs = [
        SourceSpec("result.json", "telegram", paths.chat_export_path, "telegram", "Telegram chat export with the historical discussion."),
        SourceSpec("Draft2.docx", "telegram", paths.draft_docx_path, "telegram", "Draft presentation used as an additional evidence source."),
    ]
    return histones_specs + my_test_specs + grant_specs + telegram_specs


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def count_fasta_sites(records: Sequence[SeqRecord]) -> str:
    if not records:
        return ""
    lengths = {len(record.seq) for record in records}
    return str(lengths.pop()) if len(lengths) == 1 else ""


def summarize_file(path: Path) -> tuple[str, str]:
    suffix = path.suffix.lower()
    if suffix in {".fasta", ".fa", ".faa", ".fna"}:
        try:
            records = list(SeqIO.parse(str(path), "fasta"))
        except Exception:
            return "", ""
        if not records:
            return "", ""
        return str(len(records)), count_fasta_sites(records)

    if suffix == ".phy":
        try:
            first_line = path.read_text(encoding="utf-8", errors="replace").splitlines()[0]
        except Exception:
            return "", ""
        match = re.match(r"\s*(\d+)\s+(\d+)", first_line)
        if match:
            return match.group(1), match.group(2)
        return "", ""

    return "", ""


def mtime_iso(path: Path) -> str:
    return datetime.fromtimestamp(path.stat().st_mtime).isoformat(timespec="seconds")


def write_tsv(path: Path, fieldnames: Sequence[str], rows: Sequence[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in fieldnames})


def copy_if_needed(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and sha256_file(destination) == sha256_file(source):
        return
    shutil.copy2(source, destination)


def reuse_codons_path(source_name: str) -> Path:
    return DEFAULT_CODONS_ROOT / source_name


def archive_sources(source_paths: SourcePaths, external_root: Path) -> dict[str, ArchiveRecord]:
    external_root.mkdir(parents=True, exist_ok=True)
    specs = build_source_specs(source_paths)
    records: dict[str, ArchiveRecord] = {}
    manifest_rows: list[dict[str, str]] = []

    for spec in specs:
        if not spec.source_path.exists():
            if spec.required:
                raise FileNotFoundError(f"Missing required source: {spec.source_path}")
            continue

        note = spec.note
        archived_path = external_root / spec.destination_subdir / spec.source_path.name
        source_sha = sha256_file(spec.source_path)

        if spec.reuse_codons_archive:
            codons_path = reuse_codons_path(spec.source_path.name)
            if codons_path.exists() and sha256_file(codons_path) == source_sha:
                archived_path = codons_path
                note = f"{note} Reused existing raw/codons archive copy."
            else:
                copy_if_needed(spec.source_path, archived_path)
        else:
            copy_if_needed(spec.source_path, archived_path)

        item_records, item_sites = summarize_file(archived_path)
        record = ArchiveRecord(
            logical_name=spec.logical_name,
            category=spec.category,
            source_path=spec.source_path,
            archived_path=archived_path,
            sha256=source_sha,
            mtime=mtime_iso(spec.source_path),
            records=item_records,
            sites=item_sites,
            note=note,
        )
        records[spec.logical_name] = record
        manifest_rows.append(
            {
                "category": record.category,
                "logical_name": record.logical_name,
                "source_path": str(record.source_path),
                "copied_to": str(record.archived_path),
                "sha256": record.sha256,
                "mtime": record.mtime,
                "records": record.records,
                "sites": record.sites,
                "note": record.note,
            }
        )

    manifest_rows.sort(key=lambda row: (row["category"], row["logical_name"]))
    write_tsv(external_root / "tree_archive_manifest.tsv", MANIFEST_COLUMNS, manifest_rows)
    return records


def selected_profiles(profile: str) -> list[str]:
    if profile == "both":
        return ["historical", "clean"]
    return [profile]


def flatten_notebook_output(output: dict[str, Any]) -> str:
    chunks: list[str] = []
    text_value = output.get("text")
    if isinstance(text_value, str):
        chunks.append(text_value)
    elif isinstance(text_value, list):
        chunks.extend(str(item) for item in text_value)

    for key in ["data", "traceback"]:
        value = output.get(key)
        if isinstance(value, dict):
            for nested in value.values():
                if isinstance(nested, str):
                    chunks.append(nested)
                elif isinstance(nested, list):
                    chunks.extend(str(item) for item in nested)
        elif isinstance(value, list):
            chunks.extend(str(item) for item in value)

    if output.get("ename") and output.get("evalue"):
        chunks.append(f"{output['ename']}: {output['evalue']}")
    return "\n".join(chunks)


def load_notebook(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def collapse_whitespace(text: str) -> str:
    return re.sub(r"\s+", " ", text).strip()


def truncate(text: str, *, limit: int = 320) -> str:
    clean = collapse_whitespace(text)
    if len(clean) <= limit:
        return clean
    return clean[: limit - 3].rstrip() + "..."


def collect_notebook_evidence(
    archive_records: dict[str, ArchiveRecord],
    evidence_dir: Path,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    notebook_names = [
        "02_06 File_for_alignment_outgroup.ipynb",
        "15_05 BLAST N .ipynb",
        "07_07 Фильтруем нуклеотиды из BLAST.ipynb",
        "Копия Фильтрация SHK_SQK.ipynb",
        "Важные скрипты.ipynb",
        "Phylogenetic tree.ipynb",
        "02_06 H2A.J + H2A + platypus .ipynb",
    ]

    for notebook_name in notebook_names:
        record = archive_records.get(notebook_name)
        if record is None:
            continue
        notebook = load_notebook(record.archived_path)
        for cell_index, cell in enumerate(notebook.get("cells", [])):
            cell_source = "".join(cell.get("source", []))
            output_text = "\n".join(
                flatten_notebook_output(output)
                for output in cell.get("outputs", [])
            )
            combined = "\n".join(part for part in [cell_source, output_text] if part)
            if not combined.strip():
                continue
            for match_tag, pattern in NOTEBOOK_PATTERNS:
                if re.search(pattern, combined, flags=re.IGNORECASE):
                    rows.append(
                        {
                            "source_name": notebook_name,
                            "source_path": str(record.archived_path),
                            "cell_index": str(cell_index),
                            "cell_type": str(cell.get("cell_type", "")),
                            "match_tag": match_tag,
                            "snippet": truncate(combined),
                        }
                    )

    deduped: list[dict[str, str]] = []
    seen: set[tuple[str, str, str]] = set()
    for row in rows:
        key = (row["source_name"], row["cell_index"], row["match_tag"])
        if key in seen:
            continue
        seen.add(key)
        deduped.append(row)

    write_tsv(
        evidence_dir / "notebook_evidence.tsv",
        NOTEBOOK_EVIDENCE_COLUMNS,
        deduped,
    )
    return deduped


def flatten_telegram_text(value: Any) -> str:
    if isinstance(value, str):
        return value
    if isinstance(value, list):
        return "".join(flatten_telegram_text(item) for item in value)
    if isinstance(value, dict):
        if "text" in value:
            return flatten_telegram_text(value["text"])
        return "".join(flatten_telegram_text(v) for v in value.values())
    return str(value)


def collect_chat_evidence(
    archive_records: dict[str, ArchiveRecord],
    evidence_dir: Path,
) -> list[dict[str, str]]:
    record = archive_records["result.json"]
    payload = json.loads(record.archived_path.read_text(encoding="utf-8"))
    rows: list[dict[str, str]] = []
    for message in payload.get("messages", []):
        if message.get("type") != "message":
            continue
        text = collapse_whitespace(flatten_telegram_text(message.get("text", "")))
        if not text:
            continue
        for match_tag, pattern in CHAT_PATTERNS:
            if re.search(pattern, text, flags=re.IGNORECASE):
                rows.append(
                    {
                        "date": str(message.get("date", "")),
                        "match_tag": match_tag,
                        "text": truncate(text, limit=420),
                    }
                )

    important_dates = {
        "2025-06-03",
        "2025-06-11",
        "2025-06-16",
        "2025-06-17",
        "2025-07-07",
        "2025-07-08",
        "2025-07-18",
        "2025-07-23",
        "2025-07-30",
    }
    dated_rows: list[dict[str, str]] = []
    seen_dates: set[tuple[str, str]] = set()
    for row in rows:
        row_date = row["date"][:10]
        if row_date in important_dates:
            key = (row["date"], row["text"])
            if key in seen_dates:
                continue
            seen_dates.add(key)
            dated_rows.append(row)

    write_tsv(evidence_dir / "chat_evidence.tsv", CHAT_EVIDENCE_COLUMNS, dated_rows)
    return dated_rows


def extract_docx_paragraphs(path: Path) -> list[str]:
    namespace = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
    with zipfile.ZipFile(path) as archive:
        document_xml = archive.read("word/document.xml")
    root = ET.fromstring(document_xml)
    paragraphs: list[str] = []
    for paragraph in root.findall(".//w:p", namespace):
        text_chunks = [node.text or "" for node in paragraph.findall(".//w:t", namespace)]
        text = "".join(text_chunks).strip()
        if text:
            paragraphs.append(text)
    return paragraphs


def parse_draft2_slides(paragraphs: Sequence[str]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    current_slide = ""
    current_text: list[str] = []

    def flush() -> None:
        if current_slide:
            rows.append(
                {
                    "slide": current_slide,
                    "text": collapse_whitespace(" ".join(current_text)),
                }
            )

    for paragraph in paragraphs:
        match = re.fullmatch(r"Slide\s+(\d+):", paragraph)
        if match:
            flush()
            current_slide = match.group(1)
            current_text = []
        else:
            current_text.append(paragraph)
    flush()
    return rows


def collect_draft2_evidence(
    archive_records: dict[str, ArchiveRecord],
    evidence_dir: Path,
) -> list[dict[str, str]]:
    record = archive_records["Draft2.docx"]
    paragraphs = extract_docx_paragraphs(record.archived_path)
    rows = parse_draft2_slides(paragraphs)
    write_tsv(evidence_dir / "draft2_slides.tsv", DRAFT2_COLUMNS, rows)
    (evidence_dir / "draft2_extracted.txt").write_text(
        "\n".join(paragraphs),
        encoding="utf-8",
    )
    return rows


def parse_phyml_bundle(zip_path: Path) -> dict[str, str]:
    with zipfile.ZipFile(zip_path) as archive:
        stats_text = archive.read("all_aa_0206_1__phy_phyml_stats.txt").decode(
            "utf-8", errors="replace"
        )
        stdout_text = archive.read("all_aa_0206_1__phy_stdout.txt").decode(
            "utf-8", errors="replace"
        )

    def capture(pattern: str, text: str) -> str:
        match = re.search(pattern, text, flags=re.MULTILINE)
        return match.group(1).strip() if match else ""

    return {
        "sequence_filename": capture(r"Sequence filename:\s+(.+)", stats_text),
        "data_type": capture(r"Data type:\s+(.+)", stdout_text),
        "initial_tree": capture(r"Starting tree:\s+(.+)", stdout_text),
        "model": capture(r"Model(?: of amino acids substitution| name)?:\s+(.+)", stats_text + "\n" + stdout_text),
        "taxa": capture(r"Number of taxa:\s+(.+)", stats_text),
        "log_likelihood": capture(r"Log-likelihood:\s+(.+)", stats_text),
        "rate_model": capture(r"RAS model:\s+(.+)", stdout_text),
        "rate_classes": capture(r"Number of subst\. rate catgs:\s+(.+)", stdout_text),
        "supports": capture(r"Compute approximate likelihood ratio test:\s+(.+)", stdout_text),
        "runtime": capture(r"Time used:\s+(.+)", stats_text),
        "command_line": capture(r"Command line:\s+(.+)", stdout_text),
    }


def write_phyml_runbook(
    archive_records: dict[str, ArchiveRecord],
    evidence_dir: Path,
) -> dict[str, str]:
    bundle_path = archive_records["all_aa_0206_1__phy_phyml.zip"].archived_path
    parsed = parse_phyml_bundle(bundle_path)
    lines = [
        "# Confirmed PhyML AA Run (03 June 2025)",
        "",
        f"- Archive bundle: `{bundle_path}`",
        f"- Sequence file on server: `{parsed['sequence_filename']}`",
        f"- Data type: `{parsed['data_type']}`",
        f"- Initial tree: `{parsed['initial_tree']}`",
        f"- Model: `{parsed['model']}`",
        f"- Taxa: `{parsed['taxa']}`",
        f"- Rate model: `{parsed['rate_model']}`",
        f"- Rate classes: `{parsed['rate_classes']}`",
        f"- Branch support mode: `{parsed['supports']}`",
        f"- Log-likelihood: `{parsed['log_likelihood']}`",
        f"- Runtime: `{parsed['runtime']}`",
        f"- Command line: `{parsed['command_line']}`",
        "",
        "Interpretation:",
        "- This is the strongest direct evidence for the confirmed successful web run.",
        "- It proves the archived protein tree was run as amino-acid data with `LG`, `BioNJ`, and `FreeRate(3)`.",
        "- It conflicts with the earlier chat-memory mention of `HKY85`; the bundle wins because it is the original server output.",
    ]
    runbook_path = evidence_dir / "phyml_aa_0206_runbook.md"
    runbook_path.write_text("\n".join(lines), encoding="utf-8")
    return parsed


def write_sequence_checkpoints(
    archive_records: dict[str, ArchiveRecord],
    evidence_dir: Path,
) -> list[dict[str, str]]:
    names = [
        "H2AJ.fasta",
        "cH2A.fasta",
        "platypus.fasta",
        "All_AA_0206.fasta",
        "All_AA_0206.phy",
        "All_AA_1006.fasta",
        "All_AA_1006.phy",
        "All_NUC_1606.fasta",
        "All_NUC_1606.phy",
        "SQK_nuc.fasta",
        "SQK_nuc(without short).fasta",
        "protein_from_SQK_nuc.fasta",
        "protein_from_SQK_nuc(without short).fasta",
        "cH2A_nuc + platypus_nuc + SQK_nuc.fasta",
        "All_NUC_0907.fasta",
        "All_NUC_0907.phy",
    ]
    rows: list[dict[str, str]] = []
    for name in names:
        record = archive_records.get(name)
        if record is None:
            continue
        rows.append(
            {
                "logical_name": record.logical_name,
                "archived_path": str(record.archived_path),
                "record_count": record.records,
                "site_count": record.sites,
                "source_mtime": record.mtime,
                "note": record.note,
            }
        )
    write_tsv(evidence_dir / "historical_checkpoints.tsv", CHECKPOINT_COLUMNS, rows)
    return rows


def build_evidence_summary(
    evidence_dir: Path,
    phyml_info: dict[str, str],
    checkpoints: Sequence[dict[str, str]],
) -> None:
    checkpoint_map = {row["logical_name"]: row for row in checkpoints}
    aa_0206_sites = checkpoint_map.get("All_AA_0206.phy", {}).get("site_count", "?")
    aa_1006_sites = checkpoint_map.get("All_AA_1006.phy", {}).get("site_count", "?")
    nuc_1606_sites = checkpoint_map.get("All_NUC_1606.phy", {}).get("site_count", "?")
    lines = [
        "# Tree Reconstruction Evidence Summary",
        "",
        "## Key quantitative checkpoints",
        "",
        f"- `H2AJ.fasta`: {checkpoint_map.get('H2AJ.fasta', {}).get('record_count', '?')} records",
        f"- `cH2A.fasta`: {checkpoint_map.get('cH2A.fasta', {}).get('record_count', '?')} records",
        f"- `platypus.fasta`: {checkpoint_map.get('platypus.fasta', {}).get('record_count', '?')} records",
        f"- `All_AA_0206`: {checkpoint_map.get('All_AA_0206.phy', {}).get('record_count', '?')} records, {aa_0206_sites} sites",
        f"- `All_AA_1006`: {checkpoint_map.get('All_AA_1006.phy', {}).get('record_count', '?')} records, {aa_1006_sites} sites",
        f"- `All_NUC_1606`: {checkpoint_map.get('All_NUC_1606.phy', {}).get('record_count', '?')} records, {nuc_1606_sites} sites",
        f"- `SQK_nuc.fasta`: {checkpoint_map.get('SQK_nuc.fasta', {}).get('record_count', '?')} records",
        f"- `SQK_nuc(without short).fasta`: {checkpoint_map.get('SQK_nuc(without short).fasta', {}).get('record_count', '?')} records",
        f"- `cH2A_nuc + platypus_nuc + SQK_nuc.fasta`: {checkpoint_map.get('cH2A_nuc + platypus_nuc + SQK_nuc.fasta', {}).get('record_count', '?')} records",
        "",
        "## Current clean-profile override",
        "",
        "- The `clean` outputs intentionally drop three short H2A.J/SQK records:",
        "- `Myotis-lucifugus|XM_006084274`",
        "- `Homo-sapiens|AK303301`",
        "- `Homo-sapiens|AL133626`",
        "- July clean SQK outputs now treat `SQK_nuc(without short)` and `protein_from_SQK_nuc(without short)` as the canonical H2A.J source.",
        "",
        "## Confirmed AA web run",
        "",
        f"- Data type: `{phyml_info.get('data_type', '')}`",
        f"- Model: `{phyml_info.get('model', '')}`",
        f"- Initial tree: `{phyml_info.get('initial_tree', '')}`",
        f"- Rate model: `{phyml_info.get('rate_model', '')}` with `{phyml_info.get('rate_classes', '')}` classes",
        f"- Support mode: `{phyml_info.get('supports', '')}`",
        f"- Runtime: `{phyml_info.get('runtime', '')}`",
        "",
        "## Historical interpretation",
        "",
        "- `All_AA_0206` is the early broad-outgroup protein run and not the later mammal/nonplacental-only set.",
        "- `All_AA_1006 = 338` is the file whose count matches `H2AJ (230) + cH2A (57) + nonplacental labelled as platypus (51)`.",
        "- The nucleotide failure is upstream of PhyML: BLAST HSP fragments, truncated tails, mixed modality, and short sequences all damaged the input before tree building.",
    ]
    (evidence_dir / "evidence_summary.md").write_text("\n".join(lines), encoding="utf-8")


def run_evidence_phase(
    archive_records: dict[str, ArchiveRecord],
    output_root: Path,
    external_root: Path,
) -> None:
    evidence_dir = output_root / "evidence"
    evidence_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(external_root / "tree_archive_manifest.tsv", evidence_dir / "archive_manifest.tsv")
    shutil.copy2(archive_records["result.json"].archived_path, evidence_dir / "telegram_result.json")
    collect_notebook_evidence(archive_records, evidence_dir)
    collect_chat_evidence(archive_records, evidence_dir)
    collect_draft2_evidence(archive_records, evidence_dir)
    checkpoints = write_sequence_checkpoints(archive_records, evidence_dir)
    phyml_info = write_phyml_runbook(archive_records, evidence_dir)
    build_evidence_summary(evidence_dir, phyml_info, checkpoints)


def normalize_identifier(identifier: str) -> str:
    normalized = identifier.strip()
    normalized = normalized.replace(f"{CYRILLIC_CAPITAL_ES}H2A_nuc", "cH2A_nuc")
    normalized = normalized.replace(f"{CYRILLIC_SMALL_ES}H2A_nuc", "cH2A_nuc")
    normalized = normalized.replace("|platypus-nuc", "|nonplacental_nuc")
    normalized = normalized.replace("|platypus_nuc", "|nonplacental_nuc")
    normalized = normalized.replace("|platypus", "|nonplacental")
    normalized = normalized.replace(" ", "_")
    return normalized


def should_drop_clean_record(identifier: str, drop_prefixes: Sequence[str]) -> bool:
    clean_id = normalize_identifier(identifier)
    return any(clean_id.startswith(normalize_identifier(prefix)) for prefix in drop_prefixes)


def read_fasta_records(path: Path) -> list[SeqRecord]:
    return list(SeqIO.parse(str(path), "fasta"))


def write_fasta_records(path: Path, records: Sequence[SeqRecord]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(f">{record.id}\n{str(record.seq)}\n")


def write_relaxed_phylip(path: Path, records: Sequence[SeqRecord]) -> None:
    if not records:
        raise ValueError(f"No records available to write PHYLIP: {path}")
    lengths = {len(record.seq) for record in records}
    if len(lengths) != 1:
        raise ValueError(f"Sequences must all have the same length for PHYLIP export: {path}")
    site_count = lengths.pop()
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write(f"{len(records)} {site_count}\n")
        for record in records:
            handle.write(f"{record.id} {str(record.seq)}\n")


def read_relaxed_phylip_records(path: Path) -> list[SeqRecord]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    records: list[SeqRecord] = []
    for line in lines[1:]:
        if not line.strip():
            continue
        parts = line.split(maxsplit=1)
        if len(parts) != 2:
            continue
        identifier, sequence = parts
        records.append(SeqRecord(Seq(sequence.strip()), id=identifier.strip(), description=""))
    return records


def normalize_fasta_records(
    path: Path,
    drop_prefixes: Sequence[str] | None = None,
) -> list[SeqRecord]:
    active_drop_prefixes = list(drop_prefixes or [])
    normalized: list[SeqRecord] = []
    for record in read_fasta_records(path):
        if active_drop_prefixes and should_drop_clean_record(record.id, active_drop_prefixes):
            continue
        clean_id = normalize_identifier(record.id)
        normalized.append(SeqRecord(Seq(str(record.seq)), id=clean_id, description=""))
    return normalized


def trim_all_gap_columns(records: Sequence[SeqRecord]) -> list[SeqRecord]:
    if not records:
        return []
    lengths = {len(record.seq) for record in records}
    if len(lengths) != 1:
        return [SeqRecord(Seq(str(record.seq)), id=record.id, description="") for record in records]

    site_count = lengths.pop()
    keep_indices: list[int] = []
    for index in range(site_count):
        column = [str(record.seq)[index] for record in records]
        if any(base != "-" for base in column):
            keep_indices.append(index)

    trimmed: list[SeqRecord] = []
    for record in records:
        sequence = "".join(str(record.seq)[index] for index in keep_indices)
        trimmed.append(SeqRecord(Seq(sequence), id=record.id, description=""))
    return trimmed


def write_historical_copy(
    archive_records: dict[str, ArchiveRecord],
    output_path: Path,
    logical_name: str,
) -> None:
    source_path = archive_records[logical_name].archived_path
    copy_if_needed(source_path, output_path)


def write_clean_fasta_and_phylip(
    archive_records: dict[str, ArchiveRecord],
    logical_name: str,
    fasta_path: Path,
    phylip_path: Path | None = None,
    *,
    trim_gap_only_columns: bool = False,
    drop_prefixes: Sequence[str] | None = None,
) -> None:
    records = normalize_fasta_records(archive_records[logical_name].archived_path, drop_prefixes)
    if trim_gap_only_columns:
        records = trim_all_gap_columns(records)
    write_fasta_records(fasta_path, records)
    if phylip_path is not None:
        write_relaxed_phylip(phylip_path, records)


def write_clean_from_relaxed_phylip(
    archive_records: dict[str, ArchiveRecord],
    logical_name: str,
    fasta_path: Path,
    phylip_path: Path,
    *,
    drop_prefixes: Sequence[str] | None = None,
) -> None:
    records = read_relaxed_phylip_records(archive_records[logical_name].archived_path)
    active_drop_prefixes = list(drop_prefixes or [])
    normalized: list[SeqRecord] = []
    for record in records:
        if active_drop_prefixes and should_drop_clean_record(record.id, active_drop_prefixes):
            continue
        normalized.append(
            SeqRecord(Seq(str(record.seq)), id=normalize_identifier(record.id), description="")
        )
    write_fasta_records(fasta_path, normalized)
    write_relaxed_phylip(phylip_path, normalized)


def write_clean_note(path: Path, lines: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def filtered_merged_july_nuc_records(
    archive_records: dict[str, ArchiveRecord],
    drop_prefixes: Sequence[str],
) -> list[SeqRecord]:
    merged_records = normalize_fasta_records(
        archive_records["cH2A_nuc + platypus_nuc + SQK_nuc.fasta"].archived_path
    )
    sqk_ids = {
        record.id
        for record in normalize_fasta_records(
            archive_records["SQK_nuc(without short).fasta"].archived_path,
            drop_prefixes,
        )
    }

    filtered: list[SeqRecord] = []
    for record in merged_records:
        if should_drop_clean_record(record.id, drop_prefixes):
            continue
        if record.id.endswith("|SQK_nuc") and record.id not in sqk_ids:
            continue
        filtered.append(SeqRecord(Seq(str(record.seq)), id=record.id, description=""))
    return filtered


def run_aa_phase(
    archive_records: dict[str, ArchiveRecord],
    output_root: Path,
    profile: str,
    drop_prefixes: Sequence[str],
) -> None:
    for selected in selected_profiles(profile):
        aa_dir = output_root / selected / "aa"
        aa_dir.mkdir(parents=True, exist_ok=True)
        if selected == "historical":
            historical_names = [
                "H2AJ.fasta",
                "cH2A.fasta",
                "platypus.fasta",
                "All_AA_0206.fasta",
                "All_AA_0206.phy",
                "All_AA_1006.fasta",
                "All_AA_1006.phy",
                "all_aa_0206_1__phy_phyml.zip",
            ]
            for name in historical_names:
                write_historical_copy(archive_records, aa_dir / name, name)
        else:
            write_clean_fasta_and_phylip(
                archive_records,
                "H2AJ.fasta",
                aa_dir / "h2aj_aligned_historical.fasta",
                drop_prefixes=drop_prefixes,
            )
            write_clean_fasta_and_phylip(archive_records, "cH2A.fasta", aa_dir / "canonical_h2a_outgroup.fasta")
            write_clean_fasta_and_phylip(archive_records, "platypus.fasta", aa_dir / "nonplacental_h2a_historical_labelled_platypus.fasta")
            write_clean_fasta_and_phylip(
                archive_records,
                "All_AA_0206.fasta",
                aa_dir / "aa_0206_broad_outgroup_alignment.fasta",
                aa_dir / "aa_0206_broad_outgroup_alignment.phy",
                drop_prefixes=drop_prefixes,
            )
            write_clean_from_relaxed_phylip(
                archive_records,
                "All_AA_1006.phy",
                aa_dir / "aa_1006_mammalian_cH2A_plus_nonplacental_alignment.fasta",
                aa_dir / "aa_1006_mammalian_cH2A_plus_nonplacental_alignment.phy",
                drop_prefixes=drop_prefixes,
            )


def run_nuc_phase(
    archive_records: dict[str, ArchiveRecord],
    output_root: Path,
    profile: str,
    drop_prefixes: Sequence[str],
) -> None:
    for selected in selected_profiles(profile):
        nuc_dir = output_root / selected / "nuc"
        nuc_dir.mkdir(parents=True, exist_ok=True)
        if selected == "historical":
            historical_names = [
                "Nuc_merged_sequences.fasta",
                "All_NUC_1606.fasta",
                "All_NUC_1606.phy",
                "All_NUC_0907.fasta",
                "All_NUC_0907.phy",
                "SQK_nuc.fasta",
                "SQK_nuc(without short).fasta",
                "protein_from_SQK_nuc.fasta",
                "protein_from_SQK_nuc(without short).fasta",
                "cH2A_nuc + platypus_nuc + SQK_nuc.fasta",
            ]
            for name in historical_names:
                write_historical_copy(archive_records, nuc_dir / name, name)
        else:
            write_clean_fasta_and_phylip(
                archive_records,
                "Nuc_merged_sequences.fasta",
                nuc_dir / "nuc_1606_merged_source.fasta",
                drop_prefixes=drop_prefixes,
            )
            write_clean_from_relaxed_phylip(
                archive_records,
                "All_NUC_1606.phy",
                nuc_dir / "nuc_1606_alignment.fasta",
                nuc_dir / "nuc_1606_alignment.phy",
                drop_prefixes=drop_prefixes,
            )
            write_clean_from_relaxed_phylip(
                archive_records,
                "All_NUC_0907.phy",
                nuc_dir / "nuc_0907_alignment.fasta",
                nuc_dir / "nuc_0907_alignment.phy",
                drop_prefixes=drop_prefixes,
            )
            write_clean_fasta_and_phylip(
                archive_records,
                "SQK_nuc(without short).fasta",
                nuc_dir / "sqk_nuc_full.fasta",
                drop_prefixes=drop_prefixes,
            )
            write_clean_fasta_and_phylip(
                archive_records,
                "SQK_nuc(without short).fasta",
                nuc_dir / "sqk_nuc_without_short.fasta",
                drop_prefixes=drop_prefixes,
            )
            write_clean_fasta_and_phylip(
                archive_records,
                "protein_from_SQK_nuc(without short).fasta",
                nuc_dir / "sqk_protein_from_nuc_full.fasta",
                drop_prefixes=drop_prefixes,
            )
            write_clean_fasta_and_phylip(
                archive_records,
                "protein_from_SQK_nuc(without short).fasta",
                nuc_dir / "sqk_protein_from_nuc_without_short.fasta",
                drop_prefixes=drop_prefixes,
            )
            write_fasta_records(
                nuc_dir / "nuc_0907_ch2a_nonplacental_sqk_merged.fasta",
                filtered_merged_july_nuc_records(archive_records, drop_prefixes),
            )


def find_ch2a_terminals(tree: Any) -> list[Any]:
    terminals = []
    for terminal in tree.get_terminals():
        name = terminal.name or ""
        if "H2A_nuc" in name and "H2AJ" not in name:
            terminals.append(terminal)
    return terminals


def normalize_newick_text(text: str) -> str:
    return re.sub(r"\s+", "", text)


def root_nuc_tree_with_failure_report(
    input_path: Path,
    failure_path: Path,
    first_root_path: Path,
    human_root_path: Path,
    historical_root_path: Path | None = None,
) -> None:
    failure_path.parent.mkdir(parents=True, exist_ok=True)
    base_tree = Phylo.read(str(input_path), "newick")
    outgroup_clades = find_ch2a_terminals(base_tree)
    if not outgroup_clades:
        failure_path.write_text("No terminals containing cH2A_nuc were found in the nucleotide tree.\n", encoding="utf-8")
        return

    lines: list[str] = []
    monophyletic = base_tree.is_monophyletic(outgroup_clades)
    if monophyletic:
        lines.append("The full cH2A_nuc set is monophyletic in this tree.")
    else:
        lines.append("The full cH2A_nuc set is not monophyletic.")
        mrca = base_tree.common_ancestor(outgroup_clades)
        descendants = [clade.name or "" for clade in mrca.get_terminals()]
        foreign = [name for name in descendants if not ("H2A_nuc" in name and "H2AJ" not in name)]
        lines.append(f"cH2A_nuc terminals inspected: {len(outgroup_clades)}")
        lines.append(f"MRCA descendant count: {len(descendants)}")
        lines.append("Foreign descendants inside the MRCA subtree:")
        lines.extend(f"- {name}" for name in foreign[:80])
        if len(foreign) > 80:
            lines.append(f"- ... and {len(foreign) - 80} more")
    failure_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    first_tree = Phylo.read(str(input_path), "newick")
    first_target = find_ch2a_terminals(first_tree)[0]
    first_tree.root_with_outgroup(first_target)
    Phylo.write(first_tree, str(first_root_path), "newick")

    human_tree = Phylo.read(str(input_path), "newick")
    human_target = next(
        (
            terminal
            for terminal in human_tree.get_terminals()
            if "Homo-sapiens|NM_175065.3|" in (terminal.name or "")
            and "H2A_nuc" in (terminal.name or "")
        ),
        None,
    )
    if human_target is not None:
        human_tree.root_with_outgroup(human_target)
        Phylo.write(human_tree, str(human_root_path), "newick")
    else:
        human_root_path.write_text("Human NM_175065.3 cH2A nucleotide sequence was not found in this tree.\n", encoding="utf-8")

    if historical_root_path and historical_root_path.exists():
        archived = normalize_newick_text(historical_root_path.read_text(encoding="utf-8", errors="replace"))
        generated = normalize_newick_text(first_root_path.read_text(encoding="utf-8", errors="replace"))
        comparison = (
            "Historical rooted_tree.nwk matches the regenerated first-cH2A rooting."
            if archived == generated
            else "Historical rooted_tree.nwk does not exactly match the regenerated first-cH2A rooting."
        )
        with failure_path.open("a", encoding="utf-8") as handle:
            handle.write("\n")
            handle.write(comparison + "\n")


def run_postprocess_phase(
    archive_records: dict[str, ArchiveRecord],
    output_root: Path,
    profile: str,
    drop_prefixes: Sequence[str],
) -> None:
    for selected in selected_profiles(profile):
        post_dir = output_root / selected / "postprocess"
        post_dir.mkdir(parents=True, exist_ok=True)
        write_historical_copy(archive_records, post_dir / "nuc_tree_bootstrap.nwk", "nuc_tree_bootstrap.nwk")
        if selected == "historical":
            write_historical_copy(archive_records, post_dir / "rooted_tree.nwk", "rooted_tree.nwk")
        else:
            root_nuc_tree_with_failure_report(
                input_path=archive_records["nuc_tree_bootstrap.nwk"].archived_path,
                failure_path=post_dir / "nuc_tree_all_cH2A_outgroup_failure.txt",
                first_root_path=post_dir / "nuc_tree_first_cH2A_rooted.nwk",
                human_root_path=post_dir / "nuc_tree_human_NM_175065_rooted.nwk",
                historical_root_path=archive_records["rooted_tree.nwk"].archived_path,
            )
            write_clean_note(
                post_dir / "README_filtered_clean_note.txt",
                [
                    "These NWK files are historical tree outputs kept for reference.",
                    "The clean alignments intentionally drop the following short H2A.J/SQK records:",
                    *[f"- {normalize_identifier(prefix)}" for prefix in drop_prefixes],
                    "Re-run PhyML (or another tree inference tool) if you need a tree that matches the filtered clean alignments exactly.",
                ],
            )


def main() -> None:
    args = parse_args()
    source_paths = SourcePaths(
        histones_root=Path(args.histones_root),
        my_test_root=Path(args.my_test_root),
        grant_root=Path(args.grant_root),
        chat_export_path=Path(args.chat_export),
        draft_docx_path=Path(args.draft_docx),
    )
    external_root = Path(args.external_root)
    output_root = Path(args.output_root)
    archive_records = archive_sources(source_paths, external_root)

    if args.phase == "archive":
        print(f"[archive] wrote manifest to {external_root / 'tree_archive_manifest.tsv'}")
        return
    if args.phase in {"evidence", "all"}:
        run_evidence_phase(archive_records, output_root, external_root)
        print(f"[evidence] wrote evidence outputs under {output_root / 'evidence'}")
    if args.phase in {"aa", "all"}:
        run_aa_phase(archive_records, output_root, args.profile, args.drop_clean_id)
        print(f"[aa] wrote outputs under {output_root}")
    if args.phase in {"nuc", "all"}:
        run_nuc_phase(archive_records, output_root, args.profile, args.drop_clean_id)
        print(f"[nuc] wrote outputs under {output_root}")
    if args.phase in {"postprocess", "all"}:
        run_postprocess_phase(archive_records, output_root, args.profile, args.drop_clean_id)
        print(f"[postprocess] wrote outputs under {output_root}")


if __name__ == "__main__":
    main()
