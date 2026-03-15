"""
Build symbol-style gene_name values for mammalia_H2A_merged_with_taxonomy_v4.csv.

Rules:
- Homo sapiens -> human_histones.csv hgnc_gene_name (match by accession)
- Other species -> local VGNC symbol (match by accession, then Ensembl, then NCBI)
- Fallback -> NCBI Gene-ref_locus symbol
"""

from __future__ import annotations

import csv
import time
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


ROOT = Path("CURATED_SET")
BIO_DIR = ROOT / "BioAnalyze"
DATA_DIR = BIO_DIR / "data"
MERGED_DIR = DATA_DIR / "merged"
AUDIT_DIR = BIO_DIR / "audits"
MAMMALIA_DIR = ROOT / "mammalia_genes"
HUMAN_CSV = ROOT / "human_histones.csv"

INPUT_CSV = MERGED_DIR / "mammalia_H2A_merged_with_taxonomy_v3.csv"
OUTPUT_CSV = MERGED_DIR / "mammalia_H2A_merged_with_taxonomy_v4.csv"
AUDIT_CSV = AUDIT_DIR / "audit_gene_name_unresolved_v4.csv"

NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
NCBI_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
NCBI_ELINK = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"

REQUEST_SLEEP_SECONDS = 0.15
TOP_N_REASONS = 8


def norm(value) -> str:
    if value is None:
        return ""
    return str(value).strip()


def normalize_numeric_id(value) -> str:
    text = norm(value)
    if text.endswith(".0") and text[:-2].isdigit():
        return text[:-2]
    return text


def read_csv(path: Path) -> Tuple[List[Dict[str, str]], List[str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh)
        rows = [{k: norm(v) for k, v in row.items()} for row in reader]
        return rows, list(reader.fieldnames or [])


def write_csv(path: Path, rows: Iterable[Dict[str, str]], columns: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({col: norm(row.get(col)) for col in columns})


def build_output_columns(input_cols: List[str]) -> List[str]:
    if "gene_name" in input_cols:
        return list(input_cols)

    out: List[str] = []
    inserted = False
    for col in input_cols:
        out.append(col)
        if col == "hgnc_id":
            out.append("gene_name")
            inserted = True
    if not inserted:
        out.append("gene_name")
    return out


def build_session() -> requests.Session:
    session = requests.Session()
    retry = Retry(
        total=3,
        connect=3,
        read=3,
        backoff_factor=0.5,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"],
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("http://", adapter)
    session.mount("https://", adapter)
    session.headers.update(
        {
            "Accept": "application/json",
            "User-Agent": "histonedb-bioanalyze/1.0 (+gene_symbol_v4)",
        }
    )
    return session


def build_human_lookup() -> Dict[str, str]:
    by_accession: Dict[str, str] = {}
    with HUMAN_CSV.open("r", encoding="utf-8-sig", newline="") as fh:
        for row in csv.DictReader(fh):
            if norm(row.get("type")) != "H2A":
                continue
            symbol = norm(row.get("hgnc_gene_name"))
            if not symbol:
                continue
            for key in (norm(row.get("accession")), norm(row.get("RefSeq peptide ID"))):
                if key and key not in by_accession:
                    by_accession[key] = symbol
    return by_accession


def build_vgnc_lookups() -> Tuple[Dict[str, Dict[str, str]], Dict[str, Dict[str, str]], Dict[str, Dict[str, str]]]:
    by_accession: Dict[str, Dict[str, str]] = {}
    by_ensembl: Dict[str, Dict[str, str]] = {}
    by_ncbi: Dict[str, Dict[str, str]] = {}

    for path in sorted(MAMMALIA_DIR.glob("*_genes_vgnc.csv")):
        with path.open("r", encoding="utf-8-sig", newline="") as fh:
            for row in csv.DictReader(fh):
                if norm(row.get("Histone type")) != "H2A":
                    continue
                symbol = norm(row.get("symbol"))
                if not symbol:
                    continue
                record = {
                    "symbol": symbol,
                    "vgnc_id": norm(row.get("vgnc_id")),
                    "ensembl_gene_id": norm(row.get("ensembl_gene_id")),
                    "ncbi_id": normalize_numeric_id(row.get("ncbi_id")),
                }
                accession = norm(row.get("accession")) or norm(row.get("refseq_id"))
                if accession and accession not in by_accession:
                    by_accession[accession] = record
                if record["ensembl_gene_id"] and record["ensembl_gene_id"] not in by_ensembl:
                    by_ensembl[record["ensembl_gene_id"]] = record
                if record["ncbi_id"] and record["ncbi_id"] not in by_ncbi:
                    by_ncbi[record["ncbi_id"]] = record

    return by_accession, by_ensembl, by_ncbi


def apply_local_vgnc_symbol(
    row: Dict[str, str],
    by_accession: Dict[str, Dict[str, str]],
    by_ensembl: Dict[str, Dict[str, str]],
    by_ncbi: Dict[str, Dict[str, str]],
) -> str:
    accession = norm(row.get("accession"))
    ensembl_gene_id = norm(row.get("ensembl_gene_id"))
    ncbi_id = normalize_numeric_id(row.get("ncbi_id"))

    match = None
    if accession:
        match = by_accession.get(accession)
    if match is None and ensembl_gene_id:
        match = by_ensembl.get(ensembl_gene_id)
    if match is None and ncbi_id:
        match = by_ncbi.get(ncbi_id)

    if match is None:
        return ""

    row["gene_name"] = match["symbol"]
    if not norm(row.get("vgnc_id")) and match.get("vgnc_id"):
        row["vgnc_id"] = match["vgnc_id"]
    if not norm(row.get("ensembl_gene_id")) and match.get("ensembl_gene_id"):
        row["ensembl_gene_id"] = match["ensembl_gene_id"]
    if not normalize_numeric_id(row.get("ncbi_id")) and match.get("ncbi_id"):
        row["ncbi_id"] = match["ncbi_id"]
    else:
        row["ncbi_id"] = normalize_numeric_id(row.get("ncbi_id"))
    return "local_vgnc_symbol"


def fetch_protein_uid(accession: str, session: requests.Session, cache: Dict[str, str]) -> str:
    if accession in cache:
        return cache[accession]

    uid = ""
    try:
        response = session.get(
            NCBI_ESEARCH,
            params={"db": "protein", "term": f"{accession}[Accession]", "retmode": "xml"},
            timeout=30,
        )
        response.raise_for_status()
        root = ET.fromstring(response.text)
        ids = [node.text for node in root.findall("./IdList/Id") if node.text]
        uid = ids[0] if ids else ""
    except Exception:
        uid = ""

    cache[accession] = uid
    time.sleep(REQUEST_SLEEP_SECONDS)
    return uid


def fetch_gene_id_from_protein_uid(uid: str, session: requests.Session, cache: Dict[str, str]) -> str:
    if uid in cache:
        return cache[uid]

    gene_id = ""
    try:
        response = session.get(
            NCBI_ELINK,
            params={"dbfrom": "protein", "db": "gene", "id": uid, "retmode": "xml"},
            timeout=30,
        )
        response.raise_for_status()
        root = ET.fromstring(response.text)
        ids = [node.text for node in root.findall(".//LinkSetDb/Link/Id") if node.text]
        gene_id = ids[0] if ids else ""
    except Exception:
        gene_id = ""

    cache[uid] = gene_id
    time.sleep(REQUEST_SLEEP_SECONDS)
    return gene_id


def fetch_gene_symbol(gene_id: str, session: requests.Session, cache: Dict[str, str]) -> str:
    if gene_id in cache:
        return cache[gene_id]

    symbol = ""
    try:
        response = session.get(
            NCBI_EFETCH,
            params={"db": "gene", "id": gene_id, "retmode": "xml"},
            timeout=30,
        )
        response.raise_for_status()
        root = ET.fromstring(response.text)
        symbol = norm(root.findtext(".//Gene-ref_locus"))
    except Exception:
        symbol = ""

    cache[gene_id] = symbol
    time.sleep(REQUEST_SLEEP_SECONDS)
    return symbol


def unresolved_row(row: Dict[str, str], source: str, reason: str) -> Dict[str, str]:
    return {
        "accession": norm(row.get("accession")),
        "species_name": norm(row.get("species_name")),
        "hgnc_id": norm(row.get("hgnc_id")),
        "vgnc_id": norm(row.get("vgnc_id")),
        "ncbi_id": normalize_numeric_id(row.get("ncbi_id")),
        "ensembl_gene_id": norm(row.get("ensembl_gene_id")),
        "source": source,
        "reason": reason,
        "gene_name": norm(row.get("gene_name")),
    }


def main() -> None:
    if not INPUT_CSV.exists():
        raise FileNotFoundError(INPUT_CSV)

    rows, input_columns = read_csv(INPUT_CSV)
    output_columns = build_output_columns(input_columns)
    human_lookup = build_human_lookup()
    by_accession, by_ensembl, by_ncbi = build_vgnc_lookups()

    total_rows = 0
    resolved = 0
    unresolved = 0
    reason_counts = Counter()
    unresolved_rows: List[Dict[str, str]] = []

    protein_cache: Dict[str, str] = {}
    gene_link_cache: Dict[str, str] = {}
    gene_symbol_cache: Dict[str, str] = {}

    with build_session() as session:
        for row in rows:
            total_rows += 1
            row["ncbi_id"] = normalize_numeric_id(row.get("ncbi_id"))
            row["gene_name"] = ""

            accession = norm(row.get("accession"))
            species = norm(row.get("species_name"))
            source = ""
            reason = ""

            if species == "Homo sapiens":
                symbol = human_lookup.get(accession, "")
                if symbol:
                    row["gene_name"] = symbol
                    source = "human_histones"
                else:
                    source = "human_histones"
                    reason = "accession_not_found"
            else:
                source = apply_local_vgnc_symbol(row, by_accession, by_ensembl, by_ncbi)
                if not source:
                    source = "local_vgnc_symbol"
                    reason = "no_local_symbol_match"

            if not norm(row.get("gene_name")):
                ncbi_id = normalize_numeric_id(row.get("ncbi_id"))
                if ncbi_id:
                    symbol = fetch_gene_symbol(ncbi_id, session, gene_symbol_cache)
                    if symbol:
                        row["gene_name"] = symbol
                        source = "ncbi_symbol"
                        reason = ""
                    else:
                        source = "ncbi_symbol"
                        reason = "no_symbol_in_gene_record"
                else:
                    if not accession:
                        source = "ncbi_accession"
                        reason = "missing_accession"
                    else:
                        protein_uid = fetch_protein_uid(accession, session, protein_cache)
                        if not protein_uid:
                            source = "ncbi_accession"
                            reason = "protein_not_found"
                        else:
                            gene_id = fetch_gene_id_from_protein_uid(protein_uid, session, gene_link_cache)
                            if not gene_id:
                                source = "ncbi_accession"
                                reason = "gene_link_not_found"
                            else:
                                row["ncbi_id"] = gene_id
                                symbol = fetch_gene_symbol(gene_id, session, gene_symbol_cache)
                                if symbol:
                                    row["gene_name"] = symbol
                                    source = "ncbi_symbol"
                                    reason = ""
                                else:
                                    source = "ncbi_symbol"
                                    reason = "no_symbol_in_gene_record"

            if norm(row.get("gene_name")):
                resolved += 1
            else:
                unresolved += 1
                reason_counts[reason or "unknown"] += 1
                unresolved_rows.append(unresolved_row(row, source, reason or "unknown"))

    write_csv(OUTPUT_CSV, rows, output_columns)
    write_csv(
        AUDIT_CSV,
        unresolved_rows,
        [
            "accession",
            "species_name",
            "hgnc_id",
            "vgnc_id",
            "ncbi_id",
            "ensembl_gene_id",
            "source",
            "reason",
            "gene_name",
        ],
    )

    print(f"TOTAL_ROWS={total_rows}")
    print(f"RESOLVED={resolved}")
    print(f"UNRESOLVED={unresolved}")
    print(f"TOP_{TOP_N_REASONS}_UNRESOLVED_REASONS")
    for reason, count in reason_counts.most_common(TOP_N_REASONS):
        print(f"UNRESOLVED_REASON={reason} COUNT={count}")
    print(f"OUTPUT_CSV={OUTPUT_CSV}")
    print(f"AUDIT_CSV={AUDIT_CSV}")


if __name__ == "__main__":
    main()
