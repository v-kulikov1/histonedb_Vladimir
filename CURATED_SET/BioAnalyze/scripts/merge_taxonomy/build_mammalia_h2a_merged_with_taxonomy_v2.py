import csv
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Tuple
from urllib.parse import quote

import requests


ROOT = Path("CURATED_SET")
MAMMALIA_DIR = ROOT / "mammalia_genes"
BIO_DIR = ROOT / "BioAnalyze"
DATA_DIR = BIO_DIR / "data"
MERGED_DIR = DATA_DIR / "merged"
AUDIT_DIR = BIO_DIR / "audits"

HUMAN_CSV = ROOT / "human_histones.csv"
HISTONES_CSV = ROOT / "histones.csv"

OUTPUT_CSV = MERGED_DIR / "mammalia_H2A_merged_with_taxonomy_v2.csv"
AUDIT_HGNC = AUDIT_DIR / "audit_unresolved_hgnc_symbols.csv"
AUDIT_VGNC = AUDIT_DIR / "audit_unresolved_vgnc_ensembl.csv"
AUDIT_DEDUP = AUDIT_DIR / "audit_dedup_dropped_rows.csv"

NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
HGNC_FETCH_SYMBOL = "https://rest.genenames.org/fetch/symbol/"
VGNC_SYMBOL_REPORT = "https://vertebrate.genenames.org/cgi-bin/gene/symbol-report"

OUTPUT_COLUMNS = [
    "accession",
    "type",
    "variant",
    "isoform",
    "protein_len",
    "taxonomy id",
    "organism",
    "vgnc_id",
    "hgnc_id",
    "ncbi_id",
    "ensembl_gene_id",
    "species_name",
    "order",
]


def norm(value) -> str:
    if value is None:
        return ""
    return str(value).strip()


def read_csv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        return [{k: norm(v) for k, v in row.items()} for row in csv.DictReader(fh)]


def fetch_taxonomy(taxid: str, session: requests.Session, cache: Dict[str, Dict[str, str]]) -> Dict[str, str]:
    if taxid in cache:
        return cache[taxid]

    result = {"species_name": "", "order": "", "lineage": ""}
    if not taxid:
        cache[taxid] = result
        return result

    try:
        response = session.get(
            NCBI_EFETCH,
            params={"db": "taxonomy", "id": taxid, "retmode": "xml"},
            timeout=30,
        )
        response.raise_for_status()
        root = ET.fromstring(response.text)
        taxon = root.find("./Taxon")
        if taxon is not None:
            result["species_name"] = norm(taxon.findtext("ScientificName"))
            result["lineage"] = norm(taxon.findtext("Lineage"))
            for in_lineage in taxon.findall("./LineageEx/Taxon"):
                rank = norm(in_lineage.findtext("Rank")).lower()
                if rank == "order":
                    result["order"] = norm(in_lineage.findtext("ScientificName"))
                    break
    except Exception:
        pass

    cache[taxid] = result
    time.sleep(0.34)
    return result


def is_placental_taxon(tax_data: Dict[str, str]) -> bool:
    lineage = norm(tax_data.get("lineage")).lower()
    return ("eutheria" in lineage) or ("placentalia" in lineage)


def fetch_hgnc_id_for_symbol(
    symbol: str,
    session: requests.Session,
    cache: Dict[str, str],
    unresolved: List[Dict[str, str]],
) -> str:
    if symbol in cache:
        return cache[symbol]

    hgnc_id = ""
    try:
        url = HGNC_FETCH_SYMBOL + quote(symbol, safe="")
        response = session.get(url, headers={"Accept": "application/json"}, timeout=30)
        if response.status_code == 200:
            payload = response.json()
            docs = payload.get("response", {}).get("docs", []) or []
            # Prefer exact symbol match when multiple docs are returned.
            exact = [d for d in docs if norm(d.get("symbol")) == symbol]
            selected = exact[0] if exact else (docs[0] if docs else None)
            if selected:
                hgnc_id = norm(selected.get("hgnc_id"))
            if not hgnc_id:
                unresolved.append(
                    {
                        "hgnc_symbol": symbol,
                        "reason": "no_hgnc_id_in_response",
                        "status_code": str(response.status_code),
                    }
                )
        else:
            unresolved.append(
                {
                    "hgnc_symbol": symbol,
                    "reason": "http_error",
                    "status_code": str(response.status_code),
                }
            )
    except Exception as exc:
        unresolved.append(
            {
                "hgnc_symbol": symbol,
                "reason": "request_exception",
                "status_code": str(exc),
            }
        )

    cache[symbol] = hgnc_id
    time.sleep(0.1)
    return hgnc_id


def extract_ensembl_from_vgnc_json(payload: Dict) -> List[str]:
    out = []
    ext = payload.get("external_resources", {}) or {}
    genes = ext.get("genes", []) or []
    for group in genes:
        if not isinstance(group, list):
            continue
        for item in group:
            if not isinstance(item, dict):
                continue
            gid = norm(item.get("id"))
            if gid.startswith("ENS"):
                out.append(gid)
    # preserve order while dropping duplicates
    seen = set()
    uniq = []
    for gid in out:
        if gid not in seen:
            seen.add(gid)
            uniq.append(gid)
    return uniq


def fetch_ensembl_for_vgnc(
    vgnc_id: str,
    session: requests.Session,
    cache: Dict[str, str],
    unresolved: List[Dict[str, str]],
) -> str:
    if vgnc_id in cache:
        return cache[vgnc_id]

    ensembl_id = ""
    try:
        response = session.get(VGNC_SYMBOL_REPORT, params={"id": vgnc_id}, timeout=30)
        if response.status_code == 200:
            payload = response.json()
            candidates = extract_ensembl_from_vgnc_json(payload)
            if candidates:
                ensembl_id = candidates[0]
            else:
                unresolved.append(
                    {
                        "vgnc_id": vgnc_id,
                        "reason": "no_ensembl_in_vgnc_response",
                        "status_code": str(response.status_code),
                    }
                )
        else:
            unresolved.append(
                {
                    "vgnc_id": vgnc_id,
                    "reason": "http_error",
                    "status_code": str(response.status_code),
                }
            )
    except Exception as exc:
        unresolved.append(
            {
                "vgnc_id": vgnc_id,
                "reason": "request_exception",
                "status_code": str(exc),
            }
        )

    cache[vgnc_id] = ensembl_id
    time.sleep(0.2)
    return ensembl_id


def row_completeness_score(row: Dict[str, str]) -> int:
    return sum(1 for key in OUTPUT_COLUMNS if norm(row.get(key)))


def output_signature(row: Dict[str, str]) -> Tuple[str, ...]:
    return tuple(norm(row.get(col)) for col in OUTPUT_COLUMNS)


def dedup_rows(rows: List[Dict[str, str]]) -> Tuple[List[Dict[str, str]], List[Dict[str, str]]]:
    dropped: List[Dict[str, str]] = []

    # 1) Drop exact duplicates by final output signature.
    unique_exact: List[Dict[str, str]] = []
    seen_exact = set()
    for row in rows:
        sig = output_signature(row)
        if sig in seen_exact:
            dropped.append(
                {
                    "reason": "exact_duplicate",
                    **{k: norm(row.get(k)) for k in OUTPUT_COLUMNS},
                }
            )
        else:
            seen_exact.add(sig)
            unique_exact.append(row)

    # 2) Drop empty accession.
    no_empty: List[Dict[str, str]] = []
    for row in unique_exact:
        if not norm(row.get("accession")):
            dropped.append(
                {
                    "reason": "empty_accession",
                    **{k: norm(row.get(k)) for k in OUTPUT_COLUMNS},
                }
            )
        else:
            no_empty.append(row)

    # 3) Dedup by (accession, gene_id) where gene_id=vgnc_id|hgnc_id|""
    by_key: Dict[Tuple[str, str], Dict] = {}
    order: List[Tuple[str, str]] = []

    for row in no_empty:
        accession = norm(row.get("accession"))
        gene_id = norm(row.get("vgnc_id")) or norm(row.get("hgnc_id")) or ""
        key = (accession, gene_id)
        score = row_completeness_score(row)

        if key not in by_key:
            by_key[key] = {"row": row, "score": score}
            order.append(key)
            continue

        current = by_key[key]
        if score > current["score"]:
            dropped.append(
                {
                    "reason": "dedup_key_replaced_by_more_complete",
                    **{k: norm(current["row"].get(k)) for k in OUTPUT_COLUMNS},
                }
            )
            by_key[key] = {"row": row, "score": score}
        else:
            dropped.append(
                {
                    "reason": "dedup_key_dropped",
                    **{k: norm(row.get(k)) for k in OUTPUT_COLUMNS},
                }
            )

    final_rows = [by_key[key]["row"] for key in order]
    return final_rows, dropped


def write_csv(path: Path, rows: List[Dict[str, str]], fieldnames: List[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: norm(row.get(k)) for k in fieldnames})


def main() -> None:
    MERGED_DIR.mkdir(parents=True, exist_ok=True)
    AUDIT_DIR.mkdir(parents=True, exist_ok=True)

    mammalia_files = sorted(MAMMALIA_DIR.glob("*_genes_vgnc.csv"))
    base_rows: List[Dict[str, str]] = []

    # Mammalia VGNC files
    for path in mammalia_files:
        for src in read_csv(path):
            if norm(src.get("Histone type")) != "H2A":
                continue
            protein_len = norm(src.get("sequence_length"))
            if not protein_len and norm(src.get("sequence")):
                protein_len = str(len(norm(src.get("sequence"))))

            base_rows.append(
                {
                    "accession": norm(src.get("accession")),
                    "type": "H2A",
                    "variant": norm(src.get("Histone variant")),
                    "isoform": norm(src.get("Clustered (canonical) isoform")),
                    "protein_len": protein_len,
                    "taxonomy id": norm(src.get("taxon_id")),
                    "organism": norm(src.get("Species")),
                    "vgnc_id": norm(src.get("vgnc_id")),
                    "hgnc_id": "",
                    "ncbi_id": norm(src.get("ncbi_id")),
                    "ensembl_gene_id": norm(src.get("ensembl_gene_id")),
                    "species_name": "",
                    "order": "",
                    "_hgnc_symbol": "",
                    "_source": path.name,
                }
            )

    # Human file
    for src in read_csv(HUMAN_CSV):
        if norm(src.get("type")) != "H2A":
            continue
        base_rows.append(
            {
                "accession": norm(src.get("accession")),
                "type": "H2A",
                "variant": norm(src.get("Histone variant")) or norm(src.get("variant")),
                "isoform": norm(src.get("Clustered (canonical) isoform")),
                "protein_len": norm(src.get("Protein length")),
                "taxonomy id": norm(src.get("taxonomy_id")),
                "organism": norm(src.get("organism")),
                "vgnc_id": "",
                "hgnc_id": "",
                "ncbi_id": norm(src.get("ncbi_gene_id")),
                "ensembl_gene_id": norm(src.get("Ensembl gene ID")),
                "species_name": "",
                "order": "",
                "_hgnc_symbol": norm(src.get("hgnc_gene_name")),
                "_source": HUMAN_CSV.name,
            }
        )

    taxonomy_cache: Dict[str, Dict[str, str]] = {}
    unresolved_hgnc: List[Dict[str, str]] = []
    unresolved_vgnc: List[Dict[str, str]] = []

    with requests.Session() as session:
        # histones.csv: only placental H2A
        for src in read_csv(HISTONES_CSV):
            if norm(src.get("type")) != "H2A":
                continue
            taxid = norm(src.get("taxonomy_id"))
            if not taxid:
                continue
            tax_data = fetch_taxonomy(taxid, session, taxonomy_cache)
            if not is_placental_taxon(tax_data):
                continue

            seq = norm(src.get("sequence"))
            protein_len = str(len(seq)) if seq else ""
            base_rows.append(
                {
                    "accession": norm(src.get("accession")),
                    "type": "H2A",
                    "variant": norm(src.get("variant")),
                    "isoform": "",
                    "protein_len": protein_len,
                    "taxonomy id": taxid,
                    "organism": norm(src.get("organism")),
                    "vgnc_id": "",
                    "hgnc_id": "",
                    "ncbi_id": norm(src.get("ncbi_gene_id")),
                    "ensembl_gene_id": "",
                    "species_name": "",
                    "order": "",
                    "_hgnc_symbol": norm(src.get("hgnc_gene_name")),
                    "_source": HISTONES_CSV.name,
                }
            )

        # Step 4: HGNC ID resolution by symbol
        hgnc_cache: Dict[str, str] = {}
        symbols = sorted({norm(r.get("_hgnc_symbol")) for r in base_rows if norm(r.get("_hgnc_symbol"))})
        for sym in symbols:
            fetch_hgnc_id_for_symbol(sym, session, hgnc_cache, unresolved_hgnc)
        for row in base_rows:
            sym = norm(row.get("_hgnc_symbol"))
            if sym and not norm(row.get("hgnc_id")):
                row["hgnc_id"] = hgnc_cache.get(sym, "")

        # Step 5: Ensembl enrichment via VGNC
        vgnc_to_ens_local: Dict[str, str] = {}
        for row in base_rows:
            vgnc_id = norm(row.get("vgnc_id"))
            ens = norm(row.get("ensembl_gene_id"))
            if vgnc_id and ens and vgnc_id not in vgnc_to_ens_local:
                vgnc_to_ens_local[vgnc_id] = ens
        for row in base_rows:
            vgnc_id = norm(row.get("vgnc_id"))
            if vgnc_id and not norm(row.get("ensembl_gene_id")) and vgnc_id in vgnc_to_ens_local:
                row["ensembl_gene_id"] = vgnc_to_ens_local[vgnc_id]

        vgnc_cache: Dict[str, str] = {}
        missing_vgnc_ids = sorted(
            {norm(r.get("vgnc_id")) for r in base_rows if norm(r.get("vgnc_id")) and not norm(r.get("ensembl_gene_id"))}
        )
        for vgnc_id in missing_vgnc_ids:
            fetch_ensembl_for_vgnc(vgnc_id, session, vgnc_cache, unresolved_vgnc)
        for row in base_rows:
            vgnc_id = norm(row.get("vgnc_id"))
            if vgnc_id and not norm(row.get("ensembl_gene_id")):
                row["ensembl_gene_id"] = vgnc_cache.get(vgnc_id, "")

        # Step 6: taxonomy fields for all rows
        unique_taxids = sorted({norm(r.get("taxonomy id")) for r in base_rows if norm(r.get("taxonomy id"))})
        for taxid in unique_taxids:
            fetch_taxonomy(taxid, session, taxonomy_cache)
        for row in base_rows:
            taxid = norm(row.get("taxonomy id"))
            tax_data = taxonomy_cache.get(taxid, {})
            row["species_name"] = norm(tax_data.get("species_name"))
            row["order"] = norm(tax_data.get("order"))

    # Step 7: dedup
    deduped_rows, dedup_dropped = dedup_rows(base_rows)

    # Remove helper fields before writing.
    final_rows = [{k: norm(v.get(k)) for k in OUTPUT_COLUMNS} for v in deduped_rows]

    # Step 8: write outputs
    write_csv(OUTPUT_CSV, final_rows, OUTPUT_COLUMNS)
    write_csv(AUDIT_HGNC, unresolved_hgnc, ["hgnc_symbol", "reason", "status_code"])
    write_csv(AUDIT_VGNC, unresolved_vgnc, ["vgnc_id", "reason", "status_code"])
    write_csv(AUDIT_DEDUP, dedup_dropped, ["reason"] + OUTPUT_COLUMNS)

    print(f"WROTE {OUTPUT_CSV}")
    print(f"WROTE {AUDIT_HGNC}")
    print(f"WROTE {AUDIT_VGNC}")
    print(f"WROTE {AUDIT_DEDUP}")
    print(f"INPUT_ROWS={len(base_rows)}")
    print(f"FINAL_ROWS={len(final_rows)}")
    print(f"DEDUP_DROPPED_ROWS={len(dedup_dropped)}")
    print(f"UNRESOLVED_HGNC={len(unresolved_hgnc)}")
    print(f"UNRESOLVED_VGNC_ENSEMBL={len(unresolved_vgnc)}")


if __name__ == "__main__":
    main()
