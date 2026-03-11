import csv
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import requests


ROOT = Path("CURATED_SET")
BIO_DIR = ROOT / "BioAnalyze"
MERGED_DIR = BIO_DIR / "data" / "merged"

INPUT_CSV = MERGED_DIR / "mammalia_H2A_merged.csv"
OUTPUT_CSV = MERGED_DIR / "mammalia_H2A_merged_with_taxonomy.csv"
NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def fetch_taxonomy(taxid: str, session: requests.Session) -> tuple[str, str]:
    params = {"db": "taxonomy", "id": taxid, "retmode": "xml"}
    response = session.get(NCBI_EFETCH, params=params, timeout=30)
    response.raise_for_status()

    root = ET.fromstring(response.text)
    taxon = root.find("./Taxon")
    if taxon is None:
        return "", ""

    species_name = (taxon.findtext("ScientificName") or "").strip()
    order_name = ""
    for taxon_in_lineage in taxon.findall("./LineageEx/Taxon"):
        rank = (taxon_in_lineage.findtext("Rank") or "").strip().lower()
        if rank == "order":
            order_name = (taxon_in_lineage.findtext("ScientificName") or "").strip()
            break

    return species_name, order_name


def main() -> None:
    MERGED_DIR.mkdir(parents=True, exist_ok=True)

    with INPUT_CSV.open("r", newline="", encoding="utf-8-sig") as fh:
        rows = list(csv.DictReader(fh))
        fieldnames = list(rows[0].keys()) if rows else []

    if not rows:
        raise ValueError(f"Input file has no rows: {INPUT_CSV}")

    taxid_col = "taxonomy id"
    if taxid_col not in fieldnames:
        raise ValueError(f"Missing required column: {taxid_col}")

    unique_taxids = sorted({(row.get(taxid_col) or "").strip() for row in rows if (row.get(taxid_col) or "").strip()})
    taxonomy_cache: dict[str, tuple[str, str]] = {}

    with requests.Session() as session:
        for taxid in unique_taxids:
            taxonomy_cache[taxid] = fetch_taxonomy(taxid, session)
            time.sleep(0.34)

    out_fieldnames = fieldnames + ["species_name", "order"]
    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=out_fieldnames)
        writer.writeheader()
        for row in rows:
            taxid = (row.get(taxid_col) or "").strip()
            species_name, order_name = taxonomy_cache.get(taxid, ("", ""))
            row["species_name"] = species_name
            row["order"] = order_name
            writer.writerow(row)

    print(f"Wrote: {OUTPUT_CSV}")
    print(f"Rows: {len(rows)}")
    print(f"Unique taxonomy ids: {len(unique_taxids)}")


if __name__ == "__main__":
    main()
