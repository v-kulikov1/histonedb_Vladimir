# BioAnalyze Archive

Purpose
- Stores legacy outputs and audit files that are no longer active inputs.
- Keeps working folders clean while preserving historical artifacts.

Contents
- merged_legacy/
  Legacy merged CSV outputs from earlier pipeline versions.
- audits/
  Audit CSVs from earlier runs (dedup, HGNC/VGNC unresolved, accession enrichment).

How To Use
- If you need to inspect historical results, use files under this archive.
- New runs of the pipeline write audits to CURATED_SET/BioAnalyze/audits/.
