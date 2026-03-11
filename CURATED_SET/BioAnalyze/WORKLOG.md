# BioAnalyze H2A Merge Worklog

Last updated: 2026-03-11

Goal
- Build a single H2A dataset for mammalia + human with consistent taxonomy fields and gene identifiers.
- Make the pipeline reproducible and keep legacy outputs archived.

Canonical Output (current)
- CURATED_SET/BioAnalyze/data/merged/mammalia_H2A_merged_with_taxonomy_v3.csv

Key Inputs
- CURATED_SET/mammalia_genes/*_genes_vgnc.csv
  Source per-species VGNC tables. Used for mammalia rows.
- CURATED_SET/human_histones.csv
  Human histone source (includes H2A entries).
- CURATED_SET/histones.csv
  General histone source. Only H2A and only placental taxa are added in the pipeline.

Pipeline Summary (v3)
1) Read all mammalia VGNC files and keep only rows where Histone type == H2A.
2) Normalize fields and map to final columns:
   accession, type, variant, isoform, protein_len, taxonomy id, organism,
   vgnc_id, hgnc_id, ncbi_id, ensembl_gene_id, species_name, order
3) Append H2A rows from CURATED_SET/human_histones.csv.
4) Append H2A rows from CURATED_SET/histones.csv but only for placental taxa.
5) Accession enrichment step (v3-only):
   For rows missing all IDs (vgnc_id/hgnc_id/ncbi_id/ensembl_gene_id):
   - accession -> NCBI protein UID (esearch, db=protein)
   - protein UID -> NCBI gene id (elink, db=gene)
   - gene id -> Ensembl gene id (efetch, db=gene, Dbtag=Ensembl)
   Audit is written for resolved/unresolved accessions.
6) Resolve HGNC ID by symbol (HGNC API).
7) Fill missing Ensembl gene IDs for VGNC rows via VGNC API (symbol-report).
8) Fetch taxonomy per taxon_id from NCBI taxonomy:
   species_name and order are added.
9) Deduplicate by (accession, gene_id) with a completeness score.
10) Write output CSV and audit files.

Scripts
- CURATED_SET/BioAnalyze/scripts/merge_taxonomy/build_mammalia_h2a_merged_with_taxonomy_v2.py
  Baseline pipeline with taxonomy, HGNC/VGNC enrichment, and dedup.
- CURATED_SET/BioAnalyze/scripts/merge_taxonomy/build_mammalia_h2a_merged_with_taxonomy_v3.py
  Same as v2 plus accession enrichment (NCBI protein->gene->Ensembl).

Version Notes
- v2 output: mammalia_H2A_merged_with_taxonomy_v2.csv
  Produced by v2 script. Includes taxonomy and HGNC/VGNC enrichment.
- v3 output: mammalia_H2A_merged_with_taxonomy_v3.csv
  Produced by v3 script with accession enrichment.
  This typically increases ncbi_id and ensembl_gene_id coverage for rows
  that had no IDs in v2.

Audits (v3)
- audit_accession_enrichment_resolved.csv
- audit_accession_enrichment_unresolved.csv
- audit_dedup_dropped_rows.csv
- audit_unresolved_hgnc_symbols.csv
- audit_unresolved_vgnc_ensembl.csv

Archive Actions (2026-03-11)
Moved legacy outputs and audits into archive:
- Legacy merged CSVs moved to:
  CURATED_SET/BioAnalyze/archive/merged_legacy/
    mammalia_H2A_merged.csv
    mammalia_H2A_merged_with_taxonomy.csv
    mammalia_H2A_merged_with_taxonomy_v2.csv
- Audit CSVs moved to:
  CURATED_SET/BioAnalyze/archive/audits/
    audit_accession_enrichment_resolved.csv
    audit_accession_enrichment_unresolved.csv
    audit_dedup_dropped_rows.csv
    audit_unresolved_hgnc_symbols.csv
    audit_unresolved_vgnc_ensembl.csv

Current Snapshot (from v3 file timestamp 2026-03-03 12:09:57)
- Total rows: 464
- Filled IDs:
  vgnc_id: 328
  hgnc_id: 68
  ncbi_id: 444
  ensembl_gene_id: 329
- Taxonomy coverage:
  species_name: 464
  order: 464
- Accession enrichment audit (archived):
  resolved: 48
  unresolved: 20
- Other audits (archived):
  dedup dropped rows: 274
  unresolved hgnc symbols: 0
  unresolved vgnc -> ensembl: 64

Downstream Analyses
- Built order-level statistics based on the v3 dataset (unique genes vs unique proteins/accessions),
  using the "order" taxonomy field as the grouping key.

Expression Work (2026-03-03)
Goal
- Build human H2A expression heatmaps from Bgee advanced data using the v3 gene set.

Inputs
- H2A gene set: CURATED_SET/BioAnalyze/data/merged/mammalia_H2A_merged_with_taxonomy_v3.csv
- Bgee advanced TSV (raw, external): listed in CURATED_SET/BioAnalyze/data/raw/README.md

Process Summary
1) Join H2A human genes to Bgee by ensembl_gene_id (Gene ID).
2) Filter Bgee rows to Expression == present and Call quality == gold quality.
3) Aggregate duplicates by mean Expression score.
4) Build heatmap with log10(Expression score + 1), masking NaN.
5) Create a reduced H2A-only TSV to avoid memory issues when re-running.

Current Expression Artifacts
- CURATED_SET/BioAnalyze/data/processed/Homo_sapiens_expr_advanced_H2A_only_all_conditions.tsv
  H2A-only slice from advanced Bgee (all conditions).
- CURATED_SET/BioAnalyze/data/processed/Homo_sapiens_expr_advanced_H2A_present_gold.tsv
  H2A-only + present + gold quality subset used for plotting.
- CURATED_SET/BioAnalyze/figures/heatmaps/h2a_hs_bgee_advanced_present_gold_heatmap_ensembl_rows_gene_hgnc.svg
  Canonical human H2A heatmap (Ensembl rows, HGNC gene labels).

Legacy / Variant Expression Artifacts
- CURATED_SET/BioAnalyze/figures/heatmaps/h2a_hs_bgee_advanced_present_gold_heatmap_gene_name.png
- CURATED_SET/BioAnalyze/figures/heatmaps/h2a_hs_bgee_advanced_present_gold_heatmap_gene_name.svg
  Gene name from Bgee on the X-axis (more collapsed labels).
- CURATED_SET/BioAnalyze/figures/heatmaps/h2a_hs_bgee_advanced_present_gold_heatmap.png
- CURATED_SET/BioAnalyze/figures/heatmaps/h2a_hs_bgee_advanced_present_gold_heatmap.svg
  Earlier heatmap variant using gene labels (HGNC-based) instead of Bgee Gene name.
- CURATED_SET/BioAnalyze/figures/heatmaps/h2a_hs_bgee_advanced_present_gold_heatmap_ensembl_rows_gene_hgnc.png
- CURATED_SET/BioAnalyze/figures/heatmaps/h2a_hs_bgee_advanced_present_gold_heatmap_ensembl_rows_gene_hgnc_sorted.png
- CURATED_SET/BioAnalyze/figures/heatmaps/h2a_hs_bgee_advanced_present_gold_heatmap_ensembl_rows_gene_hgnc_sorted.svg
- CURATED_SET/BioAnalyze/figures/heatmaps/h2a_hs_bgee_advanced_present_gold_heatmap_ensembl_rows_gene_hgnc_sorted_all_xticks.png
- CURATED_SET/BioAnalyze/figures/heatmaps/h2a_hs_bgee_advanced_present_gold_heatmap_ensembl_rows_gene_hgnc_sorted_all_xticks.svg
  Iterative variants (sorting and tick-density experiments). Kept for reference.

How To Regenerate v3
1) Ensure inputs exist:
   - CURATED_SET/mammalia_genes/*_genes_vgnc.csv
   - CURATED_SET/human_histones.csv
   - CURATED_SET/histones.csv
2) Run:
   python CURATED_SET/BioAnalyze/scripts/merge_taxonomy/build_mammalia_h2a_merged_with_taxonomy_v3.py
3) Output:
   - CURATED_SET/BioAnalyze/data/merged/mammalia_H2A_merged_with_taxonomy_v3.csv
   - Audits written to CURATED_SET/BioAnalyze/audits/

Notes
- All pipeline scripts must live under CURATED_SET/BioAnalyze/scripts.
- Raw large expression files are stored outside the repo:
  see CURATED_SET/BioAnalyze/data/raw/README.md
