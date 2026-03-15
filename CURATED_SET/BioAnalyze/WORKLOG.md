# BioAnalyze H2A Merge Worklog

Last updated: 2026-03-15

Goal
- Build a single H2A dataset for mammalia + human with consistent taxonomy fields, gene identifiers, and symbol-style gene names.
- Make the pipeline reproducible and keep outputs organized.

Canonical Output (current)
- CURATED_SET/BioAnalyze/data/merged/mammalia_H2A_merged_with_taxonomy_v4.csv

Key Inputs
- CURATED_SET/mammalia_genes/*_genes_vgnc.csv
  Source per-species VGNC tables. Used for mammalia rows.
- CURATED_SET/human_histones.csv
  Human histone source (includes H2A entries and HGNC-style symbols).
- CURATED_SET/histones.csv
  General histone source. Only H2A and only placental taxa are added in the pipeline.

Pipeline Summary (v4)
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
10) Write v3 CSV and audit files.
11) Post-process v3 -> v4:
   - Homo sapiens rows: fill gene_name from human_histones.csv hgnc_gene_name by accession
   - Non-human rows: fill gene_name from local VGNC symbol by accession, then Ensembl, then NCBI
   - Remaining unresolved rows: use NCBI Gene-ref_locus symbol, optionally recovering ncbi_id via accession
   - Write unresolved rows to a dedicated audit CSV

Scripts
- CURATED_SET/BioAnalyze/scripts/merge_taxonomy/build_mammalia_h2a_merged_with_taxonomy_v2.py
  Baseline pipeline with taxonomy, HGNC/VGNC enrichment, and dedup.
- CURATED_SET/BioAnalyze/scripts/merge_taxonomy/build_mammalia_h2a_merged_with_taxonomy_v3.py
  Same as v2 plus accession enrichment (NCBI protein->gene->Ensembl).
- CURATED_SET/BioAnalyze/scripts/merge_taxonomy/add_gene_symbol_as_gene_name_v4.py
  Post-processes v3 into symbol-based v4 by filling gene_name with short symbols.

Version Notes
- v2 output: mammalia_H2A_merged_with_taxonomy_v2.csv
  Produced by v2 script. Includes taxonomy and HGNC/VGNC enrichment.
- v3 output: mammalia_H2A_merged_with_taxonomy_v3.csv
  Produced by v3 script with accession enrichment.
  This typically increases ncbi_id and ensembl_gene_id coverage for rows
  that had no IDs in v2.
- v4 output: mammalia_H2A_merged_with_taxonomy_v4.csv
  Produced by the symbol post-process script using v3 as input.
  Adds gene_name immediately after hgnc_id while preserving the rest of the
  v3 column order and row count.

Audits
- audit_accession_enrichment_resolved.csv
- audit_accession_enrichment_unresolved.csv
- audit_dedup_dropped_rows.csv
- audit_unresolved_hgnc_symbols.csv
- audit_unresolved_vgnc_ensembl.csv
- audit_gene_name_unresolved_v4.csv

Current Snapshot (from v4 run on 2026-03-15)
- Total rows: 464
- v4 change:
  `gene_name` was added/backfilled in the merged output and is now the canonical
  display name source for downstream figures.
- Filled IDs:
  vgnc_id: 344
  hgnc_id: 68
  ncbi_id: 444
  ensembl_gene_id: 331
- Filled gene_name: 444
- Taxonomy coverage:
  species_name: 464
  order: 464
- gene_name unresolved summary:
  unresolved: 20
  protein_not_found: 15
  gene_link_not_found: 5

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
- CURATED_SET/BioAnalyze/data/processed/homo_sapiens/Homo_sapiens_expr_advanced_H2A_only_all_conditions.tsv
  H2A-only slice from advanced Bgee (all conditions).
- CURATED_SET/BioAnalyze/data/processed/homo_sapiens/Homo_sapiens_expr_advanced_H2A_present_gold.tsv
  H2A-only + present + gold quality subset used for plotting.

Expression Split (2026-03-11)
- Rule: canonical == "clustered H2A" in v3 `variant` column (Homo sapiens rows).
- Script:
  - CURATED_SET/BioAnalyze/scripts/expression/build_bgee_h2a_heatmaps_split_clustered_variants.py
- Outputs:
  - CURATED_SET/BioAnalyze/data/processed/homo_sapiens/h2a_hs_canonical_variant_map.tsv
  - CURATED_SET/BioAnalyze/figures/heatmaps/human/h2a_clustered.svg
  - CURATED_SET/BioAnalyze/figures/heatmaps/human/h2a_clustered.png
  - CURATED_SET/BioAnalyze/figures/heatmaps/human/h2a_variants.svg
  - CURATED_SET/BioAnalyze/figures/heatmaps/human/h2a_variants.png

Multi-Species Heatmaps (2026-03-11)
- New script for any species (separate file to keep prior work intact):
  - CURATED_SET/BioAnalyze/scripts/expression/build_bgee_h2a_heatmaps_any_species.py
- Inputs:
  - Bgee advanced TSV for target species (all conditions)
  - CURATED_SET/BioAnalyze/data/merged/mammalia_H2A_merged_with_taxonomy_v3.csv
- Processing:
  - Filter to H2A Ensembl IDs for the species
  - Expression == present and Call quality == gold quality
  - Build labels as GeneName:HGNC for Homo sapiens, otherwise GeneName:VGNC
    (fallback GeneName:ENSG if ID missing)
  - canonical == variant "clustered H2A", variants == everything else
- Outputs (species slugged, in figures/heatmaps and data/processed/<species_slug>):
  - h2a_{species}_all.(svg|png)
  - h2a_{species}_clustered.(svg|png)
  - h2a_{species}_variants.(svg|png)
  - {species}_expr_advanced_H2A_present_gold.tsv
  - {species}_h2a_canonical_variant_map.tsv

Pan troglodytes Heatmaps (2026-03-11)
- Inputs:
  - External Bgee advanced TSV:
    C:\Users\USER\Documents\GitHub\histonedb_external_storage\BioAnalyze\raw\Pan_troglodytes_expr_advanced_all_conditions.tsv
  - CURATED_SET/BioAnalyze/data/merged/mammalia_H2A_merged_with_taxonomy_v3.csv
- Processing:
  - Filtered to H2A + present + gold.
  - Labels use GeneName:VGNC (auto ID selection for non-human species).
- Outputs:
  - CURATED_SET/BioAnalyze/data/processed/pan_troglodytes/pan_troglodytes_expr_advanced_H2A_present_gold.tsv
  - CURATED_SET/BioAnalyze/data/processed/pan_troglodytes/pan_troglodytes_h2a_canonical_variant_map.tsv
  - CURATED_SET/BioAnalyze/figures/heatmaps/pan_troglodytes/h2a_pan_troglodytes_all.(svg|png)
  - CURATED_SET/BioAnalyze/figures/heatmaps/pan_troglodytes/h2a_pan_troglodytes_clustered.(svg|png)
  - CURATED_SET/BioAnalyze/figures/heatmaps/pan_troglodytes/h2a_pan_troglodytes_variants.(svg|png)
- Summary stats (from script run):
  - Rows after filter: 746
  - Heatmap matrix sizes: all=14x26, clustered=7x26, variants=7x26
2026-03-12 update:
- Rebuilt Pan troglodytes heatmaps with square-cell sizing for easier comparison:
  - Used build_bgee_h2a_heatmaps_any_species.py with --square-cells and cell-size/min-size args.

Bos taurus Heatmaps (2026-03-12)
- Inputs:
  - External Bgee advanced TSV:
    C:\Users\USER\Documents\GitHub\histonedb_external_storage\BioAnalyze\raw\Bos_taurus_expr_advanced_all_conditions.tsv
  - CURATED_SET/BioAnalyze/data/merged/mammalia_H2A_merged_with_taxonomy_v3.csv
- Processing:
  - Filtered to H2A + present + gold.
  - Labels use GeneName:VGNC (auto ID selection for non-human species).
  - Square-cell sizing enabled.
- Outputs:
  - CURATED_SET/BioAnalyze/data/processed/bos_taurus/bos_taurus_expr_advanced_H2A_present_gold.tsv
  - CURATED_SET/BioAnalyze/data/processed/bos_taurus/bos_taurus_h2a_canonical_variant_map.tsv
  - CURATED_SET/BioAnalyze/figures/heatmaps/bos_taurus/h2a_bos_taurus_all.(svg|png)
  - CURATED_SET/BioAnalyze/figures/heatmaps/bos_taurus/h2a_bos_taurus_clustered.(svg|png)
  - CURATED_SET/BioAnalyze/figures/heatmaps/bos_taurus/h2a_bos_taurus_variants.(svg|png)
- Summary stats (from script run):
  - Rows after filter: 4644
  - Heatmap matrix sizes: all=23x127, clustered=16x127, variants=7x127

Human -> Chimp Intersection Heatmaps (2026-03-12)
- Goal: build heatmaps for the dataset with more rows, keeping only exact intersections between human and chimp:
  - X axis: Anatomical entity name (exact match).
  - Y axis: Gene name (exact match).
- Rule: no normalization applied; `multicellular organism` is kept if it exists in both.
- Primary dataset is selected automatically by larger present+gold row count.
- Inputs:
  - CURATED_SET/BioAnalyze/data/processed/homo_sapiens/Homo_sapiens_expr_advanced_H2A_present_gold.tsv
  - CURATED_SET/BioAnalyze/data/processed/pan_troglodytes/pan_troglodytes_expr_advanced_H2A_present_gold.tsv
- Intersection heatmaps script (universal):
  - CURATED_SET/BioAnalyze/scripts/expression/build_hs_pan_troglodytes_aligned_heatmaps.py
  - Output TSV:
    - CURATED_SET/BioAnalyze/data/processed/intersections/homo_sapiens_expr_advanced_H2A_present_gold_intersection.tsv
  - Heatmaps:
    - CURATED_SET/BioAnalyze/figures/heatmaps/alligned_human_pan/hs_aligned_all.(svg|png)
    - CURATED_SET/BioAnalyze/figures/heatmaps/alligned_human_pan/hs_aligned_clustered.(svg|png)
    - CURATED_SET/BioAnalyze/figures/heatmaps/alligned_human_pan/hs_aligned_variants.(svg|png)

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

How To Regenerate v4
1) Ensure v3 exists:
   - CURATED_SET/BioAnalyze/data/merged/mammalia_H2A_merged_with_taxonomy_v3.csv
2) Run:
   python CURATED_SET/BioAnalyze/scripts/merge_taxonomy/add_gene_symbol_as_gene_name_v4.py
3) Output:
   - CURATED_SET/BioAnalyze/data/merged/mammalia_H2A_merged_with_taxonomy_v4.csv
   - CURATED_SET/BioAnalyze/audits/audit_gene_name_unresolved_v4.csv

Notes
- All pipeline scripts must live under CURATED_SET/BioAnalyze/scripts.
- Raw large expression files are stored outside the repo:
  see CURATED_SET/BioAnalyze/data/raw/README.md

Heatmap Label Refresh (2026-03-15)
- Goal:
  - Rebuild all H2A heatmaps so displayed `gene_name` comes from
    `CURATED_SET/BioAnalyze/data/merged/mammalia_H2A_merged_with_taxonomy_v4.csv`
    instead of Bgee `Gene name`.
- Source-data update:
  - `gene_name` was explicitly added/updated in merged `v4`, and that refreshed
    symbol set is now treated as the source of truth for displayed labels.
- Scripts updated:
  - CURATED_SET/BioAnalyze/scripts/expression/build_bgee_h2a_heatmaps_split_clustered_variants.py
  - CURATED_SET/BioAnalyze/scripts/expression/build_bgee_h2a_heatmaps_any_species.py
  - CURATED_SET/BioAnalyze/scripts/expression/build_hs_pan_troglodytes_aligned_heatmaps.py
- Behavior changes:
  - Human labels remain `GeneName:HGNC`.
  - Non-human labels remain `GeneName:VGNC`.
  - If the chosen ID is missing, labels fall back to `GeneName:<Ensembl gene ID>`.
  - Canonical map TSVs now write `gene_name` from merged v4.
  - Human/pan aligned heatmaps now intersect genes by merged v4 `gene_name`;
    the aligned processed TSV keeps raw Bgee columns and adds `merged_gene_name`.
- Rebuilt outputs:
  - Human heatmaps and map TSV were regenerated from `v4` names.
  - `pan_troglodytes` heatmaps and map TSV were regenerated from `v4` names.
  - `bos_taurus` heatmaps and map TSV were regenerated from `v4` names.
  - Human/pan aligned heatmaps and the aligned TSV were regenerated from `v4`
    names.
- Verified examples after rebuild:
  - Human labels now include `H2AC25:HGNC:20507`, `H2AL3:HGNC:53960`,
    `H2AL1Q:HGNC:53959`.
  - Pan labels now include `H2AC25:VGNC:83685`.
  - Bos fallback still works for missing VGNC, e.g.
    `LOC528006:ENSBTAG00070035097`.
