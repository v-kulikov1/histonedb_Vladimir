# BioAnalyze H2A Merge Worklog

Last updated: 2026-03-17

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

Cross-Species Gene*Tissue Comparison (2026-03-17)
- Goal:
  - Move from manual heatmap inspection to a reproducible comparison layer that
    ranks strong cross-species expression differences at the `gene × tissue`
    level.
- Shared helper layer:
  - Added `CURATED_SET/BioAnalyze/scripts/expression/gene_compare_common.py`
    to centralize:
    - shared-gene index loading
    - present+gold loading
    - per-gene long-table construction
    - matrix assembly
    - canonical output layout constants
- Cross-species comparison scripts:
  - `CURATED_SET/BioAnalyze/scripts/expression/build_gene_compare_heatmap.py`
    refactored to use the shared helper without changing CLI behavior.
  - `CURATED_SET/BioAnalyze/scripts/expression/build_gene_compare_heatmaps_batch.py`
    added to batch-build per-gene cross-species heatmaps.
  - `CURATED_SET/BioAnalyze/scripts/expression/rank_cross_species_h2a_differences.py`
    added to compute:
    - `gene_tissue_summary.csv`
    - `conservative_candidates.csv` (global p90)
    - `high_confidence_candidates.csv` (global p95)
    - `pairwise_contrasts.csv`
    - `manuscript_shortlist.md`
  - `CURATED_SET/BioAnalyze/scripts/expression/plot_cross_species_candidate_panels.py`
    added to build:
    - focused candidate panels
    - candidate overview scatter
    - class-specific p95 pages for `clustered` and `variant`
- Comparison rules:
  - exact string match on `Anatomical entity name`
  - aggregation by mean expression score per `(gene, species, tissue)`
  - conservative candidate filter uses `species_n >= 4`
  - generic tissues removed:
    - `material anatomical entity`
    - `anatomical system`
    - `multicellular organism`
- Observed summary:
  - Global p95 candidates are dominated by clustered genes.
  - Variant genes do not enter global p95, but become informative when ranked
    within the variant class.
  - `H2AZ1` and `H2AZ2` remain among the most conserved variant genes by median
    cross-species range.

Stats Layout Reorganization (2026-03-17)
- Reorganized `CURATED_SET/BioAnalyze/stats/` into subfolders:
  - `shared_genes/`
  - `ranking/tables/`
  - `ranking/reports/`
  - `ranking/plots/`
  - `accession_stats/`
- Updated output defaults in dependent scripts so new runs write into this
  layout automatically.

Heatmap Layout Reorganization (2026-03-17)
- Problem:
  - `gene_compare/` and per-species heatmap folders were previously mixed at
    the same level under `figures/heatmaps/`.
- Change:
  - moved all per-species heatmaps into:
    `CURATED_SET/BioAnalyze/figures/heatmaps/species/`
  - kept:
    - `CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/`
    - `CURATED_SET/BioAnalyze/figures/heatmaps/alligned_human_pan/`
    at the top level as separate comparison outputs
- Updated defaults in:
  - `build_bgee_h2a_heatmaps_any_species.py`
  - `build_bgee_h2a_heatmaps_remaining_species_batch.py`
  - `build_bgee_h2a_heatmaps_split_clustered_variants.py`
  - shared `DEFAULT_HEATMAP_DIR` in `gene_compare_common.py`

Remaining Species Batch Heatmaps (2026-03-16)
- Goal:
  - Extend the any-species H2A expression pipeline so the remaining species can
    be processed reproducibly from merged `v4`, while allowing partial outputs
    for species that do not have both split classes after the Bgee join.
- Script updates:
  - Refactored:
    `CURATED_SET/BioAnalyze/scripts/expression/build_bgee_h2a_heatmaps_any_species.py`
  - New batch runner:
    `CURATED_SET/BioAnalyze/scripts/expression/build_bgee_h2a_heatmaps_remaining_species_batch.py`
- Behavior changes in the refactored any-species script:
  - Added `--canonical-rule` with:
    - `legacy`: canonical only when `variant == clustered H2A`
    - `canonical_like`: canonical when `variant == clustered H2A` or starts
      with `cH2A`
  - Added `--allow-partial-splits`:
    - writes `all` plus whichever split is available
    - returns `skipped` for species with no Ensembl IDs in merged `v4`
    - returns `skipped` for species with zero `present+gold` rows after join
  - The core logic now returns structured per-species status metadata, which is
    reused by the batch runner.
- Batch inputs:
  - merged labels:
    `CURATED_SET/BioAnalyze/data/merged/mammalia_H2A_merged_with_taxonomy_v4.csv`
  - raw Bgee TSVs:
    `C:\Users\USER\Documents\GitHub\histonedb_external_storage\BioAnalyze\raw\*_expr_advanced_all_conditions.tsv`
  - batch audit:
    `CURATED_SET/BioAnalyze/audits/audit_h2a_remaining_species_batch_v4.tsv`
- Batch run settings:
  - canonical rule: `canonical_like`
  - partial outputs allowed
  - square-cell sizing enabled
- Batch outcome summary:
  - `success`: 7 species
    - `Canis lupus familiaris`: rows_after_filter=1509, all=23x58,
      clustered=15x58, variants=8x58
    - `Cavia porcellus`: rows_after_filter=29, all=2x15, clustered=1x15,
      variants=1x1
    - `Equus caballus`: rows_after_filter=950, all=21x26, clustered=14x26,
      variants=7x26
    - `Felis catus`: rows_after_filter=264, all=23x12, clustered=13x12,
      variants=10x12
    - `Macaca mulatta`: rows_after_filter=839, all=20x28, clustered=12x28,
      variants=8x28
    - `Mus musculus`: rows_after_filter=2263, all=9x323, clustered=2x65,
      variants=7x323
    - `Sus scrofa`: rows_after_filter=1242, all=12x54, clustered=6x54,
      variants=6x54
  - `partial`: 2 species
    - `Heterocephalus glaber`: reason=`missing_variants_split`,
      rows_after_filter=2, outputs=`all` + `clustered`
    - `Oryctolagus cuniculus`: reason=`missing_clustered_split`,
      rows_after_filter=13, outputs=`all` + `variants`
  - `skipped`: 2 species
    - `Callithrix jacchus`: reason=`no_ensembl_ids_in_merged_v4`
    - `Rattus norvegicus`: reason=`no_present_gold_rows_after_join`
- Generated processed outputs:
  - New species folders under `CURATED_SET/BioAnalyze/data/processed/`:
    `canis_lupus_familiaris`, `cavia_porcellus`, `equus_caballus`,
    `felis_catus`, `heterocephalus_glaber`, `macaca_mulatta`,
    `mus_musculus`, `oryctolagus_cuniculus`, `sus_scrofa`
  - Each non-skipped species writes:
    - `{species}_expr_advanced_H2A_present_gold.tsv`
    - `{species}_h2a_canonical_variant_map.tsv`
- Generated heatmaps:
  - New species folders under `CURATED_SET/BioAnalyze/figures/heatmaps/` for
    each non-skipped species listed above
  - Full-output species write `all`, `clustered`, and `variants`
  - Partial-output species write `all` plus the available split only
- Verification notes:
  - Legacy regression check passed for `Bos taurus` with unchanged output shape:
    rows_after_filter=4644, all=23x127, clustered=16x127, variants=7x127
  - Canonical-like normalization is reflected in map TSVs; for example,
    `CURATED_SET/BioAnalyze/data/processed/mus_musculus/mus_musculus_h2a_canonical_variant_map.tsv`
    now classifies `H2ac1` / `H2ac13` (`cH2A*`) as `clustered`.

Shared Gene Summary + Gene Compare (2026-03-17)
- Goal:
  - Build a reusable cross-species summary of canonical H2A `gene_name` usage
    across the species heatmaps and use it as a fast lookup source for
    per-gene comparison heatmaps.
- Summary/index script:
  - `CURATED_SET/BioAnalyze/scripts/expression/summarize_shared_h2a_gene_names_across_species.py`
- Behavior:
  - Source species are discovered from `CURATED_SET/BioAnalyze/figures/heatmaps/`
    while skipping non-species folders such as `alligned_human_pan` and
    `gene_compare`.
  - Presence is determined from each species `*_expr_advanced_H2A_present_gold.tsv`
    by `Gene ID`.
  - Canonical display name is taken from the corresponding
    `*_canonical_variant_map.tsv` by `ensembl_gene_id -> gene_name`.
  - The script now writes both:
    - summary stats in `CURATED_SET/BioAnalyze/stats/`
    - a reusable long-form detail index in
      `CURATED_SET/BioAnalyze/data/gene_compare/index/`
- Summary outputs:
  - `CURATED_SET/BioAnalyze/stats/shared_h2a_gene_names_across_species.csv`
  - `CURATED_SET/BioAnalyze/stats/shared_h2a_gene_names_across_species.png`
  - `CURATED_SET/BioAnalyze/stats/shared_h2a_gene_names_across_species.svg`
- Detail index output:
  - `CURATED_SET/BioAnalyze/data/gene_compare/index/shared_h2a_gene_names_across_species_detail.csv`
- Summary CSV columns now include:
  - `canonical_gene_name`
  - `species_with_gene`
  - `species_ensembl_ids`
  - `detail_index_path`
- Current summary snapshot:
  - detail index rows: 179
  - shared genes (`species_count > 1`): 33
  - top shared genes include `H2AJ`, `H2AZ1`, `H2AZ2`, `MACROH2A1`,
    `MACROH2A2` with 8 species each

Gene Compare Heatmaps (2026-03-17)
- Goal:
  - Build a universal cross-species heatmap for any canonical `gene_name`,
    with:
    - X axis = species
    - Y axis = tissues (`Anatomical entity name`)
    - union of tissues across matching species
- New universal script:
  - `CURATED_SET/BioAnalyze/scripts/expression/build_gene_compare_heatmap.py`
- CLI defaults:
  - `--heatmap-dir` -> `CURATED_SET/BioAnalyze/figures/heatmaps`
  - `--processed-dir` -> `CURATED_SET/BioAnalyze/data/processed`
  - `--out-fig-root` -> `CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare`
  - `--out-data-root` -> `CURATED_SET/BioAnalyze/data/gene_compare`
  - `--detail-index` -> `CURATED_SET/BioAnalyze/data/gene_compare/index/shared_h2a_gene_names_across_species_detail.csv`
  - `--tissue-mode` -> `union`
  - `--aggregate` -> `mean`
- Heatmap behavior:
  - Species are selected from the detail index by canonical `gene_name`.
  - Values are aggregated as mean `Expression score` per `(species, tissue)`.
  - Plot values are transformed as `log10(Expression score + 1)`.
  - Tissues are ordered by the number of species where they are present, then
    alphabetically.
  - Missing values remain masked in the heatmap.
- First built example:
  - target gene: `H2AJ`
  - species included:
    `bos_taurus`, `canis_lupus_familiaris`, `equus_caballus`, `felis_catus`,
    `human`, `macaca_mulatta`, `pan_troglodytes`, `sus_scrofa`
  - union tissues: 346
- H2AJ outputs:
  - figures:
    - `CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AJ/H2AJ_gene_compare_heatmap.png`
    - `CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AJ/H2AJ_gene_compare_heatmap.svg`
  - data:
    - `CURATED_SET/BioAnalyze/data/gene_compare/H2AJ/H2AJ_gene_compare_long.csv`
    - `CURATED_SET/BioAnalyze/data/gene_compare/H2AJ/H2AJ_gene_compare_matrix.csv`
    - `CURATED_SET/BioAnalyze/data/gene_compare/H2AJ/H2AJ_gene_compare_metadata.json`
- Verification:
  - `H2AJ` matrix was generated with 8 species columns and 346 tissues.
  - Unknown `gene_name` produces a clear runtime error.
  - Multi-Ensembl aggregation was validated synthetically:
    multiple Ensembl IDs for the same `(species, tissue)` collapse to one mean
    value in the matrix.

Detailed Gene*Tissue Report + Low-Tail Panels (2026-03-18)
- Goal:
  - Extend the cross-species `gene*tissue` workflow from shortlist-style output
    to a full report layer with reproducible summary tables, low-variability
    candidate extraction, species-level tendencies, and direct links to the
    supporting heatmaps/panel pages.
- Ranking/report script updates:
  - `CURATED_SET/BioAnalyze/scripts/expression/rank_cross_species_h2a_differences.py`
    now writes, in addition to the existing shortlist outputs:
    - `gene_variability_summary.csv`
    - `gene_expression_overall_summary.csv`
    - `class_variability_summary.csv`
    - `species_extrema_summary.csv`
    - `species_pair_summary.csv`
    - `low_variability_candidates_p05.csv`
    - `low_variability_candidates_p10.csv`
    - `class_low_variability_candidates.csv`
    - detailed markdown report
      `stats/ranking/reports/cross_species_expression_report_ru.md`
- Report content:
  - The new report includes:
    - overview thresholds (`p05/p10/p90/p95`)
    - class-level comparison (`clustered` vs `variant`)
    - most variable genes
    - most conserved genes
    - most variable `gene*tissue` combinations
    - most conserved `gene*tissue` combinations
    - species-level tendency summaries
    - interpretation-ready findings for manuscript drafting
    - direct links to per-gene heatmaps and panel pages where available
  - For each key `gene*tissue` case, the report now includes:
    - `min_species` / `max_species`
    - `min_score` / `max_score`
    - `range`
    - tissue-level mean and median across species
    - links to `gene_compare` heatmaps
    - links to ranking panels via `panel_membership.csv`
- Plotting updates:
  - `CURATED_SET/BioAnalyze/scripts/expression/plot_cross_species_candidate_panels.py`
    now supports:
    - class-specific high-tail pages for both `p90` and `p95`
    - global low-tail pages for `p05` and `p10`
    - class-specific low-tail pages for `p05` and `p10`
    - `panel_membership.csv` index for report linking
- New plots written to `CURATED_SET/BioAnalyze/stats/ranking/plots/`:
  - high-tail:
    - `clustered_p90_panels_page*.png/.svg`
    - `variant_p90_panels_page*.png/.svg`
    - refreshed `clustered_p95_panels_page*.png/.svg`
    - refreshed `variant_p95_panels_page*.png/.svg`
  - low-tail:
    - `global_p05_low_panels_page*.png/.svg`
    - `global_p10_low_panels_page*.png/.svg`
    - `clustered_p05_low_panels_page*.png/.svg`
    - `clustered_p10_low_panels_page*.png/.svg`
    - `variant_p05_low_panels_page*.png/.svg`
    - `variant_p10_low_panels_page*.png/.svg`
  - index:
    - `panel_membership.csv`
- Current data-level findings captured by the new summaries:
  - `clustered` genes remain much more variable than `variant` genes in the
    conservative `species_n >= 4` set:
    - clustered median range `36.63`, p95 `62.03`
    - variant median range `10.01`, p95 `37.50`
  - `H2AC14` is currently the top gene by observed max range (`81.27`).
  - The most conserved genes by median range remain `H2AZ1`, `H2AZ2`,
    and `MACROH2A1`.
  - Global low-tail (`p05/p10`) is strongly enriched for variant genes, while
    class-specific low-tail extraction is needed to surface comparatively stable
    clustered cases.
