# BioAnalyze H2A Merge Worklog

Last updated: 2026-03-31

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
    ranks strong cross-species expression differences at the `gene Г— tissue`
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

H2A.J Synteny Plot Integration (2026-03-24)
- Goal:
  - Move the ad hoc H2A.J synteny notebook workflow into the standard
    `BioAnalyze` layout so it can be rebuilt from CLI without a notebook.
- Structure changes:
  - Added reproducible builder script:
    - `CURATED_SET/BioAnalyze/scripts/h2aj_synteny/build_h2aj_synteny_plot.py`
  - Canonical figure outputs now live in:
    - `CURATED_SET/BioAnalyze/figures/h2aj_synteny/h2aj_synteny.png`
    - `CURATED_SET/BioAnalyze/figures/h2aj_synteny/h2aj_synteny.svg`
  - The repo-local placeholder folder `CURATED_SET/BioAnalyze/synteny` was
    retired in favor of the standard `scripts/` + `figures/` layout.
- Inputs:
  - External raw synteny TSV files remain at:
    - `C:\Users\USER\Documents\GitHub\histonedb_external_storage\BioAnalyze\synteny`
  - The script uses the current notebook species set and phylogenetic order:
    - Homo sapiens
    - Pan troglodytes
    - Mus musculus
    - Oryctolagus cuniculus
    - Erinaceus europaeus
    - Eptesicus fuscus
    - Canis lupus
    - Sus scrofa
    - Equus caballus
    - Manis pentadactyla
    - Loxodonta africana
    - Choloepus didactylus
    - Dasypus novemcinctus
- Behavior:
  - Preserves the latest notebook data-prep logic:
    - detect H2A.J from `Name`, with fallback to `Symbol`
    - operate within the H2A.J chromosome only
    - normalize minus-strand rows into a single rightward H2A.J layout
    - filter neighboring genes by cross-species recurrence
  - Presentation-focused styling changes:
    - large right-side legend for neighboring genes
    - no `bp` distance labels
    - no repeated per-row `H2A.J` labels
    - single `H2A.J` label under the bottom central arrow
    - species labels remain bold and black
    - order labels are placed under species names with color coding
    - extra left margin is reserved so label blocks do not overlap the arrows
    - widened rectangular canvas and enlarged arrows for slide readability
    - presentation-enlarged variant increases left label fonts, legend text, and
      bottom `H2A.J` annotation substantially for slide use
- Verification:
  - Builder script is intended for direct CLI execution and writes both PNG and
    SVG outputs for presentation use.

Codon Heatmap Integration (2026-03-24)
- Goal:
  - Move the SQK codon notebook workflow into the standard `BioAnalyze`
    layout so codon heatmaps can be rebuilt from CLI and raw FASTA inputs can
    live in external storage instead of an ad hoc local folder.
- Structure changes:
  - Added reproducible builder script:
    - `CURATED_SET/BioAnalyze/scripts/codons/build_codon_heatmaps.py`
  - Canonical figure outputs now live in:
    - `CURATED_SET/BioAnalyze/figures/codons/sqk_nuc_without_short_codon_entropy_annotated.png`
    - `CURATED_SET/BioAnalyze/figures/codons/sqk_nuc_without_short_codon_entropy_annotated.svg`
    - `CURATED_SET/BioAnalyze/figures/codons/sqk_nuc_full_codon_entropy_annotated.png`
    - `CURATED_SET/BioAnalyze/figures/codons/sqk_nuc_full_codon_entropy_annotated.svg`
- Inputs:
  - Raw SQK FASTA files were copied to external storage:
    - `C:\Users\USER\Documents\GitHub\histonedb_external_storage\BioAnalyze\raw\codons`
  - The current builder expects these four files:
    - `protein_from_SQK_nuc(without short).fasta`
    - `SQK_nuc(without short).fasta`
    - `protein_from_SQK_nuc.fasta`
    - `SQK_nuc.fasta`
- Behavior:
  - Preserves the final notebook plotting workflow:
    - normalized Shannon entropy by protein position and mammalian order
    - order filtering to `count > 3`
    - non-synonymous overlays marked with `*`
    - left-side structural region annotations
    - highlighted reference positions with color-coded amino-acid labels
    - top order-color strip and separate legend export
- Notebook usability:
  - The repo `.venv` was extended with `ipykernel` and `notebook` so the codon
    notebook can be opened from the IDE on the project environment instead of
    relying on transient Colab-only setup.

SQK Full Codon Figure Simplification (2026-03-24)
- Goal:
  - Make the `full` codon heatmap easier to assemble on slides by enlarging and
    shortening the left structural labels and exporting the supporting legends
    as standalone files.
- Output changes:
  - Refined main figure:
    - `CURATED_SET/BioAnalyze/figures/codons/sqk_nuc_full_codon_entropy_annotated.png`
    - `CURATED_SET/BioAnalyze/figures/codons/sqk_nuc_full_codon_entropy_annotated.svg`
  - New standalone support figures:
    - `CURATED_SET/BioAnalyze/figures/codons/sqk_nuc_full_order_legend.png`
    - `CURATED_SET/BioAnalyze/figures/codons/sqk_nuc_full_order_legend.svg`
    - `CURATED_SET/BioAnalyze/figures/codons/sqk_nuc_full_shannon_scale.png`
    - `CURATED_SET/BioAnalyze/figures/codons/sqk_nuc_full_shannon_scale.svg`
- Behavior:
  - The `full` heatmap now uses large shortened region labels:
    - `О±1e`, `О±1`, `L1`, `О±2`, `L2`, `О±3`, `ОІ3`
  - After presentation tuning, the left region brackets and labels were moved
    back closer to the heatmap body while keeping the larger slide-friendly
    font size.
  - The main `full` heatmap keeps the top order-color strip but removes the
    embedded Shannon colorbar.
  - A standalone legend now maps the top-strip colors to mammalian orders, and
    a separate vertical Shannon scale is exported with a large caption below
    the bar.
  - The entropy calculation, order filtering, highlighted positions, and
    non-synonymous `*` overlay remain unchanged.
  - Data review showed that the `full` Chiroptera signal is driven almost
    entirely by one added entry absent from `without-short`:
    - `Myotis-lucifugus|XM_006084274`
    - length `128 aa` versus the reference `129 aa`
    - removing this one record returns Chiroptera from `69` nonsynonymous
      positions in `full` back to `2`, matching `without-short`

SQK Without-Short Presentation Match (2026-03-25)
- Goal:
  - Make the `without-short` codon heatmap use the same enlarged presentation
    layout as the accepted `full` figure so the two outputs can be used
    interchangeably on slides.
- Output changes:
  - Refreshed:
    - `CURATED_SET/BioAnalyze/figures/codons/sqk_nuc_without_short_codon_entropy_annotated.png`
    - `CURATED_SET/BioAnalyze/figures/codons/sqk_nuc_without_short_codon_entropy_annotated.svg`
- Behavior:
  - `without-short` now uses the same presentation mode as `full`, including:
    - enlarged shortened left labels
    - closer left brackets/labels
    - the same overall slide-friendly scale and spacing
    - hidden embedded Shannon colorbar in the main heatmap

Human H2A Boxplot Builder (2026-03-25)
- Goal:
  - Add a reproducible CLI builder for a human H2A boxplot aligned to the
    normalized species heatmap data instead of manually extracting values from
    the rendered SVG.
- Script:
  - Added:
    - `CURATED_SET/BioAnalyze/scripts/expression/build_human_h2a_boxplot.py`
- Inputs:
  - Uses the existing heatmap-aligned processed files:
    - `CURATED_SET/BioAnalyze/data/processed/homo_sapiens/homo_sapiens_expr_advanced_H2A_present_gold.tsv`
    - `CURATED_SET/BioAnalyze/data/processed/homo_sapiens/h2a_hs_canonical_variant_map.tsv`
- Outputs:
  - Writes a dedicated human boxplot bundle to:
    - `CURATED_SET/BioAnalyze/figures/boxplot/human/h2a_human_boxplot.png`
    - `CURATED_SET/BioAnalyze/figures/boxplot/human/h2a_human_boxplot.svg`
    - `CURATED_SET/BioAnalyze/figures/boxplot/human/h2a_human_boxplot_source_values.tsv`
    - `CURATED_SET/BioAnalyze/figures/boxplot/human/h2a_human_boxplot_stats.tsv`
- Behavior:
  - Reuses the displayable human H2A gene selection logic from the heatmap
    pipeline so the plot covers the same `29` displayed genes.
  - Builds one distribution per gene from tissue-level `cell_mean_score`
    values, keeps `observed_zero` rows as `0`, and excludes generic tissues:
    - `multicellular organism`
    - `anatomical system`
    - `material anatomical entity`
  - Exports per-gene summary statistics needed for a boxplot, including:
    - `mean`, `std`, `q1`, `median`, `q3`, `iqr`, whiskers, and outlier count
- Verification:
  - CLI run completed successfully and produced all four expected outputs.
  - The source table contains no excluded generic tissues.
  - The stats table contains all `29` genes with nonzero tissue counts.
  - Manual spot checks for `H2AZ1` and `MACROH2A2` confirmed that `q1`,
    `median`, `q3`, and `std` match recalculation from the exported
    `source_values.tsv`.

Human H2A Presentation Boxplot Refinement (2026-03-25)
- Goal:
  - Make the human H2A boxplot slide-friendly by replacing left-side gene
    labels with numbers, enforcing a curated `clustered -> variants` order,
    increasing the X-axis label size, and exporting a separate legend image.
- Output changes:
  - Added presentation-specific outputs under:
    - `CURATED_SET/BioAnalyze/figures/boxplot/human/h2a_human_boxplot_presentation.png`
    - `CURATED_SET/BioAnalyze/figures/boxplot/human/h2a_human_boxplot_presentation.svg`
    - `CURATED_SET/BioAnalyze/figures/boxplot/human/h2a_human_boxplot_number_legend.png`
  - Updated analytical tables:
    - `CURATED_SET/BioAnalyze/figures/boxplot/human/h2a_human_boxplot_source_values.tsv`
    - `CURATED_SET/BioAnalyze/figures/boxplot/human/h2a_human_boxplot_stats.tsv`
- Behavior:
  - The builder now assigns `presentation_number` values `1..29` from top to
    bottom using an explicit human gene order:
    - clustered genes first
    - variant genes second
  - The presentation plot uses only these numbers on the Y axis and removes the
    Y-axis title entirely.
  - The X-axis label `Expression score` is rendered at double the default size
    used by the analytical plot.
  - The presentation plot keeps vertical X-grid lines and adds light horizontal
    row guides so numbered rows are easier to align on a slide.
  - A standalone PNG legend now maps each presentation number to its gene name
    in the same top-to-bottom order as the presentation plot.
- Verification:
  - The regenerated tables contain `presentation_number` and cover all `29`
    human genes.
  - The presentation SVG contains numbered Y tick labels and no embedded
    gene-name Y labels.
  - The fixed presentation order was confirmed as `clustered -> variants`.

Human H2A Presentation Styling Polish (2026-03-25)
- Goal:
  - Improve the slide-readability of the presentation boxplot by enlarging
    numeric tick labels, removing the top title, and making the standalone gene
    legend feel closer to the scientific support-legend style used elsewhere in
    `BioAnalyze`.
- Styling changes:
  - The presentation plot now uses:
    - doubled Y-side gene-number labels
    - doubled bottom X tick labels
    - extra left/bottom layout space so the larger labels do not crowd the plot
    - no top title above the presentation panel
  - The standalone legend PNG was restyled to use:
    - right-aligned gene numbers
    - fixed spacing between the number column and the gene-name column
    - small colored class swatches for clustered vs variant rows
    - cleaner manuscript-like alignment inspired by the codon order legend
- Outputs refreshed:
  - `CURATED_SET/BioAnalyze/figures/boxplot/human/h2a_human_boxplot_presentation.png`
  - `CURATED_SET/BioAnalyze/figures/boxplot/human/h2a_human_boxplot_presentation.svg`
  - `CURATED_SET/BioAnalyze/figures/boxplot/human/h2a_human_boxplot_number_legend.png`

Graphics Batch 2503 (2026-03-25)
- Goal:
  - Consolidate the current March 25 graphics refresh across BioAnalyze into a
    reproducible batch covering normalized-expression infrastructure, source
    auditing, refreshed H2A heatmaps, presentation-oriented human/pan outputs,
    candidate-panel refreshes, and the new human H2A boxplot workflow.
- Expression / audit layer:
  - Added normalized helper layer:
    - `CURATED_SET/BioAnalyze/scripts/expression/normalized_expression_common.py`
  - Added Bgee source audit builder:
    - `CURATED_SET/BioAnalyze/scripts/expression/audit_bgee_h2a_expression_sources.py`
  - Exported audit outputs:
    - `CURATED_SET/BioAnalyze/audits/bgee_h2a_expression_source_audit_3species.tsv`
    - `CURATED_SET/BioAnalyze/stats/ranking/reports/bgee_h2a_expression_source_audit_3species.md`
- Heatmaps / presentation figures:
  - Refreshed species heatmaps from the normalized expression layer, including:
    - the restored `all` human H2A heatmap
    - new `coverage_ge70` filtered panels for human and mouse
    - updated remaining-species canonical / variant / all panels
  - Refreshed cross-species gene-compare exports and added the new:
    - `H2AB3_gene_compare_*` tables and heatmap
  - Added a presentation-specific human/pan aligned figure set via:
    - `CURATED_SET/BioAnalyze/scripts/expression/build_human_pan_presentation_heatmaps.py`
    - outputs under `CURATED_SET/BioAnalyze/figures/heatmaps/presentation_human_pan/`
- Boxplots / ranking graphics:
  - Added the reproducible human H2A boxplot workflow:
    - `CURATED_SET/BioAnalyze/scripts/expression/build_human_h2a_boxplot.py`
    - analytical outputs under `CURATED_SET/BioAnalyze/figures/boxplot/human/`
    - presentation outputs with numeric labels and standalone legend
  - Refreshed ranking candidate panels, overview plots, and supporting tables
    after the normalized-expression and plotting updates.
- Commit intent:
  - This batch is intended to be committed together as a single graphics-focused
    checkpoint for `25.03`, combining the new plotting scripts with the
    regenerated tables, audits, and figure outputs they produced.

Tree Reconstruction Route (2026-03-29)
- Goal:
  - Reconstruct the historical H2A.J / cH2A / nonplacental tree-preparation
    workflow, archive the supporting evidence, and explain why the nucleotide
    tree failed.
- New folders:
  - Script:
    - `CURATED_SET/BioAnalyze/scripts/h2aj_tree/`
  - Audit:
    - `CURATED_SET/BioAnalyze/audits/h2aj_tree/`
  - Data:
    - `CURATED_SET/BioAnalyze/data/h2aj_tree/`
- New script:
  - `CURATED_SET/BioAnalyze/scripts/h2aj_tree/rebuild_h2aj_tree_history.py`
- Inputs archived outside git:
  - `C:\Users\USER\Documents\GitHub\histonedb_external_storage\BioAnalyze\raw\tree`
  - historical files from:
    - `C:\Users\USER\Documents\Р“РёСЃС‚РѕРЅС‹`
    - `C:\Users\USER\Documents\My_test`
    - `C:\Users\USER\Documents\Р Р°Р±РѕС‚Р° РЅР°Рґ РіСЂР°РЅС‚РѕРј HistoneJ`
    - `C:\Users\USER\Downloads\Telegram Desktop\ChatExport_2026-03-29\result.json`
    - `C:\Users\USER\Downloads\Telegram Desktop\Draft2.docx`
- Behavior:
  - Copies a minimal historical evidence set into external storage and writes a
    manifest with checksum / mtime / count metadata.
  - Rebuilds both `historical` and `clean` tree-input outputs for:
    - protein stages
    - nucleotide stages
    - post-PhyML rooting
  - Exports evidence tables from notebooks, Telegram chat, and `Draft2.docx`.
  - Reproduces the documented nucleotide rooting failure where the full
    `cH2A_nuc` set is not monophyletic.
- Key findings:
  - The confirmed successful archived web run is the AA tree from `2025-06-03`
    with:
    - `LG`
    - `BioNJ`
    - `FreeRate(3)`
    - `SH-like branch supports`
  - `All_AA_0206 = 387` is the early broad-outgroup protein run.
  - `All_AA_1006 = 338` matches:
    - `H2AJ = 230`
    - `cH2A = 57`
    - `nonplacental/platypus = 51`
  - The nucleotide failure is upstream of tree-building:
    - BLAST HSP fragments instead of full CDS
    - partial-codon translations
    - SQK/SHK filtering by aligned position
    - protein sequences mixed into the nucleotide tree
    - short-sequence artifacts
- Outputs:
  - `CURATED_SET/BioAnalyze/data/h2aj_tree/evidence/`
  - `CURATED_SET/BioAnalyze/data/h2aj_tree/historical/`
  - `CURATED_SET/BioAnalyze/data/h2aj_tree/clean/`
  - `CURATED_SET/BioAnalyze/audits/h2aj_tree/h2aj_tree_reconstruction_audit.md`
- Verification:
  - `py_compile` succeeded for:
    - `CURATED_SET/BioAnalyze/scripts/h2aj_tree/rebuild_h2aj_tree_history.py`
  - Full CLI run succeeded with:
    - `--phase all --profile both`
  - Rebuilt control points include:
    - `All_AA_0206 = 387 / 208`
    - `All_AA_1006 = 338 / 180`
    - `All_NUC_1606 = 315 / 396`
    - `SQK_nuc = 227`
    - `SQK_nuc(without short) = 224`
    - merged July nucleotide set = `312`

BioAnalyze Structure Cleanup (2026-03-31)
- Goal:
  - Consolidate the recovered tree workflow under a single `h2aj_tree`
    namespace and remove the abandoned nonplacental H2A.J synteny research
    branch.
- Tree changes:
  - Moved:
    - `CURATED_SET/BioAnalyze/scripts/tree/` ->
      `CURATED_SET/BioAnalyze/scripts/h2aj_tree/`
    - `CURATED_SET/BioAnalyze/audits/tree/` ->
      `CURATED_SET/BioAnalyze/audits/h2aj_tree/`
    - `CURATED_SET/BioAnalyze/stats/tree/` ->
      `CURATED_SET/BioAnalyze/data/h2aj_tree/`
  - Updated the recovered tree CLI so its default output root now points to:
    - `CURATED_SET/BioAnalyze/data/h2aj_tree/`
- Synteny cleanup:
  - Preserved the placental H2A.J synteny route:
    - `CURATED_SET/BioAnalyze/scripts/h2aj_synteny/build_h2aj_synteny_plot.py`
    - `CURATED_SET/BioAnalyze/figures/h2aj_synteny/`
  - Removed the later nonplacental branch entirely, including:
    - NCBI survey outputs
    - BLAST-first rescue route
    - `EMP1/MGP` hand-analysis route
    - local-search probing route
    - the whole `CURATED_SET/BioAnalyze/stats/h2aj_synteny/` tree
- Documentation cleanup:
  - Updated:
    - `CURATED_SET/BioAnalyze/WORKLOG.md`
    - `CURATED_SET/BioAnalyze/AGENTS.md`
    - `CURATED_SET/BioAnalyze/audits/h2aj_tree/h2aj_tree_reconstruction_audit.md`
  - Removed stale references to:
    - `scripts/tree`
    - `audits/tree`
    - `stats/tree`
    - nonplacental `h2aj_synteny` routes
- Verification:
  - Confirmed the kept placental synteny files still exist.
  - Confirmed `CURATED_SET/BioAnalyze/stats/h2aj_synteny/` is gone.
  - `py_compile` succeeded for:
    - `CURATED_SET/BioAnalyze/scripts/h2aj_tree/rebuild_h2aj_tree_history.py`

Codon Entropy Renormalization (2026-03-30)
- Goal:
  - Change SQK codon-entropy normalization so each cell is scaled by the full
    synonymous codon space of the reference amino acid, not by the number of
    codons actually observed in that order.
- Script updated:
  - `CURATED_SET/BioAnalyze/scripts/codons/build_codon_heatmaps.py`
- Behavior:
  - Added an explicit `AA_MAX_SYNONYMOUS_CODONS` table for the standard genetic
    code and a validation helper that checks it against
    `CodonTable.unambiguous_dna_by_id[1]`.
  - Entropy is now computed only from codons translating to the reference amino
    acid at that position.
  - Nonsynonymous codons and `---` are excluded from the entropy numerator.
  - The `*` overlay remains separate and still marks cells where at least one
    observed codon translates to a non-reference amino acid.
  - Cells containing only nonsynonymous codons and/or `---` now resolve to
    `NaN` instead of contributing zero diversity.
- Verification:
  - Synthetic helper cases cover:
    - `L`: `3` observed synonymous codons normalized by `6` possible
    - `I`: `2` observed synonymous codons normalized by `3` possible
    - `M` / `W`: `0.0` despite single-codon synonymous space
    - mixed synonymous plus nonsynonymous cells preserving the `*` overlay
    - cells with only nonsynonymous codons and/or `---` returning `NaN`
  - Full CLI rebuild should regenerate both `without-short` and `full`
    heatmaps with unchanged matrix dimensions but updated entropy values.

Codon Entropy Majority-AA Basis Update (2026-03-30)
- Goal:
  - Remove the white `NaN` stripes at SQK position `126` and align both entropy
    and `*` semantics to the most frequent amino acid observed in each
    `(position, order)` cell.
- Script updated:
  - `CURATED_SET/BioAnalyze/scripts/codons/build_codon_heatmaps.py`
- Behavior:
  - The entropy basis is now the most frequent translated amino acid in the
    cell, excluding `---`.
  - If multiple amino acids are tied for first place, the reference amino acid
    is preferred when it is among the tied leaders; otherwise the tie is broken
    deterministically in code.
  - Entropy still measures synonymous codon diversity only, but now relative to
    the synonymous space of the selected majority amino acid.
  - The `*` overlay now marks codons encoding amino acids different from that
    per-cell majority amino acid rather than different from the human
    reference amino acid.
- Position `126` impact:
  - The earlier reference-AA logic created exactly `5` white `NaN` cells at
    position `126` because no synonymous codons remained for `T` in:
    - `Carnivora`
    - `Chiroptera`
    - `Lagomorpha`
    - `Perissodactyla`
    - `Rodentia`
  - Under the majority-AA basis, those `5` cells should now render with finite
    entropy values while the rest of the matrix remains effectively unchanged.
- Verification:
  - Synthetic checks should cover:
    - majority-AA cells with mixed synonymous codons
    - tie resolution when reference AA is one of the leaders
    - `---`-only cells staying `NaN`
    - uniform non-reference amino acid cells getting finite entropy and no `*`
  - A full CLI rebuild should confirm:
    - both datasets remain `129 x 8`
    - the `5` white cells at `126` disappear
    - `Eulipotyphla` at `126` resolves via the tie rule

Alternative All-61 Codon Entropy Builder (2026-03-30)
- Goal:
  - Add a separate SQK codon-entropy visualization route that normalizes every
    cell by the same theoretical denominator:
    - `log2(61)` for the `61` sense codons in the standard genetic code
  - Keep the existing main codon builder untouched while generating an
    alternative interpretation set in its own output folder.
- New script:
  - `CURATED_SET/BioAnalyze/scripts/codons/build_codon_heatmaps_all61.py`
- Output folder:
  - `CURATED_SET/BioAnalyze/figures/codons/all_coding_61/`
- Behavior:
  - Entropy numerator is the raw Shannon entropy over all observed sense codons
    in a cell.
  - `---` is excluded from the entropy numerator.
  - No per-cell amino-acid basis or majority-amino-acid logic is used.
  - The `*` overlay in this variant marks cells containing more than one
    translated amino acid.
  - The script writes:
    - `without-short` heatmap PNG/SVG
    - `full` heatmap PNG/SVG
    - standalone order legend PNG/SVG
    - standalone Shannon scale PNG/SVG labeled with `/ log2(61)`
    - `all61_entropy_summary.tsv`
- Summary TSV contents:
  - global `sense_codon_count = 61`
  - per-amino-acid theoretical maxima under `log2(61)` normalization
  - observed min / max / mean entropy for each dataset
- Key observed values:
  - Current observed max entropy in both datasets is about `0.3374`
  - Current observed mean entropy is about:
    - `0.0496` for `without-short`
    - `0.0524` for `full`
  - The standalone all-61 Shannon scale is now matched to the actual `full`
    heatmap color range:
    - `vmin = 0`
    - `vmax = observed full max в‰€ 0.3374`
    This avoids the misleading earlier mismatch where the legend used `0..1`
    while the heatmap colors were auto-stretched only to the observed maximum.
  - Position `126` no longer contains the earlier white `NaN` cells caused by
    amino-acid filtering logic in the alternative all-61 variant.
- Verification:
  - `py_compile` succeeded for:
    - `CURATED_SET/BioAnalyze/scripts/codons/build_codon_heatmaps_all61.py`
  - Synthetic checks covered:
    - one-codon cells
    - two-codon equal mixtures
    - four-codon equal mixtures
    - `---`-only cells
    - star-mask behavior for one-AA vs multi-AA cells
  - Full CLI run succeeded with:
    - `.\.venv\Scripts\python.exe CURATED_SET\BioAnalyze\scripts\codons\build_codon_heatmaps_all61.py --dataset all`
  - Rebuilt matrices remained:
    - `129 x 8` for both datasets

H2A.J Tree Clean Short-Sequence Filter (2026-04-01)
- Goal:
  - Remove three known short H2A.J/SQK records from the `clean` `data/h2aj_tree` outputs without touching `historical`.
- Dropped clean ID prefixes:
  - `Myotis-lucifugus|XM_006084274`
  - `Homo-sapiens|AK303301`
  - `Homo-sapiens|AL133626`
- Script changes:
  - Added repeatable `--drop-clean-id` to `CURATED_SET/BioAnalyze/scripts/h2aj_tree/rebuild_h2aj_tree_history.py`.
  - Applied prefix-based filtering only to `clean` FASTA/PHY outputs.
  - Switched July clean SQK sources to:
    - `SQK_nuc(without short).fasta`
    - `protein_from_SQK_nuc(without short).fasta`
  - Added `clean/postprocess/README_filtered_clean_note.txt` to mark `NWK` files as historical reference only.
- Expected rebuilt counts:
  - `nuc_0907`: `312 -> 309`
  - `nuc_1606`: `315 -> 312`
  - `aa_1006`: `338 -> 335`
  - `aa_0206`: `387 -> 385`
  - `h2aj_aligned_historical`: `230 -> 227`
- Verification note:
  - `All_AA_0206.fasta` contains the two short human records but not the short bat record, so the clean rebuild correctly removes `2` rows there rather than `3`.

Human HPA vs Bgee cH2A Compare Heatmaps (2026-04-01)
- Goal:
  - Build a human-only compare layer between the new HPA raw-`nTPM` cH2A
    workflow and the existing human Bgee cH2A processed heatmap inputs.
- New script:
  - `CURATED_SET/BioAnalyze/scripts/expression/build_human_hpa_bgee_compare_heatmaps.py`
- Inputs:
  - `CURATED_SET/BioAnalyze/data/expression_nTPM/human/h2a_human_gene_ntpm_cells.tsv`
  - `CURATED_SET/BioAnalyze/data/processed/homo_sapiens/Homo_sapiens_expr_advanced_H2A_present_gold.tsv`
  - `CURATED_SET/BioAnalyze/data/processed/homo_sapiens/h2a_hs_canonical_variant_map.tsv`
- Gene mapping behavior:
  - Internal compare identity uses canonical human cH2A `gene_name` from the
    canonical map together with `ensembl_gene_id`.
  - Bgee display labels remain source-native from the processed TSV.
  - The compare layer explicitly preserves the alias case
    `H2AC25 <-> ENSG00000181218 <-> H2AW`.
- Tissue mapping behavior:
  - Exact matches are kept as-is.
  - Safe compare-only synonym pairs are applied for:
    - `breast -> mammary gland`
    - `caudate -> caudate nucleus`
    - `cervix -> uterine cervix`
    - `heart muscle -> myocardium`
    - `prostate -> prostate gland`
    - `salivary gland -> saliva-secreting gland`
    - `skeletal muscle -> skeletal muscle tissue`
  - Unresolved HPA tissues left outside the compare axis:
    - `hippocampus`
    - `skin`
- Display behavior:
  - Heatmaps keep source-native axis labels for manual inspection.
  - Safe-synonym tissue labels are visually marked by underline on both
    heatmaps.
  - Separate color scales are exported for HPA and Bgee because the sources use
    different value ranges and semantics.
- Outputs:
  - Data:
    - `CURATED_SET/BioAnalyze/data/processed/intersections/human_hpa_bgee/hpa_vs_bgee_gene_mapping.tsv`
    - `CURATED_SET/BioAnalyze/data/processed/intersections/human_hpa_bgee/hpa_vs_bgee_tissue_mapping.tsv`
    - `CURATED_SET/BioAnalyze/data/processed/intersections/human_hpa_bgee/hpa_vs_bgee_hpa_aligned_long.tsv`
    - `CURATED_SET/BioAnalyze/data/processed/intersections/human_hpa_bgee/hpa_vs_bgee_bgee_aligned_long.tsv`
    - `CURATED_SET/BioAnalyze/data/processed/intersections/human_hpa_bgee/hpa_vs_bgee_metadata.json`
  - Figures:
    - `CURATED_SET/BioAnalyze/figures/heatmaps/compare_nTPM_bgee/hpa_vs_bgee_hpa_heatmap.(png|svg)`
    - `CURATED_SET/BioAnalyze/figures/heatmaps/compare_nTPM_bgee/hpa_vs_bgee_bgee_heatmap.(png|svg)`
    - `CURATED_SET/BioAnalyze/figures/heatmaps/compare_nTPM_bgee/hpa_vs_bgee_hpa_colorbar.(png|svg)`
    - `CURATED_SET/BioAnalyze/figures/heatmaps/compare_nTPM_bgee/hpa_vs_bgee_bgee_colorbar.(png|svg)`
- Presentation update:
  - Compare heatmaps now use compact `OY` gene numbers and shared `OX`
    tissue letter codes instead of long axis labels.
  - Added presentation-style legend assets:
    - `CURATED_SET/BioAnalyze/figures/heatmaps/compare_nTPM_bgee/hpa_vs_bgee_gene_number_legend.(png|svg)`
    - `CURATED_SET/BioAnalyze/figures/heatmaps/compare_nTPM_bgee/hpa_vs_bgee_hpa_tissue_letter_legend.(png|svg)`
    - `CURATED_SET/BioAnalyze/figures/heatmaps/compare_nTPM_bgee/hpa_vs_bgee_bgee_tissue_letter_legend.(png|svg)`
  - Added matching map TSVs in the compare figure folder:
    - `hpa_vs_bgee_gene_number_map.tsv`
    - `hpa_vs_bgee_hpa_tissue_letter_map.tsv`
    - `hpa_vs_bgee_bgee_tissue_letter_map.tsv`
  - Tissue codes extend beyond `Z` to `AA-AI`; safe-synonym slots remain
    underlined on the heatmaps and in both tissue legends.
- Mean-value compare update:
  - The HPA side of the compare workflow now reads the exported HPA gene-by-
    tissue `nTPM` value instead of `median_nTPM`.
  - Bgee remains unchanged on `cell_mean_score`.
  - Because each HPA gene-by-tissue cell is unique in the source TSV, the
    compare layer is now displayed as plain `nTPM` in titles, legends, reports,
    and compare TSV exports.
- Raw-scale compare update:
  - The compare heatmap builder now also exports raw-scale comparison figures
    without the `log10(x + 1)` transform.
  - New raw compare figures:
    - `CURATED_SET/BioAnalyze/figures/heatmaps/compare_nTPM_bgee/hpa_vs_bgee_hpa_raw_heatmap.(png|svg)`
    - `CURATED_SET/BioAnalyze/figures/heatmaps/compare_nTPM_bgee/hpa_vs_bgee_bgee_raw_heatmap.(png|svg)`
    - `CURATED_SET/BioAnalyze/figures/heatmaps/compare_nTPM_bgee/hpa_vs_bgee_hpa_raw_colorbar.(png|svg)`
    - `CURATED_SET/BioAnalyze/figures/heatmaps/compare_nTPM_bgee/hpa_vs_bgee_bgee_raw_colorbar.(png|svg)`
  - The original log-scale compare figures are kept alongside the raw-scale
    versions for side-by-side interpretation.
- Zero-cell compare styling update:
  - Exact zero-valued cells on the human HPA-vs-Bgee compare heatmaps are now
    marked with a thin red outline.
  - This visual marker is applied to all four compare heatmaps:
    HPA log, Bgee log, HPA raw, and Bgee raw.
  - The `viridis` palette and color scales are unchanged so `0` is now
    distinguishable from low nonzero values such as `0.1` without changing the
    quantitative meaning of the colors.
- Compare stats update:
  - Added a separate stats/report workflow for the human HPA-vs-Bgee compare
    layer without touching the heatmap builder.
  - New script:
    - `CURATED_SET/BioAnalyze/scripts/expression/build_human_hpa_bgee_compare_stats.py`
  - Reads the aligned compare tables from:
    - `CURATED_SET/BioAnalyze/data/processed/intersections/human_hpa_bgee/`
  - Writes stats outputs to:
    - `CURATED_SET/BioAnalyze/stats/compare_nTPM_bgee/`
  - Generated outputs:
    - `paired_cells.tsv`
    - `correlation_summary.tsv`
    - `top_differing_cells.tsv`
    - `gene_difference_summary.tsv`
    - `tissue_difference_summary.tsv`
    - `report_md.md`
    - `metadata.json`
  - Comparison policy:
    - primary comparison space is `log10(value + 1)` to match the displayed
      heatmaps
    - missing values are excluded, but zero-valued cells are kept
    - raw `nTPM` and raw `cell_mean_score` are preserved alongside the
      transformed values for interpretation
  - Current compare summary:
    - total compare slots: `595`
    - non-missing paired cells: `477`
    - log-space Pearson (`all_non_missing`): `0.438963`
    - log-space Spearman (`all_non_missing`): `0.813329`
    - strongest disagreements are dominated by kidney-associated clustered
      `H2AC` rows, with one top safe-synonym case at
      `prostate -> prostate gland`
- Per-gene correlation update:
  - Added a dedicated per-gene correlation workflow on top of
    `stats/compare_nTPM_bgee/paired_cells.tsv`.
  - New script:
    - `CURATED_SET/BioAnalyze/scripts/expression/build_human_hpa_bgee_gene_correlations.py`
  - New output folder:
    - `CURATED_SET/BioAnalyze/stats/compare_nTPM_bgee/per_gene_correlations/`
  - Generated outputs:
    - `gene_correlation_keep_zeros.tsv`
    - `gene_correlation_drop_zeros.tsv`
    - `pearson_raw_keep_zeros.(png|svg)`
    - `spearman_raw_keep_zeros.(png|svg)`
    - `pearson_raw_drop_zeros.(png|svg)`
    - `spearman_raw_drop_zeros.(png|svg)`
    - `report_md.md`
  - Correlation policy:
    - one Pearson raw and one Spearman raw value per gene across tissues
    - `keep_zeros` keeps all non-missing pairs
    - `drop_zeros` keeps only rows where both HPA and Bgee are `> 0`
  - Current summary:
    - 15 comparable cH2A genes were scored in both modes
    - strongest `keep_zeros` Pearson: `H2AC1 = 0.758929`
    - strongest `keep_zeros` Spearman: `H2AC25/H2AW = 0.873274`
    - strongest `drop_zeros` Pearson: `H2AC1 = 0.990639` on 3 tissues
    - weakest `drop_zeros` Pearson: `H2AC14 = -0.170402`
    - no genes fell into `insufficient_pairs` in the current dataset
- Expected compare matrix shape:
  - `17 x 35` for both HPA and Bgee after the fixed gene/tissue mapping is applied.
