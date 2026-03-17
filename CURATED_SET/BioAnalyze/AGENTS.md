# BioAnalyze

From now on, any script-based work must be added under CURATED_SET/BioAnalyze/scripts so it is reproducible.

Location: CURATED_SET/BioAnalyze

Description:
BioAnalyze is a local data-processing pipeline for histone H2A datasets. It standardizes raw inputs, merges taxonomy, adds canonical gene names in `v4`, produces audits, and renders summary figures and stats used by the HistoneDB curation workflow.

Structure:
- audits/ (pipeline QA and validation tables)
- data/raw (placeholder; raw files are stored externally)
- data/processed (intermediate normalized tables)
- data/gene_compare (gene-level comparison indexes and per-gene exported tables)
- data/merged (final merged outputs)
- figures/heatmaps (heatmap graphics)
- figures/heatmaps/gene_compare (cross-species per-gene heatmaps)
- figures/stats (summary stats graphics)
- scripts/annotate (annotation utilities)
- scripts/merge_taxonomy (taxonomy merge pipeline)
- scripts/expression (expression, summary, and heatmap builders)
- stats (summary tables)
- WORKLOG.md (chronological process log and decisions)

Rules:
- All pipeline scripts must live under CURATED_SET/BioAnalyze/scripts (no scripts in data/ or elsewhere).
- All work performed must be recorded in CURATED_SET/BioAnalyze/WORKLOG.md.

Canonical output:
- CURATED_SET/BioAnalyze/data/merged/mammalia_H2A_merged_with_taxonomy_v4.csv

Worklog:
- CURATED_SET/BioAnalyze/WORKLOG.md (detailed, chronological record of pipeline decisions and artifacts)

Recent Outputs (2026-03-17)
- CURATED_SET/BioAnalyze/data/merged/mammalia_H2A_merged_with_taxonomy_v4.csv
- CURATED_SET/BioAnalyze/audits/audit_h2a_remaining_species_batch_v4.tsv
- CURATED_SET/BioAnalyze/data/gene_compare/index/shared_h2a_gene_names_across_species_detail.csv
- CURATED_SET/BioAnalyze/data/gene_compare/H2AJ/H2AJ_gene_compare_long.csv
- CURATED_SET/BioAnalyze/data/gene_compare/H2AJ/H2AJ_gene_compare_matrix.csv
- CURATED_SET/BioAnalyze/data/gene_compare/H2AJ/H2AJ_gene_compare_metadata.json
- CURATED_SET/BioAnalyze/data/processed/homo_sapiens/h2a_hs_canonical_variant_map.tsv
- CURATED_SET/BioAnalyze/data/processed/pan_troglodytes/pan_troglodytes_h2a_canonical_variant_map.tsv
- CURATED_SET/BioAnalyze/data/processed/bos_taurus/bos_taurus_h2a_canonical_variant_map.tsv
- CURATED_SET/BioAnalyze/data/processed/canis_lupus_familiaris/canis_lupus_familiaris_h2a_canonical_variant_map.tsv
- CURATED_SET/BioAnalyze/data/processed/cavia_porcellus/cavia_porcellus_h2a_canonical_variant_map.tsv
- CURATED_SET/BioAnalyze/data/processed/equus_caballus/equus_caballus_h2a_canonical_variant_map.tsv
- CURATED_SET/BioAnalyze/data/processed/felis_catus/felis_catus_h2a_canonical_variant_map.tsv
- CURATED_SET/BioAnalyze/data/processed/heterocephalus_glaber/heterocephalus_glaber_h2a_canonical_variant_map.tsv
- CURATED_SET/BioAnalyze/data/processed/macaca_mulatta/macaca_mulatta_h2a_canonical_variant_map.tsv
- CURATED_SET/BioAnalyze/data/processed/mus_musculus/mus_musculus_h2a_canonical_variant_map.tsv
- CURATED_SET/BioAnalyze/data/processed/oryctolagus_cuniculus/oryctolagus_cuniculus_h2a_canonical_variant_map.tsv
- CURATED_SET/BioAnalyze/data/processed/sus_scrofa/sus_scrofa_h2a_canonical_variant_map.tsv
- CURATED_SET/BioAnalyze/data/processed/pan_troglodytes/pan_troglodytes_expr_advanced_H2A_present_gold.tsv
- CURATED_SET/BioAnalyze/data/processed/pan_troglodytes/pan_troglodytes_h2a_canonical_variant_map.tsv
- CURATED_SET/BioAnalyze/data/processed/intersections/homo_sapiens_expr_advanced_H2A_present_gold_intersection.tsv
- CURATED_SET/BioAnalyze/figures/heatmaps/human/h2a_clustered.png
- CURATED_SET/BioAnalyze/figures/heatmaps/human/h2a_clustered.svg
- CURATED_SET/BioAnalyze/figures/heatmaps/human/h2a_variants.png
- CURATED_SET/BioAnalyze/figures/heatmaps/human/h2a_variants.svg
- CURATED_SET/BioAnalyze/figures/heatmaps/bos_taurus/h2a_bos_taurus_all.png
- CURATED_SET/BioAnalyze/figures/heatmaps/bos_taurus/h2a_bos_taurus_all.svg
- CURATED_SET/BioAnalyze/figures/heatmaps/bos_taurus/h2a_bos_taurus_clustered.png
- CURATED_SET/BioAnalyze/figures/heatmaps/bos_taurus/h2a_bos_taurus_clustered.svg
- CURATED_SET/BioAnalyze/figures/heatmaps/bos_taurus/h2a_bos_taurus_variants.png
- CURATED_SET/BioAnalyze/figures/heatmaps/bos_taurus/h2a_bos_taurus_variants.svg
- CURATED_SET/BioAnalyze/figures/heatmaps/pan_troglodytes/h2a_pan_troglodytes_all.png
- CURATED_SET/BioAnalyze/figures/heatmaps/pan_troglodytes/h2a_pan_troglodytes_all.svg
- CURATED_SET/BioAnalyze/figures/heatmaps/pan_troglodytes/h2a_pan_troglodytes_clustered.png
- CURATED_SET/BioAnalyze/figures/heatmaps/pan_troglodytes/h2a_pan_troglodytes_clustered.svg
- CURATED_SET/BioAnalyze/figures/heatmaps/pan_troglodytes/h2a_pan_troglodytes_variants.png
- CURATED_SET/BioAnalyze/figures/heatmaps/pan_troglodytes/h2a_pan_troglodytes_variants.svg
- CURATED_SET/BioAnalyze/figures/heatmaps/canis_lupus_familiaris/h2a_canis_lupus_familiaris_all.svg
- CURATED_SET/BioAnalyze/figures/heatmaps/cavia_porcellus/h2a_cavia_porcellus_all.svg
- CURATED_SET/BioAnalyze/figures/heatmaps/equus_caballus/h2a_equus_caballus_all.svg
- CURATED_SET/BioAnalyze/figures/heatmaps/felis_catus/h2a_felis_catus_all.svg
- CURATED_SET/BioAnalyze/figures/heatmaps/heterocephalus_glaber/h2a_heterocephalus_glaber_all.svg
- CURATED_SET/BioAnalyze/figures/heatmaps/macaca_mulatta/h2a_macaca_mulatta_all.svg
- CURATED_SET/BioAnalyze/figures/heatmaps/mus_musculus/h2a_mus_musculus_all.svg
- CURATED_SET/BioAnalyze/figures/heatmaps/oryctolagus_cuniculus/h2a_oryctolagus_cuniculus_all.svg
- CURATED_SET/BioAnalyze/figures/heatmaps/sus_scrofa/h2a_sus_scrofa_all.svg
- CURATED_SET/BioAnalyze/figures/heatmaps/alligned_human_pan/hs_aligned_all.png
- CURATED_SET/BioAnalyze/figures/heatmaps/alligned_human_pan/hs_aligned_all.svg
- CURATED_SET/BioAnalyze/figures/heatmaps/alligned_human_pan/hs_aligned_clustered.png
- CURATED_SET/BioAnalyze/figures/heatmaps/alligned_human_pan/hs_aligned_clustered.svg
- CURATED_SET/BioAnalyze/figures/heatmaps/alligned_human_pan/hs_aligned_variants.png
- CURATED_SET/BioAnalyze/figures/heatmaps/alligned_human_pan/hs_aligned_variants.svg
- CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AJ/H2AJ_gene_compare_heatmap.png
- CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AJ/H2AJ_gene_compare_heatmap.svg
- CURATED_SET/BioAnalyze/stats/shared_h2a_gene_names_across_species.csv
- CURATED_SET/BioAnalyze/stats/shared_h2a_gene_names_across_species.png
- CURATED_SET/BioAnalyze/stats/shared_h2a_gene_names_across_species.svg

Recent Script Additions (2026-03-17)
- CURATED_SET/BioAnalyze/scripts/expression/summarize_shared_h2a_gene_names_across_species.py
  Rebuilds the shared-gene summary in `stats/` and writes a reusable long-form
  gene lookup index to `data/gene_compare/index/`.
- CURATED_SET/BioAnalyze/scripts/expression/build_gene_compare_heatmap.py
  Builds a cross-species per-gene heatmap from canonical map labels and
  `*_present_gold.tsv` inputs, with:
  - X axis = species
  - Y axis = tissues
  - default `union` tissue mode
  - mean aggregation per `(species, tissue)`
  - masked missing values

External raw data:
- stored at C:\Users\USER\Documents\GitHub\histonedb_external_storage\BioAnalyze\raw
See data/raw/README.md for the raw file list.

