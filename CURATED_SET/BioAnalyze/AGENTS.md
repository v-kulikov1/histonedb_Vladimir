# BioAnalyze

From now on, any script-based work must be added under CURATED_SET/BioAnalyze/scripts so it is reproducible.

Location: CURATED_SET/BioAnalyze

Description:
BioAnalyze is a local data-processing pipeline for histone H2A datasets. It standardizes raw inputs, merges taxonomy, adds canonical gene names in `v4`, produces audits, and renders summary figures and stats used by the HistoneDB curation workflow.

Structure:
- audits/ (pipeline QA and validation tables)
- data/raw (placeholder; raw files are stored externally)
- data/processed (intermediate normalized tables)
- data/merged (final merged outputs)
- figures/heatmaps (heatmap graphics)
- figures/stats (summary stats graphics)
- scripts/annotate (annotation utilities)
- scripts/merge_taxonomy (taxonomy merge pipeline)
- stats (summary tables)
- WORKLOG.md (chronological process log and decisions)

Rules:
- All pipeline scripts must live under CURATED_SET/BioAnalyze/scripts (no scripts in data/ or elsewhere).
- All work performed must be recorded in CURATED_SET/BioAnalyze/WORKLOG.md.

Canonical output:
- CURATED_SET/BioAnalyze/data/merged/mammalia_H2A_merged_with_taxonomy_v4.csv

Worklog:
- CURATED_SET/BioAnalyze/WORKLOG.md (detailed, chronological record of pipeline decisions and artifacts)

Recent Outputs (2026-03-15)
- CURATED_SET/BioAnalyze/data/merged/mammalia_H2A_merged_with_taxonomy_v4.csv
- CURATED_SET/BioAnalyze/data/processed/homo_sapiens/h2a_hs_canonical_variant_map.tsv
- CURATED_SET/BioAnalyze/data/processed/pan_troglodytes/pan_troglodytes_h2a_canonical_variant_map.tsv
- CURATED_SET/BioAnalyze/data/processed/bos_taurus/bos_taurus_h2a_canonical_variant_map.tsv
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
- CURATED_SET/BioAnalyze/figures/heatmaps/alligned_human_pan/hs_aligned_all.png
- CURATED_SET/BioAnalyze/figures/heatmaps/alligned_human_pan/hs_aligned_all.svg
- CURATED_SET/BioAnalyze/figures/heatmaps/alligned_human_pan/hs_aligned_clustered.png
- CURATED_SET/BioAnalyze/figures/heatmaps/alligned_human_pan/hs_aligned_clustered.svg
- CURATED_SET/BioAnalyze/figures/heatmaps/alligned_human_pan/hs_aligned_variants.png
- CURATED_SET/BioAnalyze/figures/heatmaps/alligned_human_pan/hs_aligned_variants.svg

External raw data:
- stored at C:\Users\USER\Documents\GitHub\histonedb_external_storage\BioAnalyze\raw
See data/raw/README.md for the raw file list.

