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
- figures/heatmaps (heatmap graphics root)
- figures/heatmaps/species (per-species H2A heatmaps)
- figures/heatmaps/gene_compare (cross-species per-gene heatmaps)
- figures/stats (summary stats graphics)
- scripts/annotate (annotation utilities)
- scripts/merge_taxonomy (taxonomy merge pipeline)
- scripts/expression (expression, summary, and heatmap builders)
- stats/shared_genes (shared-gene summary tables and plots)
- stats/ranking/tables (gene*tissue comparison tables)
- stats/ranking/reports (manuscript-oriented text outputs)
- stats/ranking/plots (gene*tissue comparison plots and p95 panels)
- stats/accession_stats (legacy accession summary tables)
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
- CURATED_SET/BioAnalyze/data/gene_compare/H2AX/H2AX_gene_compare_long.csv
- CURATED_SET/BioAnalyze/data/gene_compare/MACROH2A2/MACROH2A2_gene_compare_long.csv
- CURATED_SET/BioAnalyze/stats/shared_genes/shared_h2a_gene_names_across_species.csv
- CURATED_SET/BioAnalyze/stats/ranking/tables/gene_tissue_summary.csv
- CURATED_SET/BioAnalyze/stats/ranking/tables/conservative_candidates.csv
- CURATED_SET/BioAnalyze/stats/ranking/tables/high_confidence_candidates.csv
- CURATED_SET/BioAnalyze/stats/ranking/tables/pairwise_contrasts.csv
- CURATED_SET/BioAnalyze/stats/ranking/reports/manuscript_shortlist.md
- CURATED_SET/BioAnalyze/stats/ranking/plots/candidate_focus_panels.png
- CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p95_panels_page1.png
- CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p95_panels_page1.png
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
- CURATED_SET/BioAnalyze/figures/heatmaps/species/human/h2a_clustered.png
- CURATED_SET/BioAnalyze/figures/heatmaps/species/human/h2a_clustered.svg
- CURATED_SET/BioAnalyze/figures/heatmaps/species/human/h2a_variants.png
- CURATED_SET/BioAnalyze/figures/heatmaps/species/human/h2a_variants.svg
- CURATED_SET/BioAnalyze/figures/heatmaps/species/bos_taurus/h2a_bos_taurus_all.png
- CURATED_SET/BioAnalyze/figures/heatmaps/species/bos_taurus/h2a_bos_taurus_all.svg
- CURATED_SET/BioAnalyze/figures/heatmaps/species/pan_troglodytes/h2a_pan_troglodytes_all.png
- CURATED_SET/BioAnalyze/figures/heatmaps/species/pan_troglodytes/h2a_pan_troglodytes_all.svg
- CURATED_SET/BioAnalyze/figures/heatmaps/alligned_human_pan/hs_aligned_all.png
- CURATED_SET/BioAnalyze/figures/heatmaps/alligned_human_pan/hs_aligned_all.svg
- CURATED_SET/BioAnalyze/figures/heatmaps/alligned_human_pan/hs_aligned_clustered.png
- CURATED_SET/BioAnalyze/figures/heatmaps/alligned_human_pan/hs_aligned_clustered.svg
- CURATED_SET/BioAnalyze/figures/heatmaps/alligned_human_pan/hs_aligned_variants.png
- CURATED_SET/BioAnalyze/figures/heatmaps/alligned_human_pan/hs_aligned_variants.svg
- CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AJ/H2AJ_gene_compare_heatmap.png
- CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AJ/H2AJ_gene_compare_heatmap.svg
Recent Script Additions (2026-03-17)
- CURATED_SET/BioAnalyze/scripts/expression/summarize_shared_h2a_gene_names_across_species.py
  Rebuilds the shared-gene summary in `stats/shared_genes/` and writes a reusable long-form
  gene lookup index to `data/gene_compare/index/`.
- CURATED_SET/BioAnalyze/scripts/expression/build_gene_compare_heatmap.py
  Builds a cross-species per-gene heatmap from canonical map labels and
  `*_present_gold.tsv` inputs, with:
  - X axis = species
  - Y axis = tissues
  - default `union` tissue mode
  - mean aggregation per `(species, tissue)`
  - masked missing values
- CURATED_SET/BioAnalyze/scripts/expression/gene_compare_common.py
  Shared constants and helpers for detail index loading, per-gene long tables,
  and output layout across shared summaries, ranking, and plotting.
- CURATED_SET/BioAnalyze/scripts/expression/build_gene_compare_heatmaps_batch.py
  Batch-builds cross-species `gene_compare` heatmaps for all shared genes above
  a chosen `species_count` threshold.
- CURATED_SET/BioAnalyze/scripts/expression/rank_cross_species_h2a_differences.py
  Produces the gene*tissue comparison layer:
  - full `gene_tissue_summary`
  - global `p90` conservative candidates
  - global `p95` high-confidence candidates
  - pairwise contrasts for shortlisted candidates
  - manuscript-ready shortlist text
- CURATED_SET/BioAnalyze/scripts/expression/plot_cross_species_candidate_panels.py
  Builds:
  - focused candidate panels
  - overview scatter
  - class-specific `p95` panel pages for clustered vs variant genes

External raw data:
- stored at C:\Users\USER\Documents\GitHub\histonedb_external_storage\BioAnalyze\raw
See data/raw/README.md for the raw file list.

