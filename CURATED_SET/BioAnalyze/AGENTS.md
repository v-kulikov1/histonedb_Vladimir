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
Recent Outputs (2026-03-24)
- CURATED_SET/BioAnalyze/scripts/h2aj_synteny/build_h2aj_synteny_plot.py
- CURATED_SET/BioAnalyze/figures/h2aj_synteny/h2aj_synteny.png
- CURATED_SET/BioAnalyze/figures/h2aj_synteny/h2aj_synteny.svg
- CURATED_SET/BioAnalyze/scripts/codons/build_codon_heatmaps.py
- CURATED_SET/BioAnalyze/figures/codons/sqk_nuc_without_short_codon_entropy_annotated.png
- CURATED_SET/BioAnalyze/figures/codons/sqk_nuc_without_short_codon_entropy_annotated.svg
- CURATED_SET/BioAnalyze/figures/codons/sqk_nuc_full_codon_entropy_annotated.png
- CURATED_SET/BioAnalyze/figures/codons/sqk_nuc_full_codon_entropy_annotated.svg
- CURATED_SET/BioAnalyze/figures/codons/sqk_nuc_full_order_legend.png
- CURATED_SET/BioAnalyze/figures/codons/sqk_nuc_full_order_legend.svg
- CURATED_SET/BioAnalyze/figures/codons/sqk_nuc_full_shannon_scale.png
- CURATED_SET/BioAnalyze/figures/codons/sqk_nuc_full_shannon_scale.svg
Recent Outputs (2026-03-25)
- CURATED_SET/BioAnalyze/audits/bgee_h2a_expression_source_audit_3species.tsv
- CURATED_SET/BioAnalyze/stats/ranking/reports/bgee_h2a_expression_source_audit_3species.md
- CURATED_SET/BioAnalyze/figures/heatmaps/species/human/h2a_all.png
- CURATED_SET/BioAnalyze/figures/heatmaps/species/human/h2a_all.svg
- CURATED_SET/BioAnalyze/figures/heatmaps/species/human/coverage_ge70/h2a_all.png
- CURATED_SET/BioAnalyze/figures/heatmaps/species/human/coverage_ge70/h2a_all.svg
- CURATED_SET/BioAnalyze/figures/heatmaps/species/mus_musculus/coverage_ge70/h2a_all.png
- CURATED_SET/BioAnalyze/figures/heatmaps/species/mus_musculus/coverage_ge70/h2a_all.svg
- CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AB3/H2AB3_gene_compare_heatmap.png
- CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AB3/H2AB3_gene_compare_heatmap.svg
- CURATED_SET/BioAnalyze/figures/heatmaps/presentation_human_pan/hs_aligned_all_presentation.png
- CURATED_SET/BioAnalyze/figures/heatmaps/presentation_human_pan/hs_aligned_all_presentation.svg
- CURATED_SET/BioAnalyze/figures/heatmaps/presentation_human_pan/h2a_pan_troglodytes_all_presentation.png
- CURATED_SET/BioAnalyze/figures/heatmaps/presentation_human_pan/h2a_pan_troglodytes_all_presentation.svg
- CURATED_SET/BioAnalyze/figures/boxplot/human/h2a_human_boxplot.png
- CURATED_SET/BioAnalyze/figures/boxplot/human/h2a_human_boxplot.svg
- CURATED_SET/BioAnalyze/figures/boxplot/human/h2a_human_boxplot_presentation.png
- CURATED_SET/BioAnalyze/figures/boxplot/human/h2a_human_boxplot_presentation.svg
- CURATED_SET/BioAnalyze/figures/boxplot/human/h2a_human_boxplot_number_legend.png
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
Recent Script Additions (2026-03-24)
- CURATED_SET/BioAnalyze/scripts/h2aj_synteny/build_h2aj_synteny_plot.py
  Builds the presentation-oriented H2A.J synteny panel from external TSV inputs
  and writes PNG/SVG outputs under `figures/h2aj_synteny/`.
- CURATED_SET/BioAnalyze/scripts/codons/build_codon_heatmaps.py
  Rebuilds SQK codon entropy heatmaps from external FASTA pairs, writes the
  `without-short` and `full` heatmaps, and exports standalone order/Shannon
  support graphics for the `full` presentation layout.
Recent Script Additions (2026-03-25)
- CURATED_SET/BioAnalyze/scripts/expression/normalized_expression_common.py
  Centralizes normalized Gene x tissue helpers used by species heatmaps,
  coverage filtering, aligned human/pan views, and downstream expression tables.
- CURATED_SET/BioAnalyze/scripts/expression/audit_bgee_h2a_expression_sources.py
  Audits which Bgee source-specific tracks underlie heatmap-relevant H2A
  expression rows and exports both TSV and Markdown summaries.
- CURATED_SET/BioAnalyze/scripts/expression/build_human_pan_presentation_heatmaps.py
  Builds presentation-oriented aligned human/pan H2A heatmaps together with
  shared legends and number-mapping support files.
- CURATED_SET/BioAnalyze/scripts/expression/build_human_h2a_boxplot.py
  Builds the analytical and presentation versions of the human H2A boxplot and
  exports the number legend plus supporting TSV tables.

External raw data:
- stored at C:\Users\USER\Documents\GitHub\histonedb_external_storage\BioAnalyze\raw
See data/raw/README.md for the raw file list.

