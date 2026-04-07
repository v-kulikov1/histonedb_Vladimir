# BioAnalyze

From now on, any script-based work must be added under CURATED_SET/BioAnalyze/scripts so it is reproducible.

Location: CURATED_SET/BioAnalyze

Description:
BioAnalyze is a local data-processing pipeline for histone H2A datasets. It standardizes raw inputs, merges taxonomy, adds canonical gene names in `v4`, builds human HPA `nTPM` summaries and HPA-vs-Bgee compare layers, and renders audits, figures, and stats used by the HistoneDB curation workflow.

Structure:
- audits/ (pipeline QA and validation tables)
- audits/h2aj_tree (tree-reconstruction audit narrative and evidence interpretation)
- data/raw (placeholder; raw files are stored externally)
- data/expression_nTPM (human HPA-derived H2A expression tables and metadata)
- data/processed (intermediate normalized tables)
- data/gene_compare (gene-level comparison indexes and per-gene exported tables)
- data/merged (final merged outputs)
- data/h2aj_tree (historical and clean H2A.J tree-reconstruction inputs plus evidence)
- figures/heatmaps (heatmap graphics root)
- figures/heatmaps/species (per-species H2A heatmaps)
- figures/heatmaps/gene_compare (cross-species per-gene heatmaps)
- figures/heatmaps/compare_nTPM_bgee (human HPA-vs-Bgee compare heatmaps and legends)
- figures/expression_nTPM (human HPA `nTPM` summary, per-gene, and per-tissue plots)
- figures/stats (summary stats graphics)
- scripts/annotate (annotation utilities)
- scripts/merge_taxonomy (taxonomy merge pipeline)
- scripts/expression (expression, summary, and heatmap builders)
- scripts/h2aj_tree (historical H2A.J tree reconstruction workflow)
- stats/shared_genes (shared-gene summary tables and plots)
- stats/gene_tissue/ranking/tables (gene*tissue comparison tables)
- stats/gene_tissue/ranking/reports (manuscript-oriented text outputs)
- stats/gene_tissue/ranking/plots (gene*tissue comparison plots and panel pages)
- stats/gene_tissue/barplot (cross-species gene*anatomical_entity mean+std barplots, legends, and support TSVs)
- stats/accession_stats (legacy accession summary tables)
- stats/compare_nTPM_bgee (human HPA-vs-Bgee paired-cell stats and per-gene correlations)
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
- CURATED_SET/BioAnalyze/stats/gene_tissue/ranking/tables/gene_tissue_summary.csv
- CURATED_SET/BioAnalyze/stats/gene_tissue/ranking/tables/conservative_candidates.csv
- CURATED_SET/BioAnalyze/stats/gene_tissue/ranking/tables/high_confidence_candidates.csv
- CURATED_SET/BioAnalyze/stats/gene_tissue/ranking/tables/pairwise_contrasts.csv
- CURATED_SET/BioAnalyze/stats/gene_tissue/ranking/reports/manuscript_shortlist.md
- CURATED_SET/BioAnalyze/stats/gene_tissue/ranking/plots/candidate_focus_panels.png
- CURATED_SET/BioAnalyze/stats/gene_tissue/ranking/plots/clustered_p95_panels_page1.png
- CURATED_SET/BioAnalyze/stats/gene_tissue/ranking/plots/variant_p95_panels_page1.png
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
Recent Outputs (2026-03-29)
- CURATED_SET/BioAnalyze/data/h2aj_tree/evidence/evidence_summary.md
- CURATED_SET/BioAnalyze/data/h2aj_tree/historical/
- CURATED_SET/BioAnalyze/data/h2aj_tree/clean/
- CURATED_SET/BioAnalyze/audits/h2aj_tree/h2aj_tree_reconstruction_audit.md
Recent Outputs (2026-04-01)
- CURATED_SET/BioAnalyze/data/h2aj_tree/clean/nuc/nuc_0907_alignment.phy
- CURATED_SET/BioAnalyze/data/h2aj_tree/clean/aa/aa_1006_mammalian_cH2A_plus_nonplacental_alignment.phy
- CURATED_SET/BioAnalyze/data/h2aj_tree/clean/postprocess/README_filtered_clean_note.txt
- CURATED_SET/BioAnalyze/audits/h2aj_tree/h2aj_tree_reconstruction_audit.md
Recent Outputs (2026-04-01)
- CURATED_SET/BioAnalyze/data/expression_nTPM/human/h2a_human_gene_ntpm_cells.tsv
- CURATED_SET/BioAnalyze/data/expression_nTPM/human/h2a_human_gene_ntpm_summary.tsv
- CURATED_SET/BioAnalyze/data/expression_nTPM/human/h2a_human_gene_ntpm_metadata.json
- CURATED_SET/BioAnalyze/data/processed/intersections/human_hpa_bgee/hpa_vs_bgee_gene_mapping.tsv
- CURATED_SET/BioAnalyze/data/processed/intersections/human_hpa_bgee/hpa_vs_bgee_tissue_mapping.tsv
- CURATED_SET/BioAnalyze/data/processed/intersections/human_hpa_bgee/hpa_vs_bgee_hpa_aligned_long.tsv
- CURATED_SET/BioAnalyze/data/processed/intersections/human_hpa_bgee/hpa_vs_bgee_bgee_aligned_long.tsv
- CURATED_SET/BioAnalyze/data/processed/intersections/human_hpa_bgee/hpa_vs_bgee_metadata.json
- CURATED_SET/BioAnalyze/figures/expression_nTPM/human/summary/h2a_human_gene_count_by_tissue.png
- CURATED_SET/BioAnalyze/figures/expression_nTPM/human/summary/h2a_human_gene_ntpm_mean_ranking.svg
- CURATED_SET/BioAnalyze/figures/expression_nTPM/human/per_gene/H2AC25_tissue_ntpm_barplot.png
- CURATED_SET/BioAnalyze/figures/expression_nTPM/human/per_tissue/kidney_gene_ntpm_ranked_barplot.svg
- CURATED_SET/BioAnalyze/figures/expression_nTPM/human/special_cases/h2a5_split_gene_boxplot.png
- CURATED_SET/BioAnalyze/figures/heatmaps/compare_nTPM_bgee/hpa_vs_bgee_hpa_heatmap.png
- CURATED_SET/BioAnalyze/figures/heatmaps/compare_nTPM_bgee/hpa_vs_bgee_bgee_heatmap.png
- CURATED_SET/BioAnalyze/figures/heatmaps/compare_nTPM_bgee/hpa_vs_bgee_hpa_raw_heatmap.svg
- CURATED_SET/BioAnalyze/figures/heatmaps/compare_nTPM_bgee/hpa_vs_bgee_bgee_raw_heatmap.svg
- CURATED_SET/BioAnalyze/figures/heatmaps/compare_nTPM_bgee/hpa_vs_bgee_gene_number_legend.png
- CURATED_SET/BioAnalyze/figures/heatmaps/compare_nTPM_bgee/hpa_vs_bgee_hpa_tissue_letter_map.tsv
- CURATED_SET/BioAnalyze/figures/heatmaps/compare_nTPM_bgee/hpa_vs_bgee_bgee_tissue_letter_map.tsv
- CURATED_SET/BioAnalyze/stats/compare_nTPM_bgee/paired_cells.tsv
- CURATED_SET/BioAnalyze/stats/compare_nTPM_bgee/correlation_summary.tsv
- CURATED_SET/BioAnalyze/stats/compare_nTPM_bgee/report_md.md
- CURATED_SET/BioAnalyze/stats/compare_nTPM_bgee/per_gene_correlations/gene_correlation_keep_zeros.tsv
- CURATED_SET/BioAnalyze/stats/compare_nTPM_bgee/per_gene_correlations/pearson_raw_keep_zeros.png
- CURATED_SET/BioAnalyze/stats/compare_nTPM_bgee/per_gene_correlations/report_md.md
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
Recent Script Additions (2026-03-29)
- CURATED_SET/BioAnalyze/scripts/h2aj_tree/rebuild_h2aj_tree_history.py
  Reconstructs the historical H2A.J tree workflow, archives source evidence,
  and rebuilds the `historical`, `clean`, and `evidence` outputs under
  `data/h2aj_tree/`.
Recent Script Changes (2026-04-01)
- CURATED_SET/BioAnalyze/scripts/h2aj_tree/rebuild_h2aj_tree_history.py
  Adds repeatable `--drop-clean-id` filtering for `clean` outputs, switches the
  July clean SQK source to the `without short` FASTA pairs, and leaves an
  explicit note that `clean/postprocess/*.nwk` remain historical reference
  trees rather than trees inferred from the newly filtered alignments.
Recent Script Additions (2026-04-01)
- CURATED_SET/BioAnalyze/scripts/expression/build_human_h2a_ntpm_gene_expression.py
  Builds the human HPA-derived H2A `nTPM` gene-by-tissue tables, summary metadata,
  ranked tissue plots, per-gene barplots, and special-case outputs under
  `data/expression_nTPM/` and `figures/expression_nTPM/`.
- CURATED_SET/BioAnalyze/scripts/expression/build_human_hpa_bgee_compare_heatmaps.py
  Aligns the human HPA `nTPM` layer with processed human Bgee cH2A data,
  exports mapping tables and aligned long tables, and renders both log-scale
  and raw-scale compare heatmaps with shared legends under
  `figures/heatmaps/compare_nTPM_bgee/`.
- CURATED_SET/BioAnalyze/scripts/expression/build_human_hpa_bgee_compare_stats.py
  Consumes the aligned HPA-vs-Bgee tables and writes paired-cell summaries,
  correlation reports, and top-difference tables to `stats/compare_nTPM_bgee/`.
- CURATED_SET/BioAnalyze/scripts/expression/build_human_hpa_bgee_gene_correlations.py
  Builds per-gene Pearson/Spearman correlation tables and plots from
  `stats/compare_nTPM_bgee/paired_cells.tsv`, including keep-zero and
  drop-zero variants.
- CURATED_SET/BioAnalyze/scripts/bioanalyze_paths.py
  Centralizes repo-local and BioAnalyze-root path helpers so expression scripts
  can resolve defaults without hard-coding `CURATED_SET/BioAnalyze/...`.
Recent Script Changes (2026-04-06)
- CURATED_SET/BioAnalyze/scripts/expression/build_cross_species_gene_tissue_boxplots.py
  Now builds fixed-species cross-species `gene:anatomical_entity` mean+std
  barplots for `k=4/5/6/7`, keeps one combined figure plus split
  `by_class/clustered` and `by_class/variant` figures per `k`, writes separate
  gene/anatomical-entity legends, and saves outputs under
  `stats/gene_tissue/barplot/`.
- CURATED_SET/BioAnalyze/scripts/expression/gene_compare_common.py
  Centralizes the cross-species class palette and default `stats/gene_tissue`
  output roots, including the current `barplot` destination.
- CURATED_SET/BioAnalyze/scripts/bioanalyze_paths.py
  Provides repo-local BioAnalyze path resolution so active scripts no longer
  rely on hard-coded `CURATED_SET/BioAnalyze/...` defaults.
Recent Outputs (2026-04-06)
- CURATED_SET/BioAnalyze/stats/gene_tissue/ranking/tables/gene_tissue_summary.csv
- CURATED_SET/BioAnalyze/stats/gene_tissue/ranking/tables/conservative_candidates.csv
- CURATED_SET/BioAnalyze/stats/gene_tissue/ranking/tables/high_confidence_candidates.csv
- CURATED_SET/BioAnalyze/stats/gene_tissue/ranking/plots/panel_membership.csv
- CURATED_SET/BioAnalyze/stats/gene_tissue/barplot/tables/cross_species_gene_tissue_barplot_coverage_summary.tsv
- CURATED_SET/BioAnalyze/stats/gene_tissue/barplot/k4/cross_species_gene_tissue_barplot_k4.png
- CURATED_SET/BioAnalyze/stats/gene_tissue/barplot/k5/cross_species_gene_tissue_barplot_k5.png
- CURATED_SET/BioAnalyze/stats/gene_tissue/barplot/k5/by_class/clustered/cross_species_gene_tissue_barplot_clustered_k5.png
- CURATED_SET/BioAnalyze/stats/gene_tissue/barplot/k5/by_class/variant/cross_species_gene_tissue_barplot_variant_k5.png
- CURATED_SET/BioAnalyze/stats/gene_tissue/barplot/k6/cross_species_gene_tissue_barplot_k6.png
- CURATED_SET/BioAnalyze/stats/gene_tissue/barplot/k7/cross_species_gene_tissue_barplot_k7.png
Recent Outputs (2026-04-07)
- CURATED_SET/BioAnalyze/stats/gene_tissue/barplot/k5/by_class/variant/without_genes_6_7/cross_species_gene_tissue_barplot_variant_k5_without_genes_6_7.png
- CURATED_SET/BioAnalyze/stats/gene_tissue/barplot/k5/by_class/variant/without_genes_6_7/cross_species_gene_tissue_barplot_variant_k5_without_genes_6_7_gene_legend.png
- CURATED_SET/BioAnalyze/stats/gene_tissue/barplot/k5/by_class/variant/without_genes_6_7/cross_species_gene_tissue_barplot_variant_k5_without_genes_6_7_stats.tsv

External raw data:
- resolved from `HISTONEDB_EXTERNAL_STORAGE\BioAnalyze\raw` when the env var is set
- otherwise resolved from the sibling folder next to the repo:
  `..\histonedb_external_storage\BioAnalyze\raw`
See data/raw/README.md for the raw file list.

