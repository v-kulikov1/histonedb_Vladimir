# BioAnalyze

From now on, any script-based work must be added under CURATED_SET/BioAnalyze/scripts so it is reproducible.

Location: CURATED_SET/BioAnalyze

Description:
BioAnalyze is a local data-processing pipeline for histone H2A datasets. It standardizes raw inputs, merges taxonomy, produces audits, and renders summary figures and stats used by the HistoneDB curation workflow.

Structure:
- archive/ (archived or legacy outputs and audits)
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

Canonical output:
- CURATED_SET/BioAnalyze/data/merged/mammalia_H2A_merged_with_taxonomy_v3.csv

Worklog:
- CURATED_SET/BioAnalyze/WORKLOG.md (detailed, chronological record of pipeline decisions and artifacts)

Archive:
- Root: CURATED_SET/BioAnalyze/archive
- Legacy merged CSVs: CURATED_SET/BioAnalyze/archive/merged_legacy/
- Audit CSVs: CURATED_SET/BioAnalyze/archive/audits/
- Full process log: CURATED_SET/BioAnalyze/WORKLOG.md

External raw data:
- stored at C:\Users\USER\Documents\GitHub\histonedb_external_storage\BioAnalyze\raw
See data/raw/README.md for the raw file list.

