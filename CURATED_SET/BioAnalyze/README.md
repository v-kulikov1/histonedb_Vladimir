# BioAnalyze

Location: CURATED_SET/BioAnalyze

Structure:
- scripts/annotate
- scripts/merge_taxonomy
- data/processed
- data/merged
- data/raw (placeholder; raw files are stored externally)
- audits
- stats
- figures/heatmaps
- figures/stats

Rules:
- All pipeline scripts must live under CURATED_SET/BioAnalyze/scripts (no scripts in data/ or elsewhere).

Canonical output:
- CURATED_SET/BioAnalyze/data/merged/mammalia_H2A_merged_with_taxonomy_v3.csv

Worklog:
- CURATED_SET/BioAnalyze/WORKLOG.md (detailed, chronological record of pipeline decisions and artifacts)

Archive:
- Legacy merged CSVs: CURATED_SET/BioAnalyze/archive/merged_legacy/
- Audit CSVs: CURATED_SET/BioAnalyze/archive/audits/
- Full process log: CURATED_SET/BioAnalyze/WORKLOG.md

External raw data:
- stored at C:\Users\USER\Documents\GitHub\histonedb_external_storage\BioAnalyze\raw
See data/raw/README.md for the raw file list.
