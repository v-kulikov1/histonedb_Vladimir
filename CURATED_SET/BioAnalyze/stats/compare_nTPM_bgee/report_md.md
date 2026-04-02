# HPA vs Bgee Compare Statistics

Generated: 2026-04-02 13:49:13

- Compare figures: `CURATED_SET/BioAnalyze/figures/heatmaps/compare_nTPM_bgee`
- HPA compare table: `CURATED_SET/BioAnalyze/data/processed/intersections/human_hpa_bgee/hpa_vs_bgee_hpa_aligned_long.tsv`
- Bgee compare table: `CURATED_SET/BioAnalyze/data/processed/intersections/human_hpa_bgee/hpa_vs_bgee_bgee_aligned_long.tsv`
- Paired stats table: `CURATED_SET\BioAnalyze\stats\compare_nTPM_bgee\paired_cells.tsv`

## Overview

- Total compare slots: 595
- Non-missing paired cells: 477
- Nonzero-any paired cells: 417
- Primary comparison space: `log10(value + 1)`
- HPA compare value: `nTPM`
- Bgee compare value: `cell_mean_score`

## Correlation Summary

| Subset | Pairs | Pearson log | Spearman log | Pearson raw | Spearman raw |
| --- | --- | --- | --- | --- | --- |
| all_non_missing | 477 | 0.438963 | 0.813329 | 0.538538 | 0.813329 |
| nonzero_any | 417 | 0.456579 | 0.761460 | 0.600948 | 0.761460 |

## Most Different Cells

| Gene | HPA tissue | Bgee tissue | Match | HPA nTPM | Bgee score | Abs log diff | Dominance |
| --- | --- | --- | --- | --- | --- | --- | --- |
| H2AC14 | kidney | kidney | exact | 0.0000 | 81.1500 | 1.914608 | bgee_higher |
| H2AC12 | kidney | kidney | exact | 0.0000 | 79.3925 | 1.905216 | bgee_higher |
| H2AC21 | kidney | kidney | exact | 0.1000 | 83.6400 | 1.886183 | bgee_higher |
| H2AC13 | kidney | kidney | exact | 0.1000 | 82.9350 | 1.882550 | bgee_higher |
| H2AC16 | kidney | kidney | exact | 0.1000 | 79.4550 | 1.864160 | bgee_higher |
| H2AC16 | prostate | prostate gland | safe_synonym | 0.0000 | 71.2290 | 1.858712 | bgee_higher |
| H2AC17 | kidney | kidney | exact | 0.2000 | 82.8900 | 1.844529 | bgee_higher |
| H2AC11 | ovary | ovary | exact | 0.1000 | 75.2060 | 1.840596 | bgee_higher |
| H2AC4 | kidney | kidney | exact | 0.0000 | 64.1633 | 1.814003 | bgee_higher |
| H2AC4 | ovary | ovary | exact | 0.0000 | 59.9900 | 1.785259 | bgee_higher |

Safe-synonym cases among the strongest differences:

- `H2AC16`: `prostate` vs `prostate gland` (`abs_diff_log=1.858712`)

## Most Different HPA-Higher Cells

| Gene | HPA tissue | Bgee tissue | HPA nTPM | Bgee score | Abs log diff |
| --- | --- | --- | --- | --- | --- |
| H2AC8 | heart muscle | myocardium | 4.0000 | 0.0000 | 0.698970 |
| H2AC15 | skeletal muscle | skeletal muscle tissue | 0.8000 | 0.0000 | 0.255273 |
| H2AC11 | heart muscle | myocardium | 0.2000 | 0.0000 | 0.079181 |
| H2AC8 | skeletal muscle | skeletal muscle tissue | 0.2000 | 0.0000 | 0.079181 |
| H2AC11 | salivary gland | saliva-secreting gland | 0.1000 | 0.0000 | 0.041393 |
| H2AC11 | skeletal muscle | skeletal muscle tissue | 0.1000 | 0.0000 | 0.041393 |
| H2AC12 | amygdala | amygdala | 0.1000 | 0.0000 | 0.041393 |
| H2AC12 | salivary gland | saliva-secreting gland | 0.1000 | 0.0000 | 0.041393 |

## Most Different Bgee-Higher Cells

| Gene | HPA tissue | Bgee tissue | HPA nTPM | Bgee score | Abs log diff |
| --- | --- | --- | --- | --- | --- |
| H2AC14 | kidney | kidney | 0.0000 | 81.1500 | 1.914608 |
| H2AC12 | kidney | kidney | 0.0000 | 79.3925 | 1.905216 |
| H2AC21 | kidney | kidney | 0.1000 | 83.6400 | 1.886183 |
| H2AC13 | kidney | kidney | 0.1000 | 82.9350 | 1.882550 |
| H2AC16 | kidney | kidney | 0.1000 | 79.4550 | 1.864160 |
| H2AC16 | prostate | prostate gland | 0.0000 | 71.2290 | 1.858712 |
| H2AC17 | kidney | kidney | 0.2000 | 82.8900 | 1.844529 |
| H2AC11 | ovary | ovary | 0.1000 | 75.2060 | 1.840596 |

## Genes With The Strongest Overall Disagreement

| Gene | Cells | Mean abs log diff | Median abs log diff | Max abs log diff | Top HPA tissue | Top Bgee tissue |
| --- | --- | --- | --- | --- | --- | --- |
| H2AC21 | 31 | 1.519783 | 1.532489 | 1.886183 | kidney | kidney |
| H2AC17 | 31 | 1.446667 | 1.515681 | 1.844529 | kidney | kidney |
| H2AC13 | 32 | 1.441386 | 1.514521 | 1.882550 | kidney | kidney |
| H2AC15 | 32 | 1.437339 | 1.479854 | 1.649682 | stomach | stomach |
| H2AC11 | 34 | 1.408089 | 1.577136 | 1.840596 | ovary | ovary |
| H2AC7 | 31 | 1.336441 | 1.433770 | 1.726994 | kidney | kidney |
| H2AC16 | 31 | 1.309800 | 1.488879 | 1.864160 | kidney | kidney |
| H2AC20 | 31 | 1.301625 | 1.295039 | 1.527478 | skeletal muscle | skeletal muscle tissue |
| H2AC12 | 31 | 1.291816 | 1.472491 | 1.905216 | kidney | kidney |
| H2AC14 | 31 | 1.133629 | 1.445760 | 1.914608 | kidney | kidney |

## Tissues With The Strongest Overall Disagreement

| Compare tissue | HPA tissue | Bgee tissue | Match | Mean abs log diff | Median abs log diff | Max abs log diff | Top gene |
| --- | --- | --- | --- | --- | --- | --- | --- |
| kidney | kidney | kidney | exact | 1.571468 | 1.730513 | 1.914608 | H2AC14 |
| stomach | stomach | stomach | exact | 1.442605 | 1.686721 | 1.772305 | H2AC21 |
| lung | lung | lung | exact | 1.435840 | 1.620164 | 1.736852 | H2AC14 |
| liver | liver | liver | exact | 1.421963 | 1.535006 | 1.640775 | H2AC16 |
| colon | colon | colon | exact | 1.352172 | 1.579767 | 1.698278 | H2AC21 |
| substantia nigra | substantia nigra | substantia nigra | exact | 1.328380 | 1.525955 | 1.602792 | H2AC16 |
| adrenal gland | adrenal gland | adrenal gland | exact | 1.317950 | 1.555094 | 1.666433 | H2AC13 |
| cerebral cortex | cerebral cortex | cerebral cortex | exact | 1.295358 | 1.439476 | 1.667990 | H2AC11 |
| caudate nucleus | caudate | caudate nucleus | safe_synonym | 1.293776 | 1.470851 | 1.591482 | H2AC21 |
| uterine cervix | cervix | uterine cervix | safe_synonym | 1.255276 | 1.439575 | 1.625032 | H2AC12 |

## Interpretation Note

- Correlation and difference ranking are computed in the same space shown on the heatmaps: `log10(value + 1)`.
- Raw HPA `nTPM` and raw Bgee `cell_mean_score` are reported next to the log-scale differences so large disagreements remain biologically interpretable.
