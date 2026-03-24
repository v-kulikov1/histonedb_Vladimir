# Подробный отчёт по cross-species сравнению экспрессии H2A на уровне gene*tissue

## Overview

- В анализ включены 21 shared H2A genes с `species_count >= 4`.
- Основной фильтр для интерпретации gene*tissue: `species_n >= 4`.
- После удаления generic tissues осталось 308 сравнимых gene*tissue сочетаний.
- Исключённые generic tissues: anatomical system, material anatomical entity, multicellular organism.
- Глобальные пороги: p05 = 3.58, p10 = 4.71, p90 = 48.02, p95 = 54.14.

## Class-level patterns

- `clustered`: 14 genes, 133 comparable gene*tissue rows, median range 36.63, mean range 36.83, p90 54.80, p95 62.03, max range 81.27.
- `variant`: 6 genes, 175 comparable gene*tissue rows, median range 10.01, mean range 14.23, p90 30.88, p95 37.50, max range 48.96.

| Class | Genes | Rows | Median range | Mean range | p90 | p95 | Max | Global p90 hits | Global p95 hits | Global p10 hits | Global p05 hits |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| clustered | 14 | 133 | 36.63 | 36.83 | 54.80 | 62.03 | 81.27 | 29 | 16 | 0 | 0 |
| variant | 6 | 175 | 10.01 | 14.23 | 30.88 | 37.50 | 48.96 | 2 | 0 | 31 | 16 |

## Most variable genes

- Рейтинг по максимальному наблюдаемому range:

| Gene | Class | Tissues | Max range | Median range | Mean range | p90 hits | p95 hits | Overall mean expr | Overall median expr | Heatmap |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| H2AC14 | clustered | 7 | 81.27 | 53.86 | 50.07 | 4 | 3 | 49.89 | 45.87 | [H2AC14 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC14/H2AC14_gene_compare_heatmap.png) |
| H2AC6 | clustered | 23 | 68.39 | 33.46 | 36.25 | 5 | 4 | 77.02 | 81.02 | [H2AC6 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC6/H2AC6_gene_compare_heatmap.png) |
| H2AC19 | clustered | 4 | 68.00 | 45.30 | 48.87 | 2 | 1 | 73.64 | 77.81 | [H2AC19 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC19/H2AC19_gene_compare_heatmap.png) |
| H2AC21 | clustered | 18 | 65.45 | 37.33 | 38.43 | 4 | 3 | 53.33 | 47.36 | [H2AC21 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC21/H2AC21_gene_compare_heatmap.png) |
| H2AC17 | clustered | 5 | 58.79 | 38.54 | 38.73 | 1 | 1 | 51.34 | 48.56 | [H2AC17 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC17/H2AC17_gene_compare_heatmap.png) |
| H2AC1 | clustered | 5 | 58.35 | 28.69 | 33.47 | 1 | 1 | 52.40 | 51.03 | [H2AC1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC1/H2AC1_gene_compare_heatmap.png) |
| H2AC20 | clustered | 4 | 56.28 | 51.86 | 51.94 | 3 | 2 | 63.64 | 63.72 | [H2AC20 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC20/H2AC20_gene_compare_heatmap.png) |
| H2AC11 | clustered | 10 | 54.67 | 26.20 | 27.20 | 1 | 1 | 51.66 | 48.90 | [H2AC11 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC11/H2AC11_gene_compare_heatmap.png) |
| H2AC8 | clustered | 12 | 54.02 | 35.77 | 34.15 | 1 | 0 | 59.96 | 59.42 | [H2AC8 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC8/H2AC8_gene_compare_heatmap.png) |
| H2AC13 | clustered | 15 | 53.12 | 34.11 | 33.55 | 2 | 0 | 54.12 | 53.02 | [H2AC13 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC13/H2AC13_gene_compare_heatmap.png) |

- Рейтинг по медианному range across tissues:

| Gene | Class | Tissues | Median range | Max range | Mean range | p90 hits | p95 hits | Overall mean expr | Overall median expr | Heatmap |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| H2AC14 | clustered | 7 | 53.86 | 81.27 | 50.07 | 4 | 3 | 49.89 | 45.87 | [H2AC14 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC14/H2AC14_gene_compare_heatmap.png) |
| H2AC20 | clustered | 4 | 51.86 | 56.28 | 51.94 | 3 | 2 | 63.64 | 63.72 | [H2AC20 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC20/H2AC20_gene_compare_heatmap.png) |
| H2AC19 | clustered | 4 | 45.30 | 68.00 | 48.87 | 2 | 1 | 73.64 | 77.81 | [H2AC19 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC19/H2AC19_gene_compare_heatmap.png) |
| H2AC25 | clustered | 18 | 43.34 | 53.07 | 40.75 | 4 | 0 | 64.11 | 63.50 | [H2AC25 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC25/H2AC25_gene_compare_heatmap.png) |
| H2AC17 | clustered | 5 | 38.54 | 58.79 | 38.73 | 1 | 1 | 51.34 | 48.56 | [H2AC17 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC17/H2AC17_gene_compare_heatmap.png) |
| H2AC21 | clustered | 18 | 37.33 | 65.45 | 38.43 | 4 | 3 | 53.33 | 47.36 | [H2AC21 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC21/H2AC21_gene_compare_heatmap.png) |
| H2AC12 | clustered | 3 | 35.81 | 48.00 | 30.77 | 0 | 0 | 48.24 | 48.32 | [H2AC12 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC12/H2AC12_gene_compare_heatmap.png) |
| H2AC8 | clustered | 12 | 35.77 | 54.02 | 34.15 | 1 | 0 | 59.96 | 59.42 | [H2AC8 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC8/H2AC8_gene_compare_heatmap.png) |
| H2AC13 | clustered | 15 | 34.11 | 53.12 | 33.55 | 2 | 0 | 54.12 | 53.02 | [H2AC13 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC13/H2AC13_gene_compare_heatmap.png) |
| H2AC6 | clustered | 23 | 33.46 | 68.39 | 36.25 | 5 | 4 | 77.02 | 81.02 | [H2AC6 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC6/H2AC6_gene_compare_heatmap.png) |

## Most conserved genes

| Gene | Class | Tissues | Median range | Max range | Mean range | p10 hits | p05 hits | Overall mean expr | Overall median expr | Heatmap |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| H2AZ1 | variant | 29 | 4.59 | 12.63 | 5.12 | 16 | 10 | 96.13 | 97.23 | [H2AZ1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AZ1/H2AZ1_gene_compare_heatmap.png) |
| H2AZ2 | variant | 29 | 7.21 | 16.93 | 8.43 | 4 | 1 | 94.55 | 95.96 | [H2AZ2 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AZ2/H2AZ2_gene_compare_heatmap.png) |
| MACROH2A1 | variant | 29 | 7.33 | 22.38 | 7.90 | 6 | 2 | 93.56 | 95.09 | [MACROH2A1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/MACROH2A1/MACROH2A1_gene_compare_heatmap.png) |
| H2AJ | variant | 28 | 12.88 | 41.82 | 15.27 | 3 | 2 | 89.98 | 93.00 | [H2AJ heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AJ/H2AJ_gene_compare_heatmap.png) |
| H2AC7 | clustered | 4 | 17.35 | 50.40 | 25.19 | 0 | 0 | 49.11 | 45.65 | [H2AC7 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC7/H2AC7_gene_compare_heatmap.png) |
| H2AX | variant | 36 | 19.01 | 48.96 | 20.43 | 2 | 1 | 82.28 | 82.82 | [H2AX heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AX/H2AX_gene_compare_heatmap.png) |
| H2AC11 | clustered | 10 | 26.20 | 54.67 | 27.20 | 0 | 0 | 51.66 | 48.90 | [H2AC11 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC11/H2AC11_gene_compare_heatmap.png) |
| H2AC1 | clustered | 5 | 28.69 | 58.35 | 33.47 | 0 | 0 | 52.40 | 51.03 | [H2AC1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC1/H2AC1_gene_compare_heatmap.png) |
| MACROH2A2 | variant | 24 | 29.06 | 48.07 | 29.38 | 0 | 0 | 71.16 | 74.43 | [MACROH2A2 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/MACROH2A2/MACROH2A2_gene_compare_heatmap.png) |
| H2AC15 | clustered | 5 | 30.57 | 41.95 | 29.41 | 0 | 0 | 50.91 | 48.02 | [H2AC15 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC15/H2AC15_gene_compare_heatmap.png) |

## Most variable gene*tissue combinations

- Global p95 cases:

| Gene | Class | Tissue | Species n | Range | Mean score | Median score | Min species | Max species | Heatmap | Panels |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| H2AC14 | clustered | testis | 5 | 81.27 | 51.52 | 58.61 | felis_catus 12.78 | sus_scrofa 94.05 | [H2AC14 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC14/H2AC14_gene_compare_heatmap.png) | [clustered_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p90_panels_page1.png), [clustered_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p95_panels_page1.png), [candidate_focus_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/candidate_focus_panels.png) |
| H2AC6 | clustered | prefrontal cortex | 6 | 68.39 | 57.10 | 65.97 | felis_catus 13.49 | human 81.88 | [H2AC6 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC6/H2AC6_gene_compare_heatmap.png) | [clustered_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p90_panels_page1.png), [clustered_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p95_panels_page1.png), [candidate_focus_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/candidate_focus_panels.png) |
| H2AC19 | clustered | bone marrow | 4 | 68.00 | 70.78 | 75.81 | equus_caballus 31.75 | pan_troglodytes 99.75 | [H2AC19 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC19/H2AC19_gene_compare_heatmap.png) | [clustered_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p90_panels_page1.png), [clustered_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p95_panels_page1.png) |
| H2AC21 | clustered | cerebral cortex | 5 | 65.45 | 46.88 | 37.75 | pan_troglodytes 27.24 | sus_scrofa 92.69 | [H2AC21 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC21/H2AC21_gene_compare_heatmap.png) | [clustered_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p90_panels_page1.png), [clustered_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p95_panels_page1.png) |
| H2AC6 | clustered | spleen | 6 | 65.41 | 72.26 | 75.36 | felis_catus 31.04 | pan_troglodytes 96.45 | [H2AC6 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC6/H2AC6_gene_compare_heatmap.png) | [clustered_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p90_panels_page1.png), [clustered_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p95_panels_page1.png) |
| H2AC6 | clustered | adult mammalian kidney | 7 | 64.90 | 65.32 | 67.43 | felis_catus 26.04 | pan_troglodytes 90.93 | [H2AC6 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC6/H2AC6_gene_compare_heatmap.png) | [clustered_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p90_panels_page1.png), [clustered_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p95_panels_page1.png) |
| H2AC21 | clustered | colon | 5 | 64.36 | 57.10 | 58.90 | macaca_mulatta 31.62 | sus_scrofa 95.98 | [H2AC21 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC21/H2AC21_gene_compare_heatmap.png) | [clustered_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p90_panels_page2.png), [clustered_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p95_panels_page2.png) |
| H2AC14 | clustered | brain | 5 | 60.47 | 36.29 | 32.52 | canis_lupus_familiaris 16.05 | felis_catus 76.52 | [H2AC14 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC14/H2AC14_gene_compare_heatmap.png) | [clustered_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p90_panels_page2.png) |
| H2AC14 | clustered | cerebral cortex | 4 | 59.20 | 43.82 | 34.77 | human 23.27 | sus_scrofa 82.47 | [H2AC14 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC14/H2AC14_gene_compare_heatmap.png) | [clustered_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p90_panels_page2.png) |
| H2AC17 | clustered | testis | 4 | 58.79 | 48.38 | 40.20 | felis_catus 27.17 | sus_scrofa 85.96 | [H2AC17 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC17/H2AC17_gene_compare_heatmap.png) | [clustered_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p90_panels_page2.png) |

- Variant-biased high-tail cases:

| Gene | Tissue | Species n | Range | Mean score | Median score | Min species | Max species | Heatmap | Panels |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| H2AX | adipose tissue | 5 | 48.96 | 70.53 | 74.88 | canis_lupus_familiaris 40.66 | sus_scrofa 89.62 | [H2AX heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AX/H2AX_gene_compare_heatmap.png) | [variant_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p90_panels_page1.png), [variant_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p95_panels_page1.png), [candidate_focus_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/candidate_focus_panels.png) |
| MACROH2A2 | zone of skin | 5 | 48.07 | 57.66 | 64.93 | felis_catus 30.50 | human 78.57 | [MACROH2A2 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/MACROH2A2/MACROH2A2_gene_compare_heatmap.png) | [variant_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p90_panels_page1.png), [variant_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p95_panels_page1.png), [candidate_focus_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/candidate_focus_panels.png) |
| MACROH2A2 | testis | 8 | 46.57 | 60.05 | 57.61 | felis_catus 32.73 | human 79.30 | [MACROH2A2 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/MACROH2A2/MACROH2A2_gene_compare_heatmap.png) | [variant_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p90_panels_page1.png), [variant_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p95_panels_page1.png) |
| H2AX | lymph node | 5 | 43.29 | 77.69 | 83.59 | canis_lupus_familiaris 46.83 | mus_musculus 90.12 | [H2AX heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AX/H2AX_gene_compare_heatmap.png) | [variant_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p90_panels_page1.png), [variant_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p95_panels_page1.png) |
| H2AJ | testis | 8 | 41.82 | 81.90 | 84.11 | felis_catus 57.99 | pan_troglodytes 99.81 | [H2AJ heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AJ/H2AJ_gene_compare_heatmap.png) | [variant_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p90_panels_page1.png), [variant_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p95_panels_page1.png), [candidate_focus_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/candidate_focus_panels.png) |
| MACROH2A2 | lung | 5 | 40.39 | 60.83 | 51.87 | macaca_mulatta 41.20 | canis_lupus_familiaris 81.59 | [MACROH2A2 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/MACROH2A2/MACROH2A2_gene_compare_heatmap.png) | [variant_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p90_panels_page1.png), [variant_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p95_panels_page1.png) |
| H2AJ | pituitary gland | 5 | 40.12 | 86.10 | 92.69 | pan_troglodytes 56.72 | canis_lupus_familiaris 96.84 | [H2AJ heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AJ/H2AJ_gene_compare_heatmap.png) | [variant_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p90_panels_page2.png), [variant_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p95_panels_page2.png), [candidate_focus_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/candidate_focus_panels.png) |
| H2AX | heart left ventricle | 4 | 39.68 | 74.45 | 72.80 | canis_lupus_familiaris 56.26 | mus_musculus 95.94 | [H2AX heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AX/H2AX_gene_compare_heatmap.png) | [variant_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p90_panels_page2.png), [variant_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p95_panels_page2.png) |
| H2AX | hypothalamus | 4 | 38.00 | 71.12 | 71.17 | canis_lupus_familiaris 52.07 | human 90.07 | [H2AX heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AX/H2AX_gene_compare_heatmap.png) | [variant_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p90_panels_page2.png), [variant_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p95_panels_page2.png) |
| MACROH2A2 | kidney | 5 | 37.29 | 65.92 | 64.16 | bos_taurus 50.04 | human 87.33 | [MACROH2A2 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/MACROH2A2/MACROH2A2_gene_compare_heatmap.png) | [variant_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p90_panels_page2.png) |

## Most conserved gene*tissue combinations

- Global low-tail `p05` cases:

| Gene | Class | Tissue | Species n | Range | Mean score | Median score | Min species | Max species | Heatmap | Panels |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| H2AZ1 | variant | thymus | 4 | 1.23 | 98.98 | 99.01 | pan_troglodytes 98.34 | human 99.57 | [H2AZ1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AZ1/H2AZ1_gene_compare_heatmap.png) | [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page1.png), [variant_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p05_low_panels_page1.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page1.png), [variant_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p10_low_panels_page1.png) |
| H2AJ | variant | metanephros cortex | 4 | 1.57 | 96.98 | 97.04 | bos_taurus 96.14 | human 97.71 | [H2AJ heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AJ/H2AJ_gene_compare_heatmap.png) | [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page1.png), [variant_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p05_low_panels_page1.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page1.png), [variant_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p10_low_panels_page1.png) |
| H2AZ1 | variant | frontal cortex | 6 | 2.04 | 96.97 | 97.26 | sus_scrofa 95.66 | pan_troglodytes 97.69 | [H2AZ1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AZ1/H2AZ1_gene_compare_heatmap.png) | [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page1.png), [variant_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p05_low_panels_page1.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page1.png), [variant_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p10_low_panels_page1.png) |
| H2AZ1 | variant | spleen | 7 | 2.46 | 97.29 | 96.93 | sus_scrofa 96.28 | human 98.74 | [H2AZ1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AZ1/H2AZ1_gene_compare_heatmap.png) | [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page1.png), [variant_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p05_low_panels_page1.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page1.png), [variant_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p10_low_panels_page1.png) |
| H2AZ1 | variant | blood | 5 | 2.49 | 96.88 | 96.47 | equus_caballus 95.70 | human 98.19 | [H2AZ1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AZ1/H2AZ1_gene_compare_heatmap.png) | [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page1.png), [variant_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p05_low_panels_page1.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page1.png), [variant_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p10_low_panels_page1.png) |
| H2AZ1 | variant | Ammon's horn | 4 | 2.62 | 96.84 | 97.26 | bos_taurus 95.12 | macaca_mulatta 97.74 | [H2AZ1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AZ1/H2AZ1_gene_compare_heatmap.png) | [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page1.png), [variant_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p05_low_panels_page1.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page1.png), [variant_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p10_low_panels_page1.png) |
| H2AZ1 | variant | cortex of kidney | 4 | 2.62 | 95.02 | 95.25 | bos_taurus 93.47 | macaca_mulatta 96.09 | [H2AZ1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AZ1/H2AZ1_gene_compare_heatmap.png) | [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page2.png), [variant_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p05_low_panels_page2.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page2.png), [variant_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p10_low_panels_page2.png) |
| MACROH2A1 | variant | retina | 4 | 2.80 | 95.52 | 95.25 | bos_taurus 94.39 | human 97.19 | [MACROH2A1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/MACROH2A1/MACROH2A1_gene_compare_heatmap.png) | [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page2.png), [variant_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p05_low_panels_page2.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page2.png), [variant_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p10_low_panels_page2.png) |
| H2AZ1 | variant | lymph node | 5 | 2.87 | 97.04 | 97.16 | pan_troglodytes 95.48 | human 98.35 | [H2AZ1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AZ1/H2AZ1_gene_compare_heatmap.png) | [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page2.png), [variant_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p05_low_panels_page2.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page2.png), [variant_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p10_low_panels_page2.png) |
| MACROH2A1 | variant | kidney | 5 | 2.90 | 94.69 | 94.65 | bos_taurus 92.97 | macaca_mulatta 95.87 | [MACROH2A1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/MACROH2A1/MACROH2A1_gene_compare_heatmap.png) | [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page2.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page2.png), [variant_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p10_low_panels_page2.png) |

- Global low-tail `p10` cases:

| Gene | Class | Tissue | Species n | Range | Mean score | Median score | Min species | Max species | Heatmap | Panels |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| H2AZ1 | variant | thymus | 4 | 1.23 | 98.98 | 99.01 | pan_troglodytes 98.34 | human 99.57 | [H2AZ1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AZ1/H2AZ1_gene_compare_heatmap.png) | [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page1.png), [variant_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p05_low_panels_page1.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page1.png), [variant_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p10_low_panels_page1.png) |
| H2AJ | variant | metanephros cortex | 4 | 1.57 | 96.98 | 97.04 | bos_taurus 96.14 | human 97.71 | [H2AJ heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AJ/H2AJ_gene_compare_heatmap.png) | [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page1.png), [variant_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p05_low_panels_page1.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page1.png), [variant_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p10_low_panels_page1.png) |
| H2AZ1 | variant | frontal cortex | 6 | 2.04 | 96.97 | 97.26 | sus_scrofa 95.66 | pan_troglodytes 97.69 | [H2AZ1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AZ1/H2AZ1_gene_compare_heatmap.png) | [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page1.png), [variant_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p05_low_panels_page1.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page1.png), [variant_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p10_low_panels_page1.png) |
| H2AZ1 | variant | spleen | 7 | 2.46 | 97.29 | 96.93 | sus_scrofa 96.28 | human 98.74 | [H2AZ1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AZ1/H2AZ1_gene_compare_heatmap.png) | [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page1.png), [variant_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p05_low_panels_page1.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page1.png), [variant_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p10_low_panels_page1.png) |
| H2AZ1 | variant | blood | 5 | 2.49 | 96.88 | 96.47 | equus_caballus 95.70 | human 98.19 | [H2AZ1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AZ1/H2AZ1_gene_compare_heatmap.png) | [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page1.png), [variant_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p05_low_panels_page1.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page1.png), [variant_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p10_low_panels_page1.png) |
| H2AZ1 | variant | Ammon's horn | 4 | 2.62 | 96.84 | 97.26 | bos_taurus 95.12 | macaca_mulatta 97.74 | [H2AZ1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AZ1/H2AZ1_gene_compare_heatmap.png) | [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page1.png), [variant_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p05_low_panels_page1.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page1.png), [variant_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p10_low_panels_page1.png) |
| H2AZ1 | variant | cortex of kidney | 4 | 2.62 | 95.02 | 95.25 | bos_taurus 93.47 | macaca_mulatta 96.09 | [H2AZ1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AZ1/H2AZ1_gene_compare_heatmap.png) | [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page2.png), [variant_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p05_low_panels_page2.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page2.png), [variant_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p10_low_panels_page2.png) |
| MACROH2A1 | variant | retina | 4 | 2.80 | 95.52 | 95.25 | bos_taurus 94.39 | human 97.19 | [MACROH2A1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/MACROH2A1/MACROH2A1_gene_compare_heatmap.png) | [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page2.png), [variant_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p05_low_panels_page2.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page2.png), [variant_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p10_low_panels_page2.png) |
| H2AZ1 | variant | lymph node | 5 | 2.87 | 97.04 | 97.16 | pan_troglodytes 95.48 | human 98.35 | [H2AZ1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AZ1/H2AZ1_gene_compare_heatmap.png) | [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page2.png), [variant_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p05_low_panels_page2.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page2.png), [variant_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p10_low_panels_page2.png) |
| MACROH2A1 | variant | kidney | 5 | 2.90 | 94.69 | 94.65 | bos_taurus 92.97 | macaca_mulatta 95.87 | [MACROH2A1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/MACROH2A1/MACROH2A1_gene_compare_heatmap.png) | [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page2.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page2.png), [variant_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p10_low_panels_page2.png) |

- Class-specific low-tail cases:

| Gene | Class | Class low | Tissue | Species n | Range | Mean score | Median score | Min species | Max species | Heatmap | Panels |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| H2AC13 | clustered | p05 | muscle tissue | 4 | 5.05 | 42.06 | 41.69 | equus_caballus 39.91 | bos_taurus 44.96 | [H2AC13 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC13/H2AC13_gene_compare_heatmap.png) | [clustered_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p05_low_panels_page1.png), [clustered_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p10_low_panels_page1.png) |
| H2AC11 | clustered | p05 | cerebral cortex | 4 | 7.19 | 45.98 | 45.34 | canis_lupus_familiaris 43.02 | human 50.21 | [H2AC11 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC11/H2AC11_gene_compare_heatmap.png) | [clustered_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p05_low_panels_page1.png), [clustered_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p10_low_panels_page1.png) |
| H2AC8 | clustered | p05 | cortex of kidney | 4 | 7.34 | 56.00 | 56.59 | macaca_mulatta 51.74 | human 59.08 | [H2AC8 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC8/H2AC8_gene_compare_heatmap.png) | [clustered_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p05_low_panels_page1.png), [clustered_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p10_low_panels_page1.png) |
| H2AC12 | clustered | p05 | liver | 4 | 8.51 | 40.92 | 40.88 | human 36.70 | macaca_mulatta 45.22 | [H2AC12 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC12/H2AC12_gene_compare_heatmap.png) | [clustered_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p05_low_panels_page1.png), [clustered_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p10_low_panels_page1.png) |
| H2AC21 | clustered | p05 | blood | 4 | 9.51 | 52.13 | 51.44 | equus_caballus 48.07 | human 57.58 | [H2AC21 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC21/H2AC21_gene_compare_heatmap.png) | [clustered_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p05_low_panels_page1.png), [clustered_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p10_low_panels_page1.png) |
| H2AC13 | clustered | p05 | colon | 4 | 12.30 | 55.85 | 57.07 | macaca_mulatta 48.48 | human 60.78 | [H2AC13 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC13/H2AC13_gene_compare_heatmap.png) | [clustered_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p05_low_panels_page1.png), [clustered_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p10_low_panels_page1.png) |
| H2AC6 | clustered | p05 | skeletal muscle tissue | 5 | 13.51 | 85.53 | 83.56 | human 81.73 | pan_troglodytes 95.24 | [H2AC6 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC6/H2AC6_gene_compare_heatmap.png) | [clustered_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p05_low_panels_page2.png), [clustered_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p10_low_panels_page2.png) |
| H2AC13 | clustered | p10 | muscle tissue | 4 | 5.05 | 42.06 | 41.69 | equus_caballus 39.91 | bos_taurus 44.96 | [H2AC13 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC13/H2AC13_gene_compare_heatmap.png) | [clustered_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p05_low_panels_page1.png), [clustered_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p10_low_panels_page1.png) |
| H2AC11 | clustered | p10 | cerebral cortex | 4 | 7.19 | 45.98 | 45.34 | canis_lupus_familiaris 43.02 | human 50.21 | [H2AC11 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC11/H2AC11_gene_compare_heatmap.png) | [clustered_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p05_low_panels_page1.png), [clustered_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p10_low_panels_page1.png) |
| H2AC8 | clustered | p10 | cortex of kidney | 4 | 7.34 | 56.00 | 56.59 | macaca_mulatta 51.74 | human 59.08 | [H2AC8 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC8/H2AC8_gene_compare_heatmap.png) | [clustered_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p05_low_panels_page1.png), [clustered_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p10_low_panels_page1.png) |

## Species tendencies

- Наиболее частые `min_species` в полном comparable set:

| Species | Hits | Share % | Subset rows |
| --- | --- | --- | --- |
| bos_taurus | 73 | 23.70 | 308 |
| canis_lupus_familiaris | 69 | 22.40 | 308 |
| felis_catus | 39 | 12.66 | 308 |
| pan_troglodytes | 35 | 11.36 | 308 |
| sus_scrofa | 26 | 8.44 | 308 |
| equus_caballus | 25 | 8.12 | 308 |
| human | 21 | 6.82 | 308 |
| macaca_mulatta | 16 | 5.19 | 308 |

- Наиболее частые `max_species` в полном comparable set:

| Species | Hits | Share % | Subset rows |
| --- | --- | --- | --- |
| human | 118 | 38.31 | 308 |
| sus_scrofa | 43 | 13.96 | 308 |
| bos_taurus | 33 | 10.71 | 308 |
| pan_troglodytes | 33 | 10.71 | 308 |
| macaca_mulatta | 24 | 7.79 | 308 |
| mus_musculus | 18 | 5.84 | 308 |
| equus_caballus | 14 | 4.55 | 308 |
| canis_lupus_familiaris | 13 | 4.22 | 308 |

- Наиболее частые `min_species` в high-tail p90:

| Species | Hits | Share % | Subset rows |
| --- | --- | --- | --- |
| canis_lupus_familiaris | 10 | 32.26 | 31 |
| felis_catus | 8 | 25.81 | 31 |
| equus_caballus | 4 | 12.90 | 31 |
| bos_taurus | 3 | 9.68 | 31 |
| human | 3 | 9.68 | 31 |
| macaca_mulatta | 2 | 6.45 | 31 |
| pan_troglodytes | 1 | 3.23 | 31 |
| mus_musculus | 0 | 0.00 | 31 |

- Наиболее частые `max_species` в high-tail p90:

| Species | Hits | Share % | Subset rows |
| --- | --- | --- | --- |
| sus_scrofa | 12 | 38.71 | 31 |
| equus_caballus | 5 | 16.13 | 31 |
| human | 5 | 16.13 | 31 |
| pan_troglodytes | 4 | 12.90 | 31 |
| macaca_mulatta | 3 | 9.68 | 31 |
| felis_catus | 2 | 6.45 | 31 |
| bos_taurus | 0 | 0.00 | 31 |
| canis_lupus_familiaris | 0 | 0.00 | 31 |

- Наиболее частые `min_species` в low-tail p10:

| Species | Hits | Share % | Subset rows |
| --- | --- | --- | --- |
| bos_taurus | 11 | 35.48 | 31 |
| sus_scrofa | 7 | 22.58 | 31 |
| pan_troglodytes | 6 | 19.35 | 31 |
| equus_caballus | 3 | 9.68 | 31 |
| mus_musculus | 2 | 6.45 | 31 |
| canis_lupus_familiaris | 1 | 3.23 | 31 |
| human | 1 | 3.23 | 31 |
| felis_catus | 0 | 0.00 | 31 |

- Наиболее частые `max_species` в low-tail p10:

| Species | Hits | Share % | Subset rows |
| --- | --- | --- | --- |
| human | 15 | 48.39 | 31 |
| macaca_mulatta | 7 | 22.58 | 31 |
| bos_taurus | 3 | 9.68 | 31 |
| pan_troglodytes | 3 | 9.68 | 31 |
| canis_lupus_familiaris | 2 | 6.45 | 31 |
| felis_catus | 1 | 3.23 | 31 |
| equus_caballus | 0 | 0.00 | 31 |
| mus_musculus | 0 | 0.00 | 31 |

- Наиболее частые пары `min_species -> max_species`:

  - `overall`:

| Pair | Count | Share % | Subset rows |
| --- | --- | --- | --- |
| bos_taurus -> human | 43 | 13.96 | 308 |
| felis_catus -> human | 19 | 6.17 | 308 |
| canis_lupus_familiaris -> human | 18 | 5.84 | 308 |
| sus_scrofa -> human | 17 | 5.52 | 308 |
| canis_lupus_familiaris -> sus_scrofa | 12 | 3.90 | 308 |
| canis_lupus_familiaris -> mus_musculus | 10 | 3.25 | 308 |
| canis_lupus_familiaris -> pan_troglodytes | 10 | 3.25 | 308 |
| equus_caballus -> human | 9 | 2.92 | 308 |
| pan_troglodytes -> human | 9 | 2.92 | 308 |
| bos_taurus -> macaca_mulatta | 8 | 2.60 | 308 |

  - `high_p90`:

| Pair | Count | Share % | Subset rows |
| --- | --- | --- | --- |
| canis_lupus_familiaris -> sus_scrofa | 5 | 16.13 | 31 |
| felis_catus -> human | 3 | 9.68 | 31 |
| canis_lupus_familiaris -> equus_caballus | 2 | 6.45 | 31 |
| equus_caballus -> macaca_mulatta | 2 | 6.45 | 31 |
| felis_catus -> pan_troglodytes | 2 | 6.45 | 31 |
| felis_catus -> sus_scrofa | 2 | 6.45 | 31 |
| human -> equus_caballus | 2 | 6.45 | 31 |
| macaca_mulatta -> sus_scrofa | 2 | 6.45 | 31 |
| bos_taurus -> equus_caballus | 1 | 3.23 | 31 |
| bos_taurus -> felis_catus | 1 | 3.23 | 31 |

  - `low_p10`:

| Pair | Count | Share % | Subset rows |
| --- | --- | --- | --- |
| bos_taurus -> human | 6 | 19.35 | 31 |
| bos_taurus -> macaca_mulatta | 5 | 16.13 | 31 |
| sus_scrofa -> human | 4 | 12.90 | 31 |
| pan_troglodytes -> human | 3 | 9.68 | 31 |
| equus_caballus -> human | 2 | 6.45 | 31 |
| mus_musculus -> pan_troglodytes | 2 | 6.45 | 31 |
| canis_lupus_familiaris -> bos_taurus | 1 | 3.23 | 31 |
| equus_caballus -> macaca_mulatta | 1 | 3.23 | 31 |
| human -> felis_catus | 1 | 3.23 | 31 |
| pan_troglodytes -> bos_taurus | 1 | 3.23 | 31 |

## Interpretation-ready findings

- В целом clustered-гены варьируют сильнее, чем variant-гены: медиана range 36.63 против 10.01, а p95 62.03 против 37.50.
- `H2AC14` имеет самый большой наблюдаемый cross-species range на уровне gene*tissue: максимальный range 81.27, медианный range по тканям 53.86.
- Среди variant-генов наиболее сильную вариабельность показывает `H2AX`: максимальный range 48.96, медианный range 19.01.
- Наиболее консервативный блок формируют `H2AZ1` (median range 4.59), `H2AZ2` (median range 7.21), `MACROH2A1` (median range 7.33); эти гены системно попадают в low-tail по межвидовому разбросу.
- Вариантный кейс `H2AX` / `adipose tissue`: range 48.96, среднее 70.53, медиана 74.88; минимум у canis_lupus_familiaris (40.66), максимум у sus_scrofa (89.62). Panels: [variant_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p90_panels_page1.png), [variant_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p95_panels_page1.png), [candidate_focus_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/candidate_focus_panels.png).
- Вариантный кейс `MACROH2A2` / `zone of skin`: range 48.07, среднее 57.66, медиана 64.93; минимум у felis_catus (30.50), максимум у human (78.57). Panels: [variant_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p90_panels_page1.png), [variant_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p95_panels_page1.png), [candidate_focus_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/candidate_focus_panels.png).
- Вариантный кейс `MACROH2A2` / `testis`: range 46.57, среднее 60.05, медиана 57.61; минимум у felis_catus (32.73), максимум у human (79.30). Panels: [variant_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p90_panels_page1.png), [variant_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p95_panels_page1.png).
- В high-tail (`p90`) нижние значения чаще всего дают canis_lupus_familiaris (10 hits), felis_catus (8 hits), тогда как верхние значения чаще наблюдаются у sus_scrofa (12 hits), equus_caballus (5 hits).
- Даже внутри clustered-класса есть сравнительно стабильные сочетания, например `H2AC13` / `muscle tissue` (class-specific p10, range 5.05).

## Caveats

- Все сопоставления тканей основаны на exact string matching по `Anatomical entity name`; синонимы и онтологическое слияние не применялись.
- Основной narrative строится только на сочетаниях с `species_n >= 4`, поэтому менее покрытые ткани в отчёт не включены как ключевые доказательства.
- Species-level выводы следует трактовать как статистические tendencies внутри текущего набора сравнимых тканей, а не как абсолютные свойства вида.
- Общие средние и медианы по гену рассчитаны по всем доступным `species × tissue` наблюдениям данного гена после удаления generic tissues.
