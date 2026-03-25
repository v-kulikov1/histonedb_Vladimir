# Подробный отчёт по cross-species сравнению экспрессии H2A на уровне gene*tissue

## Overview

- В анализ включены 22 shared H2A genes с `species_count >= 4`.
- Основной фильтр для интерпретации gene*tissue: `species_n >= 4`.
- После удаления generic tissues осталось 459 сравнимых gene*tissue сочетаний.
- Исключённые generic tissues: anatomical system, material anatomical entity, multicellular organism.
- Глобальные пороги: p05 = 3.30, p10 = 5.04, p90 = 75.40, p95 = 82.53.

## Class-level patterns

- `clustered`: 15 genes, 271 comparable gene*tissue rows, median range 53.12, mean range 52.59, p90 80.29, p95 88.06, max range 99.63.
- `variant`: 7 genes, 188 comparable gene*tissue rows, median range 10.96, mean range 18.55, p90 42.11, p95 66.29, max range 86.69.

| Class | Genes | Rows | Median range | Mean range | p90 | p95 | Max | Global p90 hits | Global p95 hits | Global p10 hits | Global p05 hits |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| clustered | 15 | 271 | 53.12 | 52.59 | 80.29 | 88.06 | 99.63 | 39 | 20 | 9 | 9 |
| variant | 7 | 188 | 10.96 | 18.55 | 42.11 | 66.29 | 86.69 | 7 | 3 | 37 | 14 |

## Most variable genes

- Рейтинг по максимальному наблюдаемому range:

| Gene | Class | Tissues | Max range | Median range | Mean range | p90 hits | p95 hits | Overall mean expr | Overall median expr | Heatmap |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| H2AC14 | clustered | 24 | 99.63 | 64.12 | 66.70 | 9 | 5 | 29.22 | 31.12 | [H2AC14 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC14/H2AC14_gene_compare_heatmap.png) |
| H2AC8 | clustered | 16 | 99.38 | 54.83 | 53.22 | 2 | 2 | 38.35 | 48.60 | [H2AC8 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC8/H2AC8_gene_compare_heatmap.png) |
| H2AC7 | clustered | 10 | 96.09 | 54.23 | 56.44 | 2 | 1 | 35.38 | 39.74 | [H2AC7 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC7/H2AC7_gene_compare_heatmap.png) |
| H2AC11 | clustered | 23 | 94.42 | 52.85 | 51.26 | 2 | 1 | 27.39 | 32.48 | [H2AC11 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC11/H2AC11_gene_compare_heatmap.png) |
| H2AC20 | clustered | 8 | 94.00 | 81.75 | 77.86 | 5 | 3 | 41.31 | 52.03 | [H2AC20 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC20/H2AC20_gene_compare_heatmap.png) |
| H2AC21 | clustered | 20 | 93.47 | 41.11 | 47.86 | 4 | 2 | 43.27 | 44.22 | [H2AC21 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC21/H2AC21_gene_compare_heatmap.png) |
| H2AC19 | clustered | 6 | 92.61 | 59.19 | 60.58 | 1 | 1 | 69.59 | 76.40 | [H2AC19 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC19/H2AC19_gene_compare_heatmap.png) |
| H2AC25 | clustered | 23 | 92.15 | 51.14 | 55.21 | 4 | 1 | 58.13 | 60.18 | [H2AC25 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC25/H2AC25_gene_compare_heatmap.png) |
| MACROH2A2 | variant | 28 | 86.69 | 32.97 | 45.09 | 6 | 3 | 62.65 | 70.55 | [MACROH2A2 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/MACROH2A2/MACROH2A2_gene_compare_heatmap.png) |
| H2AC17 | clustered | 17 | 85.96 | 61.99 | 61.06 | 2 | 2 | 38.63 | 42.88 | [H2AC17 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC17/H2AC17_gene_compare_heatmap.png) |

- Рейтинг по медианному range across tissues:

| Gene | Class | Tissues | Median range | Max range | Mean range | p90 hits | p95 hits | Overall mean expr | Overall median expr | Heatmap |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| H2AC20 | clustered | 8 | 81.75 | 94.00 | 77.86 | 5 | 3 | 41.31 | 52.03 | [H2AC20 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC20/H2AC20_gene_compare_heatmap.png) |
| H2AC14 | clustered | 24 | 64.12 | 99.63 | 66.70 | 9 | 5 | 29.22 | 31.12 | [H2AC14 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC14/H2AC14_gene_compare_heatmap.png) |
| H2AC17 | clustered | 17 | 61.99 | 85.96 | 61.06 | 2 | 2 | 38.63 | 42.88 | [H2AC17 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC17/H2AC17_gene_compare_heatmap.png) |
| H2AC13 | clustered | 26 | 60.93 | 84.28 | 58.17 | 4 | 1 | 41.67 | 46.68 | [H2AC13 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC13/H2AC13_gene_compare_heatmap.png) |
| H2AC19 | clustered | 6 | 59.19 | 92.61 | 60.58 | 1 | 1 | 69.59 | 76.40 | [H2AC19 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC19/H2AC19_gene_compare_heatmap.png) |
| H2AC12 | clustered | 25 | 54.97 | 79.39 | 55.29 | 2 | 0 | 31.31 | 36.45 | [H2AC12 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC12/H2AC12_gene_compare_heatmap.png) |
| H2AC8 | clustered | 16 | 54.83 | 99.38 | 53.22 | 2 | 2 | 38.35 | 48.60 | [H2AC8 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC8/H2AC8_gene_compare_heatmap.png) |
| H2AC7 | clustered | 10 | 54.23 | 96.09 | 56.44 | 2 | 1 | 35.38 | 39.74 | [H2AC7 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC7/H2AC7_gene_compare_heatmap.png) |
| H2AC15 | clustered | 15 | 53.56 | 64.24 | 49.01 | 0 | 0 | 31.47 | 39.09 | [H2AC15 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC15/H2AC15_gene_compare_heatmap.png) |
| H2AC11 | clustered | 23 | 52.85 | 94.42 | 51.26 | 2 | 1 | 27.39 | 32.48 | [H2AC11 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC11/H2AC11_gene_compare_heatmap.png) |

## Most conserved genes

| Gene | Class | Tissues | Median range | Max range | Mean range | p10 hits | p05 hits | Overall mean expr | Overall median expr | Heatmap |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| H2AZ1 | variant | 29 | 4.59 | 12.63 | 5.12 | 18 | 8 | 96.13 | 97.23 | [H2AZ1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AZ1/H2AZ1_gene_compare_heatmap.png) |
| H2AZ2 | variant | 29 | 7.21 | 16.93 | 8.43 | 6 | 1 | 93.79 | 95.95 | [H2AZ2 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AZ2/H2AZ2_gene_compare_heatmap.png) |
| MACROH2A1 | variant | 29 | 7.33 | 22.38 | 7.90 | 7 | 2 | 93.26 | 95.06 | [MACROH2A1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/MACROH2A1/MACROH2A1_gene_compare_heatmap.png) |
| H2AJ | variant | 28 | 12.88 | 41.82 | 15.27 | 3 | 2 | 84.67 | 92.37 | [H2AJ heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AJ/H2AJ_gene_compare_heatmap.png) |
| H2AX | variant | 36 | 19.01 | 48.96 | 20.43 | 3 | 1 | 74.75 | 81.89 | [H2AX heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AX/H2AX_gene_compare_heatmap.png) |
| MACROH2A2 | variant | 28 | 32.97 | 86.69 | 45.09 | 0 | 0 | 62.65 | 70.55 | [MACROH2A2 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/MACROH2A2/MACROH2A2_gene_compare_heatmap.png) |
| H2AC6 | clustered | 23 | 33.46 | 84.57 | 38.48 | 0 | 0 | 72.00 | 79.62 | [H2AC6 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC6/H2AC6_gene_compare_heatmap.png) |
| H2AC21 | clustered | 20 | 41.11 | 93.47 | 47.86 | 0 | 0 | 43.27 | 44.22 | [H2AC21 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC21/H2AC21_gene_compare_heatmap.png) |
| H2AC1 | clustered | 30 | 46.73 | 65.29 | 35.92 | 8 | 8 | 9.97 | 0.00 | [H2AC1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC1/H2AC1_gene_compare_heatmap.png) |
| H2AC16 | clustered | 5 | 46.96 | 50.00 | 42.56 | 0 | 0 | 23.39 | 25.56 | [H2AC16 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC16/H2AC16_gene_compare_heatmap.png) |

## Most variable gene*tissue combinations

- Global p95 cases:

| Gene | Class | Tissue | Species n | Range | Mean score | Median score | Min species | Max species | Heatmap | Panels |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| H2AC14 | clustered | bone marrow | 4 | 99.63 | 57.87 | 65.93 | canis_lupus_familiaris 0.00 | pan_troglodytes 99.63 | [H2AC14 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC14/H2AC14_gene_compare_heatmap.png) | [clustered_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p90_panels_page1.png), [clustered_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p95_panels_page1.png), [candidate_focus_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/candidate_focus_panels.png) |
| H2AC8 | clustered | bone marrow | 4 | 99.38 | 64.63 | 79.57 | bos_taurus 0.00 | pan_troglodytes 99.38 | [H2AC8 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC8/H2AC8_gene_compare_heatmap.png) | [clustered_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p90_panels_page1.png), [clustered_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p95_panels_page1.png), [candidate_focus_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/candidate_focus_panels.png) |
| H2AC14 | clustered | colon | 6 | 96.48 | 39.76 | 43.33 | canis_lupus_familiaris 0.00 | sus_scrofa 96.48 | [H2AC14 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC14/H2AC14_gene_compare_heatmap.png) | [clustered_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p90_panels_page1.png), [clustered_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p95_panels_page1.png) |
| H2AC7 | clustered | brain | 4 | 96.09 | 41.90 | 35.75 | canis_lupus_familiaris 0.00 | felis_catus 96.09 | [H2AC7 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC7/H2AC7_gene_compare_heatmap.png) | [clustered_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p90_panels_page1.png), [clustered_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p95_panels_page1.png) |
| H2AC8 | clustered | cerebellum | 5 | 94.95 | 51.80 | 53.33 | bos_taurus 0.00 | equus_caballus 94.95 | [H2AC8 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC8/H2AC8_gene_compare_heatmap.png) | [clustered_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p90_panels_page1.png), [clustered_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p95_panels_page1.png) |
| H2AC11 | clustered | bone marrow | 4 | 94.42 | 56.10 | 65.00 | bos_taurus 0.00 | equus_caballus 94.42 | [H2AC11 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC11/H2AC11_gene_compare_heatmap.png) | [clustered_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p90_panels_page1.png), [clustered_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p95_panels_page1.png) |
| H2AC14 | clustered | testis | 7 | 94.05 | 36.80 | 31.56 | canis_lupus_familiaris 0.00 | sus_scrofa 94.05 | [H2AC14 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC14/H2AC14_gene_compare_heatmap.png) | [clustered_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p90_panels_page2.png), [clustered_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p95_panels_page2.png) |
| H2AC20 | clustered | brain | 4 | 94.00 | 53.34 | 59.68 | bos_taurus 0.00 | felis_catus 94.00 | [H2AC20 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC20/H2AC20_gene_compare_heatmap.png) | [clustered_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p90_panels_page2.png), [clustered_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p95_panels_page2.png) |
| H2AC21 | clustered | cerebellum | 6 | 93.47 | 51.27 | 43.36 | bos_taurus 0.00 | equus_caballus 93.47 | [H2AC21 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC21/H2AC21_gene_compare_heatmap.png) | [clustered_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p90_panels_page2.png), [clustered_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p95_panels_page2.png) |
| H2AC19 | clustered | adult mammalian kidney | 4 | 92.61 | 53.81 | 61.32 | equus_caballus 0.00 | pan_troglodytes 92.61 | [H2AC19 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC19/H2AC19_gene_compare_heatmap.png) | [clustered_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p90_panels_page2.png), [clustered_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p95_panels_page2.png) |

- Variant-biased high-tail cases:

| Gene | Tissue | Species n | Range | Mean score | Median score | Min species | Max species | Heatmap | Panels |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| MACROH2A2 | thymus | 4 | 86.69 | 54.36 | 65.37 | pan_troglodytes 0.00 | human 86.69 | [MACROH2A2 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/MACROH2A2/MACROH2A2_gene_compare_heatmap.png) | [variant_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p90_panels_page1.png), [variant_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p95_panels_page1.png) |
| MACROH2A2 | metanephros cortex | 4 | 83.97 | 48.51 | 55.05 | sus_scrofa 0.00 | human 83.97 | [MACROH2A2 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/MACROH2A2/MACROH2A2_gene_compare_heatmap.png) | [variant_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p90_panels_page1.png), [variant_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p95_panels_page1.png) |
| MACROH2A2 | liver | 8 | 82.71 | 44.09 | 49.96 | equus_caballus 0.00 | human 82.71 | [MACROH2A2 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/MACROH2A2/MACROH2A2_gene_compare_heatmap.png) | [variant_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p90_panels_page1.png), [variant_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p95_panels_page1.png) |
| MACROH2A2 | colon | 6 | 82.14 | 50.29 | 55.12 | sus_scrofa 0.00 | human 82.14 | [MACROH2A2 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/MACROH2A2/MACROH2A2_gene_compare_heatmap.png) | [variant_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p90_panels_page1.png), [variant_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p95_panels_page1.png) |
| MACROH2A2 | lung | 6 | 81.59 | 50.69 | 50.93 | pan_troglodytes 0.00 | canis_lupus_familiaris 81.59 | [MACROH2A2 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/MACROH2A2/MACROH2A2_gene_compare_heatmap.png) | [variant_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p90_panels_page1.png), [variant_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p95_panels_page1.png) |
| MACROH2A2 | lymph node | 5 | 80.63 | 37.05 | 51.94 | pan_troglodytes 0.00 | human 80.63 | [MACROH2A2 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/MACROH2A2/MACROH2A2_gene_compare_heatmap.png) | [variant_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p90_panels_page1.png), [variant_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p95_panels_page1.png) |
| H2AB3 | testis | 5 | 76.44 | 44.12 | 46.43 | canis_lupus_familiaris 0.00 | mus_musculus 76.44 | [H2AB3 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AB3/H2AB3_gene_compare_heatmap.png) | [variant_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p90_panels_page2.png), [variant_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p95_panels_page2.png) |
| MACROH2A2 | bone marrow | 5 | 74.94 | 51.33 | 60.51 | pan_troglodytes 0.00 | human 74.94 | [MACROH2A2 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/MACROH2A2/MACROH2A2_gene_compare_heatmap.png) | [variant_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p90_panels_page2.png), [variant_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p95_panels_page2.png) |
| MACROH2A2 | spleen | 7 | 71.57 | 37.17 | 46.22 | macaca_mulatta 0.00 | human 71.57 | [MACROH2A2 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/MACROH2A2/MACROH2A2_gene_compare_heatmap.png) | [variant_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p90_panels_page2.png), [variant_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p95_panels_page2.png) |
| MACROH2A2 | blood | 5 | 70.86 | 32.21 | 40.28 | canis_lupus_familiaris 0.00 | human 70.86 | [MACROH2A2 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/MACROH2A2/MACROH2A2_gene_compare_heatmap.png) | [variant_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p90_panels_page2.png), [variant_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p95_panels_page2.png) |

## Most conserved gene*tissue combinations

- Global low-tail `p05` cases:

| Gene | Class | Tissue | Species n | Range | Mean score | Median score | Min species | Max species | Heatmap | Panels |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| H2AC1 | clustered | frontal cortex | 6 | 0.00 | 0.00 | 0.00 | bos_taurus 0.00 | bos_taurus 0.00 | [H2AC1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC1/H2AC1_gene_compare_heatmap.png) | [clustered_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p05_low_panels_page1.png), [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page1.png), [clustered_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p10_low_panels_page1.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page1.png) |
| H2AC1 | clustered | skeletal muscle tissue | 6 | 0.00 | 0.00 | 0.00 | bos_taurus 0.00 | bos_taurus 0.00 | [H2AC1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC1/H2AC1_gene_compare_heatmap.png) | [clustered_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p05_low_panels_page1.png), [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page1.png), [clustered_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p10_low_panels_page1.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page1.png) |
| H2AC1 | clustered | cerebral cortex | 5 | 0.00 | 0.00 | 0.00 | bos_taurus 0.00 | bos_taurus 0.00 | [H2AC1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC1/H2AC1_gene_compare_heatmap.png) | [clustered_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p05_low_panels_page1.png), [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page1.png), [clustered_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p10_low_panels_page1.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page1.png) |
| H2AC1 | clustered | adipose tissue | 4 | 0.00 | 0.00 | 0.00 | bos_taurus 0.00 | bos_taurus 0.00 | [H2AC1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC1/H2AC1_gene_compare_heatmap.png) | [clustered_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p05_low_panels_page1.png), [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page1.png), [clustered_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p10_low_panels_page1.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page1.png) |
| H2AC1 | clustered | hypothalamus | 4 | 0.00 | 0.00 | 0.00 | bos_taurus 0.00 | bos_taurus 0.00 | [H2AC1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC1/H2AC1_gene_compare_heatmap.png) | [clustered_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p05_low_panels_page1.png), [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page1.png), [clustered_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p10_low_panels_page1.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page1.png) |
| H2AC1 | clustered | lymph node | 4 | 0.00 | 0.00 | 0.00 | bos_taurus 0.00 | bos_taurus 0.00 | [H2AC1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC1/H2AC1_gene_compare_heatmap.png) | [clustered_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p05_low_panels_page1.png), [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page1.png), [clustered_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p10_low_panels_page1.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page1.png) |
| H2AC1 | clustered | ovary | 4 | 0.00 | 0.00 | 0.00 | bos_taurus 0.00 | bos_taurus 0.00 | [H2AC1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC1/H2AC1_gene_compare_heatmap.png) | [clustered_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p05_low_panels_page2.png), [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page2.png), [clustered_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p10_low_panels_page2.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page2.png) |
| H2AC1 | clustered | pituitary gland | 4 | 0.00 | 0.00 | 0.00 | bos_taurus 0.00 | bos_taurus 0.00 | [H2AC1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC1/H2AC1_gene_compare_heatmap.png) | [clustered_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p05_low_panels_page2.png), [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page2.png), [clustered_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p10_low_panels_page2.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page2.png) |
| H2AC11 | clustered | metanephros cortex | 4 | 0.00 | 0.00 | 0.00 | bos_taurus 0.00 | bos_taurus 0.00 | [H2AC11 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC11/H2AC11_gene_compare_heatmap.png) | [clustered_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p05_low_panels_page2.png), [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page2.png), [clustered_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p10_low_panels_page2.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page2.png) |
| H2AZ1 | variant | thymus | 4 | 1.23 | 98.98 | 99.01 | pan_troglodytes 98.34 | human 99.57 | [H2AZ1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AZ1/H2AZ1_gene_compare_heatmap.png) | [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page2.png), [variant_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p05_low_panels_page1.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page2.png), [variant_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p10_low_panels_page1.png) |

- Global low-tail `p10` cases:

| Gene | Class | Tissue | Species n | Range | Mean score | Median score | Min species | Max species | Heatmap | Panels |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| H2AC1 | clustered | frontal cortex | 6 | 0.00 | 0.00 | 0.00 | bos_taurus 0.00 | bos_taurus 0.00 | [H2AC1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC1/H2AC1_gene_compare_heatmap.png) | [clustered_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p05_low_panels_page1.png), [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page1.png), [clustered_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p10_low_panels_page1.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page1.png) |
| H2AC1 | clustered | skeletal muscle tissue | 6 | 0.00 | 0.00 | 0.00 | bos_taurus 0.00 | bos_taurus 0.00 | [H2AC1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC1/H2AC1_gene_compare_heatmap.png) | [clustered_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p05_low_panels_page1.png), [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page1.png), [clustered_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p10_low_panels_page1.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page1.png) |
| H2AC1 | clustered | cerebral cortex | 5 | 0.00 | 0.00 | 0.00 | bos_taurus 0.00 | bos_taurus 0.00 | [H2AC1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC1/H2AC1_gene_compare_heatmap.png) | [clustered_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p05_low_panels_page1.png), [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page1.png), [clustered_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p10_low_panels_page1.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page1.png) |
| H2AC1 | clustered | adipose tissue | 4 | 0.00 | 0.00 | 0.00 | bos_taurus 0.00 | bos_taurus 0.00 | [H2AC1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC1/H2AC1_gene_compare_heatmap.png) | [clustered_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p05_low_panels_page1.png), [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page1.png), [clustered_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p10_low_panels_page1.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page1.png) |
| H2AC1 | clustered | hypothalamus | 4 | 0.00 | 0.00 | 0.00 | bos_taurus 0.00 | bos_taurus 0.00 | [H2AC1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC1/H2AC1_gene_compare_heatmap.png) | [clustered_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p05_low_panels_page1.png), [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page1.png), [clustered_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p10_low_panels_page1.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page1.png) |
| H2AC1 | clustered | lymph node | 4 | 0.00 | 0.00 | 0.00 | bos_taurus 0.00 | bos_taurus 0.00 | [H2AC1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC1/H2AC1_gene_compare_heatmap.png) | [clustered_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p05_low_panels_page1.png), [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page1.png), [clustered_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p10_low_panels_page1.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page1.png) |
| H2AC1 | clustered | ovary | 4 | 0.00 | 0.00 | 0.00 | bos_taurus 0.00 | bos_taurus 0.00 | [H2AC1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC1/H2AC1_gene_compare_heatmap.png) | [clustered_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p05_low_panels_page2.png), [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page2.png), [clustered_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p10_low_panels_page2.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page2.png) |
| H2AC1 | clustered | pituitary gland | 4 | 0.00 | 0.00 | 0.00 | bos_taurus 0.00 | bos_taurus 0.00 | [H2AC1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC1/H2AC1_gene_compare_heatmap.png) | [clustered_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p05_low_panels_page2.png), [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page2.png), [clustered_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p10_low_panels_page2.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page2.png) |
| H2AC11 | clustered | metanephros cortex | 4 | 0.00 | 0.00 | 0.00 | bos_taurus 0.00 | bos_taurus 0.00 | [H2AC11 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC11/H2AC11_gene_compare_heatmap.png) | [clustered_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p05_low_panels_page2.png), [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page2.png), [clustered_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p10_low_panels_page2.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page2.png) |
| H2AZ1 | variant | thymus | 4 | 1.23 | 98.98 | 99.01 | pan_troglodytes 98.34 | human 99.57 | [H2AZ1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AZ1/H2AZ1_gene_compare_heatmap.png) | [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page2.png), [variant_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p05_low_panels_page1.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page2.png), [variant_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p10_low_panels_page1.png) |

- Class-specific low-tail cases:

| Gene | Class | Class low | Tissue | Species n | Range | Mean score | Median score | Min species | Max species | Heatmap | Panels |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| H2AC1 | clustered | p05 | frontal cortex | 6 | 0.00 | 0.00 | 0.00 | bos_taurus 0.00 | bos_taurus 0.00 | [H2AC1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC1/H2AC1_gene_compare_heatmap.png) | [clustered_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p05_low_panels_page1.png), [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page1.png), [clustered_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p10_low_panels_page1.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page1.png) |
| H2AC1 | clustered | p05 | skeletal muscle tissue | 6 | 0.00 | 0.00 | 0.00 | bos_taurus 0.00 | bos_taurus 0.00 | [H2AC1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC1/H2AC1_gene_compare_heatmap.png) | [clustered_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p05_low_panels_page1.png), [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page1.png), [clustered_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p10_low_panels_page1.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page1.png) |
| H2AC1 | clustered | p05 | cerebral cortex | 5 | 0.00 | 0.00 | 0.00 | bos_taurus 0.00 | bos_taurus 0.00 | [H2AC1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC1/H2AC1_gene_compare_heatmap.png) | [clustered_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p05_low_panels_page1.png), [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page1.png), [clustered_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p10_low_panels_page1.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page1.png) |
| H2AC1 | clustered | p05 | adipose tissue | 4 | 0.00 | 0.00 | 0.00 | bos_taurus 0.00 | bos_taurus 0.00 | [H2AC1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC1/H2AC1_gene_compare_heatmap.png) | [clustered_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p05_low_panels_page1.png), [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page1.png), [clustered_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p10_low_panels_page1.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page1.png) |
| H2AC1 | clustered | p05 | hypothalamus | 4 | 0.00 | 0.00 | 0.00 | bos_taurus 0.00 | bos_taurus 0.00 | [H2AC1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC1/H2AC1_gene_compare_heatmap.png) | [clustered_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p05_low_panels_page1.png), [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page1.png), [clustered_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p10_low_panels_page1.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page1.png) |
| H2AC1 | clustered | p05 | lymph node | 4 | 0.00 | 0.00 | 0.00 | bos_taurus 0.00 | bos_taurus 0.00 | [H2AC1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC1/H2AC1_gene_compare_heatmap.png) | [clustered_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p05_low_panels_page1.png), [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page1.png), [clustered_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p10_low_panels_page1.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page1.png) |
| H2AC1 | clustered | p05 | ovary | 4 | 0.00 | 0.00 | 0.00 | bos_taurus 0.00 | bos_taurus 0.00 | [H2AC1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC1/H2AC1_gene_compare_heatmap.png) | [clustered_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p05_low_panels_page2.png), [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page2.png), [clustered_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p10_low_panels_page2.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page2.png) |
| H2AC1 | clustered | p05 | pituitary gland | 4 | 0.00 | 0.00 | 0.00 | bos_taurus 0.00 | bos_taurus 0.00 | [H2AC1 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC1/H2AC1_gene_compare_heatmap.png) | [clustered_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p05_low_panels_page2.png), [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page2.png), [clustered_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p10_low_panels_page2.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page2.png) |
| H2AC11 | clustered | p05 | metanephros cortex | 4 | 0.00 | 0.00 | 0.00 | bos_taurus 0.00 | bos_taurus 0.00 | [H2AC11 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC11/H2AC11_gene_compare_heatmap.png) | [clustered_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p05_low_panels_page2.png), [global_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p05_low_panels_page2.png), [clustered_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p10_low_panels_page2.png), [global_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/global_p10_low_panels_page2.png) |
| H2AC13 | clustered | p05 | muscle tissue | 4 | 5.05 | 42.06 | 41.69 | equus_caballus 39.91 | bos_taurus 44.96 | [H2AC13 heatmap](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/figures/heatmaps/gene_compare/H2AC13/H2AC13_gene_compare_heatmap.png) | [clustered_p05_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p05_low_panels_page2.png), [clustered_p10_low_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/clustered_p10_low_panels_page2.png) |

## Species tendencies

- Наиболее частые `min_species` в полном comparable set:

| Species | Hits | Share % | Subset rows |
| --- | --- | --- | --- |
| bos_taurus | 130 | 28.32 | 459 |
| canis_lupus_familiaris | 112 | 24.40 | 459 |
| equus_caballus | 41 | 8.93 | 459 |
| macaca_mulatta | 40 | 8.71 | 459 |
| pan_troglodytes | 39 | 8.50 | 459 |
| sus_scrofa | 34 | 7.41 | 459 |
| felis_catus | 33 | 7.19 | 459 |
| human | 25 | 5.45 | 459 |

- Наиболее частые `max_species` в полном comparable set:

| Species | Hits | Share % | Subset rows |
| --- | --- | --- | --- |
| human | 164 | 35.73 | 459 |
| bos_taurus | 67 | 14.60 | 459 |
| sus_scrofa | 55 | 11.98 | 459 |
| mus_musculus | 41 | 8.93 | 459 |
| pan_troglodytes | 39 | 8.50 | 459 |
| macaca_mulatta | 29 | 6.32 | 459 |
| canis_lupus_familiaris | 26 | 5.66 | 459 |
| equus_caballus | 22 | 4.79 | 459 |

- Наиболее частые `min_species` в high-tail p90:

| Species | Hits | Share % | Subset rows |
| --- | --- | --- | --- |
| canis_lupus_familiaris | 13 | 28.26 | 46 |
| bos_taurus | 11 | 23.91 | 46 |
| equus_caballus | 5 | 10.87 | 46 |
| macaca_mulatta | 5 | 10.87 | 46 |
| pan_troglodytes | 4 | 8.70 | 46 |
| felis_catus | 3 | 6.52 | 46 |
| sus_scrofa | 3 | 6.52 | 46 |
| human | 2 | 4.35 | 46 |

- Наиболее частые `max_species` в high-tail p90:

| Species | Hits | Share % | Subset rows |
| --- | --- | --- | --- |
| human | 13 | 28.26 | 46 |
| sus_scrofa | 12 | 26.09 | 46 |
| equus_caballus | 7 | 15.22 | 46 |
| pan_troglodytes | 4 | 8.70 | 46 |
| felis_catus | 3 | 6.52 | 46 |
| macaca_mulatta | 3 | 6.52 | 46 |
| mus_musculus | 3 | 6.52 | 46 |
| canis_lupus_familiaris | 1 | 2.17 | 46 |

- Наиболее частые `min_species` в low-tail p10:

| Species | Hits | Share % | Subset rows |
| --- | --- | --- | --- |
| bos_taurus | 23 | 50.00 | 46 |
| sus_scrofa | 7 | 15.22 | 46 |
| pan_troglodytes | 6 | 13.04 | 46 |
| canis_lupus_familiaris | 3 | 6.52 | 46 |
| equus_caballus | 3 | 6.52 | 46 |
| mus_musculus | 2 | 4.35 | 46 |
| human | 1 | 2.17 | 46 |
| macaca_mulatta | 1 | 2.17 | 46 |

- Наиболее частые `max_species` в low-tail p10:

| Species | Hits | Share % | Subset rows |
| --- | --- | --- | --- |
| human | 19 | 41.30 | 46 |
| bos_taurus | 13 | 28.26 | 46 |
| macaca_mulatta | 7 | 15.22 | 46 |
| pan_troglodytes | 4 | 8.70 | 46 |
| canis_lupus_familiaris | 2 | 4.35 | 46 |
| felis_catus | 1 | 2.17 | 46 |
| equus_caballus | 0 | 0.00 | 46 |
| mus_musculus | 0 | 0.00 | 46 |

- Наиболее частые пары `min_species -> max_species`:

  - `overall`:

| Pair | Count | Share % | Subset rows |
| --- | --- | --- | --- |
| bos_taurus -> human | 58 | 12.64 | 459 |
| canis_lupus_familiaris -> human | 28 | 6.10 | 459 |
| canis_lupus_familiaris -> bos_taurus | 25 | 5.45 | 459 |
| canis_lupus_familiaris -> sus_scrofa | 21 | 4.58 | 459 |
| sus_scrofa -> human | 21 | 4.58 | 459 |
| equus_caballus -> human | 16 | 3.49 | 459 |
| felis_catus -> human | 16 | 3.49 | 459 |
| bos_taurus -> macaca_mulatta | 15 | 3.27 | 459 |
| bos_taurus -> mus_musculus | 13 | 2.83 | 459 |
| canis_lupus_familiaris -> mus_musculus | 13 | 2.83 | 459 |

  - `high_p90`:

| Pair | Count | Share % | Subset rows |
| --- | --- | --- | --- |
| canis_lupus_familiaris -> sus_scrofa | 7 | 15.22 | 46 |
| bos_taurus -> equus_caballus | 6 | 13.04 | 46 |
| bos_taurus -> macaca_mulatta | 3 | 6.52 | 46 |
| macaca_mulatta -> human | 3 | 6.52 | 46 |
| canis_lupus_familiaris -> human | 2 | 4.35 | 46 |
| equus_caballus -> human | 2 | 4.35 | 46 |
| felis_catus -> human | 2 | 4.35 | 46 |
| pan_troglodytes -> human | 2 | 4.35 | 46 |
| sus_scrofa -> human | 2 | 4.35 | 46 |
| bos_taurus -> felis_catus | 1 | 2.17 | 46 |

  - `low_p10`:

| Pair | Count | Share % | Subset rows |
| --- | --- | --- | --- |
| bos_taurus -> bos_taurus | 9 | 19.57 | 46 |
| bos_taurus -> human | 9 | 19.57 | 46 |
| bos_taurus -> macaca_mulatta | 5 | 10.87 | 46 |
| sus_scrofa -> human | 4 | 8.70 | 46 |
| pan_troglodytes -> human | 3 | 6.52 | 46 |
| canis_lupus_familiaris -> bos_taurus | 2 | 4.35 | 46 |
| equus_caballus -> human | 2 | 4.35 | 46 |
| mus_musculus -> pan_troglodytes | 2 | 4.35 | 46 |
| canis_lupus_familiaris -> human | 1 | 2.17 | 46 |
| equus_caballus -> macaca_mulatta | 1 | 2.17 | 46 |

## Interpretation-ready findings

- В целом clustered-гены варьируют сильнее, чем variant-гены: медиана range 53.12 против 10.96, а p95 88.06 против 66.29.
- `H2AC14` имеет самый большой наблюдаемый cross-species range на уровне gene*tissue: максимальный range 99.63, медианный range по тканям 64.12.
- Среди variant-генов наиболее сильную вариабельность показывает `MACROH2A2`: максимальный range 86.69, медианный range 32.97.
- Наиболее консервативный блок формируют `H2AZ1` (median range 4.59), `H2AZ2` (median range 7.21), `MACROH2A1` (median range 7.33); эти гены системно попадают в low-tail по межвидовому разбросу.
- Вариантный кейс `MACROH2A2` / `thymus`: range 86.69, среднее 54.36, медиана 65.37; минимум у pan_troglodytes (0.00), максимум у human (86.69). Panels: [variant_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p90_panels_page1.png), [variant_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p95_panels_page1.png).
- Вариантный кейс `MACROH2A2` / `metanephros cortex`: range 83.97, среднее 48.51, медиана 55.05; минимум у sus_scrofa (0.00), максимум у human (83.97). Panels: [variant_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p90_panels_page1.png), [variant_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p95_panels_page1.png).
- Вариантный кейс `MACROH2A2` / `liver`: range 82.71, среднее 44.09, медиана 49.96; минимум у equus_caballus (0.00), максимум у human (82.71). Panels: [variant_p90_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p90_panels_page1.png), [variant_p95_panels](C:/Users/USER/Documents/GitHub/histonedb/CURATED_SET/BioAnalyze/stats/ranking/plots/variant_p95_panels_page1.png).
- В high-tail (`p90`) нижние значения чаще всего дают canis_lupus_familiaris (13 hits), bos_taurus (11 hits), тогда как верхние значения чаще наблюдаются у human (13 hits), sus_scrofa (12 hits).
- Даже внутри clustered-класса есть сравнительно стабильные сочетания, например `H2AC1` / `frontal cortex` (class-specific p10, range 0.00).

## Caveats

- Все сопоставления тканей основаны на exact string matching по `Anatomical entity name`; синонимы и онтологическое слияние не применялись.
- Основной narrative строится только на сочетаниях с `species_n >= 4`, поэтому менее покрытые ткани в отчёт не включены как ключевые доказательства.
- Species-level выводы следует трактовать как статистические tendencies внутри текущего набора сравнимых тканей, а не как абсолютные свойства вида.
- Общие средние и медианы по гену рассчитаны по всем доступным `species × tissue` наблюдениям данного гена после удаления generic tissues.
