# Per-Gene HPA vs Bgee Correlations

Generated: 2026-04-02 13:49:43

- Source paired table: `CURATED_SET/BioAnalyze/stats/compare_nTPM_bgee/paired_cells.tsv`
- Compare figures: `CURATED_SET/BioAnalyze/figures/heatmaps/compare_nTPM_bgee`

## Modes

- `keep_zeros`: all non-missing paired tissues are kept; zeroes are preserved.
- `drop_zeros`: only tissues with both `nTPM > 0` and `cell_mean_score > 0` are used.
- Correlations are raw-only here: `pearson_raw` and `spearman_raw`.

## keep_zeros

### Pearson raw: strongest

| Gene | Bgee label | Tissues | Pearson raw | Status |
| --- | --- | --- | --- | --- |
| H2AC1 | H2AC1 | 30 | 0.758929 | ok |
| H2AC25 | H2AW | 31 | 0.694110 | ok |
| H2AC20 | H2AC20 | 31 | 0.667709 | ok |
| H2AC6 | H2AC6 | 35 | 0.609999 | ok |
| H2AC7 | H2AC7 | 31 | 0.508740 | ok |

### Pearson raw: weakest

| Gene | Bgee label | Tissues | Pearson raw | Status |
| --- | --- | --- | --- | --- |
| H2AC12 | H2AC12 | 31 | 0.279630 | ok |
| H2AC15 | H2AC15 | 32 | 0.296210 | ok |
| H2AC21 | H2AC21 | 31 | 0.308403 | ok |
| H2AC16 | H2AC16 | 31 | 0.320666 | ok |
| H2AC14 | H2AC14 | 31 | 0.326711 | ok |

### Spearman raw: strongest

| Gene | Bgee label | Tissues | Spearman raw | Status |
| --- | --- | --- | --- | --- |
| H2AC25 | H2AW | 31 | 0.873274 | ok |
| H2AC1 | H2AC1 | 30 | 0.800726 | ok |
| H2AC20 | H2AC20 | 31 | 0.745534 | ok |
| H2AC7 | H2AC7 | 31 | 0.734403 | ok |
| H2AC8 | H2AC8 | 34 | 0.726592 | ok |

### Spearman raw: weakest

| Gene | Bgee label | Tissues | Spearman raw | Status |
| --- | --- | --- | --- | --- |
| H2AC16 | H2AC16 | 31 | 0.364942 | ok |
| H2AC12 | H2AC12 | 31 | 0.399651 | ok |
| H2AC21 | H2AC21 | 31 | 0.507929 | ok |
| H2AC4 | H2AC4 | 32 | 0.520618 | ok |
| H2AC6 | H2AC6 | 35 | 0.555478 | ok |

### Genes with insufficient or constant input

- None.

## drop_zeros

### Pearson raw: strongest

| Gene | Bgee label | Tissues | Pearson raw | Status |
| --- | --- | --- | --- | --- |
| H2AC1 | H2AC1 | 3 | 0.990639 | ok |
| H2AC25 | H2AW | 31 | 0.694110 | ok |
| H2AC20 | H2AC20 | 31 | 0.667709 | ok |
| H2AC6 | H2AC6 | 35 | 0.609999 | ok |
| H2AC8 | H2AC8 | 32 | 0.607787 | ok |

### Pearson raw: weakest

| Gene | Bgee label | Tissues | Pearson raw | Status |
| --- | --- | --- | --- | --- |
| H2AC14 | H2AC14 | 12 | -0.170402 | ok |
| H2AC21 | H2AC21 | 19 | -0.025297 | ok |
| H2AC16 | H2AC16 | 14 | 0.082993 | ok |
| H2AC4 | H2AC4 | 14 | 0.094588 | ok |
| H2AC17 | H2AC17 | 23 | 0.237043 | ok |

### Spearman raw: strongest

| Gene | Bgee label | Tissues | Spearman raw | Status |
| --- | --- | --- | --- | --- |
| H2AC25 | H2AW | 31 | 0.873274 | ok |
| H2AC20 | H2AC20 | 31 | 0.745534 | ok |
| H2AC8 | H2AC8 | 32 | 0.723455 | ok |
| H2AC7 | H2AC7 | 29 | 0.675062 | ok |
| H2AC15 | H2AC15 | 31 | 0.636447 | ok |

### Spearman raw: weakest

| Gene | Bgee label | Tissues | Spearman raw | Status |
| --- | --- | --- | --- | --- |
| H2AC4 | H2AC4 | 14 | 0.026489 | ok |
| H2AC14 | H2AC14 | 12 | 0.070182 | ok |
| H2AC16 | H2AC16 | 14 | 0.150677 | ok |
| H2AC21 | H2AC21 | 19 | 0.203373 | ok |
| H2AC13 | H2AC13 | 21 | 0.416580 | ok |

### Genes with insufficient or constant input

- None.

