# Bgee H2A Expression Source Audit

Generated: 2026-03-24 18:41:02

- Merged H2A file: `CURATED_SET\BioAnalyze\data\merged\mammalia_H2A_merged_with_taxonomy_v4.csv`
- Raw Bgee dir: `C:\Users\USER\Documents\GitHub\histonedb_external_storage\BioAnalyze\raw`
- Detailed audit TSV: `CURATED_SET\BioAnalyze\audits\bgee_h2a_expression_source_audit_3species.tsv`

## Overview

| Species | H2A genes scanned | Matched raw rows | Qualifying rows | RNA-Seq only | Multi-source | Non-RNA only | Unresolved |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Homo sapiens | 29 | 32568 | 21425 | 12476 | 4305 | 4644 | 0 |
| Bos taurus | 28 | 7212 | 4644 | 4644 | 0 | 0 | 0 |
| Mus musculus | 9 | 5453 | 2263 | 595 | 134 | 1534 | 0 |

## Homo sapiens

- H2A genes scanned: 29
- Matched raw rows: 32568
- Heatmap-qualifying rows (`present` + `gold quality` + non-empty overall `Expression score`): 21425
- Qualifying source classes: `rna_seq_only=12476` (58.23%), `multi_source=4305` (20.09%), `non_rna_only=4644` (21.68%), `unresolved=0` (0.00%)
- Conclusion: For qualifying H2A rows in Homo sapiens, expression score is predominantly RNA-Seq-based.

Non-`rna_seq_only` qualifying examples:

| Gene ID | Gene name | Anatomical entity | Source class | Active sources | Expression score | RNA-Seq expression score |
| --- | --- | --- | --- | --- | ---: | ---: |
| ENSG00000099284 | MACROH2A2 | oocyte | non_rna_only | Affymetrix | 79.61 | - |
| ENSG00000099284 | MACROH2A2 | male germ line stem cell (sensu Vertebrata) in testis | non_rna_only | single-cell RNA-Seq | 79.89 | - |
| ENSG00000099284 | MACROH2A2 | secondary oocyte | non_rna_only | Affymetrix | 80.51 | - |
| ENSG00000099284 | MACROH2A2 | primordial germ cell in gonad | non_rna_only | single-cell RNA-Seq | 85.60 | - |
| ENSG00000099284 | MACROH2A2 | mononuclear cell | multi_source | Affymetrix|RNA-Seq | 71.62 | 71.64 |

## Bos taurus

- H2A genes scanned: 28
- Matched raw rows: 7212
- Heatmap-qualifying rows (`present` + `gold quality` + non-empty overall `Expression score`): 4644
- Qualifying source classes: `rna_seq_only=4644` (100.00%), `multi_source=0` (0.00%), `non_rna_only=0` (0.00%), `unresolved=0` (0.00%)
- Conclusion: For qualifying H2A rows in Bos taurus, expression score is exclusively RNA-Seq-based.

No non-`rna_seq_only` qualifying examples found.

## Mus musculus

- H2A genes scanned: 9
- Matched raw rows: 5453
- Heatmap-qualifying rows (`present` + `gold quality` + non-empty overall `Expression score`): 2263
- Qualifying source classes: `rna_seq_only=595` (26.29%), `multi_source=134` (5.92%), `non_rna_only=1534` (67.79%), `unresolved=0` (0.00%)
- Conclusion: For qualifying H2A rows in Mus musculus, expression score is not exclusively RNA-Seq-based.

Non-`rna_seq_only` qualifying examples:

| Gene ID | Gene name | Anatomical entity | Source class | Active sources | Expression score | RNA-Seq expression score |
| --- | --- | --- | --- | --- | ---: | ---: |
| ENSMUSG00000040456 | H2ap | multicellular organism | multi_source | Affymetrix|RNA-Seq | 37.31 | 23.03 |
| ENSMUSG00000040456 | H2ap | testis | non_rna_only | Affymetrix | 98.10 | - |
| ENSMUSG00000040456 | H2ap | testis | non_rna_only | Affymetrix | 98.71 | - |
| ENSMUSG00000040456 | H2ap | testis | non_rna_only | Affymetrix | 98.55 | - |
| ENSMUSG00000040456 | H2ap | testis | non_rna_only | Affymetrix | 89.18 | - |
