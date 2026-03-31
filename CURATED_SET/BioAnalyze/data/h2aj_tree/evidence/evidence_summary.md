# Tree Reconstruction Evidence Summary

## Key quantitative checkpoints

- `H2AJ.fasta`: 230 records
- `cH2A.fasta`: 57 records
- `platypus.fasta`: 51 records
- `All_AA_0206`: 387 records, 208 sites
- `All_AA_1006`: 338 records, 180 sites
- `All_NUC_1606`: 315 records, 396 sites
- `SQK_nuc.fasta`: 227 records
- `SQK_nuc(without short).fasta`: 224 records
- `cH2A_nuc + platypus_nuc + SQK_nuc.fasta`: 312 records

## Current clean-profile override

- The `clean` outputs intentionally drop three short H2A.J/SQK records:
- `Myotis-lucifugus|XM_006084274`
- `Homo-sapiens|AK303301`
- `Homo-sapiens|AL133626`
- July clean SQK outputs now treat `SQK_nuc(without short)` and `protein_from_SQK_nuc(without short)` as the canonical H2A.J source.

## Confirmed AA web run

- Data type: `aa`
- Model: `LG`
- Initial tree: `BioNJ`
- Rate model: `FreeRate` with `3` classes
- Support mode: `yes (SH-like branch supports)`
- Runtime: `0h22m57s (1377 seconds)`

## Historical interpretation

- `All_AA_0206` is the early broad-outgroup protein run and not the later mammal/nonplacental-only set.
- `All_AA_1006 = 338` is the file whose count matches `H2AJ (230) + cH2A (57) + nonplacental labelled as platypus (51)`.
- The nucleotide failure is upstream of PhyML: BLAST HSP fragments, truncated tails, mixed modality, and short sequences all damaged the input before tree building.