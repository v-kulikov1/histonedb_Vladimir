# human_histones.csv revision information

## Structure of human_histones.csv
What information do we need:

1. Histone type
2. Histone variant (Maximally specific, histone variant group will be devised from classification JSON).
3. Canonical isoform (This is experimental)
4. HGNC symbol
5. NCBI gene ID
6. Ensembl gene ID
7. Transcript stable ID
8. RefSeq mRNA ID
9. RefSeq protein ID
10. Protein length
11. References
12. Info


H2A.Z.1.s2_Homo_sapiens_H2AZ2_XP_123456
cH2A.1.1_Homo_sapiens
variant_group=>

## Where the info was taken from
1. New histone gene nomenclature and paper from HGNC
2. Look at data in RefSeq, Ensembl
3. Original research papers

## What proteins were chosen

Protein coding transcripts from Ensemble release 105 that 
1. Are in CCDS
2. Have a TSL1, TSL2 or TSL3, unless literature data is available to support lower levels of TSL.
3. If transcripts give rise to identical protein sequence only one record is retained.
4. If recent literature data is available suggesting presence of currently unannotated proteins, this was added on top.

#TODO
1. add geneid to human histones

#History