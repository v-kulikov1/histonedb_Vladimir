# CURATED SET of histone sequences and their classification

This directory includes:
1. List of all curated histone sequences in csv format [histones.csv](histones.csv) - this file is manully edited and curated.
2. [classification.json](classification.json) - this file includes all info about histone classification, types, variants, sub-variants, etc. [classificationDef.json](classificationDef.json) - defines JSON Schema.
3. [features.json](features.json) - this file includes info about sequence features of different histones and variants. [featuresDef.json](featuresDef.json) - defines JSON Schema.
4. a library to analyze and vizualize histone sequence information which produces web-pages rendered with GitHub pages.
5. a library to prepare inputs for HistoneDB, namely [histones_processed.csv](histones_processed.csv), which contains actual sequences (and not only NCBI accessions).


## histones.csv

Accession.version, Type, Variant, Subvariant, Doublet, GI, TaxonomyID, Organism, Taxonomy_group, Info, Sequence (if no NCBI accession is present)

- Accession.version - this is NCBI id of the sequence, also will be used as an id in HistoneDB. If the sequence in NCBI is not present, than a custom id is used (might be NONCBI_VARIANTNAME_NUMBER), and sequence is provided in Sequence column.
- Type - according to classification, H3, H4, H2A, H2B or Archaeal.
- Variant - this is the id of the top-level variant classification in our hierarchy attributable to the sequence. If sequence has to be attributed to two variants wirte them separated with a space, e.g. "H2A.X H2A.Z"
- Subvariant - this is the id of the most specific level variant classification in our hierarchy attributable to the sequence or none if equal to Variant. If sequence has to be attributed to two subvariants wirte them separated with a space.
- Doublet - true if it is a doublet. if it is of both types - write them with a space in Type column. E.g. "H3 H4"
- GI - legacy field for GIs.
- TaxonomyID - NCBI taxonomy id of the sequence.
- Organism - NCBI human readable taxonomy name.
- Taxonomy_group - this is usually taxonomy class if available or higher order rank if not. We adhere to NCBI current taxnomy name E.g. Mammalia, Magnoliopsida, etc.
- Info - information about this particular sequence - including its function and references to literature as [PMID].
- Sequence - sequence as a string if no NCBI accession is available.

### Special cases
- If histone is a doublet - we include it ones - see above.

## classification.json

This is a very important file detailing all classification of hisone variants currently used, and all info about them.
The classification is hierarchichal.
[classificationDef.json](classificationDef.json) - defines and describes json schema.

### Here are key points about our classification:
- Classification is hierarchical.
- Top level is type: H3, H4, H2A, H2B or Archaeal
- Next level is top-level variants, e.g. canonical H2A in Metazoa
- Variants are always specify in their name the taxonomic span of this particular variant.
- Every variant has an id of the following form VARIANTNAME_(Taxa)
- Variant ids are case-insensitive in the database, but during representation case is important.
- Variants also have a full name. E.g.

