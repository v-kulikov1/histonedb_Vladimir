from ete3 import NCBITaxa

HIGHLVL_TAXA = ['Metazoa', 'Fungi',
                'Rhodophyta', 'Viridiplantae',
                'Rhizaria', 'Alveolata', 'Stramenopiles',
                'Archaea']
NCBI_TAXA = NCBITaxa()
TAXIDS_DICT = NCBI_TAXA.get_name_translator(HIGHLVL_TAXA)
TAXIDS = [TAXIDS_DICT[name][0] for name in HIGHLVL_TAXA]