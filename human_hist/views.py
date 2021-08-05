from django.shortcuts import render

from browse.models import *
from human_hist.models import * 
from djangophylocore.models import *

def human_histones(request):
    human_hist = Histone_Human_genes.objects.all()
    data = []
    for gene in human_hist:
        human_prot = gene.human_proteins.all() #Human_Histone_proteins.objects.filter(gene=gene) или без .all()
        for prot in human_prot:

                data.append(  {
                    "hgnc_symbol" : gene.hgnc_symbol,
                    "prev_hgnc_symb" : gene.prev_hgnc_symb,
                    "ncbi_gene_id" : gene.ncbi_gene_id,
                    "ensg" : gene.ensg,
                    "expr_timing" : gene.expr_timing,
                    "expr_pattern" : gene.expr_pattern,
                    "biotype" : gene.biotype,
                    "bona_fidecanonical" : gene.bona_fidecanonical,
                    "pmids" : gene.pmids,
                    "variant" : gene.variant, 
                    "hist_type" : gene.hist_type,
                    "enst" : prot.enst,
                    "refseq_transcript_id" : prot.refseq_transcript_id,
                    "refseq_protein_id" : prot.refseq_protein_id, 
                    "prot_lenght" : prot.prot_lenght,
                    "isoform" : prot.isoform,
                } )

    return render(request, 'human_hist.html', {'human_histones': data})

def browse_variant_with_highlighted_sequence(request, histone_type, variant, accession):
    return browse_variant(request, histone_type, variant, accession)

def hist_mutations(request):
    return render(request, 'hist_mutations.html')