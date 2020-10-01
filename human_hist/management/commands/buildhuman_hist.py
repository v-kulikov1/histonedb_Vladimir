from django.core.management.base import BaseCommand, CommandError
from django.core.exceptions import ObjectDoesNotExist

import os
import logging
from django.conf import settings
from browse.models import Variant, OldStyleVariant, Histone, Publication, TemplateSequence, Feature
from djangophylocore.models import Taxonomy
import json
from Bio import SeqIO
from itertools import groupby
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from human_hist.models import * # Histone_Human_genes, Histone_Human_proteins, Histone_Human_mutations, Histone_Human_cancers

import pandas as pd


class Command(BaseCommand):
    info_directory = os.path.join(settings.STATIC_ROOT_AUX, "human_hist", "info")

    # Logging info

    logging.basicConfig(filename='log/buildhuman_hist.log',
                        format='%(asctime)s %(name)s %(levelname)-8s %(message)s',
                        level=logging.INFO,
                        datefmt='%Y-%m-%d %H:%M:%S')
    log = logging.getLogger(__name__)

    def add_arguments(self, parser):
        parser.add_argument(
            "-f",
            "--force",
            default=False,
            action="store_true",
            help="Force the regeneration of HMM from seeds, HUMMER search in db_file, Test models and loading of results to database")

    def handle(self, *args, **options):
        self.log.info('=======================================================')
        self.log.info('===            buildhuman_hist START              ===')
        self.log.info('=======================================================')

        Histone_Human_genes.objects.all().delete()
        Histone_Human_proteins.objects.all().delete()
        Histone_Human_mutations.objects.all().delete()
        Histone_Human_cancers.objects.all().delete()

        # Genes
        self.log.info('===   GENES START  ===')
        human_hist_genes_info_path = os.path.join(self.info_directory, "human_hist_genes.csv")
        human_genes_table = pd.read_csv(human_hist_genes_info_path)

        for index, gene in human_genes_table.iterrows():
            self.log.info('{}'.format(index))
            try:
                variant = Variant.objects.get(id=gene['Histone variant'])
            except ObjectDoesNotExist:
                self.log.info('not in Variant model {}'.format(gene['Histone variant']))
                continue

            obj = Histone_Human_genes(hgnc_symbol=gene["HGNC Symbol"], prev_hgnc_symb=gene['Previous HGNC Symbol'], \
                                      ncbi_gene_id=gene['NCBI gene ID'], ensg=gene['Ensembl gene ID'],
                                      expr_timing=gene['Expr. timing'], \
                                      expr_pattern=gene['Expr. pattern'], biotype=gene['Biotype'],
                                      bona_fidecanonical=gene['Bona fide canonical'], \
                                      pmids=gene['PMIDs'], variant = variant )

            obj.save()


            if obj:
                self.log.info("{} was created".format(obj.hgnc_symbol))

        # Proteins
        self.log.info('===   PROTEINS START  ===')
        human_hist_proteins_info_path = os.path.join(self.info_directory, "human_hist_proteins.csv")
        human_proteins_table = pd.read_csv(human_hist_proteins_info_path)

        for index, protein in human_proteins_table.iterrows():
            self.log.info('{}'.format(index))

            genes = Histone_Human_genes.objects.filter(hgnc_symbol=protein['HGNC Symbol'])

            obj = Histone_Human_proteins(enst=protein['Transcript stable ID'],
                                         refseq_transcript_id=protein['RefSeq mRNA ID'], \
                                         refseq_protein_id=protein['RefSeq peptide ID'], prot_lenght=protein['Protein length'],\
                                       isoform = protein['canonical_isoform'] )
            obj.save()
            obj.gene.add(*genes)



            if obj:
                self.log.info("{} was created".format(obj.enst))

        # CAncer
        self.log.info('===   CANCERS START  ===')
        human_hist_cancers_info_path = os.path.join(self.info_directory, "human_hist_cancers.csv")
        human_cancers_table = pd.read_csv(human_hist_cancers_info_path)
        
        for index, patient in human_cancers_table.iterrows():
            obj = Histone_Human_cancers(case_id=patient['case_id'], \
                                      sequencing_center=patient['sequencing_center'], cancer=patient['cancer'])
            obj.save()

            if obj:
                self.log.info("{} was created".format(obj.case_id))

        # Mutations
        self.log.info('===   MUTATIONS START  ===')
        human_hist_mutations_info_path = os.path.join(self.info_directory, "human_hist_mutations.csv")
        human_mutations_table = pd.read_csv(human_hist_mutations_info_path)

        for index, mutation in human_mutations_table.iterrows():

            genes = Histone_Human_genes.objects.filter(hgnc_symbol=mutation['HGNC Symbol'])
            # without 'chr','start_position','end_position',
            cases = Histone_Human_cancers.objects.filter(case_id=mutation['case_id'])

            obj = Histone_Human_mutations(aa_change=mutation['amino_acid_change'], \
                                  mutation_type=mutation['mutation_type'], ref_allele=mutation['reference_allele'], \
                                  var_allele=mutation['variant_allele'],
                                  var_allele_freq=mutation['variant_allele_freq_tumor'])
            obj.save()

            obj.case.add(*cases)
            obj.gene.add(*genes)
            if obj:
                self.log.info("{} was created".format(obj.aa_change))














