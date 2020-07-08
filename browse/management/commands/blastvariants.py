from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from browse.models import Sequence, SequenceBlast, ScoreBlast
import subprocess
import os, sys
from tools.blast_search import make_blastp, load_blast_search

from Bio import SeqIO

import logging
import pickle

from cProfile import Profile

BLAST_PROCS=200

class Command(BaseCommand):
    help = 'Blast sequences from HistoneDB and load classification data'
    seed_directory = os.path.join(settings.STATIC_ROOT_AUX, "browse", "seeds")
    blast_file = os.path.join(settings.STATIC_ROOT_AUX, "browse", "blast", "search", "blast_search.out")
    # Logging info
    logging.basicConfig(filename='log/blastvariants.log',
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
        parser.add_argument(
            "--profile",
            default=False,
            action="store_true",
            help="Profile the command")

    def _handle(self, *args, **options):
        self.log.info('=======================================================')
        self.log.info('===               blastvariants START               ===')
        self.log.info('=======================================================')

        if options["force"]:
            # Clean the DB, removing all sequence/variants/etc
            SequenceBlast.objects.all().delete()
            ScoreBlast.objects.all().delete()

        self.load_curated()

        # Make BLASTDB for curated sequences and blast sequences extracted by HMMs
        # self.make_blastdb()
        # if options["force"] or not os.path.isfile(self.blast_file + "0"):
        if not os.path.isfile(self.blast_file + "0"):
            self.search_blast()

        self.load_in_db()

        self.log.info('=======================================================')
        self.log.info('===       blastvariants SUCCESSFULLY finished       ===')
        self.log.info('=======================================================')

    def handle(self, *args, **options):
        if options['profile']:
            profiler = Profile()
            profiler.runcall(self._handle, *args, **options)
            profiler.print_stats()
        else:
            self._handle(*args, **options)

    def load_curated(self):
        self.log.info('Loading curated...')
        for s in Sequence.objects.filter(reviewed=True):
            seq = SequenceBlast(
                accession=s.id,
                variant=s.variant,
            )
            seq.save()

    # def search_and_load_blast(self):
    #     # sequences = Sequence.objects.filter(reviewed=False).values_list('sequence')
    #     sequences_all = [seq.format(format='fasta') for seq in Sequence.objects.filter(reviewed=False)]
    #     self.log.info("Running BLASTP for {} sequences...".format(len(sequences_all)))
    #     split_count = int(len(sequences_all)/BLAST_PROCS)
    #     for i in range(BLAST_PROCS+1):
    #         saved_pickle = os.path.join(settings.STATIC_ROOT_AUX, "browse", "blast", "search_pickles", "search_result_{}.pickle".format(i))
    #         if os.path.isfile(saved_pickle):
    #             self.log.info('Loading Blast results from pickle {}'.format(saved_pickle))
    #             with open(saved_pickle, 'rb') as f:
    #                 blastp_res, error = pickle.load(f)
    #         else:
    #             sequences = sequences_all[split_count * i:split_count * i + split_count]
    #             self.log.info('Starting Blast sequences for {}/{}'.format(i, BLAST_PROCS))
    #             blastp_res, error = make_blastp(sequences, save_to=saved_pickle)
    #             self.log.info('Blast results saved to {}'.format(saved_pickle))
    #         # blastp_res, error = make_blastp(sequences, save_to=saved_pickle)
    #         self.log.info("Loading BLASTP data into HistoneDB...")
    #         load_blast_search(blastp_res)
    #         # self.log.info('Loaded {} BlastRecords'.format(len(sequences)))
    #         self.log.info('Loaded {}/{} BlastRecords'.format(i,BLAST_PROCS))

    def search_blast(self):
        # sequences = Sequence.objects.filter(reviewed=False).values_list('sequence')
        sequences_all = [seq.format(format='fasta') for seq in Sequence.objects.filter(reviewed=False)]
        self.log.info("Running BLASTP for {} sequences...".format(len(sequences_all)))
        split_count = int(len(sequences_all)/BLAST_PROCS)

        for i in range(BLAST_PROCS+1):
            sequences = sequences_all[split_count * i:split_count * i + split_count]
            save_to = self.blast_file + "%d" % i
            self.log.info('Starting Blast sequences for {}/{}'.format(i, BLAST_PROCS))
            make_blastp(sequences, save_to=save_to)

    def load_in_db(self):
        self.log.info("Loading BLASTP data into HistoneDB...")
        # for i in range(BLAST_PROCS + 1):
        for i in range(2):
            load_blast_search(self.blast_file + "%d" % i)
            # self.log.info('Loaded {} BlastRecords'.format(len(sequences)))
            self.log.info('Loaded {}/{} BlastRecords'.format(i,BLAST_PROCS))


