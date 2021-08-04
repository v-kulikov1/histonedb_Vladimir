from django.core.management.base import BaseCommand, CommandError

from tools.browse_service import *

import subprocess
import logging
from cProfile import Profile

class Command(BaseCommand):
    help = 'Build HMMs for histone types and histone variants'

    # Logging info
    logging.basicConfig(filename=os.path.join(LOG_DIRECTORY, "buildhmms.log"),
                        format='%(asctime)s %(name)s %(levelname)-8s %(message)s',
                        level=logging.INFO,
                        datefmt='%Y-%m-%d %H:%M:%S')
    log = logging.getLogger(__name__)

    def add_arguments(self, parser):
        parser.add_argument(
            "--test",
            default=False,
            action="store_true",
            help="Use this option for test running if you need a small version of specified database")

        parser.add_argument(
            "--profile",
            default=False,
            action="store_true",
            help="Profile the command")

    def _handle(self, *args, **options):
        self.log.info('=======================================================')
        self.log.info('===                 buildhmms START                 ===')
        self.log.info('=======================================================')

        #Create HMMs from seeds and compress them to one HMM file tp faster search with hmmpress.
        self.build_hmms_from_seeds()
        self.press_hmms()

        self.log.info('=======================================================')
        self.log.info('===         buildhmms SUCCESSFULLY finished         ===')
        self.log.info('=======================================================')

    def handle(self, *args, **options):
        if options['profile']:
            profiler = Profile()
            profiler.runcall(self._handle, *args, **options)
            profiler.print_stats()
        else:
            self._handle(*args, **options)

    def build_hmms_from_seeds(self):
        """Build HMMs from seed histone sequences,
        and outputing them to
        static/browse/hmms/
        to individual dirs as well as combining to pne file combined_hmm_file
        """
        self.log.info("Building HMMs...")

        with open(COMBINED_HMM_HISTTYPES_FILE, "w") as combined_hmm_histtypes, open(COMBINED_HMM_VARIANTS_FILE, "w") as combined_hmm_variants:
            for hist_type, seed in get_seeds(seeds_name=SEEDS_FOR_HMM, combined_alignments=True):
                #Build HMMs
                hmm_dir = os.path.join(HMM_DIRECTORY, hist_type)
                if not os.path.exists(hmm_dir):
                    os.makedirs(hmm_dir)
                hmm_file = os.path.join(hmm_dir, "{}.hmm".format(seed[:-6]))
                self.build_hmm(seed[:-6], hmm_file, os.path.join(FOLD_SEED_DIRECTORY, hist_type, seed))
                if hist_type=='': #combine all combined by hist_types sequences in one file (by Preety)
                    with open(hmm_file) as hmm:
                        print(hmm.read().rstrip(), file=combined_hmm_histtypes)
                    continue
                with open(hmm_file) as hmm:
                    print(hmm.read().rstrip(), file=combined_hmm_variants)

    def build_hmm(self, name, db, seqs):
        self.log.info(' '.join(["hmmbuild", "-n", name, db, seqs]))
        subprocess.call(["hmmbuild", "-n", name, db, seqs])


    def press_hmms(self):
        self.press(COMBINED_HMM_HISTTYPES_FILE)
        self.press(COMBINED_HMM_VARIANTS_FILE)

    def press(self, combined_hmm):
        """Press the HMMs into a single HMM file, overwriting if present"""
        self.log.info("Pressing HMMs...")
        # print >> self.stdout, "Pressing HMMs..."
        subprocess.call(["hmmpress", "-f", combined_hmm])