from django.core.management.base import BaseCommand, CommandError

from browse.models import Sequence
from tools.browse_service import *

from Bio import SeqIO

import subprocess
import logging
from cProfile import Profile

class Command(BaseCommand):
    help = 'Build BLASTDBs for histone types'

    # Logging info
    logging.basicConfig(filename=os.path.join(LOG_DIRECTORY, "buildblastdbs.log"),
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
        self.log.info('===               buildblastdbs START               ===')
        self.log.info('=======================================================')
        self.create_seq_dbs()
        self.build_blast_dbs()
        self.log.info('=======================================================')
        self.log.info('===       buildblastdbs SUCCESSFULLY finished       ===')
        self.log.info('=======================================================')

    def handle(self, *args, **options):
        if options['profile']:
            profiler = Profile()
            profiler.runcall(self._handle, *args, **options)
            profiler.print_stats()
        else:
            self._handle(*args, **options)

    def create_seq_dbs(self):
        if not os.path.exists(BLASTDBS_DIR):
            os.makedirs(BLASTDBS_DIR)
        for hist_type in ['H1', 'H2A', 'H2B', 'H3', 'H4']:
            with open(BLASTDB_FILE.format(hist_type), "w") as seqs:
                for s in Sequence.objects.filter(reviewed=True).filter(variant__hist_type__id=hist_type):
                    SeqIO.write(s.to_biopython(ungap=True), seqs, "fasta")

    def build_blast_dbs(self):
        child_procs = []
        for hist_type in ['H1', 'H2A', 'H2B', 'H3', 'H4']:
            self.log.info(" ".join(["makeblastdb", "-in", BLASTDB_FILE.format(hist_type), "-dbtype", "prot","-title", "HistoneDB"]))
            p = subprocess.Popen(["makeblastdb", "-in", BLASTDB_FILE.format(hist_type), "-dbtype", "prot","-title", "HistoneDB"])
            child_procs.append(p)
        for cp in child_procs:
            cp.wait()