from django.core.management.base import BaseCommand, CommandError
from django.conf import settings

from browse.models import Sequence, Score, Histone, Variant
from tools.blast_search import get_best_hsp
from tools.test_model import plot_scores, test_curated
from tools.blast_search import make_blastp, load_blast_search, load_blast_search_diagnosis, add_score

from Bio import SeqIO, AlignIO #, pairwise2
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from Bio.Emboss.Applications import NeedleCommandline

import subprocess
import os, sys, io
import re
import uuid
import logging
import pickle

from datetime import date, datetime

from cProfile import Profile

# BLAST_PROCS=25
BLAST_PROCS=5 # for small random data

class Command(BaseCommand):
    help = 'Blast sequences from HistoneDB and load classification data'
    seed_directory = os.path.join(settings.STATIC_ROOT_AUX, "browse", "seeds")
    blast_directory = os.path.join(settings.STATIC_ROOT_AUX, "browse", "blast")
    blast_file = os.path.join(blast_directory, "search", "blast_search.out")
    model_evaluation = os.path.join(blast_directory, "model_evaluation")
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
            help="Force the regeneration of BLAST from Database, BLASTP search in db_file and loading of results to database")

        parser.add_argument(
            "-t",
            "--test",
            default=False,
            action="store_true",
            help="Test algorithm of BLAST classification on curated")

        parser.add_argument(
            "--profile",
            default=False,
            action="store_true",
            help="Profile the command")

    def _handle(self, *args, **options):
        self.log.info('=======================================================')
        self.log.info('===               blastvariants START               ===')
        self.log.info('=======================================================')
        self.start_time = datetime.now()
        if options["force"]:
            # Clean the DB, removing all sequence/scores/etc
            Score.objects.all().delete()
            for s in Sequence.objects.all():
                s.variant = None
                s.save()

        if options["test"]:
            # self.load_curated()
            # self.add_scores_for_curated()
            self.add_identities_for_curated()
            # self.estimate_thresholds()
        else:
            self.load_curated()
            self.add_scores_for_curated()
            # self.estimate_thresholds()

            if options["force"] or not os.path.isfile(self.blast_file + "0"):
            # if not os.path.isfile(self.blast_file + "0"):
                self.search_blast()

            self.load_in_db()
            self.get_stats()

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

    def estimate_thresholds(self):
        """
            Estimate BLAST thresholds that we will use for variant classification.
            # Construct two sets for every variant:
            #     negative: The seed alignmnents from every other variant
            #     positive: the current seed alignment for the variant
            # And estimate params from ROC-curves.
        """
        test_curated(save_dir=self.model_evaluation)

    def load_curated(self):
        self.log.info('Loading curated...')
        for s in Sequence.objects.filter(reviewed=True):
            s.variant = s.variant_hmm
            s.save()

    def add_scores_for_curated(self):
        # Score.objects.filter(sequence__reviewed=True).delete()
        for sequence in Sequence.objects.filter(reviewed=True):
            seqs_file = os.path.join(settings.STATIC_ROOT_AUX, "browse", "blast", "HistoneDB_sequences.fa")

            blastp = os.path.join(os.path.dirname(sys.executable), "blastp")
            # output = os.path.join("/", "tmp", "{}.xml".format(uuid.uuid4()))
            blastp_cline = NcbiblastpCommandline(
                cmd=blastp,
                db=seqs_file,
                evalue=.01, outfmt=5)
            result, error = blastp_cline(stdin="\n".join([s.format("fasta") for s in [sequence]]))

            resultFile = io.BytesIO()
            resultFile.write(result.encode("utf-8"))
            resultFile.seek(0)

            for i, blast_record in enumerate(NCBIXML.parse(resultFile)):
                if len(blast_record.alignments) == 0:
                    self.log.error('No BLAST record alignments for {} during adding scores for curated sequences'.format(sequence.id))
                    continue

                for algn in blast_record.alignments:
                    algn_self = False if sequence.id != algn.hit_def.split("|")[0] else True # alignment on itself
                    for hsp in algn.hsps:
                        add_score(sequence, sequence.variant, hsp, algn.hit_def.split("|")[0], best=algn_self) ### looks like there is an error

    def add_identities_for_curated(self):
        for asequence in Sequence.objects.filter(reviewed=True):
            for bsequence in Sequence.objects.filter(reviewed=True):
                # global_align = pairwise2.align.globalxx(asequence.format(format='fasta', ungap=True),
                #                                         bsequence.format(format='fasta', ungap=True))
                # print(global_align)
                # print(global_align[0])
                with open(os.path.join(self.model_evaluation, "aseq.fasta"), 'w') as f:
                    f.write(asequence.format(format='fasta', ungap=True))
                with open(os.path.join(self.model_evaluation, "bseq.fasta"), 'w') as f:
                    f.write(bsequence.format(format='fasta', ungap=True))

                n2 = str(uuid.uuid4())
                needle_results = os.path.join(self.model_evaluation, "needle_{}.txt".format(n2))
                cmd = os.path.join(os.path.dirname(sys.executable), "needle")

                if not os.path.isfile(cmd):
                    cmd = "needle"
                needle_cline = NeedleCommandline(
                    cmd=cmd,
                    asequence=os.path.join(self.model_evaluation, "aseq.fasta"),
                    bsequence=os.path.join(self.model_evaluation, "bseq.fasta"),
                    gapopen=10,
                    gapextend=1,
                    nobrief=True,
                    # similarity=True,
                    outfile=needle_results)
                # needle_cline.similarity = True
                stdout, stderr = needle_cline()
                with open(needle_results) as f:
                    out_split = f.readlines()
                # print(needle_cline)
                # print(needle_cline.nobrief)
                # print(needle_cline.similarity)
                print(out_split[24])
                print(out_split[25])

                identity_split = re.search(r'\d+\/\d+', out_split[24]).group(0).split('/')
                similarity_split = re.search(r'\d+\/\d+', out_split[25]).group(0).split('/')
                identity = float(identity_split[0])/float(identity_split[1])
                similarity = float(similarity_split[0])/float(similarity_split[1])
                print(identity)
                print(similarity)
                self.add_identity_score(asequence, bsequence.variant, identity, similarity, bsequence.id,
                                        best=False if asequence.id != bsequence.id else True)
                break
            break

    def add_identity_score(self, seq, variant_model, identity, similarity, hit_accession, best=False):
        """Add score for a given sequence"""
        score = ScoreIdentity(
            # id                      = score_num,
            sequence=seq,
            variant=variant_model,
            score=identity,
            similarity=similarity,
            hit_accession=hit_accession,
            used_for_classification=best,
        )
        score.save()

    def search_blast(self):
        # sequences = Sequence.objects.filter(reviewed=False).values_list('sequence')
        for hist_type in ['H1', 'H2A', 'H2B', 'H3', 'H4']:
            blastdb_file = os.path.join(settings.STATIC_ROOT_AUX, "browse", "blast", hist_type,
                                        "BLASTDB_sequences_{}.fa".format(hist_type))
            self.create_blastdb(Sequence.objects.filter(reviewed=True).filter(variant__hist_type__id=hist_type), blastdb_file)

            sequences_all = [seq.format(format='fasta') for seq in Sequence.objects.filter(reviewed=False).filter(variant_hmm__hist_type__id=hist_type)]
            # sequences_all = [seq.format(format='fasta') for seq in Sequence.objects.filter(reviewed=False)]
            self.log.info("Running BLASTP for {} sequences of {}...".format(len(sequences_all), hist_type))
            split_count = int(len(sequences_all)/BLAST_PROCS)

            for i in range(BLAST_PROCS+1):
                sequences = sequences_all[split_count * i:split_count * i + split_count]
                save_to = self.blast_file + hist_type + "%d" % i
                self.log.info('Starting Blast sequences for {}/{} of {}'.format(i, BLAST_PROCS, hist_type))
                make_blastp(sequences, blastdb_file, save_to=save_to)

    def load_in_db(self):
        self.log.info("Loading BLASTP data into HistoneDB...")
        for hist_type in ['H1', 'H2A', 'H2B', 'H3', 'H4']:
            for i in range(BLAST_PROCS+1):
                with open(self.blast_file + hist_type + "%d" % i) as blastFile:
                    if re.search(r'^\s*$', blastFile.read()):
                        self.log.warning('Results of {} is empty'.format(self.blast_file + hist_type + "%d" % i))
                        continue
                load_blast_search(self.blast_file + hist_type + "%d" % i)
                self.log.info('Loaded {}/{} BlastRecords of {}'.format(i,BLAST_PROCS, hist_type))
                self.log.info('Classified {} from {}'.format(Sequence.objects.exclude(variant=None).count(),Sequence.objects.all().count()))

        non_classified = Sequence.objects.filter(variant=None)
        for s in non_classified:
            s.variant = Variant.objects.get(id='generic_{}'.format(s.variant_hmm.hist_type.id))
            s.save()
            score = Score(
                sequence=s,
                variant=s.variant,
                score=.000001,
                bitScore=None,
                evalue=.011,
                blastStart=None,
                blastEnd=None,
                seqStart=None,
                seqEnd=None,
                align_length=None,
                used_for_classification=True,
            )
            score.save()


    def create_blastdb(self, sequences, seqs_file):
        # seqs_file = os.path.join(settings.STATIC_ROOT_AUX, "browse", "blast", "HistoneDB_sequences.fa")
        with open(seqs_file, "w") as seqs:
            for s in sequences:  # here we restrict the blast DB to reviewed seqs
                SeqIO.write(s.to_biopython(ungap=True), seqs, "fasta")

        makeblastdb = os.path.join(os.path.dirname(sys.executable), "makeblastdb")
        subprocess.call(["makeblastdb", "-in", seqs_file, "-dbtype", "prot", "-title", "HistoneDB"])

    def get_stats(self, filename_suff = ''):
        self.log.info('Outputting statistics file ...')

        with open('NR_VERSION', 'r') as nrv:
            nr_version = nrv.read()

        now = datetime.now()
        dt_string = now.strftime("%Y%m%d-%H%M%S")
        with open('log/blast_db_stat_'+dt_string+filename_suff,'w') as f:
            f.write("Variant database regeneration statistics after BLAST\n")
            f.write("DB regen start time: %s \n"%self.start_time)
            f.write("DB regen end time: %s\n"%now)
            f.write("Time taken for regeneration of variants: %f hours\n"%(float((now-self.start_time).total_seconds())/3600.))
            f.write("Parallel threads used %d\n"%BLAST_PROCS)
            f.write("DB file used: %s\n" % nr_version)
            f.write('---Database statistics----\n')
            f.write('Total seqs = %d\n'%Sequence.objects.all().count())
            f.write('Reviewed seqs = %d\n'%Sequence.objects.filter(reviewed=True).count())
            f.write('Automatic seqs = %d\n'%Sequence.objects.filter(reviewed=False).count())

            f.write('\n---Histone type statistics----\n')
            f.write('Type        | Total  |Reviewed|  Auto  \n')
            for h in Histone.objects.all():
                tot=Sequence.objects.filter(variant__hist_type=h).count()
                rev=Sequence.objects.filter(variant__hist_type=h,reviewed=True).count()
                auto=Sequence.objects.filter(variant__hist_type=h,reviewed=False).count()

                f.write('%12s|%8d|%8d|%8d\n'%(h.id,tot,rev,auto))

            f.write('\n---Histone variant statistics----\n')
            f.write('Variant     | Total  |Reviewed|  Auto  \n')
            for v in Variant.objects.all():
                tot=Sequence.objects.filter(variant=v).count()
                rev=Sequence.objects.filter(variant=v,reviewed=True).count()
                auto=Sequence.objects.filter(variant=v,reviewed=False).count()

                f.write('%12s|%8d|%8d|%8d\n'%(v.id,tot,rev,auto))


