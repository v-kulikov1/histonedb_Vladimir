from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from browse.models import Sequence, SequenceBlast, ScoreBlast, Histone, Variant
from tools.blast_search import get_best_hsp
from tools.test_model import plot_scores
import subprocess
import os, sys, io
from tools.blast_search import make_blastp, load_blast_search, load_blast_search_diagnosis, add_score

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML

import logging
import pickle

from datetime import date, datetime

from cProfile import Profile

# BLAST_PROCS=100
BLAST_PROCS=25

class Command(BaseCommand):
    help = 'Blast sequences from HistoneDB and load classification data'
    seed_directory = os.path.join(settings.STATIC_ROOT_AUX, "browse", "seeds")
    blast_directory = os.path.join(settings.STATIC_ROOT_AUX, "browse", "blast")
    blast_file = os.path.join(blast_directory, "search", "blast_search.out")
    blast_curated_file = os.path.join(blast_directory, "search", "blast_search_curated.out")
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
            # Clean the DB, removing all sequence/variants/etc
            ScoreBlast.objects.all().delete()
            for s in Sequence.objects.all():
                s.variant = None
                s.save()

        if options["test"]:
            self.test_curated()
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
        thresholds_scores = {}
        for sequence in Sequence.objects.filter(reviewed=True):
            if 'generic' in sequence.variant.id: continue
            
            hist_type = sequence.variant.hist_type.id
            seqs_file = os.path.join(self.model_evaluation, hist_type, "BLASTDB_sequences_{}.fa".format(sequence.id))
            # self.create_blastdb(Sequence.objects.filter(reviewed=True).exclude(id=sequence.id), seqs_file)
            self.create_blastdb(Sequence.objects.filter(reviewed=True), seqs_file)

            blastp = os.path.join(os.path.dirname(sys.executable), "blastp")
            # output = os.path.join("/", "tmp", "{}.xml".format(uuid.uuid4()))
            blastp_cline = NcbiblastpCommandline(
                cmd=blastp,
                db=seqs_file,
                evalue=.01, outfmt=5)
            # evalue=0.004, outfmt=5)
            result, error = blastp_cline(stdin="\n".join([s.format("fasta") for s in [sequence]]))

            resultFile = io.BytesIO()
            resultFile.write(result.encode("utf-8"))
            resultFile.seek(0)

            for i, blast_record in enumerate(NCBIXML.parse(resultFile)):
                if len(blast_record.alignments) == 0:
                    continue

                best_alignment = blast_record.alignments[0]
                best_hsp = get_best_hsp(best_alignment.hsps)
                best_scores = [best_hsp.score]
                for alignment in blast_record.alignments[1:]:
                    best_alignment_hsp = get_best_hsp(alignment.hsps)
                    best_scores.append(best_alignment_hsp.score)
                    if best_alignment_hsp.score > best_hsp.score:
                        best_alignment = alignment
                        best_hsp = best_alignment_hsp
                # best_scores = [get_best_hsp(alignment.hsps).score for alignment in blast_record.alignments]
                # add_score(sequence, sequence.variant, best_hsp, best_alignment.hit_def.split("|")[0], best=True)
                best_scores.sort(reverse=True)
                self.log.info('DEBUG::{}'.format(best_scores))

                i_treshold = 1
                seq_treshold = best_scores[i_treshold]
                self.log.info('DEBUG::{}-{}-{}'.format(best_scores[i_treshold-1], seq_treshold, best_scores[i_treshold+1]))
                self.log.info('DEBUG::{}-{}'.format(best_scores[i_treshold+1]-seq_treshold, (seq_treshold-best_scores[i_treshold-1])+.00001))
                while (best_scores[i_treshold+1]-seq_treshold)/((seq_treshold-best_scores[i_treshold-1])+.00001) > .8:
                # while best_scores[i_treshold+1]/seq_treshold > .9:
                    i_treshold+=1
                    seq_treshold = best_scores[i_treshold]
                    if i_treshold+1==len(best_scores): break
                self.log.info('DEBUG::i_treshold-{}'.format(i_treshold))
                # seq_treshold = best_scores[2]
                if sequence.variant.id not in thresholds_scores:
                    thresholds_scores[sequence.variant.id] = []
                thresholds_scores[sequence.variant.id].append(seq_treshold)

                plot_scores(sequence.id, save_dir=os.path.join(self.model_evaluation, hist_type),
                            scores_sorted=best_scores, best_threshold=seq_treshold)

        for var_name in thresholds_scores:
            var = Variant.objects.get(id=var_name)
            var.blastthreshold = min(thresholds_scores[var_name])
            var.save()
            # print('{}: {}'.format(var_name, var.blastthreshold))

    def load_curated(self):
        self.log.info('Loading curated...')
        for s in Sequence.objects.filter(reviewed=True):
            s.variant = s.variant_hmm
            s.save()

    def add_scores_for_curated(self):
        # ScoreBlast.objects.filter(sequence__reviewed=True).delete()
        for sequence in Sequence.objects.filter(reviewed=True):
            hist_type = sequence.variant.hist_type.id
            seqs_file = os.path.join(settings.STATIC_ROOT_AUX, "browse", "blast", hist_type, "BLASTDB_sequence_{}.fa".format(sequence.id))
            # self.create_blastdb(Sequence.objects.filter(reviewed=True).exclude(id=sequence.id), seqs_file)
            self.create_blastdb(Sequence.objects.filter(id=sequence.id, reviewed=True), seqs_file)

            blastp = os.path.join(os.path.dirname(sys.executable), "blastp")
            # output = os.path.join("/", "tmp", "{}.xml".format(uuid.uuid4()))
            blastp_cline = NcbiblastpCommandline(
                cmd=blastp,
                db=seqs_file,
                evalue=.01, outfmt=5)
            # evalue=0.004, outfmt=5)
            result, error = blastp_cline(stdin="\n".join([s.format("fasta") for s in [sequence]]))

            resultFile = io.BytesIO()
            resultFile.write(result.encode("utf-8"))
            resultFile.seek(0)

            for i, blast_record in enumerate(NCBIXML.parse(resultFile)):
                if len(blast_record.alignments) == 0:
                    self.log.error('No BLAST record alignments for {}'.format(sequence.id))
                    continue
                if len(blast_record.alignments) > 0:
                    self.log.error('More than 1 BLAST record alignments for {}'.format(sequence.id))

                best_alignment = blast_record.alignments[0]
                add_score(sequence, sequence.variant, get_best_hsp(best_alignment.hsps), best_alignment.hit_def.split("|")[0], best=True)

    def search_blast(self):
        # sequences = Sequence.objects.filter(reviewed=False).values_list('sequence')
        for hist_type in ['H1', 'H2A', 'H2B', 'H3', 'H4']:
            blastdb_file = os.path.join(settings.STATIC_ROOT_AUX, "browse", "blast", hist_type,
                                        "BLASTDB_sequences_{}.fa".format(hist_type))
            self.create_blastdb(Sequence.objects.filter(reviewed=True).filter(variant__hist_type__id=hist_type), blastdb_file)

            sequences_all = [seq.format(format='fasta') for seq in Sequence.objects.filter(reviewed=False).filter(variant_hmm__hist_type__id=hist_type)]
            # sequences_all = [seq.format(format='fasta') for seq in Sequence.objects.filter(reviewed=False)]
            self.log.info("Running BLASTPb for {} sequences of {}...".format(len(sequences_all), hist_type))
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
                load_blast_search(self.blast_file + hist_type + "%d" % i)
                self.log.info('Loaded {}/{} BlastRecords of {}'.format(i,BLAST_PROCS, hist_type))
                self.log.info('Classified {} from {}'.format(Sequence.objects.exclude(variant=None).count(),Sequence.objects.all().count()))

        non_classified = Sequence.objects.filter(variant=None)
        for s in non_classified:
            s.variant = Variant.objects.get(id='generic_{}'.format(s.variant_hmm.hist_type.id))
            s.save()
            score = ScoreBlast(
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


    def search_blast_old(self):
        # sequences = Sequence.objects.filter(reviewed=False).values_list('sequence')
        sequences_all = [seq.format(format='fasta') for seq in Sequence.objects.filter(reviewed=False)]
        self.log.info("Running BLASTP for {} sequences...".format(len(sequences_all)))
        split_count = int(len(sequences_all) / BLAST_PROCS)

        for i in range(BLAST_PROCS + 1):
            sequences = sequences_all[split_count * i:split_count * i + split_count]
            save_to = self.blast_file + "%d" % i
            self.log.info('Starting Blast sequences for {}/{}'.format(i, BLAST_PROCS))
            blastdb_file = os.path.join(settings.STATIC_ROOT_AUX, "browse", "blast", "HistoneDB_sequences.fa")
            make_blastp(sequences, blastdb=blastdb_file, save_to=save_to)

    def load_in_db_old(self):
        self.log.info("Loading BLASTP data into HistoneDB...")
        for i in range(BLAST_PROCS+1):
            load_blast_search(self.blast_file + "%d" % i)
            self.log.info('Loaded {}/{} BlastRecords'.format(i,BLAST_PROCS))
            self.log.info('Classified {} from {}'.format(Sequence.objects.exclude(variant=None).count(),Sequence.objects.all().count()))

        non_classified = Sequence.objects.filter(variant=None)
        for s in non_classified:
            s.variant = Variant.objects.get(id='generic_{}'.format(s.variant_hmm.hist_type.id))
            s.save()
            score = ScoreBlast(
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

    def test_curated(self):
        sequences = [seq.format(format='fasta') for seq in Sequence.objects.filter(reviewed=True)]
        self.log.info("Running BLASTP for {} curated sequences...".format(len(sequences)))

        save_to = self.blast_curated_file
        self.log.info('Starting Blast sequences for {}'.format(len(sequences)))
        make_blastp(sequences, save_to=save_to)

        self.log.info("Loading BLASTP data into HistoneDB...")
        load_blast_search(self.blast_curated_file)
        self.log.info('Loaded {} BlastRecords'.format(len(sequences)))

    def create_blastdb(self, sequences, seqs_file):
        # seqs_file = os.path.join(settings.STATIC_ROOT_AUX, "browse", "blast", "HistoneDB_sequences.fa")
        with open(seqs_file, "w") as seqs:
            for s in sequences:  # here we restrict the blast DB to reviewed seqs
                SeqIO.write(s.to_biopython(ungap=True), seqs, "fasta")

        makeblastdb = os.path.join(os.path.dirname(sys.executable), "makeblastdb")
        subprocess.call(["makeblastdb", "-in", seqs_file, "-dbtype", "prot", "-title", "HistoneDB"])

    def get_stats(self, filename_suff = ''):
        self.log.info('Outputting statistics file ...')

        now = datetime.now()
        dt_string = now.strftime("%Y%m%d-%H%M%S")
        with open('log/blast_db_stat_'+dt_string+filename_suff,'w') as f:
            f.write("Variant database regeneration statistics after BLAST\n")
            f.write("DB regen start time: %s \n"%self.start_time)
            f.write("DB regen end time: %s\n"%now)
            f.write("Time taken for regeneration of variants: %f hours\n"%(float((now-self.start_time).total_seconds())/3600.))
            f.write("Parallel threads used %d\n"%BLAST_PROCS)
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


