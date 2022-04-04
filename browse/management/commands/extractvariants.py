from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from browse.models import Histone, Variant, Sequence, ScoreHmm, ScoreForHistoneType, Feature
from tools.load_hmmsearch import load_hmm_results, add_histone_score, add_generic_score, \
    get_many_prot_seqrec_by_accession, load_hmm_classification_results, load_generic_scores
from tools.test_model import test_model
import subprocess
import os, sys
import re
import io
from tools.taxonomy_from_accessions import taxonomy_from_header, easytaxonomy_from_header, fetch_taxids, update_taxonomy
from tools.blast_search import make_blastp, load_blast_search
from tools.hist_ss import get_variant_features
from tools.browse_service import *

from Bio import SearchIO
from Bio import SeqIO
from tqdm import tqdm

import logging
from datetime import date, datetime

from cProfile import Profile

# This command is the main one in creating the histone database system from seed alignments
# and by using HMMs constructed based on these alignment to classify the bigger database.
# see handle() for the workflow description.
# INPUT NEEDED: static/browse/seeds
# db_file (nr) in the main directory - the database of raw sequences.

class Command(BaseCommand):
    help = 'Build HistoneDB by first loading the seed sequences and then parsing the database file'

    # Logging info
    logging.basicConfig(filename=os.path.join(LOG_DIRECTORY, "extractvariants.log"),
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
            "--db",
            dest="db_file",
            default="nr",
            help="Specify the database file, by default will use or download nr")

        parser.add_argument(
            "--adjust_hmmer_procs",
            default=20,
            help="Adjust hmmer_procs")

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
        self.log.info('===             extractvariants START               ===')
        self.log.info('=======================================================')
        self.start_time = datetime.now()

        ## If hmmer_procs in options change the value
        if 'adjust_hmmer_procs' in options:
            HMMER_PROCS = options["adjust_hmmer_procs"]
        self.log.info('HMMER_PROCS !!!! {}'.format(HMMER_PROCS))

        ##If no nr file is present in the main dir, will download nr from the NCBI ftp.
        self.db_file = options['db_file']
        
        if self.db_file == "nr":
            if options["force"] or not os.path.isfile('nr'):
                self.get_nr()
        if self.db_file == "swissprot":
            if options["force"] or not os.path.isfile('swissprot'):
                self.get_swissprot()
        if ('http://' in self.db_file) or ('https://' in self.db_file) or ('ftp://' in self.db_file):
            self.log.info(
                'Provided db file is a link %s - Expecting a gzipped file. Attempting to download ...' % self.db_file)
            subprocess.call(["wget", self.db_file, '-O', 'db.gz'])
            subprocess.call(["gunzip", "db.gz"])
            self.db_file = 'db'
            
        self.set_nr_version()
        
        if options["force"]:
            # Clean the DB, removing all sequence/scores/etc
            Sequence.objects.all().delete()
            ScoreForHistoneType.objects.all().delete()
            # Load our curated sets taken from seed alignments into the database an run classification algorithm
            self.load_curated()
            self.load_scores_for_curated()

        self.load_from_db_parallel(force=options["force"])
        self.get_stats()

        seq_num = Sequence.objects.count()
        seqauto_num = Sequence.objects.filter(reviewed=False).count()

        self.log.info(' The database has %d sequences now !!!' % seq_num)
        self.log.info(' %d sequences came from automatic search !!!' % seqauto_num)
        self.log.info('=======================================================')
        self.log.info('===      extractvariants SUCCESSFULLY finished      ===')
        self.log.info('=======================================================')

    def handle(self, *args, **options):
        if options['profile']:
            profiler = Profile()
            profiler.runcall(self._handle, *args, **options)
            profiler.print_stats()
        else:
            self._handle(*args, **options)


    # get Database
    def get_nr(self):
        """Download nr if not present"""
        if not os.path.isfile(self.db_file):
            self.log.info("Downloading nr...")
            # print >> self.stdout, "Downloading nr..."
            with open("nr.gz", "w") as nrgz:
                subprocess.call(["curl", "-#", "ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz"], stdout=nrgz)
            subprocess.call(["gunzip", "nr.gz"])

    def get_swissprot(self):
        """Download nr if not present"""
        if not os.path.isfile(self.db_file):
            self.log.info("Downloading swissprot...")
            # print >> self.stdout, "Downloading swissprot..."
            with open("swissprot.gz", "w") as swissprotgz:
                subprocess.call(["curl", "-#", "ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/swissprot.gz"],
                                stdout=swissprotgz)
            subprocess.call(["gunzip", "swissprot.gz"])

    def set_nr_version(self):
        with open(NR_VERSION_FILE, 'w') as nrv:
            nrv.write(self.db_file)


    # Curated
    def load_curated(self):
        """
        Extracts sequences from seed alignments in static/browse/seeds
        Loads them into the database with flag reviewed=True (which means curated)
        An important fact:
        the seqs in seeds, should have a special header currently:
        >Ixodes|XP_002403551.1|macroH2A Ixodes_macroH2A
        we accept only this patterns to extract ACCESSIONs
        """
        accessions = []
        for hist_type, seed in get_seeds(generic=True):
            variant_name = seed[:-6]
            self.log.info(' '.join([variant_name, "==========="]))
            seed_aln_file = os.path.join(SEED_DIRECTORY, hist_type, seed)
            for s in SeqIO.parse(seed_aln_file, "fasta"):
                s.seq = s.seq.ungap("-")
                accession = s.id.split("|")[1]
                # if accession.startswith("NOGI"):
                #     self.log.info("NO GI detected {}".format(s.id))
                #     taxid= easytaxonomy_from_header(s.id).id
                # else:
                #     #trick to make taxid retrieval faster
                #     # taxonomy = taxonomy_from_header("", gi=gi)
                #     taxid=1
                #     accessions.append(accession)
                # trick to make taxid retrieval faster
                # taxonomy = taxonomy_from_header("", gi=gi)
                taxid = 1
                accessions.append(accession)
                self.log.info("Loading {}".format(s.id))

                variant_model = Variant.objects.get(id=variant_name)
                seq = Sequence(
                    id=accession,
                    variant=variant_model,
                    variant_hmm=variant_model,
                    histone_type = Histone.objects.get(id=hist_type),
                    gene=None,
                    splice=None,
                    taxonomy_id=taxid,
                    header="CURATED SEQUENCE: {}".format(s.description),
                    sequence=s.seq,
                    reviewed=True,
                )
                seq.save()

        # Now let's lookup taxid for those having ACCESSIONs via NCBI.
        update_taxonomy(accessions)

    def load_scores_for_curated(self):
        """
        Load the scores of histone types for curated sequences
        """
        from sklearn.metrics import classification_report

        # collecting all sequences and inserting taget histone type at the end of description
        sequences = CURATED_ALL_FASTA
        if not os.path.isfile(sequences):
            with open(sequences, "w") as f:
                for hist_type, seed in get_seeds(generic=True):
                    seed_aln_file = os.path.join(SEED_DIRECTORY, hist_type, seed)
                    for s in SeqIO.parse(seed_aln_file, "fasta"):
                        s.seq = s.seq.ungap("-")
                        s.description = s.description + '|target={}'.format(hist_type)
                        SeqIO.write(s, f, "fasta")

        # prediction
        predicted_result = self.extract_predict_histtypes(sequences=sequences, hmmout=CURATED_HISTTYPE_RESULTS_FILE,
                                                          hmmdb=COMBINED_HMM_HISTTYPES_FILE, E=10,
                                                          result_file=os.path.join(HMM_DIRECTORY, 'curated_parsed_results.csv'),
                                                          load_to_db=False)

        # adding scores
        for i, pred_res in predicted_result.iterrows():
            seqs = Sequence.objects.filter(id=pred_res['accession'].split('|')[1])
            if len(seqs) > 0:
                seq = seqs.first()
                try:
                    add_histone_score(
                        seq,
                        Histone.objects.get(id=pred_res['histone_type']),
                        pred_res['hsp'],
                        best=pred_res['best'])
                except Exception as e:
                    self.log.warning('Failed to add histone type score: {}'.format(str(e)))
                    pass
            else:
                self.log.error('There is no curated sequence in HistoneDB with accession {}'.format(pred_res['accession'].split('|')[2]))

        # saving classification report
        predicted_result = predicted_result[predicted_result['best']]
        target = [pred_res['description'].split('target=')[1] for i, pred_res in predicted_result.iterrows()]
        predicted_result['target_histone_type'] = target
        with open(TYPE_CLASSIFICATION_REPORT_FILE, 'w') as type_class_report_file:
            type_class_report_file.write(classification_report(predicted_result['target_histone_type'], predicted_result['histone_type']))


    # extract variants parallel

    def load_from_db_parallel(self, force=True):
        """Use HMMs to search the nr database if force and load data into the histone database"""
        self.log.info("Searching HMMs in parallel...")
        accessions = self.extract_predict_histtypes(sequences=self.db_file, hmmout=DB_HISTTYPE_RESULTS_FILE,
                                                    hmmdb=COMBINED_HMM_HISTTYPES_FILE, E=10,
                                                    result_file=DB_HISTTYPE_PARSED_RESULTS_FILE, do_search=force,
                                                    procs=HMMER_PROCS)
        for i, acc_list in enumerate(accessions):
            with open(IDS_FILE.format('%02d' % (i + 1)), 'w') as ids_file:
                for acc in acc_list:
                    ids_file.write('{}\n'.format(acc))
        self.extract_full_sequences_parallel()

    def extract_full_sequences_parallel(self, sequences=None):
        """Create database to extract full length sequences"""

        if sequences is None:
            sequences = self.db_file

        child_procs=[]
        for i in range(HMMER_PROCS):
            #1) Create and index of sequence file
            self.log.info("Indexing sequence database "+"db_split/split%02d"%(i+1))
            self.log.info(" ".join(["esl-sfetch", "--index", "db_split/split%02d"%(i+1)]))
            p=subprocess.Popen(["esl-sfetch", "--index", "db_split/split%02d"%(i+1)])
            child_procs.append(p)
        for cp in child_procs:
            cp.wait()


        child_procs=[]
        for i in range(HMMER_PROCS):

            #2) Extract all ids
            self.log.info("Extracting full length sequences...")
            self.log.info(" ".join(["esl-sfetch", "-o", FULL_LENGTH_SEQS_FILE.format('%02d' % (i + 1)), "-f", "db_split/split%02d"%(i+1), IDS_FILE.format('%02d' % (i + 1))]))
            p=subprocess.Popen(["esl-sfetch", "-o", FULL_LENGTH_SEQS_FILE.format('%02d' % (i + 1)), "-f", "db_split/split%02d"%(i+1), IDS_FILE.format('%02d' % (i + 1))])
            child_procs.append(p)
        for cp in child_procs:
            cp.wait()


        #3) Update sequences with full length NR sequences -- is there a faster way?
        self.log.info("Updating records with full length sequences...")
        counter=0
        counter_dne=0
        for i in range(HMMER_PROCS):
            self.log.info("Updating sequences from file {}".format(FULL_LENGTH_SEQS_FILE.format('%02d' % (i + 1))))

            for record in SeqIO.parse(FULL_LENGTH_SEQS_FILE.format('%02d' % (i + 1)), "fasta"):
                headers = record.description.split('\x01')
                for header in headers:
                    accession = header.split(" ")[0]
                    try:
                        seq = Sequence.objects.get(id=accession)
                        # self.log.info("Updating sequence: {}".format(seq.description))
                        seq.sequence = str(record.seq)
                        seq.save()
                        counter=counter+1
                    except Sequence.DoesNotExist:
                        counter_dne=counter_dne+1
                        #These seqs likely did not exceed the threshold.
                        # self.log.error("Strangely sequence %s does not exist in database - unable to update"%gi)
                        pass
        self.log.info("Updated %d sequences"%counter)
        self.log.info("%d sequences where attempded to update, but were not found in the database"%counter_dne)

    def extract_full_sequences_from_ncbi(self):
        """Exract full seq by direct call to NCBI servers"""
        self.log.info("Getting full sequences of automatically annotated proteins from NCBI====")
        accessions = Sequence.objects.filter(reviewed=False).values_list('id', flat=True)
        fasta_dict = get_many_prot_seqrec_by_accession(accessions)

        # 3) Update sequences with full length NR sequences -- is there a faster way?
        for accession, record in tqdm(fasta_dict.items()):
            # self.log.info('::DEBUG::buildvariants:: record:\n{}\n'.format(record))
            # headers = record.description.split(" >")
            # for header in headers:
            #     accession = header.split(" ", 1)[0]
            #     # self.log.info('::DEBUG::buildvariants:: accession: {}'.format(accession))
            #     try:
            #         seq = Sequence.objects.get(id=accession)
            #         seq.sequence = str(record.seq)
            #         seq.save()
            #     except Sequence.DoesNotExist:
            #         self.log.error('Sequence with accession {} does not exist in DB.'.format(accession))
            #         pass
            try:
                seq = Sequence.objects.get(id=accession)
                seq.sequence = str(record.seq)
                seq.save()
            except Sequence.DoesNotExist:
                self.log.error('Sequence with accession {} does not exist in DB.'.format(accession))
                pass


    # Algorithm
    def extract_predict_histtypes(self, sequences, hmmout, hmmdb, E=10, result_file=None, do_search=True, load_to_db=True, procs=1):
        """
        This method do a search among given sequences and extracts just histone variants classified by histone types.
        All the results saves to the result_file, e.g.:
        accession, histone_type, score, best, description
        NP_563627.1, H3, 2.14, True, Histone superfamily protein [Arabidopsis thaliana]
        DAA13058.1, H2B, 3.11, False, [Bos taurus]
        Parameters:
        ___________
        sequences: fasta file
        hmmout: hmmout file
        hmmdb: hmmdb file
        result_file: file to save results
        """

        if procs > 1:
            self.db_split(sequences, split_num=procs)
            # sequences = ["db_split/split%02d" % (i+1) for i in range(procs)]
            sequences = "db_split/split%02d"

        if do_search:
            self.search(hmms_db=hmmdb, out=hmmout, sequences=sequences, procs=procs, E=E)

        if procs == 1:
            self.log.info("Parsing search results...")
            result = self.parse_search_out(hmmout=hmmout)
            if load_to_db:
                self.load_hmm_results(result)
            result = pd.DataFrame(result)
            if result_file:
                result.to_csv(result_file, index=False)
                self.log.info("Predicted sequences saved to {}".format(result_file))
            return result

        accessions = []
        for i in range(procs):
            self.log.info("Parsing search results from {}...".format(hmmout.format('%02d' % (i + 1))))
            result = self.parse_search_out(hmmout=hmmout.format('%02d' % (i + 1)))
            accessions.append(set([r['accession'] for r in result]))
            # accessions.append([header.split(" ")[0] for r in result for header in "{} {}".format(r['accession'], r['description']).split('\x01')])
            if load_to_db:
                self.load_hmm_results(result)
            result = pd.DataFrame(result)
            if result_file:
                result.to_csv(result_file.format('%02d' % (i + 1)), index=False)
                self.log.info("Predicted sequences saved to {}".format(result_file.format('%02d' % (i + 1))))
        return accessions

    def search(self, hmms_db, out, sequences, procs, E=10):
        """Use HMMs to search the nr database"""
        self.log.info("Searching HMMs...")
        self.log.info("Launching %d processes" % procs)
        child_procs = []
        for i in range(procs):
            try:
                self.log.info(" ".join(["nice", "hmmsearch", "-o", out.format('%02d' % (i + 1)), "-E", str(E), "--notextw", hmms_db, sequences % (i + 1)]))
                p=subprocess.Popen(["nice", "hmmsearch", "-o", out.format('%02d' % (i + 1)), "-E", str(E), "--notextw", hmms_db, sequences % (i+1)])
            except TypeError: #if procs=1 we have sequences as string without formating %02d st the end
                self.log.info(" ".join(["nice", "hmmsearch", "-o", out, "-E", str(E), "--notextw", hmms_db, sequences]))
                p = subprocess.Popen(["nice", "hmmsearch", "-o", out, "-E", str(E), "--notextw", hmms_db, sequences])
            child_procs.append(p)
        for cp in child_procs:
            cp.wait()

    def parse_search_out(self, hmmout):
        """Parse hmmout file and return results as dict.
          Parameters:
          ___________
          hmmout: hmmout file
          return e.g.: [{'accession': NP_563627.1, 'histone_type': H3, 'score': 2.14, 'best': True,
                    'description': Histone superfamily protein [Arabidopsis thaliana], 'hsp': best HSPObject},
                    {'accession': DAA13058.1, 'histone_type': H2B, 'score': 3.11, 'best': False,
                    'description': [Bos taurus], 'hsp': best HSPObject},]
          """
        result = []
        count_new = 0
        for histone_query in tqdm(SearchIO.parse(hmmout, "hmmer3-text")):
            self.log.info("Loading histone: {}".format(histone_query.id))  # histone_query.id is histone type
            for hit in tqdm(histone_query):
                best_hsp = max(hit, key=lambda hsp: hsp.bitscore)
                is_new = True
                for res in result:
                    if hit.id != res['accession'] or not res['best']: continue
                    is_new = False
                    if best_hsp.bitscore > res['score']:
                        res['best'] = False
                        result.append({'accession': hit.id,  # hit.id is a first accession
                                       'histone_type': histone_query.id,
                                       'score': best_hsp.bitscore,
                                       'best': True,
                                       'description': hit.description,
                                       'hsp': best_hsp,
                                       })
                    else:
                        result.append({'accession': hit.id,  # hit.id is a first accession
                                       'histone_type': histone_query.id,
                                       'score': best_hsp.bitscore,
                                       'best': False,
                                       'description': hit.description,
                                       'hsp': best_hsp,
                                       })
                if is_new:
                    count_new += 1
                    result.append({'accession': hit.id,  # hit.id is a first accession
                                   'histone_type': histone_query.id,
                                   'score': best_hsp.bitscore,
                                   'best': True,
                                   'description': hit.description,
                                   'hsp': best_hsp,
                                   })
        self.log.info("From file %s we got %d sequences" % (hmmout, count_new))
        return result

    def load_hmm_results(self, hmm_results):
        """Save hmm_results given as list of dicts to HistoneDB.
          Parameters:
          ___________
          hmm_results: list of dists
            e.g.: [{'accession': NP_563627.1, 'histone_type': H3, 'score': 2.14, 'best': True,
                    'description': Histone superfamily protein [Arabidopsis thaliana]},
                    {'accession': DAA13058.1, 'histone_type': H2B, 'score': 3.11, 'best': False,
                    'description': [Bos taurus]},]
          """

        self.log.info("Loading prediction results to HistoneDB...")
        for hmmres in hmm_results:
            # Below we are fetching a list of headers if there are multiple headers for identical sequences
            # Technically HUMMER might put the second and on accessions in description column.
            # The format should be strictly the genbank format: XP_027215218.1 histone H1.5 [Pan troglodytes]PNI76178.1 HIST1H1B isoform 1 [Pan troglodytes]
            # print('!!!!!!!!!!!!!!!!!!!!!!!!!!')
            # print("{} {}".format(hmmres['accession'], hmmres['description']))
            # print("{} {}".format(hmmres['accession'], hmmres['description']).split('\x01'))
            for header in "{} {}".format(hmmres['accession'], hmmres['description']).split('\x01'):
                accession = header.split(" ")[0]
                seqs = Sequence.objects.filter(id=accession)
                if len(seqs) > 0:
                    seq = seqs.first()
                    if hmmres['best']:
                        seq.histone_type = Histone.objects.get(id=hmmres['histone_type'])
                        seq.save()
                else:
                    seq = add_sequence(
                        accession,
                        Histone.objects.get(id=hmmres['histone_type']),
                        taxonomy_from_header(header, accession),
                        header,
                        '')
                try:
                    add_histone_score(
                        seq,
                        Histone.objects.get(id=hmmres['histone_type']),
                        hmmres['hsp'],
                        best=hmmres['best'])
                except Exception as e:
                    self.log.warning('Failed to add histone type score: {}'.format(str(e)))
                    pass
        self.log.info("Initiating taxonomy update for %d seqs where it is not identified" % (
            Sequence.objects.filter(taxonomy__name="unidentified").count()))
        # Now let's lookup taxid for those we could not pare from header using NCBI eutils.
        update_taxonomy(Sequence.objects.filter(taxonomy__name="unidentified").values_list("id", flat=True))
        self.log.info("Total loaded %d sequences"%(Sequence.objects.all().count()))

    def db_split(self, sequences, split_num=1, path_name='db_split'):
        self.log.info("Splitting database file into %d parts" % split_num)
        # !!!!!!!!!!!!!!!
        os.system('rm -rf db_split')
        os.system('mkdir db_split')
        # This is tricky tricky to make it fast
        size = os.path.getsize(sequences)
        split_size = int(size / split_num) + 1
        os.system('split --bytes=%d --numeric-suffixes=1 %s db_split/split' % (split_size, sequences))
        # We need to heal broken fasta records
        for i in range(1, split_num + 1):
            for k in range(1, 10000):
                outsp = subprocess.check_output(['head', '-n', '%d' % k, 'db_split/split%02d' % i])
                if (outsp.split(b'\n')[-2].decode("utf-8")[0] == '>'):
                    print(k)
                    break
            if (k > 1):
                os.system('head -n %d db_split/split%02d >>db_split/split%02d ' % (k - 1, i, i - 1))
                os.system('tail -n +%d db_split/split%02d> db_split/temp' % (k, i))
                os.system('mv db_split/temp db_split/split%02d' % i)


    # Statistics
    def get_stats(self):
        get_stats(start_time=self.start_time, filename='extractvariants')

