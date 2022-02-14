import os, subprocess
import json
import logging

from Bio import SeqIO

from utils import *
from tools.path_variables import *
from tools.general_utils import get_seeds
'''
This model will:
    - create HMMs from curated set
    - extract and classify sequences from db using HMMs
    - classify sequences using BLAST
    - test model
'''

log = logging.getLogger(__name__)

class HistonedbClassifier(object):
    def __init__(self, classification_tree=None):
        self.classification_tree = classification_tree
        self.combined_hmm_file = None
        self.blast_dbs_dir = None
        self.predicted_types = None
        self.predicted_variants = None
        self.types_prediction_info = None
        self.variants_prediction_info = None

        if not self.classification_tree:
            try:
                with open(VARIANTS_JSON) as f:
                    variant_json = json.loads(f.read())
                    self.classification_tree = variant_json['tree']
            except FileNotFoundError:
                raise FileNotFoundError(f'{VARIANTS_JSON} is not found. Please contact with technical support.')

    def create_hmms(self, seed_directory=None, hmm_directory=None, force=False):
        log.info("Building HMMs...")
        if self.combined_hmm_file and not force: return self
        seed_directory = seed_directory if seed_directory else SEED_CORE_DIRECTORY
        hmm_directory = hmm_directory if hmm_directory else HMM_DIRECTORY
        self.combined_hmm_file = os.path.join(hmm_directory, 'combined_hmm', 'types_combined.hmm')
        with open(self.combined_hmm_file, "w") as combined_hmm:
            for hist_type in self.classification_tree:
                # Build HMMs
                hmm_dir = os.path.join(HMM_DIRECTORY, hist_type)
                if not os.path.exists(hmm_dir): os.makedirs(hmm_dir)
                hmm_file = os.path.join(hmm_dir, f"{hist_type}.hmm")
                log.info(' '.join(["hmmbuild", "-n", hist_type, hmm_file, os.path.join(seed_directory, f'{hist_type}.fasta')]))
                subprocess.call(["hmmbuild", "-n", hist_type, hmm_file, os.path.join(seed_directory, f'{hist_type}.fasta')])
                with open(hmm_file) as hmm: print(hmm.read().rstrip(), file=combined_hmm)
        # Pressing HMMs
        log.info("Pressing HMMs for histone types...")
        subprocess.call(["hmmpress", "-f", self.combined_hmm_file])
        return self

    def create_blastdbs(self, seed_directory=None, blast_dbs_dir=None, force=False):
        log.info("Building BLASTDBs...")
        if self.blast_dbs_dir and not force: return self
        self.blast_dbs_dir = blast_dbs_dir if blast_dbs_dir else BLASTDBS_DIR
        seed_directory = seed_directory if seed_directory else SEED_DIRECTORY
        if not os.path.exists(self.blast_dbs_dir): os.makedirs(self.blast_dbs_dir)
        for hist_type in self.classification_tree:
            with open(os.path.join(self.blast_dbs_dir, f"BLASTDB_sequences_{hist_type}.fa"), "w") as seqs:
                for s in list(SeqIO.parse(os.path.join(seed_directory, f'{hist_type}.fasta'), "fasta")):
                    SeqIO.write(s.seq.ungap("-"), seqs, "fasta")
            log.info(" ".join(["makeblastdb", "-in", os.path.join(self.blast_dbs_dir, f"{hist_type}.fa"), "-dbtype", "prot", "-title", "HistoneDB"]))
            subprocess.call(["makeblastdb", "-in", os.path.join(self.blast_dbs_dir, f"{hist_type}.fa"), "-dbtype", "prot", "-title", "HistoneDB"])
        return self

    def predict_type(self, sequences, hmmout=None, E=10, save_to=None):
        if self.predicted_types: return self.predicted_types
        hmmout = hmmout if hmmout else DB_HISTTYPE_RESULTS_FILE
        self.types_prediction_info = predict_types(sequences=sequences, hmmout=hmmout, hmmdb=self.combined_hmm_file,
                                                   E=E, result_file=save_to)
        self.predicted_types = [{k: v for k, v in d.iteritems() if k in ['accession','histone_type']} for d in list(filter(lambda d: d['best'], self.types_prediction_info))]
        return self.predicted_types

    def predict_variant(self, sequences, hmmout=None, E_hmm=10,
                        blastout=None, E_bast=10, save_to=None):
        if self.predicted_variants: return self.predicted_variants
        hmmout = hmmout if hmmout else DB_HISTTYPE_RESULTS_FILE
        self.types_prediction_info = predict_types(sequences=sequences, hmmout=hmmout, hmmdb=self.combined_hmm_file)
        #
        blastout = blastout if blastout else DB_HISTVARIANTS_BLAST_RESULTS_FILE
        self.variants_prediction_info = predict_variants(sequences,blastout=blastout,blastdb=self.blast_db_file,result_file=save_to)
        return self.predicted_variants
    

