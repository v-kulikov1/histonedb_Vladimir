import os, subprocess
import json
import logging

from Bio import SeqIO

from prediction.utils import *
from tools.path_variables import *
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
        self.blast_db_file = None
        self.predicted_types = None
        self.predicted_variants = None
        self.types_prediction_info = None
        self.variants_prediction_info = None

        self.hmmout = None
        self.blastout = None

        if not self.classification_tree:
            with open(VARIANTS_JSON) as f:
                variant_json = json.loads(f.read())
                self.classification_tree = variant_json['tree']

    def create_hmms(self, seed_directory=None):
        log.info("Building HMMs...")
        seed_directory = seed_directory if seed_directory else SEED_CORE_DIRECTORY
        if not os.path.exists(COMBINED_HMM_DIRECTORY): os.makedirs(COMBINED_HMM_DIRECTORY)
        self.combined_hmm_file = os.path.join(COMBINED_HMM_DIRECTORY, 'types_combined.hmm')
        with open(self.combined_hmm_file, "w") as combined_hmm:
            for hist_type in self.classification_tree:
                # Build HMMs
                # hmm_dir = os.path.join(HMM_DIRECTORY, hist_type)
                # if not os.path.exists(hmm_dir): os.makedirs(hmm_dir)
                hmm_file = os.path.join(HMM_DIRECTORY, f"{hist_type}.hmm")
                log.info(' '.join(["hmmbuild", "-n", hist_type, hmm_file, os.path.join(seed_directory, f'{hist_type}.fasta')]))
                subprocess.call(["hmmbuild", "-n", hist_type, hmm_file, os.path.join(seed_directory, f'{hist_type}.fasta')])
                with open(hmm_file) as hmm: print(hmm.read().rstrip(), file=combined_hmm)
        # Pressing HMMs
        log.info("Pressing HMMs for histone types...")
        subprocess.call(["hmmpress", "-f", self.combined_hmm_file])
        return self

    def create_blastdbs(self, seed_directory=None):
        log.info("Building BLASTDBs...")
        seed_directory = seed_directory if seed_directory else SEED_DIRECTORY
        if not os.path.exists(BLASTDBS_DIR): os.makedirs(BLASTDBS_DIR)
        self.blast_db_file = os.path.join(BLASTDBS_DIR, "{}.fa")
        for hist_type in self.classification_tree:
            with open(self.blast_db_file.format(hist_type), "w") as seqs:
                for s in list(SeqIO.parse(os.path.join(seed_directory, f'{hist_type}.fasta'), "fasta")):
                    SeqIO.write(s.seq.ungap("-"), seqs, "fasta")
            log.info(" ".join(["makeblastdb", "-in", self.blast_db_file.format(hist_type), "-dbtype", "prot", "-title", "HistoneDB"]))
            subprocess.call(["makeblastdb", "-in", self.blast_db_file.format(hist_type), "-dbtype", "prot", "-title", "HistoneDB"])
        return self

    def predict_type(self, sequences, E=10):
        if self.predicted_types: return self.predicted_types
        self.hmmout = os.path.join(PREDICTION_RESULTS_DIRECTORY, "hmm_search.out")
        self.types_prediction_info = predict_types(sequences=sequences, hmmout=self.hmmout, hmmdb=self.combined_hmm_file, E=E)
        self.predicted_types = [{k: v for k, v in d.iteritems() if k in ['accession','type', 'description']} for d in list(filter(lambda d: d['best'], self.types_prediction_info))]
        return self.predicted_types

    def predict_variant(self, sequences, E_hmm=10, E_blast=10):
        if self.predicted_variants: return self.predicted_variants
        self.hmmout = os.path.join(PREDICTION_RESULTS_DIRECTORY, "hmm_search.out")
        self.types_prediction_info = predict_types(sequences=sequences, hmmout=self.hmmout, hmmdb=self.combined_hmm_file, E=E_hmm)
        self.predicted_types = [{k: v for k, v in d.iteritems() if k in ['accession', 'type', 'description']}
                                for d in list(filter(lambda d: d['best'], self.types_prediction_info))]
        #
        self.blastout = os.path.join(PREDICTION_RESULTS_DIRECTORY, "blast_search_{}.out")
        self.variants_prediction_info = []
        for hist_type in self.classification_tree:
            tmpfile_ids = os.path.join(PREDICTION_RESULTS_DIRECTORY, f"extracted_sequences_{hist_type}.ids")
            full_length_seq_file = os.path.join(PREDICTION_RESULTS_DIRECTORY, f"extracted_sequences_{hist_type}.fasta")
            with open(tmpfile_ids, 'w') as ids_file:
                for p in list(filter(lambda d: d['type']==hist_type, self.predicted_types)):
                    ids_file.write('{}\n'.format(p['accession']))
            extract_full_sequences(tmpfile_ids, sequences, full_length_seq_file)
            seqrec_sequences = list(SeqIO.parse(full_length_seq_file, "fasta"))
            pred_variants = predict_variants(seqrec_sequences, blastout=self.blastout.format(hist_type),
                                             blastdb=self.blast_db_file.format(hist_type), E=E_blast)
            self.variants_prediction_info.append(list(map(lambda x: x.update({'type': hist_type}), pred_variants)))
            # os.remove(tmpfile_ids)
            # os.remove(full_length_seq_file)
        self.predicted_types = [{k: v for k, v in d.iteritems() if k in ['accession', 'type', 'variant', 'description']}
                                for d in list(filter(lambda d: d['best'], self.types_prediction_info))]
        return self.predicted_variants

    def save_types_prediction_info(self, file_name):
        import pandas as pd
        pd.DataFrame(self.types_prediction_info).fillna('').drop(columns=['hsp']).to_csv(file_name, index=False)
        log.info(f"Predicted sequences saved to {file_name}")

    def save_variants_prediction_info(self, file_name):
        import pandas as pd
        pd.DataFrame(self.variants_prediction_info).fillna('').drop(columns=['hsp']).to_csv(file_name, index=False)
        log.info(f"Predicted sequences saved to {file_name}")
    

