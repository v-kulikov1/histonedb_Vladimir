import os
from django.conf import settings

# Working paths
LOG_DIRECTORY                       = os.path.join("log")
LOG_DB_STAT_DIRECTORY               = os.path.join("log", "db_stat")
DATA_DIRECTORY                      = os.path.join("data")
SEED_DIRECTORY                      = os.path.join(DATA_DIRECTORY, "seeds")
SEED_CORE_DIRECTORY                 = os.path.join(DATA_DIRECTORY, "seeds_core")
PREDICTION_DIRECTORY                = os.path.join(DATA_DIRECTORY, "PREDICTION")
HMM_DIRECTORY                       = os.path.join(PREDICTION_DIRECTORY, "hmms")
HMM_MODEL_EVALUATION                = os.path.join(HMM_DIRECTORY, "model_evaluation")
BLAST_DIRECTORY                     = os.path.join(PREDICTION_DIRECTORY, "blast")
BLAST_MODEL_EVALUATION              = os.path.join(BLAST_DIRECTORY, "model_evaluation")
BLASTDBS_DIR                        = os.path.join(BLAST_DIRECTORY, 'blastdbs')

VARIANTS_JSON                       = os.path.join(DATA_DIRECTORY, 'classification.json')
CURATED_ALL_FASTA                   = os.path.join(HMM_DIRECTORY, 'model_evaluation',"curated.fasta")
CURATED_GENERICLESS_FASTA           = os.path.join(HMM_DIRECTORY, 'model_evaluation', "curated_genericless.fasta")
CURATED_GENERIC_FASTA               = os.path.join(HMM_DIRECTORY, 'model_evaluation', "curated_generic.fasta")

DB_HISTTYPE_RESULTS_FILE            = os.path.join(HMM_DIRECTORY, "results", "histone_type", "search{}.out")
CURATED_HISTTYPE_RESULTS_FILE       = os.path.join(HMM_DIRECTORY, "results", "histone_type", "curated_search.out")
DB_HISTTYPE_PARSED_RESULTS_FILE     = os.path.join(HMM_DIRECTORY, "results", "histone_type", "search_parsed{}.csv")
IDS_FILE                            = os.path.join(HMM_DIRECTORY, "results", "histone_type", "extracted_sequences{}.ids")
FULL_LENGTH_SEQS_FILE               = os.path.join(HMM_DIRECTORY, "results", "histone_type", "HistoneDB_sequences{}.fasta")
DB_HISTVARIANTS_HMM_RESULTS_FILE    = os.path.join(HMM_DIRECTORY, "results", "histone_variants", "search{}.out")
CURATED_HISTVAR_RESULTS_FILE        = os.path.join(HMM_DIRECTORY, "results", "histone_variants", "curated_search.out")
CURATED_GEN_HISTVAR_RESULTS_FILE    = os.path.join(HMM_DIRECTORY, "results", "histone_variants", "curated_generic_search.out")
TYPE_CLASSIFICATION_REPORT_FILE     = os.path.join(HMM_MODEL_EVALUATION, "type_classification_report.txt")
DB_HISTVARIANTS_BLAST_RESULTS_FILE  = os.path.join(BLAST_DIRECTORY, "results", "{}", "blast_search.out")
DB_HISTVARIANTS_PARSED_RESULTS_FILE = os.path.join(BLAST_DIRECTORY, "results", "{}", "blast_search_parsed.out")
