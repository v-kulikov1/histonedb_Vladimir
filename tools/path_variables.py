import os
from django.conf import settings

# Working paths
LOG_DIRECTORY                    = os.path.join("log")
LOG_DB_STAT_DIRECTORY            = os.path.join("log", "db_stat")
DATA_DIRECTORY                   = os.path.join("data")
SEED_DIRECTORY                   = os.path.join(DATA_DIRECTORY, "seeds")
SEED_CORE_DIRECTORY              = os.path.join(DATA_DIRECTORY, "seeds_core")
PREDICTION_DIRECTORY             = os.path.join("prediction")
HMM_DIRECTORY                    = os.path.join(PREDICTION_DIRECTORY, "hmms")
PREDICTION_RESULTS_DIRECTORY     = os.path.join(PREDICTION_DIRECTORY, "results")
COMBINED_HMM_DIRECTORY           = os.path.join(HMM_DIRECTORY, "combined_hmm")
HMM_MODEL_EVALUATION_DIRECTORY   = os.path.join(HMM_DIRECTORY, "model_evaluation")
BLAST_DIRECTORY                  = os.path.join(PREDICTION_DIRECTORY, "blast")
BLAST_MODEL_EVALUATION_DIRECTORY = os.path.join(BLAST_DIRECTORY, "model_evaluation")
BLASTDBS_DIR                     = os.path.join(BLAST_DIRECTORY, 'blastdbs')

VARIANTS_JSON                    = os.path.join(DATA_DIRECTORY, 'classification.json')
CURATED_ALL_FASTA                = os.path.join(HMM_MODEL_EVALUATION_DIRECTORY,"curated.fasta")
CURATED_GENERICLESS_FASTA        = os.path.join(HMM_MODEL_EVALUATION_DIRECTORY, "curated_genericless.fasta")
CURATED_GENERIC_FASTA            = os.path.join(HMM_MODEL_EVALUATION_DIRECTORY, "curated_generic.fasta")

# HISTTYPE_RESULTS_FILE            = os.path.join(PREDICTION_RESULTS_DIRECTORY, "histone_type", "search{}.out")
# CURATED_HISTTYPE_RESULTS_FILE    = os.path.join(PREDICTION_RESULTS_DIRECTORY, "histone_type", "curated_search.out")
# HISTTYPE_PARSED_RESULTS_FILE     = os.path.join(PREDICTION_RESULTS_DIRECTORY, "histone_type", "search_parsed{}.csv")
# TYPE_CLASSIFICATION_REPORT_FILE  = os.path.join(BLAST_MODEL_EVALUATION_DIRECTORY, "type_classification_report.txt")
# VARIANTS_BLAST_RESULTS_FILE      = os.path.join(PREDICTION_RESULTS_DIRECTORY, "results", "{}", "blast_search.out")
# VARIANTS_PARSED_RESULTS_FILE     = os.path.join(PREDICTION_RESULTS_DIRECTORY, "results", "{}", "blast_search_parsed.out")
