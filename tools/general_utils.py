from prediction_app.prediction.path_variables import *

# Seeds
def get_seeds(seeds_name='seeds',combined_alignments=False, generic=False):
    """
    Goes through static/browse/seeds directories and returns histone type names and fasta file name of variant (without path).
    If combined_alignments returns histone types as seed and '' as hist_type, additionally to variants
    """
    seed_dir = SEED_DIRECTORY
    if seeds_name == 'seeds_core':
        seed_dir = SEED_CORE_DIRECTORY
    for i, (root, _, files) in enumerate(os.walk(seed_dir)):
        hist_type = os.path.basename(root)
        if hist_type == seeds_name and not combined_alignments:  # means we are in top dir, we skip,
            # combinded alignmnents for hist types are their, but we do not use them in database constuction,
            # only in visualization on website
            continue
        elif hist_type == seeds_name:
            hist_type = ""
        for seed in files:
            if not seed.endswith(".fasta") or (not generic and 'generic' in seed): continue
            yield hist_type, seed