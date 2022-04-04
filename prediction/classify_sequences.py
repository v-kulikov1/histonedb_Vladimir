import subprocess, json, logging, argparse

from prediction.histonedb_classifier import HistonedbVariantClassifier
from path_variables import *
'''
This executable script will:
    + download nr or other db
    - classify sequences using classification_model
'''
LOG_DIRECTORY = os.path.join("log")
logging.basicConfig(filename=os.path.join(LOG_DIRECTORY, "classification.log"),
                        format='%(asctime)s %(name)s %(levelname)-8s %(message)s',
                        level=logging.INFO,
                        datefmt='%Y-%m-%d %H:%M:%S')
log = logging.getLogger(__name__)

parser = argparse.ArgumentParser(description='Process some classification parameters')
parser.add_argument('--db', dest='db', default='nr',
                    help='defines the database or href for download')
parser.add_argument('--procs', dest='procs', default='20',
                    help='defines the number of procs')

args = parser.parse_args()

def get_db(db):
    if db == "nr" and not os.path.isfile('nr'): download_nr()
    if db == "swissprot" and not os.path.isfile('swissprot'): download_swissprot()
    if ('http://' in db) or ('https://' in db) or ('ftp://' in db):
        log.info(f'Provided db file is a link {db} - Expecting a gzipped file. Attempting to download ...')
        subprocess.call(["wget", db, '-O', 'db.gz'])
        subprocess.call(["gunzip", "db.gz"])
        db = 'db'
    return db

def download_nr():
    """Download nr if not present"""
    if not os.path.isfile('nr'):
        log.info("Downloading nr...")
        # print >> .stdout, "Downloading nr..."
        with open("nr.gz", "w") as nrgz:
            subprocess.call(["curl", "-#", "ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz"], stdout=nrgz)
        subprocess.call(["gunzip", "nr.gz"])

def download_swissprot():
    """Download nr if not present"""
    if not os.path.isfile('swissprot'):
        log.info("Downloading swissprot...")
        # print >> .stdout, "Downloading swissprot..."
        with open("swissprot.gz", "w") as swissprotgz:
            subprocess.call(["curl", "-#", "ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/swissprot.gz"],
                            stdout=swissprotgz)
        subprocess.call(["gunzip", "swissprot.gz"])

def db_split(sequences, procs):
    log.info(f"Splitting database file into {procs} parts")
    # !!!!!!!!!!!!!!!
    os.system('rm -rf db_split')
    os.system('mkdir db_split')
    # This is tricky tricky to make it fast
    size = os.path.getsize(sequences)
    split_size = int(size / procs) + 1
    os.system(f'split --bytes={split_size} --numeric-suffixes=1 {sequences} db_split/split')
    # We need to heal broken fasta records
    for i in range(1, procs + 1):
        for k in range(1, 10000):
            outsp = subprocess.check_output(['head', '-n', str(k), f'db_split/split{i:02d}'])
            if (outsp.split(b'\n')[-2].decode("utf-8")[0] == '>'):
                print(k)
                break
        if (k > 1):
            os.system(f'head -n {k} db_split/split{i:02d} >>db_split/split{i-1:02d} ')
            os.system(f'tail -n +{k} db_split/split{i:02d}> db_split/temp')
            os.system(f'mv db_split/temp db_split/split{i:02d}')

def main():
    db_file = get_db(args.db)
    print(args.db)
    print(args.procs)

    db_split(args.db, args.procs)

    for i in range(args.procs):
        with open(VARIANTS_JSON) as f:
            variant_json = json.loads(f.read())
            classification_tree = variant_json['tree']
        # classifier = HistonedbVariantClassifier(classification_tree=classification_tree)
        # classifier.create_hmms(seed_directory=SEED_CORE_DIRECTORY)
        # classifier.create_blastdbs(seed_directory=SEED_DIRECTORY)
        classification_tree.pop('Archaeal')
        classifier = HistonedbVariantClassifier(classification_tree=classification_tree)
        classifier.create_hmms(seed_directory=os.path.join(DATA_DIRECTORY, "draft_seeds"))
        classifier.create_blastdbs(seed_directory=os.path.join(DATA_DIRECTORY, "draft_seeds"))
        res = classifier.predict(sequences=os.path.join(PREDICTION_DIRECTORY, f'db_split/split{i:02d}'))
        classifier.dump_results(file_name=DUMPS_PICKLE.format(i))


if __name__ == "__main__":
    main()
