import os, subprocess
import logging, argparse

# from tools.browse_service import *
'''
This executable script will:
    - download nr or other db
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

def main():
    db_file = get_db(args.db)
    print(args.db)
    print(args.procs)


if __name__ == "__main__":
    main()
