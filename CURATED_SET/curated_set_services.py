from Bio import Entrez, SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Align.AlignInfo import SummaryInfo

from pynucl.hist_features import hist_shf4seq
from pytexshade import ipyshade

import pandas as pd
import numpy as np
import collections
import sys, os, re, subprocess, io
from datetime import datetime

HISTONES_FILE = 'histones.csv'
HISTONES_PROCESSED_FILE = 'histones_processed.csv'
FEATURES_FILE = 'features.json'
BACKUP_DIR = 'backups'

NONCBI_IDENTIFICATOR = 'NONCBI'

#Auxillary functions

from IPython.display import display

def show_msa_jl(msa):
    """
    This requires jupyterlab-fasta extenstion and works only in jupyterlab
    """
    data=msa.format('fasta')
    bundle = {}
    bundle['application/vnd.fasta.fasta'] = data
    bundle['text/plain'] = data
    display(bundle, raw=True)

# def show_msa_with_features(self, hist_type=None, variant=None, subvariant=None, taxonomyid=None, organism=None, phylum=None, taxonomy_class=None,
#                              shading_modes=['charge_functional'], legend=False, title='', logo=False, hideseqs=False, splitN=20, setends=[],
#                              ruler=False, show_seq_names=True, funcgroups=None, show_seq_length=False, debug=False):
#         df = self.data
#         if hist_type: df = df.loc[self.data['type'] == hist_type]
#         if variant: df = df.loc[self.data['variant'] == variant]
#         if subvariant: df = df.loc[self.data['subvariant'] == subvariant]
#         if taxonomyid: df = df.loc[self.data['taxonomyid'] == taxonomyid]
#         if organism: df = df.loc[self.data['organism'] == organism]
#         if phylum: df = df.loc[self.data['phylum'] == phylum]
#         if taxonomy_class: df = df.loc[self.data['class'] == taxonomy_class]
#         if df.shape[0] < 1: raise AssertionError(f'No such sequences')
#         return self.show_aln(list(df['accession']), shading_modes=shading_modes, legend=legend, title=title, logo=logo, hideseqs=hideseqs, splitN=splitN,
#                              setends=setends, ruler=ruler, show_seq_names=show_seq_names, funcgroups=funcgroups, show_seq_length=show_seq_length, debug=debug)

# def show_aln(self, accessions,
#                  shading_modes=['charge_functional'], legend=False, title='', logo=False, hideseqs=False, splitN=20, setends=[],
#                  ruler=False, show_seq_names=True, funcgroups=None, show_seq_length=False, debug=True):
#         msa = self.muscle_aln(accessions=accessions)
#         a = SummaryInfo(msa)
#         cons = a.dumb_consensus(threshold=0.1, ambiguous='X')
#         features = hist_shf4seq(cons)
#         return ipyshade.shadedmsa(msa,
#                                   shading_modes=shading_modes,
#                                   legend=legend,
#                                   features=features,
#                                   title=title,
#                                   logo=logo,
#                                   hideseqs=hideseqs,
#                                   splitN=splitN,
#                                   setends=setends,
#                                   ruler=ruler,
#                                   show_seq_names=show_seq_names,
#                                   funcgroups=funcgroups,
#                                   show_seq_length=show_seq_length,
#                                   debug=debug
#                                  )


#Main class    

class CuratedSet(object):
    def __init__(self):
        self.data = self.read_data(HISTONES_FILE)
        self.fasta_seqrec = {} # keys - accession, values SeqRec Object
        self.msa_variant = {} # keys - variant, values MultipleSeqAlignment Object
        self.msa_type = {} # keys - type, values MultipleSeqAlignment Object
        self.create_fasta_seqrec(self.data[self.data['sequence']!=''])

    def read_data(self, filename):
        df = pd.read_csv(filename, sep=',|;', engine='python').fillna('')
        df['taxonomyid'] = df['taxonomyid'].astype(str)
        df.index = list(df['accession'])
        return df

    def has_duplicates(self):
        '''
        This method allows us to check that our curated set does not have identifiers that are assigned different variants.
        :return: None (if no such identifiers) or iterator of accessions (that are assigned different variants)
        '''
        acc_list = [acc.split('.')[0] for acc in self.get_accessions_list()] #retrieving accessions without version
        if len(set(acc_list)) == self.get_count(): return
        c = collections.Counter()
        for acc in acc_list:
            c[acc] += 1
        for acc_count in c.most_common():
            if acc_count[1] == 1: continue
            else: yield acc_count[0]

    def create_fasta_seqrec(self, data):
        for i, row in data.iterrows():
            self.fasta_seqrec[row['accession']] = SeqRecord(Seq(row['sequence']), id=row['variant']+'_'+row['organism'].replace(' ','_')+'_'+row['accession'],description=f"{row['accession']} histone: {row['type']} variant: {row['variant']} organism: {row['organism']}")

    def get_count(self): return self.data.shape[0]

    def get_accessions_list(self): return list(self.data['accession'])

    def get_gis_list(self): return list(self.data['gi'])

    def get_variant(self, accession): return self.data.loc[self.data['accession'] == accession]['variant'].iloc[0]

    def get_type(self, accession): return self.data.loc[self.data['accession'] == accession]['type'].iloc[0]

    def get_gi(self, accession): return self.data.loc[self.data['accession'] == accession]['gi'].iloc[0]

    def get_taxid_genus(self, accession):
        if not 'taxonomyid' in self.data:
            print('No taxonomyid loaded yet. Update update_taxids fisrt.')
            return
        return (self.data.loc[self.data['accession'] == accession]['taxonomyid'].iloc[0],
                self.data.loc[self.data['accession'] == accession]['organism'].iloc[0])

    def get_ncbi_data(self): return self.data[self.data['accession'].str.startswith(NONCBI_IDENTIFICATOR, na=False) == False]

    def get_noncbi_data(self): return self.data[self.data['accession'].str.startswith(NONCBI_IDENTIFICATOR, na=False)]

    def get_blank_data(self, columns):
        data = self.get_ncbi_data()
        return data.iloc[list(np.where(data[columns]=='')[0])]

    def save(self, filename=None, processed=False):
        # compare current data with data loaded from file
        data_old = self.read_data(HISTONES_FILE)
        for i, row in data_old.compare(self.data).iterrows():
            print(row)

        # backup file first to history
        backup_file = os.path.join(BACKUP_DIR, f'{HISTONES_FILE}-{datetime.now().strftime("%b%d%y%H%M%S")}')
        # data_old.to_csv(backup_file, mode='w', index=False)
        print(f'cp {HISTONES_FILE} {backup_file}')
        os.system(f'cp {HISTONES_FILE} {backup_file}')
        print(f'Previous data backuped to {backup_file}')

        # saving
        if not filename: filename = HISTONES_FILE
        if processed: filename = HISTONES_PROCESSED_FILE
        self.data.to_csv(filename, mode='w', index=False)
        print(f'Results saved to {filename}')

    def update_accession_version(self):
        updating_data = self.get_ncbi_data()
        for i in range(10):
            try:
                handle = Entrez.efetch(db="protein", id=",".join([i.split('.')[0] for i in list(updating_data['accession'])]), rettype="gb", retmode="text")
                sequences = list(SeqIO.parse(handle, "gb"))
                if (len(updating_data['accession']) == len(sequences)):
                    new_accessions = [s.id for s in sequences]
                    for new_acc, acc in zip(new_accessions, updating_data['accession']):
                        if acc!=new_acc: print(f'{acc} changes to {new_acc}')
                    updating_data['accession'] = new_accessions
                    self.data = updating_data.append(self.get_noncbi_data())
                    break
                else:
                    print("Mismatch:", len(updating_data['accession']), " ", len(sequences))
            except:
                print("Unexpected error: {}, Retrying, attempt {}".format(sys.exc_info()[0], i))
                if i == 9:
                    print(
                        "FATAL ERROR could not get seqs from NCBI after 10 attempts for %s. Will return empty list!" % (
                            ",".join(list(updating_data))))
                else:
                    continue

    def update_taxids(self, blank_data=False, accessions=None, exclude=None):
        if blank_data and (accessions or exclude): raise AssertionError(f'blank_data cannot be true when accessions list is given. Please, use only one of the options.')
        updating_data = self.get_ncbi_data()
        if blank_data: updating_data = self.get_blank_data(['taxonomyid', 'organism', 'taxonomy_group'])
        if accessions: updating_data = self.data.loc[accessions]
        if exclude: updating_data = updating_data.drop(exclude)

        sequences = []
        # if len(self.data['accession'].values) == 0:
        #     sequences = []
        # Bug in eutils
        # E.g. 6A5L_FF cannot be retrieved, 6A5L_f cannot
        # This is likely a bad fix!!! but nothing else  can be done at the moment
        # Here is an addhock fix.
        accessions = []
        p1 = re.compile("(\w{4})_([a-zA-Z]{2})")
        p2 = re.compile("(\w{4})_([a-z]{1})")

        for ac in updating_data['accession'].values:
            m1 = p1.match(ac)
            m2 = p2.match(ac)

            if m1: accessions.append(m1.group(1) + '_' + m1.group(2)[0].upper())
            elif m2: accessions.append(m2.group(1) + '_' + m2.group(2)[0].upper())
            else: accessions.append(ac)

        for i in range(10):
            try:
                handle = Entrez.efetch(db="protein", id=",".join(accessions), rettype="gb", retmode="text")
                sequences = list(SeqIO.parse(handle, "gb"))
                if (len(accessions) == len(sequences)): break
            except:
                print("Unexpected error: {}, Retrying, attempt {}".format(sys.exc_info()[0], i))
                if i == 9:
                    print("FATAL ERROR could not get seqs from NCBI after 10 attempts for %s. Will return empty list!" % (",".join(accessions)))
                else: continue

        taxids, genus, phylum, taxonomy_class = [], [], [], []
        for s in sequences:
            genus.append(s.annotations["organism"])
            try:
                for a in s.features[0].qualifiers['db_xref']:
                    text = re.search('(\S+):(\S+)', a).group(1)
                    id = re.search('(\S+):(\S+)', a).group(2)
                    if (text == "taxon"):
                        print("Fetched taxid from NCBI {}".format(id))
                        taxids.append(id)
                    else: continue
            except:
                print("!!!!!!Unable to get TAXID for \n {} setting it to 1".format(s))
                taxids.append(1)  # unable to identify
            handle = Entrez.efetch(id=taxids[-1], db="taxonomy", retmode="xml")
            tax_data = Entrez.read(handle)
            lineage = {d['Rank']: d['ScientificName'] for d in
                       tax_data[0]['LineageEx'] if d['Rank'] in ['class', 'phylum']}
            phylum.append(lineage.get('phylum', 'None'))
            taxonomy_class.append(lineage.get('class', 'None'))

        for t_new, t, g_new, g, ph_new, ph, tc_new, tc in zip(taxids, updating_data['taxonomyid'],
                                      genus, updating_data['organism'],
                                      phylum, updating_data['phylum'],
                                      taxonomy_class, updating_data['class']):
            if t != t_new: print(f'{t} changes to {t_new}')
            if g != g_new: print(f'{g} changes to {g_new}')
            if ph != ph_new: print(f'{ph} changes to {ph_new}')
            if tc != tc_new: print(f'{tc} changes to {tc_new}')
        updating_data['taxonomyid'] = taxids
        updating_data['organism'] = genus
        updating_data['phylum'] = phylum
        updating_data['class'] = taxonomy_class
        self.data = updating_data.append(self.data.drop(list(updating_data['accession'])))

    def update_sequence(self, noncbi=False, blank_data=False):
        """
        Download a dictionary of fasta SeqsRec from NCBI given a list of ACCESSIONs.
        :return: dict of fasta SeqRecords
        """

        print("Downloading FASTA SeqRecords by ACCESSIONs from NCBI")
        updating_data = self.get_ncbi_data()
        if blank_data: updating_data = self.get_blank_data(['sequence'])

        num = updating_data.shape[0]
        for i in range(10):
            fasta_seqrec_with_acc = dict()
            try:
                print("Fetching %d seqs" % (num))
                strn = ",".join(list(updating_data['accession']))
                request = Entrez.epost(db="protein", id=strn)
                result = Entrez.read(request)
                webEnv = result["WebEnv"]
                queryKey = result["QueryKey"]
                handle = Entrez.efetch(db="protein", rettype='fasta', retmode='text', webenv=webEnv, query_key=queryKey)
                for r in SeqIO.parse(handle, 'fasta'):
                    id = r.id
                    if '|' in r.id: id = r.id.split('|')[1]
                    if id not in updating_data['accession']:
                        raise AssertionError(f'{id} is not in accessions list')
                    fasta_seqrec_with_acc[id] = r
                    if self.data.loc[id, 'sequence'] != str(r.seq):
                        print(f'Sequence for {id} changes from {self.data.loc[id, "sequence"]} to {str(r.seq)}')
                        self.data.at[id, 'sequence'] = str(r.seq)
            except AssertionError: raise
            except Exception as e: continue
            if (len(fasta_seqrec_with_acc) == num):
                break
            else:
                print("Mismatch:", num, " ", len(fasta_seqrec_with_acc))
        if (len(fasta_seqrec_with_acc) != num):
            print('FASTA Records could not be fetched')
            return
        print(f"FASTA Records downloaded: {len(fasta_seqrec_with_acc)}")
        print(f"Added SeqRecords for {fasta_seqrec_with_acc.keys()}")

        self.fasta_seqrec.update(fasta_seqrec_with_acc)

        if noncbi:
            print("Creating FASTA Records for NONCBI sequences...")
            self.create_fasta_seqrec(self.get_noncbi_data())

        return self.fasta_seqrec # maybe not need

    def muscle_aln(self, accessions, options=[],debug = False):
        """
        Align with muscle all sequences from defined accessions
        accessions: list of accessions
        :return: MultipleSeqAlignment object
        """
        # if not self.fasta_seqrec: self.update_sequence()
        if not set(accessions).issubset(self.fasta_seqrec.keys()): self.update_sequence()

        muscle = os.path.join(os.path.dirname(sys.executable), "muscle")
        process = subprocess.Popen([muscle]+options, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        sequences = "\n".join([s.format("fasta") for key, s in self.fasta_seqrec.items() if key in accessions])

        aln, error = process.communicate(sequences.encode('utf-8'))
        if debug:
            print(sequences)
            print()
            print("Stderr:")
            print(error.decode('utf-8')) 
            print("Stdout:")
            print(aln.decode('utf-8')) 

        seqFile = io.StringIO()
        seqFile.write(aln.decode('utf-8'))
        seqFile.seek(0)
        sequences = list(SeqIO.parse(seqFile, "fasta"))  # Not in same order, but does it matter?
        msa = MultipleSeqAlignment(sequences)
        return msa



    def muscle_aln_for_variant(self):
        """
        Align with muscle sequences grouped by histone variant
        :return: dict MultipleSeqAlignment for variants object
        """
        self.msa_variant = {variant: self.muscle_aln(self.data[self.data['variant'] == variant]['accession']) for variant in set(self.data['variant'])}
        return self.msa_variant

    def muscle_aln_for_type(self):
        """
        Align with muscle sequences grouped by histone types
        :return: dict MultipleSeqAlignment for types object
        """
        self.msa_type = {hist_type: self.muscle_aln(self.data[self.data['type'] == hist_type]['accession']) for hist_type in set(self.data['type'])}
        return self.msa_type

    def generate_seeds(self, directory='draft_seeds'):
        if not self.msa_variant: self.muscle_aln_for_variant()

        if not os.path.exists(directory):
            os.makedirs(directory)
        for hist_type in ['H2A', 'H2B', 'H3', 'H4', 'H1']:
            for hist_var in set(self.data[self.data['type'] == hist_type]['variant']):
                print("##########Starting", hist_var, hist_type, f'{hist_var}.fasta')
                if not os.path.exists(os.path.join(directory, hist_type)):
                    os.makedirs(os.path.join(directory, hist_type))
                print(self.msa_variant[hist_var])
                # generate new msa with refactored title
                msa_r = MultipleSeqAlignment([])
                for i in self.msa_variant[hist_var]:
                    # print(f'-------------------{i.id}')
                    accession = i.id
                    if '|' in accession: accession = accession.split('|')[1]
                    try:
                        genus = re.search(r"\[(\S+)\s+.+\S+\]", i.description).group(1)
                    except:
                        genus = self.get_taxid_genus(accession)[1]

                    i.id = genus + "|" + accession + "|" + hist_var
                    i.description = genus + "_" + hist_var + "_" + accession
                    msa_r.append(i)
                msa_r.sort()
                print(msa_r)
                AlignIO.write(msa_r, os.path.join("draft_seeds", hist_type, f'{hist_var}.fasta'), 'fasta')

            # combines MSA
            if not self.msa_type: self.muscle_aln_for_type()
            # generate new msa with refactored title
            msa_r = MultipleSeqAlignment([])
            for i in self.msa_type[hist_type]:
                print(i.description)
                accession = i.id
                if '|' in accession: accession = accession.split('|')[1]
                try:
                    genus=re.search(r"\[(\S+)\s+.+\S+\]",i.description).group(1)
                except:
                    genus = self.get_taxid_genus(accession)[1]
                # text = re.search(r"(\S+)\|(\d+)\|(\S+)", i.id)
                print(f'-------------------{i.id}')
                i.id = genus + "|" + accession + "|" + hist_type
                # i.description=genus+"_"+variant+"_"+gi
                msa_r.append(i)
            msa_r.sort()
            AlignIO.write(msa_r, os.path.join("draft_seeds", f'{hist_type}.fasta'), 'fasta')
