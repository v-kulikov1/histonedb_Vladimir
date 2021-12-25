from Bio import Entrez, SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

import pandas as pd
import collections
import sys, os, re, subprocess, io
from datetime import datetime

HISTONES_FILE = 'histones.csv'
HISTONES_PROCESSED_FILE = 'histones_processed.csv'
FEATURES_FILE = 'features.json'
BACKUP_DIR = 'backups'

NONCBI_IDENTIFICATOR = 'NONCBI'


class CuratedSet(object):
    def __init__(self):
        self.data = pd.read_csv(HISTONES_FILE).fillna('') # needs sprting
        self.data.index = list(self.data['accession'])
        self.fasta_seqrec = {} # keys - accession, values SeqRec Object
        self.msa_variant = {} # keys - variant, values MultipleSeqAlignment Object
        self.msa_type = {} # keys - type, values MultipleSeqAlignment Object
        self.features_variant = {} # keys - variant, values Feature Object
        self.features_type = {} # keys - type, values Feature Object

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

    def get_count(self): return self.data.shape[0]

    def get_accessions_list(self): return list(self.data['accession'])

    def get_gis_list(self): return list(self.data['gi'])

    def get_variant(self, accession): return self.data.loc[self.data['accession'] == accession]['variant'].iloc[0]

    def get_type(self, accession): return self.data.loc[self.data['accession'] == accession]['type'].iloc[0]

    def get_gi(self, accession): return self.data.loc[self.data['accession'] == accession]['gi'].iloc[0]

    def get_taxid_genus(self, accession):
        if not 'taxonomyid' in self.data:
            print('No taxid loaded yet. Update update_taxids fisrt.')
            return
        return (self.data.loc[self.data['accession'] == accession]['taxonomyid'].iloc[0],
                self.data.loc[self.data['accession'] == accession]['organism'].iloc[0])

    def save(self, filename=None, processed=False):
        #backup file first to history
        data_old = pd.read_csv(HISTONES_FILE).fillna('')
        backup_file = os.path.join(BACKUP_DIR, f'{HISTONES_FILE}-{datetime.now().strftime("%b%d%y%H%M%S")}')
        data_old.to_csv(backup_file, mode='w', index=False)
        print(f'Previous data backuped to {backup_file}')
        if not filename: filename = HISTONES_FILE
        if processed: filename = HISTONES_PROCESSED_FILE
        self.data.to_csv(filename, mode='w', index=False)
        print(f'Results saved to {filename}')

    def add_curated_data(self, filemane, save=False):
        add_curated_data = pd.read_csv(filemane).fillna('')
        self.data = self.data.append(add_curated_data)
        f = self.has_duplicates()
        if f: print(f'Duplicates: {list(f)}')
        elif save: self.save()

    def update_accession_version(self):
        curated_data_with_acc = self.data[self.data['accession'].str.startswith(NONCBI_IDENTIFICATOR, na=False) == False]
        curated_data_noncbi = self.data[self.data['accession'].str.startswith(NONCBI_IDENTIFICATOR, na=False)]
        for i in range(10):
            try:
                handle = Entrez.efetch(db="protein", id=",".join(list(curated_data_with_acc['accession'])), rettype="gb", retmode="text")
                sequences = list(SeqIO.parse(handle, "gb"))
                if (len(curated_data_with_acc['accession']) == len(sequences)):
                    new_accessions = [s.id for s in sequences]
                    for new_acc, acc in zip(new_accessions, curated_data_with_acc['accession']):
                        if acc!=new_acc: print(f'{acc} changes to {new_acc}')
                    curated_data_with_acc['accession'] = new_accessions
                    self.data = curated_data_with_acc.append(curated_data_noncbi)
                    break
                else:
                    print("Mismatch:", len(curated_data_with_acc['accession']), " ", len(sequences))
            except:
                print("Unexpected error: {}, Retrying, attempt {}".format(sys.exc_info()[0], i))
                if i == 9:
                    print(
                        "FATAL ERROR could not get seqs from NCBI after 10 attempts for %s. Will return empty list!" % (
                            ",".join(list(curated_data_with_acc))))
                else:
                    continue

    def update_taxonomy_group(self):
        curated_data_with_acc = self.data[self.data['accession'].str.startswith(NONCBI_IDENTIFICATOR, na=False) == False]
        curated_data_noncbi = self.data[self.data['accession'].str.startswith(NONCBI_IDENTIFICATOR, na=False)]

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

        for ac in curated_data_with_acc['accession'].values:
            m1 = p1.match(ac)
            m2 = p2.match(ac)

            if m1: accessions.append(m1.group(1) + '_' + m1.group(2)[0].upper())
            elif m2: accessions.append(m2.group(1) + '_' + m2.group(2)[0].upper())
            else: accessions.append(ac)
        else:
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

        taxids, genus = [], []
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

        # taxids for NONCBI set as root. We will change them manually!
        taxids.extend([1]*curated_data_noncbi.shape[0])
        genus.extend(['root']*curated_data_noncbi.shape[0])
        self.data['taxonomyid'] = taxids
        self.data['organism'] = genus

    def update_taxids(self):
        curated_data_with_acc = self.data[self.data['accession'].str.startswith(NONCBI_IDENTIFICATOR, na=False) == False]
        curated_data_noncbi = self.data[self.data['accession'].str.startswith(NONCBI_IDENTIFICATOR, na=False)]

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

        for ac in curated_data_with_acc['accession'].values:
            m1 = p1.match(ac)
            m2 = p2.match(ac)

            if m1: accessions.append(m1.group(1) + '_' + m1.group(2)[0].upper())
            elif m2: accessions.append(m2.group(1) + '_' + m2.group(2)[0].upper())
            else: accessions.append(ac)
        else:
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

        taxids, genus = [], []
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

        # taxids for NONCBI set as root. We will change them manually!
        taxids.extend([1]*curated_data_noncbi.shape[0])
        genus.extend(['root']*curated_data_noncbi.shape[0])
        self.data['taxonomyid'] = taxids
        self.data['organism'] = genus

    def update_fasta_seqrec(self):
        """
        Download a dictionary of fasta SeqsRec from NCBI given a list of ACCESSIONs.
        :return: dict of fasta SeqRecords
        """

        print("Downloading FASTA SeqRecords by ACCESSIONs from NCBI")
        curated_data_with_acc = self.data[self.data['accession'].str.startswith(NONCBI_IDENTIFICATOR, na=False) == False]
        curated_data_noncbi = self.data[self.data['accession'].str.startswith(NONCBI_IDENTIFICATOR, na=False)]

        num = curated_data_with_acc.shape[0]
        for i in range(10):
            self.fasta_seqrec = dict()
            try:
                print("Fetching %d seqs" % (num))
                strn = ",".join(list(curated_data_with_acc['accession']))
                request = Entrez.epost(db="protein", id=strn)
                result = Entrez.read(request)
                webEnv = result["WebEnv"]
                queryKey = result["QueryKey"]
                handle = Entrez.efetch(db="protein", rettype='fasta', retmode='text', webenv=webEnv, query_key=queryKey)
                for r in SeqIO.parse(handle, 'fasta'):
                    if r.id not in curated_data_with_acc['accession']:
                        raise Exception(f'{r.id} is not in accessions list')
                    if '|' in r.id:
                        self.fasta_seqrec[r.id.split('|')[1]] = r
                        self.data.at[r.id.split('|')[1], 'sequence'] = str(r.seq)
                    else:
                        self.fasta_seqrec[r.id] = r
                        self.data.at[r.id, 'sequence'] = str(r.seq)
            except Exception as e:
                continue
            if (len(self.fasta_seqrec) == num):
                break
            else:
                print("Mismatch:", num, " ", len(self.fasta_seqrec))
        if (len(self.fasta_seqrec) != num):
            print('FASTA Records could not be fetched')
            return
        print(f"FASTA Records downloaded: {len(self.fasta_seqrec)}")

        print("Creating FASTA Records for NONCBI sequences...")
        for i, row in curated_data_noncbi.iterrows():
            self.fasta_seqrec[row['accession']] = SeqRecord(Seq(row['sequence']), id=row['accession'],
                                                            description=f"{row['accession']} histone {row['type']}")

        return self.fasta_seqrec

    def muscle_aln(self, accessions):
        """
        Align with muscle all sequences from defined accessions
        accessions: list of accessions
        :return: MultipleSeqAlignment object
        """
        if not self.fasta_seqrec: self.update_fasta_seqrec()

        muscle = os.path.join(os.path.dirname(sys.executable), "muscle")
        process = subprocess.Popen([muscle], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        sequences = "\n".join([s.format("fasta") for key, s in self.fasta_seqrec.items() if key in [accessions]])
        print(sequences)
        aln, error = process.communicate(sequences.encode('utf-8'))
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
        if not self.fasta_seqrec: self.update_fasta_seqrec()
        self.msa_variant = {variant: self.muscle_aln(self.data[self.data['variant'] == variant]['accession']) for variant in set(self.data['variant'])}
        return self.msa_variant

    def muscle_aln_for_type(self):
        """
        Align with muscle sequences grouped by histone types
        :return: dict MultipleSeqAlignment for types object
        """
        if not self.fasta_seqrec: self.update_fasta_seqrec()
        self.msa_type = {hist_type: self.muscle_aln(self.data[self.data['type'] == hist_type]['accession']) for hist_type in set(self.data['type'])}
        return self.msa_type

    def generate_seeds(self, directory='draft_seeds'):
        if not self.fasta_seqrec: self.update_fasta_seqrec()
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


class Feature(object):
    def __init__(self):
        id          = None
        template    = None
        start       = None
        end         = None
        name        = None
        description = None
        color       = None
        objects     = None

    def __str__(self):
        """Returns Jalview GFF format"""
        return self.gff(str(self.template))

    def gff(self, sequence_label=None, featureType="{}"):
        tmp = ""
        if sequence_label is None:
            tmp += "color1\t{}\n".format(self.color)
            sequence_label = str(self.template)
            featureType = "color1"

        tmp += "\t".join((self.name, sequence_label, "-1", str(self.start), str(self.end), featureType))
        tmp += "\n"
        return tmp

    def get_variant_features(self, sequence, variants=None, save_dir="", save_not_found=False, save_gff=True,
                             only_general=False):
        """Get the features of a sequence based on its variant.

        Parameters:
        -----------
        sequence: Sequence django model
            The seuqence to add get features for with identified variant
        variants: List of Variant models
            Anntate others variants. Optional.
        save_dir: str
            Path to save temp files.
        save_not_found: bool
            Add Features even if they weren't found. Indices will be (-1, -1)

        Return:
        -------
        A string containing the gff file of all features
        """
        # Save query fasta to a file to EMBOSS needle can read it
        n2 = str(uuid.uuid4())
        test_record = sequence.to_biopython()
        query_file = os.path.join(save_dir, "query_{}.fa".format(n2))
        SeqIO.write(test_record, query_file, 'fasta')

        # A list of updated Features for the query
        variant_features = set()

        if not variants:
            variants = [sequence.variant]

        for variant in variants:
            templates = [variant.id, "General{}".format(variant.hist_type.id)] if not only_general else [
                "General{}".format(variant.hist_type.id)]
            for template_variant in templates:
                try:
                    features = Feature.objects.filter(template__variant=template_variant)
                except:
                    continue
                # Find features with the same taxonomy
                tax_features = features.filter(template__taxonomy=sequence.taxonomy)
                if len(tax_features) == 0:
                    # Find features with closest taxonomy => rank class
                    tax_features = features.filter(
                        template__taxonomy__parent__parent__parent=sequence.taxonomy.parent.parent.parent)
                if len(tax_features) == 0:
                    # Nothing, use unidentified which is the standard
                    tax_features = features.filter(template__taxonomy__name="undefined")
                features = tax_features
                for updated_feature in transfer_features_from_template_to_query(features, query_file, save_dir=save_dir,
                                                                                save_not_found=save_not_found):
                    variant_features.add(updated_feature)

        os.remove(query_file)
        if save_gff:
            return Feature.objects.gff(sequence.id, variant_features)
        return variant_features



# curated_set = CuratedSet()
# curated_set.rename('macroH2A', '')
# # print(list(curated_set.has_duplicates()))
# curated_set.update_accession_version()
# curated_set.update_fasta_seqrec()
# curated_set.muscle_aln_for_variant()
# print(curated_set.msa_variant)
# curated_set.muscle_aln_for_type()
# print(curated_set.msa_type)
# curated_set.generate_seeds()
# print(curated_set.get_taxid_genus('XP_846259.1'))
# for i in curated_set.get_taxid_genus('XP_846259.1'):
#     print(i)