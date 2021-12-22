from Bio import Entrez, SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

import pandas as pd
import collections
import sys, os, re, subprocess, io

class CuratedSet(object):
    def __init__(self):
        self.curated_filename = 'curated_histones.csv'
        self.curated_data = pd.read_csv(self.curated_filename).fillna('')
        self.curated_data.index = list(self.curated_data['accession'])
        self.fasta_seqrec = {} # keys - accession, values SeqRec Object
        self.msa_variant = {} # keys - variant, values MultipleSeqAlignment Object
        self.msa_type = {} # keys - type, values MultipleSeqAlignment Object
        self.msa_full = None # MultipleSeqAlignment object

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

    def get_count(self): return self.curated_data.shape[0]

    def get_accessions_list(self): return list(self.curated_data['accession'])

    def get_gis_list(self): return list(self.curated_data['gi'])

    def get_variant(self, accession):
        if not hasattr(accession, '__iter__'):
            return self.curated_data.loc[self.curated_data['accession']==accession]['variant'].iloc[0]
        for acc in accession:
            yield self.curated_data.loc[self.curated_data['accession']==acc]['variant'].iloc[0]

    def get_type(self, accession):
        if not hasattr(accession, '__iter__'):
            return self.curated_data.loc[self.curated_data['accession']==accession]['type'].iloc[0]
        for acc in accession:
            yield self.curated_data.loc[self.curated_data['accession']==acc]['type'].iloc[0]

    def get_gi(self, accession):
        if not hasattr(accession, '__iter__'):
            return self.curated_data.loc[self.curated_data['accession']==accession]['gi'].iloc[0]
        for acc in accession:
            yield self.curated_data.loc[self.curated_data['accession']==acc]['gi'].iloc[0]

    def get_taxid_genus(self, accession):
        if not 'taxid' in self.curated_data:
            print('No taxid loaded yet. Update update_taxids fisrt.')
            return
        return (self.curated_data.loc[self.curated_data['accession']==accession]['taxid'].iloc[0],
                self.curated_data.loc[self.curated_data['accession'] == accession]['organism'].iloc[0])

    def save(self):
        self.curated_data.to_csv(self.curated_filename, mode='w', index=False)
        print(f'Results saved to {self.curated_filename}')

    def add_curated_data(self, filemane, save=False):
        add_curated_data = pd.read_csv(filemane).fillna('')
        self.curated_data = self.curated_data.append(add_curated_data)
        f = self.has_duplicates()
        if f: print(f'Duplicates: {list(f)}')
        elif save: self.save()

    def update_accession_version(self, save=False):
        curated_data_with_acc = self.curated_data[self.curated_data['accession'].str.startswith('NOGI', na=False) == False]
        curated_data_nogi = self.curated_data[self.curated_data['accession'].str.startswith('NOGI', na=False)]
        for i in range(10):
            try:
                handle = Entrez.efetch(db="protein", id=",".join(list(curated_data_with_acc['accession'])), rettype="gb", retmode="text")
                sequences = list(SeqIO.parse(handle, "gb"))
                if (len(curated_data_with_acc['accession']) == len(sequences)):
                    new_accessions = [s.id for s in sequences]
                    new_accessions.extend(list(curated_data_nogi['accession']))
                    self.curated_data['accession'] = new_accessions
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
        if save: self.save()

    def update_taxonomy_group(self, save=False):
        curated_data_with_acc = self.curated_data[self.curated_data['accession'].str.startswith('NOGI', na=False)==False]
        curated_data_nogi = self.curated_data[self.curated_data['accession'].str.startswith('NOGI', na=False)]

        sequences = []
        # if len(self.curated_data['accession'].values) == 0:
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

        # taxids for NOGI set as root. We will change them manually!
        taxids.extend([1]*curated_data_nogi.shape[0])
        genus.extend(['root']*curated_data_nogi.shape[0])
        self.curated_data['taxid'] = taxids
        self.curated_data['organism'] = genus
        if save: self.save()

    def update_taxids(self, save=False):
        curated_data_with_acc = self.curated_data[self.curated_data['accession'].str.startswith('NOGI', na=False)==False]
        curated_data_nogi = self.curated_data[self.curated_data['accession'].str.startswith('NOGI', na=False)]

        sequences = []
        # if len(self.curated_data['accession'].values) == 0:
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

        # taxids for NOGI set as root. We will change them manually!
        taxids.extend([1]*curated_data_nogi.shape[0])
        genus.extend(['root']*curated_data_nogi.shape[0])
        self.curated_data['taxid'] = taxids
        self.curated_data['organism'] = genus
        if save: self.save()

    def update_fasta_seqrec(self):
        """
        Download a dictionary of fasta SeqsRec from NCBI given a list of ACCESSIONs.
        :return: dict of fasta SeqRecords
        """

        print("Downloading FASTA SeqRecords by ACCESSIONs from NCBI")
        curated_data_with_acc = self.curated_data[self.curated_data['accession'].str.startswith('NOGI', na=False)==False]
        curated_data_nogi = self.curated_data[self.curated_data['accession'].str.startswith('NOGI', na=False)]

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
                    if '|' in r.id:
                        self.fasta_seqrec[r.id.split('|')[1]] = r
                        self.curated_data.at[r.id.split('|')[1],'sequence'] = str(r.seq)
                    else:
                        self.fasta_seqrec[r.id] = r
                        self.curated_data.at[r.id,'sequence'] = str(r.seq)
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

        print("Creating FASTA Records for NOGI sequences...")
        for i, row in curated_data_nogi.iterrows():
            self.fasta_seqrec[row['accession']] = SeqRecord(Seq(row['sequence']), id=row['accession'],
                                                            description=f"{row['accession']} histone {row['type']}")

        return self.fasta_seqrec

    def muscle_aln(self):
        """
        Align with muscle all sequences if only you have updated it
        :return: MultipleSeqAlignment object if only you have updated it else None
        """
        if not self.fasta_seqrec:
            print('Empty data. Update fasta_seqrec fisrt.')
            return

        muscle = os.path.join(os.path.dirname(sys.executable), "muscle")
        process = subprocess.Popen([muscle], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        sequences = "\n".join([s.format("fasta") for key, s in self.fasta_seqrec.items()])
        print(sequences)
        aln, error = process.communicate(sequences)
        seqFile = io.BytesIO()
        seqFile.write(aln)
        seqFile.seek(0)
        sequences = list(SeqIO.parse(seqFile, "fasta"))  # Not in same order, but does it matter?
        self.msa_full = MultipleSeqAlignment(sequences)
        return self.msa_full

    def muscle_aln_for_variant(self):
        """
        Align with muscle sequences grouped by histone variant if only you have updated it
        :return: dict MultipleSeqAlignment for variants object if only you have updated it else None
        """
        self.msa_variant = {}
        if not self.fasta_seqrec:
            print('Empty data. Update fasta_seqrec fisrt.')
            return
        for variant in set(self.curated_data['variant']):
            variant_fasta_seqrec = {acc: self.fasta_seqrec[acc] for acc in self.curated_data[self.curated_data['variant']==variant]['accession']}
            muscle = os.path.join(os.path.dirname(sys.executable), "muscle")
            process = subprocess.Popen([muscle], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
            sequences = "\n".join([s.format("fasta") for key, s in variant_fasta_seqrec.items()])
            print(sequences)
            aln, error = process.communicate(sequences.encode('utf-8'))
            seqFile = io.StringIO()
            seqFile.write(aln.decode('utf-8'))
            seqFile.seek(0)
            sequences = list(SeqIO.parse(seqFile, "fasta"))
            self.msa_variant[variant] = MultipleSeqAlignment(sequences)
        return self.msa_variant

    def muscle_aln_for_type(self):
        """
        Align with muscle sequences grouped by histone types if only you have updated it
        :return: dict MultipleSeqAlignment for types object if only you have updated it else None
        """
        self.msa_type = {}
        if not self.fasta_seqrec:
            print('Empty data. Update fasta_seqrec fisrt.')
            return
        for type in set(self.curated_data['type']):
            variant_fasta_seqrec = {acc: self.fasta_seqrec[acc] for acc in self.curated_data[self.curated_data['type']==type]['accession']}
            muscle = os.path.join(os.path.dirname(sys.executable), "muscle")
            process = subprocess.Popen([muscle], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
            sequences = "\n".join([s.format("fasta") for key, s in variant_fasta_seqrec.items()])
            print(sequences)
            aln, error = process.communicate(sequences.encode('utf-8'))
            seqFile = io.StringIO()
            seqFile.write(aln.decode('utf-8'))
            seqFile.seek(0)
            sequences = list(SeqIO.parse(seqFile, "fasta"))  # Not in same order, but does it matter?
            self.msa_type[type] = MultipleSeqAlignment(sequences)
        return self.msa_type

    def generate_seeds(self, directory='draft_seeds'):
        if not self.fasta_seqrec: self.update_fasta_seqrec()
        if not self.msa_variant: self.muscle_aln_for_variant()

        if not os.path.exists(directory):
            os.makedirs(directory)
        for hist_type in ['H2A', 'H2B', 'H3', 'H4', 'H1']:
            for hist_var in set(self.curated_data[self.curated_data['type']==hist_type]['variant']):
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