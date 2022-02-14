import os, sys, subprocess
import logging
from tqdm import tqdm

import pandas as pd
from Bio import SearchIO
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML

log = logging.getLogger(__name__)

class TypePredictionInfo(object):
    def __init__(self):
        self.accession = None
        self.histone_type = None
        self.score = None
        self.best = None
        self.description = None
        self.hsp = None

def predict_types(sequences, hmmout, hmmdb, E=10, result_file=None):
    """
    This method do a search among given sequences and extracts just histone variants classified by histone types.
    All the results saves to the result_file, e.g.:
    accession, histone_type, score, best, description
    NP_563627.1, H3, 2.14, True, Histone superfamily protein [Arabidopsis thaliana]
    DAA13058.1, H2B, 3.11, False, [Bos taurus]
    ___________
    Parameters:
    ___________
    sequences: fasta file
    hmmout: hmmout file
    hmmdb: hmmdb file
    result_file: file to save results
    """

    # Use HMMs to search the nr database
    log.info("Searching HMMs...")
    log.info(" ".join(["nice", "hmmsearch", "-o", hmmout, "-E", str(E), "--notextw", hmmdb, sequences]))
    subprocess.call(["hmmsearch", "-o", hmmout, "-E", str(E), "--notextw", hmmdb, sequences])

    log.info("Parsing search results...")
    return parse_hmm_search_out(hmmout=hmmout, save_to=result_file)

def predict_variants(sequences, blastout, blastdb, result_file=None):
    """
    This method do a search among given sequences and classify by histone variants.
    All the results saves to the result_file, e.g.: needs correction
    accession, histone_type, score, best, description
    NP_563627.1, H3, 2.14, True, Histone superfamily protein [Arabidopsis thaliana]
    DAA13058.1, H2B, 3.11, False, [Bos taurus]
    ___________
    Parameters:
    ___________
    sequences: list of fasta str
    blastout: hmmout file
    blastdb: hmmdb file
    result_file: file to save results
    """

    # log.info("Predicting variants via BLAST")
    log.info("Running BLASTP for {} sequences...".format(len(sequences)))
    make_blastp(sequences, blastdb, save_to=blastout)
    return parse_blast_search_out(blast_file=blastout, save_to=result_file)
    # checkH2AX()

#HMM
def parse_hmm_search_out(hmmout, save_to=None):
    """Parse hmmout file and return results as dict.
      Parameters:
      ___________
      hmmout: hmmout file
      return e.g.: [{'accession': NP_563627.1, 'histone_type': H3, 'score': 2.14, 'best': True,
                'description': Histone superfamily protein [Arabidopsis thaliana], 'hsp': best HSPObject},
                {'accession': DAA13058.1, 'histone_type': H2B, 'score': 3.11, 'best': False,
                'description': [Bos taurus], 'hsp': best HSPObject},]
      """
    result = []
    for histone_query in tqdm(SearchIO.parse(hmmout, "hmmer3-text")):
        log.info("Loading histone: {}".format(histone_query.id))  # histone_query.id is histone type
        for hit in tqdm(histone_query):
            best_hsp = max(hit, key=lambda hsp: hsp.bitscore)
            best_f = True
            existing_best_result = list(filter(lambda d: d['accession'] == hit.id and d['best'], result))[0]
            if best_hsp.bitscore > existing_best_result['score']:
                result.remove(existing_best_result)
                existing_best_result['best'] = False
                result.append(existing_best_result)
            else: best_f = False
            result.append({'accession': hit.id,  # hit.id is a first accession
                           'histone_type': histone_query.id,
                           'score': best_hsp.bitscore,
                           'best': best_f,
                           'description': hit.description,
                           'hsp': best_hsp,
                           })
    if save_to:
        pd.DataFrame(result).fillna('').drop(columns=['hsp']).to_csv(save_to, index=False)
        log.info("Predicted sequences saved to {}".format(save_to))
    return result

# BLAST
def make_blastp(sequences, blastdb, save_to):
    # log.error('Error:: sequences {}'.format(len(sequences)))
    # log.error('Error:: hasattr {}'.format(hasattr(sequences, '__iter__')))
    if not hasattr(sequences, '__iter__'):
        sequences = [sequences]

    blastp = os.path.join(os.path.dirname(sys.executable), "blastp")
    # output = os.path.join("/", "tmp", "{}.xml".format(uuid.uuid4()))
    blastp_cline = NcbiblastpCommandline(
        cmd=blastp,
        # db=os.path.join(settings.STATIC_ROOT_AUX, "browse", "blast", "HistoneDB_sequences.fa"),
        db=blastdb,
        evalue=.01, outfmt=5)
    # evalue=0.004, outfmt=5)
    result, error = blastp_cline(stdin="\n".join([s.format("fasta") for s in sequences]))

    with open(save_to, 'w') as f:
        f.write(result)

    log.info('Blast results saved to {}'.format(save_to))

def get_best_hsp(hsps, align_longer=0):
    best_alignment_hsp = hsps[0]
    for hsp in hsps[1:]:
        if hsp.score > best_alignment_hsp.score and hsp.align_length > align_longer:
            best_alignment_hsp = hsp
    return best_alignment_hsp

def check_features_macroH2A(query_accession, hsp_hit_start, hsp_hit_end):
    feature = Feature.objects.filter(template__variant='macroH2A', name='Macro domain').first()
    ratio = (min(hsp_hit_end, feature.end) - max(hsp_hit_start, feature.start)) / (feature.end - feature.start)
    log.info('{} expected as macroH2A with ratio={} of macro domain contained in hsp'.format(query_accession, ratio))
    if ratio < .8:
        log.info('{} expected as macroH2A cannot pass 0.8 ratio_treshhold'.format(query_accession))
        return False
    return True

def parse_blast_search_out(blast_file, save_to=None):
    """Parse blastFile file and return results as dict.
        Parameters:
        ___________
        blastFile: blastFile
        return e.g.: [{'accession': NP_563627.1, 'histone_variant': cenH3, 'score': 2.14, 'best': True,
                        'description': Histone superfamily protein [Arabidopsis thaliana], 'hsp': best HSPObject, 'hit_accession': DAA13058.1},
                        {'accession': DAA13058.1, 'histone_variant': H2B.W, 'score': 3.11, 'best': False,
                        'description': [Bos taurus], 'hsp': best HSPObject, 'hit_accession': DAA13058.1},]
    """
    blast_file_handle = open(blast_file)
    result = []
    for i, blast_record in enumerate(NCBIXML.parse(blast_file_handle)): # <Iteration>
        query_split = blast_record.query.split('|')
        accession, description = query_split[0], query_split[-1]
        if accession == 'pir' or accession == 'prf':
            accession = '|'.join(query_split[:3])
            log.info('Non-standard accession {} got from {}'.format(accession, blast_record.query))

        if len(blast_record.alignments) == 0: # <Hit> count = 0
            # log.info("No blast hits for {} with e-value {}".format(blast_record.query, blast_record.descriptions[0]))
            log.info("No blast hits for {}".format(blast_record.query))
            # raise InvalidFASTA("No blast hits for {}.".format(blast_record.query))
            result.append({'accession': accession, 'variant': 'generic',
                           'score': .000001, 'best': True, 'description': description,
                           'hsp': None, 'hit_accession': None})
            continue

        best_alignments = []
        for alignment in blast_record.alignments: # <Hit>
            best_algn_hsp = get_best_hsp(alignment.hsps, align_longer=.25 * blast_record.query_letters)
            best_alignments.append({'hit_accession': alignment.hit_def.split("|")[0],
                                    'hit_variant': alignment.hit_def.split("|")[2],
                                    'best_hsp': best_algn_hsp,
                                    'score': best_algn_hsp.score})
        best_alignments = sorted(best_alignments, key=lambda algn: algn['score'], reverse=True)

        # If hsp contains macro domain?
        if best_alignments[0]['hit_variant'] == 'macroH2A':
            f = check_features_macroH2A(query_accession=accession,
                                        hsp_hit_start=best_alignments[0]['best_hsp'].sbjct_start, # get start and end of hit HSP
                                        hsp_hit_end=best_alignments[0]['best_hsp'].sbjct_end)
            if not f: continue

        histone_variant = best_alignments[0]['hit_variant']
        # if best_alignments[0]['best_hsp'].score < variant_model.blastthreshold:
        #   histone_variant = 'generic'

        result.append({'accession': accession, 'variant': histone_variant,
                       'score': best_alignments[0]['score'], 'best': True, 'description': description,
                       'hsp': best_alignments[0]['best_hsp'], 'hit_accession': best_alignments[0]['hit_accession']})

        for best_algn in best_alignments[1:]:
            result.append({'accession': accession, 'variant': best_algn['hit_variant'],
                           'score': best_algn['score'], 'best': True, 'description': description,
                           'hsp': best_algn['best_hsp'], 'hit_accession': best_algn['hit_accession']})

    blast_file_handle.close()
    if save_to:
        pd.DataFrame(result).fillna('').drop(columns=['hsp', 'hit_accession']).to_csv(save_to, index=False)
        log.info(f"Predicted sequences saved to {save_to}")
    return result