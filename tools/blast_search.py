import os
import sys
import subprocess
import io
import logging

# BioPython
from Bio import SeqIO, SearchIO
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio import Entrez

# Django libraires
from django.db.models import Max
from django.db.utils import IntegrityError
from django.db import close_old_connections
from django.conf import settings
from django.db import transaction

# Custom librairies
from tools.taxonomy_from_accessions import taxonomy_from_header, easytaxonomy_from_header, taxonomy_from_accessions, \
    update_taxonomy

from tqdm import tqdm
import pickle

from browse.models import *
from djangophylocore.models import Taxonomy


class InvalidFASTA(Exception):
    pass


log = logging.getLogger(__name__)


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


def add_score(seq, variant_model, hsp, hit_accession, best=False):
    """Add score for a given sequence"""
    score = ScoreBlast(
        # id                      = score_num,
        sequence=seq,
        variant=variant_model,
        score=hsp.score,
        bitScore=hsp.bits,
        evalue=hsp.expect,
        blastStart=hsp.query_start,
        blastEnd=hsp.query_end,
        seqStart=hsp.sbjct_start,
        seqEnd=hsp.sbjct_end,
        align_length=hsp.align_length,
        match=hsp.match,
        hit_accession=hit_accession,
        used_for_classification=best,
    )
    score.save()


def load_blast_search(blastFile_name):
    blastFile = open(blastFile_name)

    for i, blast_record in enumerate(NCBIXML.parse(blastFile)):
        query_split = blast_record.query.split('|')
        accession = query_split[0]
        if accession == 'pir' or accession == 'prf':
            accession = '|'.join(query_split[:3])
            log.info('Non-standard accession {} got from {}'.format(accession, blast_record.query))

        seq = Sequence.objects.get(id=accession)

        if len(blast_record.alignments) == 0:
            # log.info("No blast hits for {} with e-value {}".format(blast_record.query, blast_record.descriptions[0]))
            log.info("No blast hits for {}".format(blast_record.query))
            # raise InvalidFASTA("No blast hits for {}.".format(blast_record.query))
            seq.variant = Variant.objects.get(id='generic_{}'.format(seq.variant_hmm.hist_type.id))
            seq.save()
            score = ScoreBlast(
                sequence=seq,
                variant=seq.variant,
                score=.000001,
                bitScore=None,
                evalue=.011,
                blastStart=None,
                blastEnd=None,
                seqStart=None,
                seqEnd=None,
                align_length=None,
                used_for_classification=True,
            )
            score.save()
            continue

        # best_alignment = blast_record.alignments[0]
        # best_hsp = get_best_hsp(best_alignment.hsps, align_longer=.25*seq_len)
        # for alignment in blast_record.alignments[1:]:
        #   best_alignment_hsp = get_best_hsp(alignment.hsps, align_longer=.25*seq_len)
        #   if best_alignment_hsp.score > best_hsp.score:
        #     best_alignment = alignment
        #     best_hsp = best_alignment_hsp

        # variant_model = Variant.objects.get(id=best_alignment.hit_def.split("|")[2])
        # seq.variant_hmm = seq.variant
        # # print('{}: {}'.format(variant_model.id, variant_model.blastthreshold))
        # if best_hsp.score < variant_model.blastthreshold:
        #   seq.variant = Variant.objects.get(id='generic_{}'.format(variant_model.hist_type.id))
        # else:
        #   seq.variant = variant_model
        # seq.save()
        # # log.info('Hit accession {}'.format(best_alignment.hit_def))
        # add_score(seq, variant_model, best_hsp, best_alignment.hit_def.split("|")[0], best=True)

        best_alignments = []
        for alignment in blast_record.alignments:
            best_algn_hsp = get_best_hsp(alignment.hsps, align_longer=.25 * len(seq.sequence))
            best_alignments.append({'accession': alignment.hit_def.split("|")[0],
                                    'variant': alignment.hit_def.split("|")[2],
                                    'best_hsp': best_algn_hsp,
                                    'score': best_algn_hsp.score})
        best_alignments = sorted(best_alignments, key=lambda algn: algn['score'], reverse=True)

        # If hsp contains macro domain?
        if best_alignments[0]['variant'] == 'macroH2A':
            feature = Feature.objects.filter(template__variant=best_alignments[0]['variant'],
                                             name='Macro domain').first()
            hsp_start = best_alignments[0]['best_hsp'].sbjct_start  # get start and end of hit HSP
            hsp_end = best_alignments[0]['best_hsp'].sbjct_end
            ratio = (min(hsp_end, feature.end) - max(hsp_start, feature.start)) / (feature.end - feature.start)
            log.info('{} expected as macroH2A with ratio={} of macro domain contained in hsp'.format(seq.id, ratio))
            if ratio < .8:
                log.info('{} expected as macroH2A cannot pass 0.8 ratio_treshhold'.format(seq.id))
                continue

        variant_model = Variant.objects.get(id=best_alignments[0]['variant'])
        # print('DEBUG::{}: {}'.format(variant_model.id, variant_model.blastthreshold))
        # print('DEBUG::best_score: {}'.format(best_alignments[0]['best_hsp'].score))
        # print('--------------------------------------------------------------')
        # break
        # if best_alignments[0]['best_hsp'].score < variant_model.blastthreshold:
        #   seq.variant = Variant.objects.get(id='generic_{}'.format(variant_model.hist_type.id))
        # else:
        #   seq.variant = variant_model
        seq.variant = variant_model
        seq.save()
        # log.info('Hit accession {}'.format(best_alignment.hit_def))
        add_score(seq, variant_model, best_alignments[0]['best_hsp'], best_alignments[0]['accession'], best=True)

        # for best_algn in best_alignments[1:]:
        #   add_score(seq, Variant.objects.get(id=best_algn['variant']), best_algn['best_hsp'], best_algn['accession'], best=False)

    blastFile.close()


def get_best_hsp(hsps, align_longer=0):
    best_alignment_hsp = hsps[0]
    for hsp in hsps[1:]:
        if hsp.score > best_alignment_hsp.score and hsp.align_length > align_longer:
            best_alignment_hsp = hsp
    return best_alignment_hsp


def load_blast_search_diagnosis(blastFile_name):
    # blastFile = io.StringIO()
    # blastFile.write(blastp_result)
    # blastFile.seek(0)
    blastFile = open(blastFile_name)
    # log.info('Loaded Blast results from {}'.format(blastFile_name))

    for i, blast_record in enumerate(NCBIXML.parse(blastFile)):
        query_split = blast_record.query.split('|')
        if len(query_split) == 3:
            accession = query_split[0]
        else:
            accession = '|'.join(query_split[:3])
            log.info('Non-standard accession {} got from {}'.format(accession, blast_record.query))
        # accession = blast_record.query.split('|')[0]
        # log.info('Loading {}: {}'.format(accession, blast_record.query))
        if len(blast_record.alignments) == 0:
            log.info("No blast hits for {}".format(blast_record.query))
            # raise InvalidFASTA("No blast hits for {}.".format(blast_record.query))
            continue
        for alignment in blast_record.alignments:
            variant_model = Variant.objects.get(id=alignment.hit_def.split("|")[2])
            # log.info('Alignment of variant {}'.format(alignment.hit_def))
            load_blasthsps_diagnosis(accession, blast_record.query, alignment.hsps, variant_model,
                                     alignment.hit_def.split("|")[0])

        # if i%10000==0:
        #   log.info('Loaded {} BlastRecords'.format(i))


def load_blasthsps_diagnosis(accession, header, hsps, variant_model, hit_accession):
    ###Iterate through high scoring fragments.
    for hsp in hsps:
        seqs = Sequence.objects.filter(id=accession)
        if len(seqs) < 1:
            log.error("New sequence is found: {}. This is strange.".format(accession))
            continue

        seq = seqs.first()
        best_scores = seq.all_model_blast_scores.filter(used_for_classification=True)
        if len(best_scores) > 0:
            ##Sequence have passed the threshold for one of previous models.
            best_score = best_scores.first()
            if hsp.score > best_score.score and hsp.align_length > len(seq.sequence):
                # best scoring
                seq.variant_hmm = seq.variant
                seq.variant = variant_model
                best_score_2 = ScoreBlast.objects.get(id=best_score.id)
                best_score_2.used_for_classification = False
                best_score_2.save()
                seq.save()
                # print "UPDATED VARIANT"
                add_score(seq, variant_model, hsp, hit_accession, best=True)
            else:
                add_score(seq, variant_model, hsp, hit_accession, best=False)
        elif hsp.align_length > len(seq.sequence):
            # No previous model passed the threshold, it is the first
            seq.variant_hmm = seq.variant
            seq.variant = variant_model
            seq.save()
            add_score(seq, variant_model, hsp, hit_accession, best=True)
        else:
            add_score(seq, variant_model, hsp, hit_accession, best=False)
