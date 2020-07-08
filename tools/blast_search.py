import os
import sys
import subprocess
import io
import logging
import uuid

#BioPython
from Bio import SeqIO, SearchIO
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio import Entrez

#Django libraires
from django.db.models import Max
from django.db.utils import IntegrityError
from django.db import close_old_connections
from django.conf import settings
from django.db import transaction

#Custom librairies
from tools.taxonomy_from_accessions import taxonomy_from_header, easytaxonomy_from_header, taxonomy_from_accessions, update_taxonomy


from tqdm import tqdm
import pickle

from browse.models import *
from djangophylocore.models import Taxonomy

class InvalidFASTA(Exception):
    pass

log = logging.getLogger(__name__)

def make_blastp(sequences, save_to):

  if not isinstance(sequences, list):
    sequences = [sequences]

  blastp =  os.path.join(os.path.dirname(sys.executable), "blastp")
  output = os.path.join("/", "tmp", "{}.xml".format(uuid.uuid4()))
  blastp_cline = NcbiblastpCommandline(
    cmd=blastp,
    db=os.path.join(settings.STATIC_ROOT_AUX, "browse", "blast", "HistoneDB_sequences.fa"),
    evalue=0.005, outfmt=5)
  result, error = blastp_cline(stdin="\n".join([s.format("fasta") for s in sequences]))

  with open(save_to, 'w') as f:
    f.write(result)

  log.info('Blast results saved to {}'.format(save_to))
  # if save_to:
  #   with open(save_to, 'wb') as f:
  #     pickle.dump(result, f)
  # return result

def load_blast_search(blastFile_name):
  # blastFile = io.StringIO()
  # blastFile.write(blastp_result)
  # blastFile.seek(0)
  blastFile = open(blastFile_name)
  # log.info('Loaded Blast results from {}'.format(blastFile_name))

  for i, blast_record in enumerate(NCBIXML.parse(blastFile)):
    accession = blast_record.query.split('|')[0]
    # accession = blast_record.query.split("|")[1]
    # log.info('Loading {}: {}'.format(accession, blast_record.query))
    if len(blast_record.alignments) == 0:
      raise InvalidFASTA("No blast hits for {}.".format(blast_record.query))
    for alignment in blast_record.alignments:
      variant_model = Variant.objects.get(id=alignment.hit_def.split("|")[2])
      # log.info('Alignment of variant {}'.format(alignment.hit_def))
      load_blasthsps(accession, blast_record.query, alignment.hsps, variant_model)

    # if i%10000==0:
    #   log.info('Loaded {} BlastRecords'.format(i))

def add_sequence(accession, variant_model):
  """Add sequence into the database, autfilling empty Parameters"""
  seq = SequenceBlast(
    accession      = accession,
    variant        = variant_model,
    )
  seq.save()
  return seq

def add_score(seq_blast, variant_model, hsp, best=False):
  """Add score for a given sequence"""
  score = ScoreBlast(
    # id                      = score_num,
    sequence                = seq_blast,
    variant                 = variant_model,
    score                   = hsp.score,
    bitScore                = hsp.bits,
    evalue                  = hsp.expect,
    blastStart              = hsp.query_start,
    blastEnd                = hsp.query_end,
    seqStart                = hsp.sbjct_start,
    seqEnd                  = hsp.sbjct_end,
    align_length            = hsp.align_length,
    used_for_classification = best,
    )
  score.save()

def load_blasthsps(accession, header, hsps, variant_model):
  ###Iterate through high scoring fragments.
  for hsp in hsps:
    # hmmthreshold_passed = hsp.score >= 100.0

    try:
      # seq = Sequence.objects.filter(id=accession).first()
      seqs_blast = SequenceBlast.objects.filter(accession=accession)
    except AttributeError:
      log.info('AttributeError: BlastRecord for {} with full header {}'.format(accession, header))
      accession = '|'.join(header.split('|')[:3])
      log.info('New accession {} got'.format(accession))
      # seq = Sequence.objects.filter(id=accession).first()
      seqs_blast = SequenceBlast.objects.filter(accession=accession)

    if len(seqs_blast) > 0:
      seq_blast = seqs_blast.first()
      best_scores = seq_blast.all_model_scores.filter(used_for_classification=True)
      if len(best_scores) > 0:
        ##Sequence have passed the threshold for one of previous models.
        best_score = best_scores.first()  # Why we extract first?
        if hsp.score > best_score.score:
          # best scoring
          seq_blast.variant = variant_model
          best_score_2 = ScoreBlast.objects.get(id=best_score.id)
          best_score_2.used_for_classification = False
          best_score_2.save()
          seq_blast.save()
          # print "UPDATED VARIANT"
          add_score(seq_blast, variant_model, hsp, best=True)
        else:
          add_score(seq_blast, variant_model, hsp, best=False)
      else:
        # No previous model passed the threshold, it is the first
        seq_blast.variant = variant_model
        seq_blast.save()
        add_score(seq_blast, variant_model, hsp, best=True)
    else:
      # log.error("New Sequence found for {}".format(header))
      ##A new sequence is found that passed treshold.
      try:
        seq_blast = add_sequence(
          accession,
          variant_model)
        add_score(seq_blast, variant_model, hsp, best=True)
      except IntegrityError as e:
        log.error("Error adding sequence {}".format(seq_blast))
