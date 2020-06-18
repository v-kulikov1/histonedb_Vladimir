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

from browse.models import *
from djangophylocore.models import Taxonomy

class InvalidFASTA(Exception):
    pass

log = logging.getLogger(__name__)

def make_blastp(sequences):

  if not isinstance(sequences, list):
    sequences = [sequences]

  blastp =  os.path.join(os.path.dirname(sys.executable), "blastp")
  output = os.path.join("/", "tmp", "{}.xml".format(uuid.uuid4()))
  blastp_cline = NcbiblastpCommandline(
    cmd=blastp,
    db=os.path.join(settings.STATIC_ROOT_AUX, "browse", "blast", "HistoneDB_sequences_0.fa"),
    evalue=0.005, outfmt=5)
  return blastp_cline(stdin="\n".join([s.format("fasta") for s in sequences]))

def load_blast_search(blastp_result):
  blastFile = io.StringIO()
  blastFile.write(blastp_result)
  blastFile.seek(0)

  for i, blast_record in enumerate(NCBIXML.parse(blastFile)):
    if len(blast_record.alignments) == 0:
      raise InvalidFASTA("No blast hits for {}.".format(blast_record.query))
    for alignment in blast_record.alignments:
      variant_model = Variant.objects.get(id=alignment.hit_def.split("|")[2])
      load_blasthsps(blast_record.query, alignment.hsps, variant_model)

def add_sequence(accession, variant_model, taxonomy, header, sequence):
  """Add sequence into the database, autfilling empty Parameters"""
  seq = Sequence(
    id             = accession,
    variant        = None,
    variant_blast  = variant_model,
    gene           = None,
    splice         = None,
    taxonomy       = taxonomy,
    header         = header[:250],
    sequence       = str(sequence).replace("-", "").upper(),
    reviewed       = False,
    )
  seq.save()
  return seq

def add_score(seq, variant_model, hsp, best=False):
  """Add score for a given sequence"""
  # score_num = Score.objects.count()+1
  score = ScoreBlast(
    # id                      = score_num,
    sequence                = seq,
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

def load_blasthsps(header, hsps, variant_model):
  ###Iterate through high scoring fragments.
  for hsp in hsps:
    # hmmthreshold_passed = hsp.score >= 100.0
    accession = header.split('|')[2]
    # accession = header.split("|")[1]
    seqs = Sequence.objects.filter(id=accession)
    if len(seqs) > 0:
      # Sequence already exists. Compare bit scores, if loaded bit score is
      # greater than current, reassign variant and update scores. Else, append score
      seq = seqs.first()
      if (seq.reviewed == True):
        continue  # we do not want to alter a reviewed sequence!
      # print "Already in database", seq

      best_scores = seq.all_model_scores.filter(used_for_classification=True)
      if len(best_scores) > 0:
        ##Sequence have passed the threshold for one of previous models.
        best_score = best_scores.first()  # Why we extract first?
        if hsp.bitscore > best_score.score:
          # best scoring
          seq.variant_blast = variant_model
          seq.sequence = str(hsp.query)
          best_score_2 = Score.objects.get(id=best_score.id)
          best_score_2.used_for_classification = False
          best_score_2.save()
          seq.save()
          # print "UPDATED VARIANT"
          add_score(seq, variant_model, hsp, best=True)
        else:
          add_score(seq, variant_model, hsp, best=False)
      else:
        # No previous model passed the threshold, it is the first
        seq.variant_blast = variant_model
        seq.sequence = str(hsp.query)
        seq.save()
        add_score(seq, variant_model, hsp, best=True)
    else:
      ##A new sequence is found that passed treshold.
      taxonomy = taxonomy_from_header(header, accession)
      sequence = Seq(str(hsp.query))
      try:
        seq = add_sequence(
          accession,
          variant_model,
          taxonomy,
          header,
          sequence)
        add_score(seq, variant_model, hsp, best=True)
      except IntegrityError as e:
        log.error("Error adding sequence {}".format(seq))
        global already_exists
        already_exists.append(accession)
        continue
