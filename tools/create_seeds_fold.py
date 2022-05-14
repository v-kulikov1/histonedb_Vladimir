"""
@author: preety
Goes through static/browse/seeds directories and creates similar structure in static/browse/seeds_fold with clipped sequences.
Directory static/browse/seeds_fold contains similar seed sequences clipped to contain only histone folds.
"""

from tools.browse_service import *
from tools.hist_ss import get_variant_features

from Bio.Align import MultipleSeqAlignment
from Bio.Align.AlignInfo import SummaryInfo
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

import os

if not os.path.exists(FOLD_SEED_DIRECTORY):
    os.makedirs(FOLD_SEED_DIRECTORY)
for i, (root, _, files) in enumerate(os.walk(SEED_DIRECTORY)):
    for seed in files:
        if not seed.endswith(".fasta"): continue
        hist_type = os.path.basename(root) if os.path.basename(root) != 'seeds' else seed[:-6]
        create_file = os.path.join(FOLD_SEED_DIRECTORY, hist_type, seed) if os.path.basename(
            root) != 'seeds' else os.path.join(FOLD_SEED_DIRECTORY, seed)
        variant = Variant.objects.get(
            id='canonical_{}'.format(hist_type)) if hist_type != 'H1' else Variant.objects.get(id='generic_H1')
        seqFile = os.path.join(root, seed)
        sequences = list(SeqIO.parse(seqFile, "fasta"))
        try:
            msa = MultipleSeqAlignment(sequences)
            a = SummaryInfo(msa)
            cons = Sequence(id="Consensus", variant_id=variant.id, taxonomy_id=1,
                            sequence=str(a.dumb_consensus(threshold=0.1, ambiguous='X')))
            features = get_variant_features(cons, variants=[variant], save_gff=False, only_general=True)
            cutting_params = [next(
                filter(lambda x: x.id.split('_', 2)[-1].strip() == 'General{}_root_alpha1'.format(hist_type),
                       features)).start,
                              next(filter(
                                  lambda x: x.id.split('_', 2)[-1].strip() == 'General{}_root_alpha3'.format(hist_type),
                                  features)).end]
            # cutting_params = list(filter(None, [feature.start if feature.id=='General{}_root_alpha1'.format(seed[:-6]) else feature.end if feature.id=='General{}_root_alpha3'.format(seed[:-6]) else None for feature in features]))
            try:
                fd = open(create_file, 'w')
            except FileNotFoundError:
                os.makedirs(os.path.join(FOLD_SEED_DIRECTORY, hist_type))
                fd = open(create_file, 'w')
            for s in SeqIO.parse(seqFile, "fasta"):
                s.seq = Seq(str(s.seq)[cutting_params[0]:cutting_params[1] + 1], IUPAC.protein)
                SeqIO.write(s, fd, "fasta")
            fd.close()
        except StopIteration:
            self.log.warning('StopIteration: hist_type {}, seed {}'.format(hist_type, seed))
            try:
                fd = open(create_file, 'w')
            except FileNotFoundError:
                os.makedirs(os.path.join(FOLD_SEED_DIRECTORY, hist_type))
                fd = open(create_file, 'w')
            for s in SeqIO.parse(seqFile, "fasta"):
                SeqIO.write(s, fd, "fasta")
            fd.close()
self.log.info('Sequences clipped. See static/browse/seeds_fold directory.')
