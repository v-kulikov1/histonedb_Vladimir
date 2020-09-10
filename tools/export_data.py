from browse.models import Sequence, Score, ScoreBlast, SequenceBlast
from django.conf import settings
# from django.db.models import Q

#This script is used to export tables from database for futher use in research

import os

# with open(os.path.join(settings.STATIC_ROOT_AUX, "browse", "dumps", "{}.txt".format('seqs')),'w') as f:
#     f.write("accession,hist_type,hist_var,taxid,curated\n")
#     for seq in Sequence.objects.all():
#         f.write("%s,%s,%s,%s,%s\n"%(seq.id,seq.variant.hist_type,seq.variant,seq.taxonomy_id,seq.reviewed))
#
# with open(os.path.join(settings.STATIC_ROOT_AUX, "browse", "dumps", "{}.txt".format('scores')),'w') as f:
#     f.write("accession,hmm_model,score,used_for_classification\n")
#     for s in Score.objects.all():
#         f.write("%s,%s,%s,%s\n"%(s.sequence.id,s.variant,s.score,s.used_for_classification))
#
with open(os.path.join(settings.STATIC_ROOT_AUX, "browse", "dumps", "{}.txt".format('scores_blast')),'w') as f:
    f.write("accession,blast_model,score,bit_score,evalue,hsp_length,used_for_classification,hit_accession,sequence,hit_sequence,match,blastStart,blastEnd,seqStart,seqEnd\n")
    for s in ScoreBlast.objects.all():
        sequence_obj = Sequence.objects.get(id=s.sequence.id)
        print(s.hit_accession)
        hit_seq = Sequence.objects.get(id=s.hit_accession)
        f.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n"%(s.sequence.id,s.variant,s.score,s.bitScore,s.evalue,s.align_length,
                                                                  s.used_for_classification,s.hit_accession,sequence_obj.sequence,hit_seq.sequence,s.match,
                                                                  s.blastStart,s.blastEnd,s.seqStart,s.seqEnd))

with open(os.path.join(settings.STATIC_ROOT_AUX, "browse", "dumps", "{}.txt".format('seqs_with_blast')),'w') as f:
    f.write("accession,hist_type,hist_var_hmm,hist_var_blast,taxid,curated,score_hmm,score_blast,bitscore_blast\n")
    for seq in Sequence.objects.all():
        try:
            seq_blast = SequenceBlast.objects.get(id=seq.id)
            hist_var_blast = seq_blast.variant.id
            score_obj = seq_blast.all_model_scores.filter(used_for_classification=True).first()
            score_blast = score_obj.score
            bitscore_blast = score_obj.bitScore
        except:
            hist_var_blast = ''
            score_blast = ''
            bitscore_blast = ''
        f.write("%s,%s,%s,%s,%s,%s,%s,%s,%s\n"%(seq.id,seq.variant.hist_type,seq.variant,hist_var_blast,seq.taxonomy_id,seq.reviewed,
                                             seq.all_model_scores.filter(used_for_classification=True).first().score,
                                             score_blast, bitscore_blast))

