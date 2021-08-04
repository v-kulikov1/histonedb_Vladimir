from browse.models import Sequence, ScoreHmm, Score, SequenceBlast, Feature
from django.conf import settings
# from django.db.models import Q

#This script is used to export tables from database for futher use in research

import os
from datetime import date, datetime

from django.core.exceptions import ObjectDoesNotExist

now = datetime.now()
dt_string = now.strftime("%Y%m%d-%H%M%S")

with open(os.path.join(settings.STATIC_ROOT_AUX, "browse", "dumps", "{}_{}.txt".format('seqs', dt_string)),'w') as f:
    f.write("accession,hist_type,hist_var, hist_var_hmm,taxid,curated\n")
    for seq in Sequence.objects.all():
        if seq.variant is None: continue
        f.write("%s,%s,%s,%s,%s,%s\n"%(seq.id,seq.variant.hist_type,seq.variant,seq.variant_hmm,seq.taxonomy_id,seq.reviewed))

with open(os.path.join(settings.STATIC_ROOT_AUX, "browse", "dumps", "{}_{}.txt".format('scores', dt_string)),'w') as f:
    f.write("accession,hmm_model,score,used_for_classification\n")
    for s in ScoreHmm.objects.all():
        f.write("%s,%s,%s,%s\n"%(s.sequence.id,s.variant_hmm,s.score,s.used_for_classification))

with open(os.path.join(settings.STATIC_ROOT_AUX, "browse", "dumps", "{}_{}.txt".format('scores_blast', dt_string)),'w') as f:
    f.write("accession,blast_model,score,bit_score,evalue,hsp_length,used_for_classification,hit_accession,sequence,hit_sequence,match,blastStart,blastEnd,seqStart,seqEnd\n")
    for s in Score.objects.all():
        # sequence_obj = Sequence.objects.get(id=s.sequence.id)
        # print(s.hit_accession)
        if s.hit_accession=='':
            hit_seq = ''
        else:
            hit_seq = Sequence.objects.get(id=s.hit_accession).sequence
        f.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n"%(s.sequence.id,s.variant,s.score,s.bitScore,s.evalue,s.align_length,
                                                                  s.used_for_classification,s.hit_accession,s.sequence.sequence,hit_seq,s.match,
                                                                  s.blastStart,s.blastEnd,s.seqStart,s.seqEnd))

with open(os.path.join(settings.STATIC_ROOT_AUX, "browse", "dumps", "{}_{}..txt".format('features', dt_string)),'w') as f:
    f.write("id,variant,name,start,end\n")
    for feature in Feature.objects.all():
        f.write("%s,%s,%s,%s,%s\n"%(feature.id,feature.template.variant,feature.name,feature.start,feature.end))


