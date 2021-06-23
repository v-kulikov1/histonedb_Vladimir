from browse.models import Sequence, Score
from django.conf import settings
# from django.db.models import Q

#This script is used to export tables from database for futher use in research

import os
from human_hist.models import Human_variants

with open(os.path.join(settings.STATIC_ROOT_AUX, "browse", "dumps", "{}.txt".format('seqs')),'w') as f:
    f.write("accession,hist_type,hist_var,taxid,curated\n")
    for seq in Sequence.objects.all():
        f.write("%s,%s,%s,%s,%s\n"%(seq.id,seq.variant.hist_type,seq.variant,seq.taxonomy_id,seq.reviewed))

with open(os.path.join(settings.STATIC_ROOT_AUX, "browse", "dumps", "{}.txt".format('scores')),'w') as f:
    f.write("accession,hmm_model,score,used_for_classification\n")
    for s in Score.objects.all():
        f.write("%s,%s,%s,%s\n"%(s.sequence.id,s.variant,s.score,s.used_for_classification))

with open(os.path.join(settings.STATIC_ROOT_AUX, "browse", "dumps", "{}.txt".format('human_hist')),'w') as f:
    f.write("variant,id\n")
    for s in Human_variants.objects.all():
        f.write("%s,%s\n"%(s.variant,s.id))

