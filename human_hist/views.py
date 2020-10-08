from django.shortcuts import render

from browse.models import *
from human_hist.models import Histone_Human_genes
from djangophylocore.models import *

def human_histones(request):
    return render(request, 'human_hist.html', {'human_histones': Histone_Human_genes.objects.all()})


