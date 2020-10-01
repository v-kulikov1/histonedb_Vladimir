from django.shortcuts import render

from browse.models import *
from human_hist.models import *

def human_histones(request):
    return render(request, 'human_hist.html', {'human_histones': Human_histones_genes.objects.all()})
