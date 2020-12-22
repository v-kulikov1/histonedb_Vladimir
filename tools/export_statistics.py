from browse.models import Sequence, Score, ScoreBlast, SequenceBlast, Histone, Variant
from django.conf import settings
# from django.db.models import Q

#This script is used to export tables from database for futher use in research

import os
from datetime import date, datetime
from pathlib import Path

import numpy as np
import pandas as pd

from ete3 import NCBITaxa
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

NOW = datetime.now()
DT_STRING = NOW.strftime("%Y%m%d-%H%M%S")

STAT_DIR = os.path.join(settings.STATIC_ROOT_AUX, "browse", "statistics")
CURR_STAT_DIR = os.path.join(STAT_DIR, DT_STRING)
CURR_STAT_DIR = os.path.join(STAT_DIR, '20201219-195135')
DATA_DIR = os.path.join(CURR_STAT_DIR, 'data')

SEQS_FILE = os.path.join(DATA_DIR, 'seqs.csv')
SCORES_FILE = os.path.join(DATA_DIR, 'scores.csv')
SCORES_HMM_FILE = os.path.join(DATA_DIR, 'scores_hmm.csv')

HIGHLVL_TAXA = ['Metazoa', 'Fungi',
                'Rhodophyta', 'Viridiplantae',
                'Rhizaria', 'Alveolata', 'Stramenopiles']
NCBI_TAXA = NCBITaxa()
TAXIDS_DICT = NCBI_TAXA.get_name_translator(HIGHLVL_TAXA)
TAXIDS = [TAXIDS_DICT[name][0] for name in HIGHLVL_TAXA]

def create_directories():
    Path(CURR_STAT_DIR).mkdir(parents=True, exist_ok=True)
    Path(DATA_DIR).mkdir(parents=True, exist_ok=True)
    for hlvl in HIGHLVL_TAXA:
        Path(os.path.join(CURR_STAT_DIR, hlvl)).mkdir(parents=True, exist_ok=True)
    for h in ['H1', 'H2A', 'H2B', 'H3', 'H4']:
        Path(os.path.join(CURR_STAT_DIR, h)).mkdir(parents=True, exist_ok=True)
        for hlvl in HIGHLVL_TAXA:
            Path(os.path.join(CURR_STAT_DIR, h, hlvl)).mkdir(parents=True, exist_ok=True)

def get_seqs_data(rewrite = False, curated=None, hist_type=None, hist_var=None, hist_var_hmm=None, taxid=None, highlvl_taxid=None):
    if Path(SEQS_FILE).exists() and not rewrite:
        # return pd.read_csv(SEQS_FILE)
        return filter_seqs_data(pd.read_csv(SEQS_FILE), curated, hist_type, hist_var, hist_var_hmm, taxid, highlvl_taxid)
    data = []
    for seq in Sequence.objects.all():
        # print(seq.id)
        if seq.variant is None: continue
        try:
            lineage = NCBI_TAXA.get_lineage(seq.taxonomy_id)
            intersection = set(TAXIDS).intersection(set(lineage))
            if len(intersection) == 0:
                highlvl_taxid = ''
            elif len(intersection) == 1:
                highlvl_taxid = intersection.pop()
        except ValueError as e:
            highlvl_taxid = ''
        data.append([seq.id, seq.variant.hist_type, seq.variant, seq.variant_hmm, seq.taxonomy_id, highlvl_taxid, seq.reviewed])
    data = pd.DataFrame(np.array(data),
                        columns=['accession', 'hist_type', 'hist_var', 'hist_var_hmm', 'taxid', 'highlvl_taxid', 'curated'])
    data.to_csv(SEQS_FILE)
    return filter_seqs_data(data, curated, hist_type, hist_var, hist_var_hmm, taxid, highlvl_taxid)

def filter_seqs_data(data, curated=None, hist_type=None, hist_var=None, hist_var_hmm=None, taxid=None, highlvl_taxid=None):
    data = data[data['curated']] if curated else data[data['curated'] == False] if curated == False else data
    data = data[data['hist_type'] == hist_type] if hist_type else data
    data = data[data['hist_var'] == hist_var] if hist_var else data
    data = data[data['hist_var_hmm'] == hist_var_hmm] if hist_var_hmm else data
    data = data[data['taxid'] == taxid] if taxid else data
    data = data[data['highlvl_taxid'] == highlvl_taxid] if highlvl_taxid else data
    return data

def get_scores_data(rewrite = False, curated=None, hist_type=None, blast_model=None,
                       score_1=None, score_2=None, hsp_length_1=None, hsp_length_2=None,
                       seq_taxid=None, seq_highlvl_taxid=None):
    if Path(SCORES_FILE).exists() and not rewrite:
        # return pd.read_csv(SCORES_FILE)
        return filter_scores_data(pd.read_csv(SCORES_FILE), curated, hist_type, blast_model,
                                  score_1, score_2, hsp_length_1, hsp_length_2,
                                  seq_taxid, seq_highlvl_taxid)
    data = []
    for s in ScoreBlast.objects.all():
        print(s.id)
        if s.hit_accession == '':
            hit_seq = ''
        else:
            hit_seq = Sequence.objects.get(id=s.hit_accession).sequence
        try:
            lineage = NCBI_TAXA.get_lineage(s.sequence.taxonomy_id)
            intersection = set(TAXIDS).intersection(set(lineage))
            if len(intersection) == 0:
                highlvl_taxid = ''
            elif len(intersection) == 1:
                highlvl_taxid = intersection.pop()
        except ValueError as e:
            highlvl_taxid = ''
        data.append([s.sequence.id, s.variant.hist_type, s.variant,
                     s.score, s.bitScore, s.evalue, s.align_length, s.used_for_classification,
                     s.hit_accession, s.sequence.sequence, hit_seq, s.match,
                     s.blastStart, s.blastEnd, s.seqStart, s.seqEnd,
                     s.sequence.taxonomy_id, highlvl_taxid, s.sequence.reviewed])
    data = pd.DataFrame(np.array(data),
                        columns=['accession', 'hist_type', 'blast_model',
                                 'score', 'bit_score', 'evalue', 'hsp_length', 'used_for_classification',
                                 'hit_accession', 'sequence', 'hit_sequence', 'match',
                                 'blastStart', 'blastEnd', 'seqStart', 'seqEnd',
                                 'seq_taxid', 'seq_highlvl_taxid', 'curated'])
    data.to_csv(SCORES_FILE)
    # return data
    return filter_scores_data(pd.read_csv(SCORES_FILE), curated, hist_type, blast_model,
                              score_1, score_2, hsp_length_1, hsp_length_2,
                              seq_taxid, seq_highlvl_taxid)

def filter_scores_data(data, curated=None, hist_type=None, blast_model=None,
                       score_1=None, score_2=None, hsp_length_1=None, hsp_length_2=None,
                       seq_taxid=None, seq_highlvl_taxid=None):
    data = data[data['curated']] if curated else data[data['curated'] == False] if curated == False else data
    data = data[data['hist_type'] == hist_type] if hist_type else data
    data = data[data['blast_model'] == blast_model] if blast_model else data
    data = data[data['score'] > score_1] if score_1 else data
    data = data[data['score'] <= score_2] if score_2 else data
    data = data[data['hsp_length'] > hsp_length_1] if hsp_length_1 else data
    data = data[data['hsp_length'] <= hsp_length_2] if hsp_length_2 else data
    data = data[data['seq_taxid'] == taxid] if seq_taxid else data
    data = data[data['highlvl_taxid'] == seq_highlvl_taxid] if seq_highlvl_taxid else data
    return data

def general_statistics():
    with open(os.path.join(CURR_STAT_DIR, "general_stat"), 'w') as f:
        f.write('---Database statistics----\n')
        f.write('Total seqs = %d\n' % Sequence.objects.all().count())
        f.write('Reviewed seqs = %d\n' % Sequence.objects.filter(reviewed=True).count())
        f.write('Automatic seqs = %d\n' % Sequence.objects.filter(reviewed=False).count())
        f.write('\n---Histone type statistics----\n')
        f.write('Type        | Total  |Reviewed|  Auto  \n')
        for h in Histone.objects.all():
            tot = Sequence.objects.filter(variant__hist_type=h).count()
            rev = Sequence.objects.filter(variant__hist_type=h, reviewed=True).count()
            auto = Sequence.objects.filter(variant__hist_type=h, reviewed=False).count()
            f.write('%12s|%8d|%8d|%8d\n' % (h.id, tot, rev, auto))
        f.write('\n---Histone variant statistics----\n')
        f.write('Variant     | Total  |Reviewed|  Auto  \n')
        for v in Variant.objects.all():
            tot = Sequence.objects.filter(variant=v).count()
            rev = Sequence.objects.filter(variant=v, reviewed=True).count()
            auto = Sequence.objects.filter(variant=v, reviewed=False).count()
            f.write('%12s|%8d|%8d|%8d\n' % (v.id, tot, rev, auto))
        # f.write('\n---Matching with HMM algorithm----\n')
        # f.write('Variant     | Total  |Reviewed|  Auto  \n')
    pp = PdfPages(os.path.join(CURR_STAT_DIR, "Variants_hist.pdf"))
    data = get_seqs_data(curated=False)
    plt.figure(figsize=(15, 16))
    ax = sns.countplot(y='hist_var', data=data)
    for p in ax.patches:
        x = p.get_bbox().get_points()[1, 0]
        y = p.get_bbox().get_points()[:, 1]
        ax.annotate('{}'.format(int(x)), (x + 1, y.mean()), ha='left', va='center')
    plt.title('Distribution of variants')
    pp.savefig()
    pp.close()
    pp = PdfPages(os.path.join(CURR_STAT_DIR, "Variants_hist_curated.pdf"))
    data = get_seqs_data(curated=True)
    plt.figure(figsize=(15, 16))
    ax = sns.countplot(y='hist_var', data=data)
    for p in ax.patches:
        x = p.get_bbox().get_points()[1, 0]
        y = p.get_bbox().get_points()[:, 1]
        ax.annotate('{}'.format(int(x)), (x + 1, y.mean()), ha='left', va='center')
    plt.title('Distribution of curated variants')
    pp.savefig()
    pp.close()

def var_distrib_within_histtypes():
    for histone_type in ['H1', 'H2A', 'H2B', 'H3', 'H4']:
        data = get_seqs_data(curated=False, hist_type=histone_type)
        pp = PdfPages(os.path.join(os.path.join(CURR_STAT_DIR, histone_type), "Variants_distribution_{}.pdf".format(histone_type)))
        plt.figure(figsize=(15, 16))
        ax = sns.countplot(y='hist_var', data=data)
        for p in ax.patches:
            x = p.get_bbox().get_points()[1, 0]
            y = p.get_bbox().get_points()[:, 1]
            ax.annotate('{}'.format(int(x)), (x + 1, y.mean()), ha='left', va='center')
        plt.title('Distribution of variants within {}'.format(histone_type))
        pp.savefig()
        pp.close()

def var_distrib_for_highlvltaxa():
    for hlvl_taxa in HIGHLVL_TAXA:
        data = get_seqs_data(curated=False, highlvl_taxid=TAXIDS_DICT[hlvl_taxa][0])
        pp = PdfPages(os.path.join(os.path.join(CURR_STAT_DIR, hlvl_taxa), "Variants_distribution_{}.pdf".format(hlvl_taxa)))
        plt.figure(figsize=(15, 16))
        ax = sns.countplot(y='hist_var', data=data)
        for p in ax.patches:
            x = p.get_bbox().get_points()[1, 0]
            y = p.get_bbox().get_points()[:, 1]
            ax.annotate('{}'.format(int(x)), (x + 1, y.mean()), ha='left', va='center')
        plt.title('Distribution of variants for {}'.format(hlvl_taxa))
        pp.savefig()
        pp.close()

def var_distrib_within_histtypes_for_highlvltaxa():
    for histone_type in ['H1', 'H2A', 'H2B', 'H3', 'H4']:
        for hlvl_taxa in HIGHLVL_TAXA:
            data = get_seqs_data(curated=False, hist_type=histone_type, highlvl_taxid=TAXIDS_DICT[hlvl_taxa][0])
            pp = PdfPages(os.path.join(os.path.join(CURR_STAT_DIR, histone_type, hlvl_taxa),
                                       "Variants_distribution_{}_{}.pdf".format(histone_type, hlvl_taxa)))
            plt.figure(figsize=(15, 16))
            ax = sns.countplot(y='hist_var', data=data)
            for p in ax.patches:
                x = p.get_bbox().get_points()[1, 0]
                y = p.get_bbox().get_points()[:, 1]
                ax.annotate('{}'.format(int(x)), (x + 1, y.mean()), ha='left', va='center')
            plt.title('Distribution of variants within {} for {}'.format(histone_type, hlvl_taxa))
            pp.savefig()
            pp.close()

def scores_boxplot():
    data = get_scores_data(curated=False)
    pp = PdfPages(
        os.path.join(os.path.join(CURR_STAT_DIR), "Scores_boxplot.pdf"))
    plt.figure(figsize=(20, 5))
    ax = sns.boxplot(data['score'])
    for p in ax.patches:
        x = p.get_bbox().get_points()[1, 0]
        y = p.get_bbox().get_points()[:, 1]
        ax.annotate('{}'.format(int(x)), (x + 1, y.mean()), ha='left', va='center')
    plt.title('Alignment scores')
    pp.savefig()
    pp.close()

def scores_boxplot_within_variants():
    variants = Variant.objects.all()
    pp = PdfPages(
        os.path.join(os.path.join(CURR_STAT_DIR), "Scores_boxplot_within_variants.pdf"))
    f, axes = plt.subplots(len(variants), figsize=(10, 5*len(variants)))
    for i, var in enumerate(variants):
        data = get_scores_data(curated=False, blast_model=var.id)
        sns.boxplot(data['score'], ax=axes[i])
        axes[i].set_title('Alignment scores for {}'.format(var.id))
    pp.savefig()
    pp.close()

def curated_scores_boxplot_within_variants():
    variants = Variant.objects.all()
    pp = PdfPages(
        os.path.join(os.path.join(CURR_STAT_DIR), "Curated_scores_boxplot_within_variants.pdf"))
    f, axes = plt.subplots(len(variants), figsize=(10, 5*len(variants)))
    for i, var in enumerate(variants):
        data = get_scores_data(curated=True, blast_model=var.id)
        sns.boxplot(data['score'], ax=axes[i])
        axes[i].set_title('Alignment scores for curated {}'.format(var.id))
    pp.savefig()
    pp.close()

def compare_scores_boxplot_within_variants():
    variants = Variant.objects.all()
    pp = PdfPages(
        os.path.join(os.path.join(CURR_STAT_DIR), "Compare_scores_boxplot_within_variants.pdf"))
    f, axes = plt.subplots(len(variants), figsize=(10, 5*len(variants)))
    for i, var in enumerate(variants):
        data = get_scores_data(curated=False, blast_model=var.id)
        data_curated = get_scores_data(curated=True, blast_model=var.id)
        sns.boxplot(data['score'], ax=axes[i], color=sns.xkcd_rgb["denim blue"])
        sns.boxplot(data_curated['score'], ax=axes[i], color=sns.xkcd_rgb["pale red"])
        # sns.boxplot(data['score'], ax=axes[i], color=sns.xkcd_rgb["denim blue"], label='Automatically extracted sequences')
        # sns.boxplot(data_curated['score'], ax=axes[i], color=sns.xkcd_rgb["pale red"], label='Curated sequences')
        axes[i].legend(loc="upper left")
        axes[i].set_title('Alignment scores for {}'.format(var.id))
    pp.savefig()
    pp.close()

# create_directories()
# general_statistics()
# var_distrib_within_histtypes()
# var_distrib_for_highlvltaxa()
# var_distrib_within_histtypes_for_highlvltaxa()
# scores_boxplot()
# scores_boxplot_within_variants()
# curated_scores_boxplot_within_variants()
compare_scores_boxplot_within_variants()