from browse.models import Sequence, Score, ScoreBlast, SequenceBlast, Histone, Variant
from django.conf import settings
# from django.db.models import Q

from tools.stat_taxa import *

#This script is used to export tables from database for futher use in research

import os
from datetime import date, datetime
from pathlib import Path

import numpy as np
import pandas as pd
import math

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

NOW = datetime.now()
DT_STRING = NOW.strftime("%Y%m%d-%H%M%S")

def get_nr_version():
    with open('NR_VERSION', 'r') as nrv:
        return nrv.read()

STAT_DIR = os.path.join(settings.STATIC_ROOT_AUX, "browse", "statistics")
# CURR_STAT_DIR = os.path.join(STAT_DIR, '{}_{}'.format(get_nr_version(), DT_STRING))
CURR_STAT_DIR = os.path.join(STAT_DIR, 'nr_small_per10_v4_20210622-152129')
# CURR_STAT_DIR = os.path.join(STAT_DIR, DT_STRING)
# CURR_STAT_DIR = os.path.join(STAT_DIR, '20201219-195135')
# CURR_STAT_DIR = os.path.join(STAT_DIR, 'nr_small_per10_v4_20210126-183517')
DATA_DIR = os.path.join(CURR_STAT_DIR, 'data')

SEQS_FILE = os.path.join(DATA_DIR, 'seqs.csv')
SCORES_FILE = os.path.join(DATA_DIR, 'scores.csv')
SCORES_HMM_FILE = os.path.join(DATA_DIR, 'scores_hmm.csv')

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
        data.append([seq.id, seq.variant.hist_type.id, seq.variant.id, seq.variant_hmm.id, seq.taxonomy_id, highlvl_taxid, seq.reviewed])
    data = pd.DataFrame(data,
                        columns=['accession', 'hist_type', 'hist_var', 'hist_var_hmm', 'taxid', 'highlvl_taxid', 'curated'])
    data.to_csv(SEQS_FILE, index=False)
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
        data.append([s.sequence.id, s.variant.hist_type.id, s.variant.id,
                     s.score, s.bitScore, s.evalue, s.align_length, s.used_for_classification,
                     s.hit_accession, s.sequence.sequence, hit_seq, s.match,
                     s.blastStart, s.blastEnd, s.seqStart, s.seqEnd,
                     s.sequence.taxonomy_id, highlvl_taxid, s.sequence.reviewed])
    data = pd.DataFrame(data,
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
    ax = sns.countplot(y='hist_var', data=data, order = data['hist_var'].value_counts().index)
    print(ax.get_yticklabels())
    print(ax.get_legend_handles_labels())
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
    ax = sns.countplot(y='hist_var', data=data, order = data['hist_var'].value_counts().index)
    for p in ax.patches:
        x = p.get_bbox().get_points()[1, 0]
        y = p.get_bbox().get_points()[:, 1]
        ax.annotate('{}'.format(int(x)), (x + 1, y.mean()), ha='left', va='center')
    plt.title('Distribution of curated variants')
    pp.savefig()
    pp.close()

def general_statistics_pickle():
    import pickle
    import mpld3
    from mpld3 import plugins
    data = get_seqs_data(curated=False)
    fig = plt.figure(figsize=(17, 8))
    ax = sns.countplot(y='hist_var', data=data, order = data['hist_var'].value_counts().index)
    for p in ax.patches:
        x = p.get_bbox().get_points()[1, 0]
        y = p.get_bbox().get_points()[:, 1]
        ax.annotate('{}'.format(int(x)), (x + 1, y.mean()), ha='left', va='center')
        ax.annotate('{}'.format(int(x)), (x + 1, y.mean()), ha='right', va='center')
    # handles, labels = ax.get_legend_handles_labels()  # return lines and labels
    # interactive_legend = plugins.InteractiveLegendPlugin(zip(handles,
    #                                                          ax.collections),
    #                                                      labels,
    #                                                      alpha_unsel=0.5,
    #                                                      alpha_over=1.5,
    #                                                      start_visible=True)
    # plugins.connect(fig, interactive_legend)
    ax.set_xlabel('Count')
    ax.set_ylabel('Histone variant')
    ax.set_title('Distribution of variants', size=20)
    html_fig = mpld3.fig_to_html(fig, template_type='general')
    plt.close(fig)
    with open(os.path.join(CURR_STAT_DIR, 'Variants_hist.pickle'), 'wb') as f:
        pickle.dump(html_fig, f)
    data = get_seqs_data(curated=True)
    fig = plt.figure(figsize=(15, 16))
    ax = sns.countplot(y='hist_var', data=data, order = data['hist_var'].value_counts().index)
    for p in ax.patches:
        x = p.get_bbox().get_points()[1, 0]
        y = p.get_bbox().get_points()[:, 1]
        ax.annotate('{}'.format(int(x)), (x + 1, y.mean()), ha='left', va='center')
    handles, labels = ax.get_legend_handles_labels()  # return lines and labels
    interactive_legend = plugins.InteractiveLegendPlugin(zip(handles,
                                                             ax.collections),
                                                         labels,
                                                         alpha_unsel=0.5,
                                                         alpha_over=1.5,
                                                         start_visible=True)
    plugins.connect(fig, interactive_legend)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('Distribution of curated variants', size=20)
    html_fig = mpld3.fig_to_html(fig, template_type='general')
    plt.close(fig)
    with open(os.path.join(CURR_STAT_DIR, 'Variants_hist_curated.pickle'), 'wb') as f:
        pickle.dump(html_fig, f)

def general_statistics_html():
    import mpld3
    from mpld3 import plugins
    data = get_seqs_data(curated=False)
    fig = plt.figure(figsize=(10, 8))
    ax = sns.countplot(y='hist_var', data=data, order = data['hist_var'].value_counts().index)
    # ax.set_yticklabels(ax.get_yticklabels(), fontsize=7)
    # print(ax.get_yticklabels())
    # print(ax.get_legend_handles_labels())
    for p in ax.patches:
        x = p.get_bbox().get_points()[1, 0]
        y = p.get_bbox().get_points()[:, 1]
        ax.annotate('{}'.format(int(x)), (x + 1, y.mean()), ha='left', va='center')
    # handles, labels = ax.get_legend_handles_labels()  # return lines and labels
    # interactive_legend = plugins.InteractiveLegendPlugin(zip(handles,
    #                                                          ax.collections),
    #                                                      labels,
    #                                                      alpha_unsel=0.5,
    #                                                      alpha_over=1.5,
    #                                                      start_visible=True)
    # plugins.connect(fig, interactive_legend)
    ax.set_xlabel('Count')
    ax.set_ylabel('Histone variant')
    ax.set_title('Distribution of variants', size=20)
    # print(ax.get_legend_handles_labels())
    # tooltip = mpld3.plugins.PointLabelTooltip(ax.get_yticklabels(), labels=data['hist_var'].value_counts().index)
    # mpld3.plugins.connect(fig, tooltip)
    # html_fig = mpld3.fig_to_html(fig, template_type='general')
    print(mpld3.fig_to_dict(fig))
    mpld3.save_html(fig,os.path.join('templates', 'var_hist.html'), template_type='general')
    plt.close(fig)
    # with open(os.path.join(CURR_STAT_DIR, 'Variants_hist.pickle'), 'wb') as f:
    #     pickle.dump(html_fig, f)
    # data = get_seqs_data(curated=True)
    # fig = plt.figure(figsize=(15, 16))
    # ax = sns.countplot(y='hist_var', data=data, order = data['hist_var'].value_counts().index)
    # for p in ax.patches:
    #     x = p.get_bbox().get_points()[1, 0]
    #     y = p.get_bbox().get_points()[:, 1]
    #     ax.annotate('{}'.format(int(x)), (x + 1, y.mean()), ha='left', va='center')
    # handles, labels = ax.get_legend_handles_labels()  # return lines and labels
    # interactive_legend = plugins.InteractiveLegendPlugin(zip(handles,
    #                                                          ax.collections),
    #                                                      labels,
    #                                                      alpha_unsel=0.5,
    #                                                      alpha_over=1.5,
    #                                                      start_visible=True)
    # plugins.connect(fig, interactive_legend)
    # ax.set_xlabel('x')
    # ax.set_ylabel('y')
    # ax.set_title('Distribution of curated variants', size=20)
    # # html_fig = mpld3.fig_to_html(fig, template_type='general')
    # mpld3.save_html(fig,os.path.join('templates', 'Variants_hist_curated.html'), template_type='general')
    # plt.close(fig)
    # # with open(os.path.join(CURR_STAT_DIR, 'Variants_hist_curated.pickle'), 'wb') as f:
    # #     pickle.dump(html_fig, f)

def general_statistics_pickle_test():
    import pickle
    import mpld3
    from mpld3 import plugins
    np.random.seed(9615)
    N = 100
    df = pd.DataFrame((.1 * (np.random.random((N, 5)) - .5)).cumsum(0),
                      columns=['a', 'b', 'c', 'd', 'e'], )
    # plot line + confidence interval
    fig, ax = plt.subplots()
    ax.grid(True, alpha=0.3)
    for key, val in df.iteritems():
        l, = ax.plot(val.index, val.values, label=key)
        ax.fill_between(val.index,
                        val.values * .5, val.values * 1.5,
                        color=l.get_color(), alpha=.4)
    handles, labels = ax.get_legend_handles_labels()  # return lines and labels
    interactive_legend = plugins.InteractiveLegendPlugin(zip(handles,
                                                             ax.collections),
                                                         labels,
                                                         alpha_unsel=0.5,
                                                         alpha_over=1.5,
                                                         start_visible=True)
    plugins.connect(fig, interactive_legend)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('Distribution of variants', size=20)
    html_fig = mpld3.fig_to_html(fig, template_type='general')
    plt.close(fig)
    # with open(os.path.join(CURR_STAT_DIR, 'test_fig.pickle'), 'wb') as f:
    #     pickle.dump(fig, f)
    # with open(os.path.join(CURR_STAT_DIR, 'test_ax.pickle'), 'wb') as f:
    #     pickle.dump(ax, f)
    with open(os.path.join(CURR_STAT_DIR, 'test_html_fig.pickle'), 'wb') as f:
        pickle.dump(html_fig, f)

def general_statistics_2():
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
    data = get_seqs_data()
    plt.figure(figsize=(30, 40))
    ax = sns.countplot(y='hist_var', hue='curated', data=data, order = data['hist_var'].value_counts().index)
    for p in ax.patches:
        x = p.get_bbox().get_points()[1, 0]
        y = p.get_bbox().get_points()[:, 1]
        # print(x)
        # print(y)
        # print(p.get_bbox())
        import math
        if math.isnan(x):
            # print(x)
            x = 1.0
        ax.annotate('{}'.format(int(x)), (x + 1, y.mean()), ha='left', va='center')
    plt.title('Distribution of variants')
    pp.savefig()
    pp.close()

def var_distrib_within_histtypes(filename='Variants_distribution'):
    for histone_type in ['H1', 'H2A', 'H2B', 'H3', 'H4']:
        # data = get_seqs_data(curated=False, hist_type=histone_type)
        data = get_seqs_data(hist_type=histone_type)
        pp = PdfPages(os.path.join(os.path.join(CURR_STAT_DIR, histone_type), "{}_{}.pdf".format(filename, histone_type)))
        plt.figure(figsize=(15, 5))
        ax = sns.countplot(y='hist_var', hue='curated', data=data, order = data['hist_var'].value_counts().index)
        for p in ax.patches:
            x = p.get_bbox().get_points()[1, 0]
            y = p.get_bbox().get_points()[:, 1]
            if math.isnan(x): x = .0
            ax.annotate('{}'.format(int(x)), (x + 1, y.mean()), ha='left', va='center')
        plt.title('Distribution of variants within {}'.format(histone_type))
        pp.savefig()
        pp.close()

def var_distrib_for_highlvltaxa(filename='Variants_distribution'):
    for hlvl_taxa in HIGHLVL_TAXA:
        # data = get_seqs_data(curated=False, highlvl_taxid=TAXIDS_DICT[hlvl_taxa][0])
        data = get_seqs_data(highlvl_taxid=TAXIDS_DICT[hlvl_taxa][0])
        pp = PdfPages(os.path.join(os.path.join(CURR_STAT_DIR, hlvl_taxa), "{}_{}.pdf".format(filename, hlvl_taxa)))
        plt.figure(figsize=(25, 20))
        ax = sns.countplot(y='hist_var', hue='curated', data=data, linewidth=0, order = data['hist_var'].value_counts().index)
        for p in ax.patches:
            x = p.get_bbox().get_points()[1, 0]
            y = p.get_bbox().get_points()[:, 1]
            if math.isnan(x): x = .0
            ax.annotate('{}'.format(int(x)), (x + 1, y.mean()), ha='left', va='center')
        plt.title('Distribution of variants for {}'.format(hlvl_taxa))
        pp.savefig()
        pp.close()

def var_distrib_within_histtypes_for_highlvltaxa(filename='Variants_distribution'):
    for histone_type in ['H1', 'H2A', 'H2B', 'H3', 'H4']:
        for hlvl_taxa in HIGHLVL_TAXA:
            # data = get_seqs_data(curated=False, hist_type=histone_type, highlvl_taxid=TAXIDS_DICT[hlvl_taxa][0])
            data = get_seqs_data(hist_type=histone_type, highlvl_taxid=TAXIDS_DICT[hlvl_taxa][0])
            if data.empty: continue
            pp = PdfPages(os.path.join(os.path.join(CURR_STAT_DIR, histone_type, hlvl_taxa),
                                       "{}_{}_{}.pdf".format(filename, histone_type, hlvl_taxa)))
            plt.figure(figsize=(15, 16))
            ax = sns.countplot(y='hist_var', hue='curated', data=data, linewidth=0, order = data['hist_var'].value_counts().index)
            for p in ax.patches:
                x = p.get_bbox().get_points()[1, 0]
                y = p.get_bbox().get_points()[:, 1]
                if math.isnan(x): x = .0
                ax.annotate('{}'.format(int(x)), (x + 1, y.mean()), ha='left', va='center')
            plt.title('Distribution of variants within {} for {}'.format(histone_type, hlvl_taxa))
            pp.savefig()
            pp.close()

def scores_boxplot(filename='Scores_boxplot'):
    variants = Variant.objects.all()
    pp = PdfPages(
        os.path.join(os.path.join(CURR_STAT_DIR), filename+".pdf"))
    plt.figure(figsize=(25, 50))
    data = get_scores_data()
    sns.boxplot(x='score', y='blast_model', hue='curated', data=data)
    plt.title('Alignment scores')
    pp.savefig()
    pp.close()

# create_directories()
# general_statistics()
# var_distrib_within_histtypes()
# var_distrib_for_highlvltaxa()
# var_distrib_within_histtypes_for_highlvltaxa()
# scores_boxplot()

# general_statistics_pickle()
general_statistics_html()
# general_statistics_pickle_test()