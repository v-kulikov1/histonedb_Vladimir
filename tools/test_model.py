#!/usr/local/bin/python
# Author: Eli Draizen
# Date: 16-3-2014
# File: classify.py

import matplotlib
matplotlib.use('Agg')
 
#Standard Libraries
import argparse
import json
import os
import logging
from pathlib import Path

#BioPython
from Bio import SearchIO
from Bio import SeqIO
from Bio.Seq import Seq

#Required Libraries
import numpy as np
import seaborn as sns
import pandas as pd

import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, matthews_corrcoef
from matplotlib.backends.backend_pdf import PdfPages

from scipy.interpolate import interp1d

from browse.models import Sequence, ScoreBlast, Histone, Variant
from django.conf import settings
# from tools.stat_taxa import *

sns.set(style="white", context="talk")

# HMMsearch

#Thresholds decided by manualy looking at the headers in the hmmsearch output. Useful to compare.
original_thresholds = {
    "H2A": {
        "H2A.X":275, 
        "H2A.Z":230, 
        "H2A.B":115, 
        "H2A.L":145, 
        "H2A.M":80, 
        "macroH2A":270, 
        "canonical_H2A":265
    }
}

def get_model_scores(model_output):
    """Get the bit score for each hit/domain in a hmmersearch result

    Parameters:
    -----------
    model_output: str or File-like object
        Path to hmmersearch output file

    Return:
    -------
    A list of all bitscores
    """
    return [hsp.bitscore for query in SearchIO.parse(model_output, "hmmer3-text") \
        for hit in query for hsp in hit] 

def test_model(model_name, save_dir, postive_file, negative_file, measure="SPC", measure_threshold=0.95):
    """Test the model by calcuating

    Returns:
    A dictionary with containg the AUCROC and Threshold. An image is also saved
    with the ROC and score histograms. 
    """
    logging.info(model_name)
    postive_scores = get_model_scores(postive_file)
    negative_scores = get_model_scores(negative_file)
    all_scores = postive_scores+negative_scores
    # print all_scores

    if len(negative_scores) == 0:
        return {"roc_auc":0, "threshold":min(postive_scores)}

    y_true = [1]*len(postive_scores) + [0]*len(negative_scores)
    y_score = np.array(all_scores)

    fpr, tpr, thresholds = roc_curve(y_true, y_score)
    roc_auc = auc(fpr, tpr)

    best_threshold, thresholds, values = calcualte_threshold(
        postive_scores, 
        negative_scores, 
        measure=measure,
        measure_threshold=measure_threshold, 
        thresholds=reversed(thresholds))


    pp = PdfPages(os.path.join(save_dir, "{}_model_evaluation.pdf".format(model_name)))

    sns.set(style="darkgrid")
    f, axes = plt.subplots(3)
    trans = f.transFigure.inverted()
    colors = sns.color_palette("Set2", 7)

    try:
        sns.kdeplot(np.array(postive_scores), shade=True, color=sns.xkcd_rgb["denim blue"], label="Scores for postive examples", ax=axes[0])
        sns.kdeplot(np.array(negative_scores), shade=True, color=sns.xkcd_rgb["pale red"], label="Scores for negative examples",  ax=axes[0])
    except:
        pass
        # axes[0].hist(postive_scores, color=sns.xkcd_rgb["denim blue"], label="Scores for postive examples")
        # axes[0].hist(negative_scores, color=sns.xkcd_rgb["pale red"], label="Scores for negative examples")
    # except np.linalg.LinAlgError as err:
    #     if 'singular matrix' in str(err):
    #         axes[0].hist(postive_scores, color=sns.xkcd_rgb["denim blue"], label="Scores for postive examples")
    #         axes[0].hist(negative_scores, color=sns.xkcd_rgb["pale red"], label="Scores for negative examples")
    #     else:
    #         raise
    axes[0].set_xlabel("Bit score")
    axes[0].set_ylabel("Density")
    axes[0].legend(loc="upper left")
    #axes[0].set_title("Kernel Density of Scores")
    axes[1].set_xlim([0, 1.0])
    axes[1].set_ylim([0.0, 1.05])

    
    axes[1].plot(fpr,tpr, color=colors[0], lw=3., label="ROC (AUC: {})".format(roc_auc))
    axes[1].set_xlabel("False Positive Rate")
    axes[1].set_ylabel("True Positive Rate")
    axes[1].legend(loc="lower right")
    axes[1].set_xlim([-0.05, 1.0])
    axes[1].set_ylim([0.0, 1.05])
    #axes[1].set_title("ROC")
 
    for i, (measure, values) in enumerate(values.items()):
        label = "SPC: (>={})".format(best_threshold) if measure=="SPC" else measure
        axes[2].plot(list(thresholds), values, label=label, linewidth=2, color=colors[i])
    axes[2].axvline(best_threshold)

    axes[2].legend()
    #axes[2].set_title("Coosing Cutoff")
    axes[2].set_ylabel("Rate")
    axes[2].set_xlabel("Threshold")

    f.suptitle("{} Model Evaluation".format(model_name), fontsize=20)

    pp.savefig()
    pp.close()

    return {"roc_auc":roc_auc, "threshold":best_threshold}

def calcualte_threshold(positives, negatives, measure="SPC", measure_threshold=0.95, thresholds=None, attempt=0):
    """Plot the TPR the FPR vs threshold values
 
    Input:
    postives - list of scores of postive runs
    negatives - list of scores of negative runs
    measure - choose coffectiong by 95% Specificity ("SPC"), or matthews_corrcoef ("MCC")
    """
    assert measure in ["TPR", "FPR", "SPC", "MCC", "PPV", "FDR", "ACC"]
    y_true = [1]*len(positives)+[0]*len(negatives)
    values = {name:[] for name in ["TPR", "FPR", "SPC", "MCC", "PPV", "FDR", "ACC"]}
    saveThreshold = None
    saveValue = 1.0
    thresholds = list(sorted(thresholds or [i/10. for i in range(1,10000)]))

    for threshold in thresholds:
        TN = sum([1 for score in negatives if score < threshold])
        FP = sum([1 for score in negatives if score >= threshold])
        TP = sum([1 for score in positives if score >= threshold])
        FN = sum([1 for score in positives if score < threshold])

        values["FPR"].append(float(FP)/(FP+TN))
        values["TPR"].append(float(TP)/(TP+FN))
        values["SPC"].append(float(TN)/(FP+TN))

        y_pred = [int(score >= threshold) for scores in (positives, negatives) for score in scores]
        values["MCC"].append(matthews_corrcoef(y_true, y_pred))
        values["PPV"].append(float(TP)/(TP+FP) if TP+FP>0 else 0.0)
        values["FDR"].append(float(FP)/(TP+FP) if TP+FP>0 else 0.0)
        values["ACC"].append(float(TP+TN)/(len(positives)+len(negatives)))
        
    specificity_curve_inverse = interp1d(values[measure], thresholds)
    saveThreshold = specificity_curve_inverse(measure_threshold) #modified by Alexey 0.95 => measure_threshold
                
    return saveThreshold, thresholds, values


# BLASTsearch

def scores_model_to_dataframe(hist_type=None):
    """Prepares scores of curated sequences for statistical analysis, returns a DataFrame according to filter parameters"""

    if Path(os.path.join(settings.STATIC_ROOT_AUX, "browse", "blast", "model_evaluation", "score.csv")).exists():
        return filter_scores_data(pd.read_csv(os.path.join(settings.STATIC_ROOT_AUX, "browse", "blast", "model_evaluation", "score.csv")), hist_type)

    data = []
    # for s in ScoreBlast.objects.filter(sequence__reviewed=True).exclude(sequence__variant__startswith='generic'):
    for s in ScoreBlast.objects.filter(sequence__reviewed=True):
        if 'generic' in s.variant.id: continue

        # hit_seq = Sequence.objects.get(id=s.hit_accession).sequence
        # hit_var = Sequence.objects.get(id=s.hit_accession).variant
        hit = Sequence.objects.get(id=s.hit_accession)
        is_positive = 'True positive' if s.variant.id == hit.variant.id else 'False positive'
        relative_length = s.align_length/len(hit.sequence)
        # try:
        #     lineage = NCBI_TAXA.get_lineage(s.sequence.taxonomy_id)
        #     intersection = set(TAXIDS).intersection(set(lineage))
        #     if len(intersection) == 0:
        #         highlvl_taxid = ''
        #     elif len(intersection) == 1:
        #         highlvl_taxid = intersection.pop()
        # except ValueError as e:
        #     highlvl_taxid = ''
        data.append([s.sequence.id, s.variant.hist_type, s.variant, s.score, s.evalue,
                     s.align_length, relative_length,
                     s.used_for_classification,
                     s.hit_accession, hit.variant, is_positive,
                     s.sequence.sequence, hit.sequence,
                     s.blastStart, s.blastEnd, s.seqStart, s.seqEnd])
    data = pd.DataFrame(np.array(data),
                        columns=['accession', 'hist_type', 'variant', 'score', 'evalue',
                                 'hsp_length', 'relative_length',
                                 'used_for_classification',
                                 'hit_accession', 'expected_variant', 'is_positive',
                                 'sequence', 'hit_sequence',
                                 'query_start', 'query_end', 'sbjct_start', 'sbjct_end'])
    data.to_csv(os.path.join(settings.STATIC_ROOT_AUX, "browse", "blast", "model_evaluation", "score.csv"))

    return filter_scores_data(data, hist_type)

def filter_scores_data(data, hist_type=None):
    data = data[data['hist_type'] == hist_type] if hist_type else data
    return data

def test_curated(save_dir):
    """Test the model by visualizing score data"""

    visualize_general_stat(save_dir)
    # visualize_variants_stat(save_dir)

def visualize_general_stat(save_dir):
    """Visualize all scores, hsp_lengths and relative_lengths via violinplots"""

    hsps_data = scores_model_to_dataframe()
    hsps_data['score'] = hsps_data['score'].astype('float64')

    pp = PdfPages(os.path.join(save_dir, "General_statistics.pdf"))

    # Scores
    sns.set(style="darkgrid")
    colors = sns.color_palette("Set2", 7)
    ax = sns.catplot(x="score", y="variant", hue="is_positive",
                     data=hsps_data, orient="h", palette="Set2", split=True,
                     kind="violin")
    ax.fig.set_size_inches(7, 10)
    plt.title('Distribution of scores')
    pp.savefig()
    # plt.close()

    # HSP_length
    sns.set(style="darkgrid")
    colors = sns.color_palette("Set2", 7)
    ax = sns.catplot(x="hsp_length", y="variant", hue="is_positive",
                     data=hsps_data, orient="h", palette="Set2", split=True,
                     kind="violin")
    ax.fig.set_size_inches(7, 10)
    plt.title('Distribution of HSP length')
    pp.savefig()

    # Relative_length
    sns.set(style="darkgrid")
    colors = sns.color_palette("Set2", 7)
    ax = sns.catplot(x="relative_length", y="variant", hue="is_positive",
                     data=hsps_data, orient="h", palette="Set2", split=True,
                     kind="violin")
    ax.fig.set_size_inches(7, 10)
    plt.title('Distribution of relative length')
    pp.savefig()

    pp.close()

def visualize_variants_stat(save_dir):
    """Visualize scores, hsp_lengths and relative_lengths for each histone variant"""
    for histone_type in ['H1', 'H2A', 'H2B', 'H3', 'H4']:
        hsps_data = scores_model_to_dataframe(hist_type=histone_type)
        hsps_data['score'] = hsps_data['score'].astype('float64')

        for variant in set(hsps_data['variant']):
            print('{}'.format(variant))
            pp = PdfPages(os.path.join(os.path.join(save_dir, histone_type), "{}.pdf".format(variant)))

            # Scores
            sns.set(style="darkgrid")
            colors = sns.color_palette("Set2", 7)
            ax = sns.catplot(x="score", y="accession", hue="is_positive", row="variant",
                             data=hsps_data[hsps_data['variant']==variant], orient="h", palette="Set2", split=True,
                             kind="violin")
            ax.fig.set_size_inches(7, 10)
            plt.title('Distribution of scores within {}'.format(variant))
            pp.savefig()

            # HSP length
            sns.set(style="darkgrid")
            colors = sns.color_palette("Set2", 7)
            ax = sns.catplot(x="hsp_length", y="accession",
                             hue="is_positive", row="variant",
                             data=hsps_data[hsps_data['variant'] == variant],
                             orient="h", palette="Set2", split=True, kind="violin")
            ax.fig.set_size_inches(7, 10)
            plt.title('Distribution of HSP length within {}'.format(variant))
            pp.savefig()

            # Relative length
            sns.set(style="darkgrid")
            colors = sns.color_palette("Set2", 7)
            ax = sns.catplot(x="relative_length", y="accession",
                             hue="is_positive", row="variant",
                             data=hsps_data[hsps_data['variant'] == variant],
                             orient="h", palette="Set2", split=True, kind="violin")
            ax.fig.set_size_inches(7, 10)
            plt.title('Distribution of relative length within {}'.format(variant))
            pp.savefig()

            # Relation between score and HSP length
            sns.set(style="darkgrid")
            colors = sns.color_palette("Set2", 7)
            ax = sns.relplot(data=hsps_data[hsps_data['variant'] == variant],
                             x="hsp_length", y="score",
                             hue="is_positive", style="accession", palette="Set2")
            ax.fig.set_size_inches(7, 10)
            plt.title('Relationship between score and HSP length within {}'.format(variant))
            pp.savefig()

            # Relation between score and relative length
            sns.set(style="darkgrid")
            colors = sns.color_palette("Set2", 7)
            ax = sns.relplot(data=hsps_data[hsps_data['variant'] == variant],
                             x="relative_length", y="score",
                             hue="is_positive", style="accession", palette="Set2")
            ax.fig.set_size_inches(7, 10)
            plt.title('Relationship between score and relative length within {}'.format(variant))
            pp.savefig()

            pp.close()

def plot_scores(seq_accession, save_dir, scores_sorted, best_threshold=None):
    """Plot the scores
    Returns:
    A dictionary with images of dynamic score changes
    """

    pp = PdfPages(os.path.join(save_dir, "{}_model_evaluation.pdf".format(seq_accession)))

    sns.set(style="darkgrid")
    colors = sns.color_palette("Set2", 7)

    plt.plot(range(len(scores_sorted)), scores_sorted, linewidth=2)
    if best_threshold: plt.axhline(best_threshold)

    plt.legend()
    plt.ylabel("Score")
    plt.xlabel("Threshold")

    plt.title("{} Model Evaluation".format(seq_accession), fontsize=20)

    pp.savefig()
    pp.close()

