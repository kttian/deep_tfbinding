#!/usr/bin/env python
from __future__ import print_function, division

import logging
import sys
import os
import matplotlib as mpl
mpl.use('Agg')

import numpy as np
import modisco
import theano
import sys
import argparse

def debug_signal_handler(signal, frame):
    import pdb
    pdb.set_trace()
import signal
signal.signal(signal.SIGINT, debug_signal_handler)

logging.basicConfig(
        format='%(asctime)s %(levelname)-5s %(message)s',
        level=logging.DEBUG,
        datefmt='%Y-%m-%d %H:%M:%S')

logging.info(" ".join(sys.argv))


logging.debug("Theano version:" + str(theano.__version__))
logging.debug(sys.version)


# ### Functions for one-hot encoding sequences

import gzip

def one_hot_encode_along_channel_axis(sequence):
    #theano dim ordering, uses row axis for one-hot
    to_return = np.zeros((len(sequence),4), dtype=np.int8)
    seq_to_one_hot_fill_in_array(zeros_array=to_return,
                                 sequence=sequence, one_hot_axis=1)
    return to_return

def seq_to_one_hot_fill_in_array(zeros_array, sequence, one_hot_axis):
    assert one_hot_axis==0 or one_hot_axis==1
    if (one_hot_axis==0):
        assert zeros_array.shape[1] == len(sequence)
    elif (one_hot_axis==1): 
        assert zeros_array.shape[0] == len(sequence)
    #will mutate zeros_array
    for (i,char) in enumerate(sequence):
        if (char=="A" or char=="a"):
            char_idx = 0
        elif (char=="C" or char=="c"):
            char_idx = 1
        elif (char=="G" or char=="g"):
            char_idx = 2
        elif (char=="T" or char=="t"):
            char_idx = 3
        elif (char=="N" or char=="n"):
            continue #leave that pos as all 0's
        else:
            raise RuntimeError("Unsupported character: "+str(char))
        if (one_hot_axis==0):
            zeros_array[char_idx,i] = 1
        elif (one_hot_axis==1):
            zeros_array[i,char_idx] = 1

from merge_overlaps import MergeOverlaps
from merge_overlaps import merge_overlaps

def parse_args(args = None):
    parser = argparse.ArgumentParser('run_tfmodisco.py',
                                     description='run tfmodisco',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--scores', type=str, help="prefix for the hypothetical score files")
    parser.add_argument('--fasta', type=str, help="fasta input")
    parser.add_argument('--tsv', type=str, help="tsv input")
    parser.add_argument('--start-task', type=int, default=0, help="start tast")
    parser.add_argument('--end-task', type=int, default=5, help="end task")
    parser.add_argument('--fdr', type=float, default=0.01, help="target FDR")
    parser.add_argument('--min-seqlets', type=int, default=1000, help="min seqlets")
    args = parser.parse_args(args)
    return args

args = parse_args()

score_prefix = args.scores
input_name   = args.fasta
input_tsv    = args.tsv
start_task   = args.start_task
end_task     = args.end_task
target_fdr   = args.fdr

logging.debug("method file prefix is %s, input seq file is %s, input tsv is %s, start_task is %d end_task is %d, fdr is %f", 
              score_prefix, input_name, input_tsv, start_task, end_task, target_fdr)

#https://www.biostars.org/p/710/
from itertools import groupby
def fasta_iter(fasta_name):
    """
        given a fasta file, yield tuples of (header, sequence)
    """
    fh = open(fasta_name) # file handle
    # ditch the boolean (x[0]) and just keep the header or sequence since they alternate
    fa_iter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in fa_iter:
        header = header.next()[1:].strip() # drop the ">" from the header
        seq = "".join(s.strip() for s in fa_iter.next()) # join all sequence lines to one
        yield header, seq

fasta_sequences = []
fasta = fasta_iter(input_name)

for header, seq in fasta:   
    fasta_sequences.append(seq)
logging.debug("lenth of sequences = %d", len(fasta_sequences))

#onehot_data = [one_hot_encode_along_channel_axis(seq) for seq in fasta_sequences]
#logging.debug("shape of onehot" + str(onehot_data[0].shape))

# ## Prepare the data for input into TF-MoDISCo
# 
# You need a numpy array of importance scores and hypothetical importance scores for every task.

from collections import OrderedDict

task_to_scores = OrderedDict()
task_to_hyp_scores = OrderedDict()

# locations of deeplift scores
scores_loc = []
task_names = []
for i in range(start_task , end_task):
    loc_i = score_prefix + str(i) + ".npy"
    scores_loc.append(loc_i)
    task_names.append("task" + str(i))

# scores & their one-hot encodings

merged_seq_list        = []
merged_tsv_list        = []
merged_onehot_list     = []
num_tasks = end_task - start_task
for t in range(start_task, end_task):
    merged_hyp_scores_list     = []
    merged_contrib_scores_list = []

    task = task_names[t - start_task]
    hyp_scores_all = np.load(scores_loc[t - start_task])
    merge_overlaps(input_tsv, hyp_scores_all, merged_hyp_scores_list, fasta_sequences,
                   merged_seq_list = merged_seq_list if t==start_task else None,
                   merged_tsv_list = merged_tsv_list if t==start_task else None)

    for i in range(len(merged_hyp_scores_list)):
        onehot_seq = one_hot_encode_along_channel_axis(merged_seq_list[i])
        contrib_scores = merged_hyp_scores_list[i] * onehot_seq
        merged_contrib_scores_list.append(contrib_scores)
        if t == start_task:
            merged_onehot_list.append(onehot_seq)

    task_to_hyp_scores[task] = merged_hyp_scores_list
    task_to_scores[task]     = merged_contrib_scores_list

    if t == start_task:
        logging.debug("shape of hyp_score " + str(task_to_hyp_scores[task][0].shape))
        logging.debug("shape of score " + str(task_to_scores[task][0].shape))

#save the scores and tsv
'''
np.savez("merged_hyp_scores.npz", *merged_hyp_scores_list)
with open("merged_overlaps.tsv", 'w') as fh:
    for tsv in merged_tsv_list:
        fields = [str(f) for f in tsv]
        fh.write("\t".join(fields) + "\n")
'''

import h5py
import numpy as np
import modisco
#reload(modisco)
import modisco.backend
#reload(modisco.backend.theano_backend)
#reload(modisco.backend)
import modisco.nearest_neighbors
#reload(modisco.nearest_neighbors)
import modisco.affinitymat
#reload(modisco.affinitymat.core)
#reload(modisco.affinitymat.transformers)
import modisco.tfmodisco_workflow.seqlets_to_patterns
#reload(modisco.tfmodisco_workflow.seqlets_to_patterns)
import modisco.tfmodisco_workflow.workflow
#reload(modisco.tfmodisco_workflow.workflow)
import modisco.aggregator
#reload(modisco.aggregator)
import modisco.cluster
#reload(modisco.cluster.core)
#reload(modisco.cluster.phenograph.core)
#reload(modisco.cluster.phenograph.cluster)
import modisco.core
#reload(modisco.core)
import modisco.coordproducers
#reload(modisco.coordproducers)
import modisco.metaclusterers
#reload(modisco.metaclusterers)

factory = modisco.tfmodisco_workflow.seqlets_to_patterns.TfModiscoSeqletsToPatternsFactory(
              trim_to_window_size=30,
              initial_flank_to_add=10,
              kmer_len=8, num_gaps=3,
              num_mismatches=2,
              final_min_cluster_size=30)

tfmodisco_results = modisco.tfmodisco_workflow.workflow.TfModiscoWorkflow(
                            sliding_window_size=21,
                            flank_size=10,
                            target_seqlet_fdr=target_fdr,
                            seqlets_to_patterns_factory=factory,
                            min_seqlets_per_task=args.min_seqlets
                        )(
                            task_names=task_names,
                            contrib_scores        = task_to_scores,
                            hypothetical_contribs = task_to_hyp_scores,
                            one_hot=merged_onehot_list)

"""
#Slight modifications from the default settings
factory = modisco.tfmodisco_workflow.seqlets_to_patterns.TfModiscoSeqletsToPatternsFactory(
              trim_to_window_size=15,
              initial_flank_to_add=5,
              kmer_len=5, num_gaps=1,
              num_mismatches=0,
              final_min_cluster_size=60)

tfmodisco_results = modisco.tfmodisco_workflow.workflow.TfModiscoWorkflow(
                            sliding_window_size=15,
                            flank_size=5,
                            target_seqlet_fdr=0.01,
                            seqlets_to_patterns_factory=factory
                        )(
                            task_names=task_names,
                            contrib_scores        = task_to_scores,
                            hypothetical_contribs = task_to_hyp_scores,
                            one_hot=merged_onehot_list)
"""

logging.debug("**************** workflow done *********************")

# ## Save and print the results

# In[8]:

import h5py
import modisco.util
#reload(modisco.util)
os.system('rm -f results.hdf5')
grp = h5py.File("results.hdf5")
tfmodisco_results.save_hdf5(grp)

logging.debug("**************** result saved *********************")



