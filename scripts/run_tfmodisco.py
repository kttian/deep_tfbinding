#!/usr/bin/env python
from __future__ import print_function, division

import logging
logging.basicConfig(
        format='%(asctime)s %(levelname)-5s %(message)s',
        level=logging.DEBUG,
        datefmt='%Y-%m-%d %H:%M:%S')

import sys
logging.info(" ".join(sys.argv))

import os
import matplotlib as mpl
mpl.use('Agg')

import numpy as np
import modisco
import theano
import sys

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

#do_test()
#quit()

import sys
if len(sys.argv) < 5 or len(sys.argv) > 6:
    print("Syntax: ", sys.argv[0] , " <score prefix> <sequence fa file> <sequence tsv> {<number of tasks> | <start task> <end task>}")
    quit()

score_prefix = sys.argv[1]                 # "./scores/hyp_scores_task_"
input_name   = sys.argv[2]                 # subset.fa, sequences 
input_tsv    = sys.argv[3]                 # subset.tsv
if len(sys.argv) == 5:
    start_task    = 0
    end_task    = int(sys.argv[4])            # 
else:
    start_task    = int(sys.argv[4])            # 
    end_task    = int(sys.argv[5])            # 


logging.debug("method file prefix is %s, input seq file is %s, input tsv is %s, start_task is %d end_task is %d", 
              score_prefix, input_name, input_tsv, start_task, end_task)

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
merged_onehot_list     = []
num_tasks = end_task - start_task
for t in range(num_tasks):
    merged_hyp_scores_list     = []
    merged_contrib_scores_list = []

    task = task_names[t]
    hyp_scores_all = np.load(scores_loc[t])
    merge_overlaps(input_tsv, hyp_scores_all, merged_hyp_scores_list, fasta_sequences,
                   merged_seq_list = merged_seq_list if t==0 else None)

    for i in range(len(merged_hyp_scores_list)):
        onehot_seq = one_hot_encode_along_channel_axis(merged_seq_list[i])
        contrib_scores = merged_hyp_scores_list[i] * onehot_seq
        merged_contrib_scores_list.append(contrib_scores)
        if t == 0:
            merged_onehot_list.append(onehot_seq)

    task_to_hyp_scores[task] = merged_hyp_scores_list
    task_to_scores[task]     = merged_contrib_scores_list

    if t == 0:
        logging.debug("shape of hyp_score " + str(task_to_hyp_scores[task][0].shape))
        logging.debug("shape of score " + str(task_to_scores[task][0].shape))

"""
for i in range(num_tasks):
    task = 'task'+str(i)
    task_to_hyp_scores[task] = np.load(scores_loc[i])
    logging.debug("shape of hyp_score", task_to_hyp_scores['task0'].shape)
    task_to_scores[task]     = task_to_hyp_scores[task] * onehot_data
    logging.debug("shape of score",     task_to_scores['task0'].shape)

logging.debug(onehot_data[0][:5])
logging.debug("hyp_scores")
logging.debug(task_to_hyp_scores["task0"][0][:5])
logging.debug("scores")
logging.debug(task_to_scores["task0"][0][:5])

# for hdf5 input
f = h5py.File(score_file,"r")
tasks = f["contrib_scores"].keys()
logging.debug("tasks are: ", tasks)

for task in tasks:
    #Note that the sequences can be of variable lengths;
    #in this example they all have the same length (200bp) but that is
    #not necessary.
    task_to_scores[task] = [np.array(x) for x in f['contrib_scores'][task][:]]
    task_to_hyp_scores[task] = [np.array(x) for x in f['hyp_contrib_scores'][task][:]]
"""

# ## Run TF-MoDISco
# 
# TF-MoDISco first identifies seqlets, then splits the seqlets into
# "metaclusters" according to their pattern of activity across all the tasks,
# and then performs clustering within each task. Since there is just one task,
# there are only 2 possible metaclusters: +1 for the task and -1 for the task.
# The -1 metacluster does not turn up any motifs after noise filtering, but
# the +1 metacluster produces two motifs.
# 
# To demonstrate customization, the code below has slight modifications from
# default settings in the following ways: - Because the TAL and GATA motifs
# are relatively short compared to something like CTCF, it uses a sliding
# window size of 15 (rather than the default of 21) and flanks of 5 (rather
# than the default of 10). The sliding window size and flanks should be
# adjusted according to the expected length of the core motif and its flanks.
# If the window size or flank sizes are too long, you risk picking up more
# noise.  - During the seqlet clustering, motifs are trimmed to the central
# `trim_to_window_size` bp with the highest importance. `trim_to_window_size`
# is set to 10 rather than the default of 30. After the trimming is done, the
# seqlet is expanded on either side by `initial_flank_to_add`. This is set to
# 3 rather than the default of 10.  - The `final_min_cluster_size` is set to
# 60 rather than the default of 30. This is used to filter out small clusters
# with relatively weak support (in this case, fewer than 50 seqlets).  - It
# uses kmers of length 5 with 1 gap and no mismatches to compute the "quick
# and dirty" affinity matrix across all seqlets. The "quick and dirty"
# affinity matrix is used both for noise filtering and as a first pass to
# speed up computation of the continuous jaccard affinity matrix (the latter
# affinities are only computed between seqlets deemed to be close together by
# the "quick and dirty" method). I made the kmer length smaller to keep memory
# usage on the GPU down when testing on my macbook pro. The default is to use
# kmers of length 8 with 3 gaps and 2 mismatches, and this works fine on more
# modern GPUs than the one in my 4-year-old macbook.  - `target_seqlet_fdr`
# controls the noisiness of the seqelts. For a particular task, "significant"
# seqlets are identified by fitting a laplace distribution to the left and
# right tails of the values obtained after smoothing the importance scores
# with a window of size `sliding_window_size`. This laplace distribution is
# assumed to represent the null distribution of random seqlet importance
# scores. A threshold is then identified such that the false discovery rate
# (computed as the ratio of the expected number of seqlets according to the
# laplace null to the observed number of seqlets above the threshold) is less
# that `target_seqlet_fdr`. This is what is meant by "Est. FDR" printed in the
# logs below. If "Est. FDR" is above the target threshold, that means there
# was no significant increase in the number of seqlets relative to the null.
# You'll see below that "Est. FDR" for negative scores for any task is above
# this threshold, which fits with the simulation because there were no
# "negative set" motifs.

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
                            target_seqlet_fdr=0.01,
                            seqlets_to_patterns_factory=factory
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



