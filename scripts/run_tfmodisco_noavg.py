
# coding: utf-8

# In[1]:

from __future__ import print_function, division

import datetime
import os

start = datetime.datetime.now()
print ("start time " + str(start))

import matplotlib as mpl
mpl.use('Agg')

#try:
#    reload  # Python 2.7
#except NameError:
#    try:
#        from importlib import reload  # Python 3.4+
#    except ImportError:
#        from imp import reload  # Python 3.0 - 3.3


# In[2]:

import numpy as np
import modisco
import theano
print("Theano version:",theano.__version__)
import sys
print (sys.version)


# ### Functions for one-hot encoding sequences

# In[4]:

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

import sys
if len(sys.argv) != 4:
    print("Syntax: ", sys.argv[0] , " <score prefix> <sequence fa file> <number of tasks>")
    quit()

score_prefix = sys.argv[1]                 # "./scores/hyp_scores_task_"
input_file   = sys.argv[2]                 # subset.fa, sequences 
num_tasks    = int(sys.argv[3])            # 

print("method file prefix is ", score_prefix, ", input seq file is ", input_file, ", number of tasks is ", num_tasks)

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
fasta = fasta_iter(input_file)

for header, seq in fasta:   
    fasta_sequences.append(seq)
print("lenth of sequences = ", len(fasta_sequences))

onehot_data = [one_hot_encode_along_channel_axis(seq) for seq in fasta_sequences]

# ## Prepare the data for input into TF-MoDISCo
# 
# You need a numpy array of importance scores and hypothetical importance scores for every task.

# In[5]:

from collections import OrderedDict

task_to_scores = OrderedDict()
task_to_hyp_scores = OrderedDict()

# locations of deeplift scores
scores_loc = []
task_names = []
for i in range(num_tasks):
    loc_i = score_prefix + str(i) + ".npy"
    scores_loc.append(loc_i)
    task_names.append("task" + str(i))


print("shape of onehot",     onehot_data[0].shape)

# scores & their one-hot encodings
task_to_scores     = {}
task_to_hyp_scores = {}
for i in range(num_tasks):
    task = 'task'+str(i)
    task_to_hyp_scores[task] = np.load(scores_loc[i])
    print("shape of hyp_score", task_to_hyp_scores['task0'].shape)
    #task_to_scores[task]     = onehot_data * task_to_hyp_scores[task][:, :, None]
    task_to_scores[task]     = task_to_hyp_scores[task] * onehot_data
    print("shape of score",     task_to_scores['task0'].shape)

scores = task_to_scores
#scores_for_idx_print = original_onehot*scores_for_idx[:,None]
print(onehot_data[0][:5])
print("hyp_scores")
print(task_to_hyp_scores["task0"][0][:5])
print("scores")
print(task_to_scores["task0"][0][:5])




"""
f = h5py.File(score_file,"r")
tasks = f["contrib_scores"].keys()
print("tasks are: ", tasks)

for task in tasks:
    #Note that the sequences can be of variable lengths;
    #in this example they all have the same length (200bp) but that is
    #not necessary.
    task_to_scores[task] = [np.array(x) for x in f['contrib_scores'][task][:]]
    task_to_hyp_scores[task] = [np.array(x) for x in f['hyp_contrib_scores'][task][:]]
"""

# ## Run TF-MoDISco
# 
# TF-MoDISco first identifies seqlets, then splits the seqlets into "metaclusters" according to their pattern of activity across all the tasks, and then performs clustering within each task. Since there is just one task, there are only 2 possible metaclusters: +1 for the task and -1 for the task. The -1 metacluster does not turn up any motifs after noise filtering, but the +1 metacluster produces two motifs.
# 
# To demonstrate customization, the code below has slight modifications from default settings in the following ways:
# - Because the TAL and GATA motifs are relatively short compared to something like CTCF, it uses a sliding window size of 15 (rather than the default of 21) and flanks of 5 (rather than the default of 10). The sliding window size and flanks should be adjusted according to the expected length of the core motif and its flanks. If the window size or flank sizes are too long, you risk picking up more noise.
# - During the seqlet clustering, motifs are trimmed to the central `trim_to_window_size` bp with the highest importance. `trim_to_window_size` is set to 10 rather than the default of 30. After the trimming is done, the seqlet is expanded on either side by `initial_flank_to_add`. This is set to 3 rather than the default of 10.
# - The `final_min_cluster_size` is set to 60 rather than the default of 30. This is used to filter out small clusters with relatively weak support (in this case, fewer than 50 seqlets).
# - It uses kmers of length 5 with 1 gap and no mismatches to compute the "quick and dirty" affinity matrix across all seqlets. The "quick and dirty" affinity matrix is used both for noise filtering and as a first pass to speed up computation of the continuous jaccard affinity matrix (the latter affinities are only computed between seqlets deemed to be close together by the "quick and dirty" method). I made the kmer length smaller to keep memory usage on the GPU down when testing on my macbook pro. The default is to use kmers of length 8 with 3 gaps and 2 mismatches, and this works fine on more modern GPUs than the one in my 4-year-old macbook.
# - `target_seqlet_fdr` controls the noisiness of the seqelts. For a particular task, "significant" seqlets are identified by fitting a laplace distribution to the left and right tails of the values obtained after smoothing the importance scores with a window of size `sliding_window_size`. This laplace distribution is assumed to represent the null distribution of random seqlet importance scores. A threshold is then identified such that the false discovery rate (computed as the ratio of the expected number of seqlets according to the laplace null to the observed number of seqlets above the threshold) is less that `target_seqlet_fdr`. This is what is meant by "Est. FDR" printed in the logs below. If "Est. FDR" is above the target threshold, that means there was no significant increase in the number of seqlets relative to the null. You'll see below that "Est. FDR" for negative scores for any task is above this threshold, which fits with the simulation because there were no "negative set" motifs.

# In[7]:

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
                            seqlets_to_patterns_factory=factory)(

                            task_names=task_names,
                            contrib_scores=task_to_scores,
                            hypothetical_contribs=task_to_hyp_scores,
                            one_hot=onehot_data)


print("**************** workflow done *********************")

# ## Save and print the results

# In[8]:

import h5py
import modisco.util
#reload(modisco.util)
os.system('rm -f results.hdf5')
grp = h5py.File("results.hdf5")
tfmodisco_results.save_hdf5(grp)

end = datetime.datetime.now()
print ("start time " + str(start))
print ("end time " + str(end))
print("**************** result saved *********************")



