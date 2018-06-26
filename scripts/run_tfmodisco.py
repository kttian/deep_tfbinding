from __future__ import print_function, division

import datetime

start = datetime.datetime.now()
print ("start time " + str(start))

import matplotlib as mpl
mpl.use('Agg')

import os
resultDir = "../results/"
os.chdir(resultDir)

try:
    reload  # Python 2.7
except NameError:
    try:
        from importlib import reload  # Python 3.4+
    except ImportError:
        from imp import reload  # Python 3.0 - 3.3

import numpy as np
import modisco
import theano
print("Theano version:",theano.__version__)
import sys
print (sys.version)


# ### Functions for one-hot encoding sequences
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


# ## Prepare the data for input into TF-MoDISCo
# 
# You need a numpy array of importance scores and hypothetical importance scores for every task.


num_tasks = 5

# onehot_data is the one hot encoding of the input sequences used
onehot_data = np.array([one_hot_encode_along_channel_axis(x.rstrip()) 
                        for x in open("../results/sequences.txt")])
                        # generated during preprocessing


# locations of deeplift scores
scores_loc = []
for i in range(num_tasks):
    loc_i = "../results/scores/rescale_conv_revealcancel_fc_multiref_10_task_" + str(i) + ".npy"
    print(temp_loc)
    scores_loc.append(loc_i)

# scores & their one-hot encodings
orig_scores = {}
task_to_scores = {}
for i in range(num_tasks):
    orig_scores['task'+str(i)] = np.load(scores_loc[i])
    task_to_scores['task'+str(i)] = onehot_data * np.load(scores_loc[i])[:, :, None]

print("shape of onehot", onehot_data.shape)
print("shape of orig_score",  orig_scores['task0'].shape)
print("shape of score",  task_to_scores['task0'].shape)

scores = task_to_scores
#scores_for_idx_print = original_onehot*scores_for_idx[:,None]
print(onehot_data[0][:10])
print("orig_scores")
print(orig_scores["task0"][0][:10])
print("scores")
print(task_to_scores["task0"][0][:10])



import modisco.visualization
from modisco.visualization import viz_sequence

#viz_sequence.plot_weights(task_to_scores["task0"][0])
#viz_sequence.plot_weights(onehot_data[0][:100])


# ## Run TF-MoDISco
# 
# TF-MoDISco first identifies seqlets, then splits the seqlets into "metaclusters" according to their pattern of activity across all the tasks, and then performs clustering within each task. Since there is just one task, there are only 2 possible metaclusters: +1 for the task and -1 for the task. The -1 metacluster does not turn up any motifs after noise filtering, but the +1 metacluster produces two motifs.
# 
# To demonstrate customization, the code below has slight modifications from default settings in the following ways:
# - Because the TAL and GATA motifs are relatively short compared to something like CTCF, it uses a sliding window size of 15 (rather than the default of 21) and flanks of 5 (rather than the default of 10). The sliding window size and flanks should be adjusted according to the expected length of the core motif and its flanks. If the window size or flank sizes are too long, you risk picking up more noise.
# - During the seqlet clustering, motifs are trimmed to the central `trim_to_window_size` bp with the highest importance. `trim_to_window_size` is set to 10 rather than the default of 30. After the trimming is done, the seqlet is expanded on either side by `initial_flank_to_add`. This is set to 3 rather than the default of 10.
# - The `final_min_cluster_size` is set to 60 rather than the default of 30. This is used to filter out small clusters with relatively weak support (in this case, fewer than 50 seqlets).
# - It uses kmers of length 5 with 1 gap and no mismatches to compute the "quick and dirty" affinity matrix across all seqlets. The "quick and dirty" affinity matrix is used both for noise filtering and as a first pass to speed up computation of the continuous jaccard affinity matrix (the latter affinities are only computed between seqlets deemed to be close together by the "quick and dirty" method). I made the kmer length smaller to keep memory usage on the GPU down when testing on my macbook pro. The default is to use kmers of length 8 with 3 gaps and 2 mismatches, and this works fine on more modern GPUs than the one in my 4-year-old macbook.


import h5py
import numpy as np
#get_ipython().run_line_magic('matplotlib', 'inline')
import modisco
reload(modisco)
import modisco.backend
reload(modisco.backend.theano_backend)
reload(modisco.backend)
import modisco.nearest_neighbors
reload(modisco.nearest_neighbors)
import modisco.affinitymat
reload(modisco.affinitymat.core)
reload(modisco.affinitymat.transformers)
import modisco.tfmodisco_workflow.seqlets_to_patterns
reload(modisco.tfmodisco_workflow.seqlets_to_patterns)
import modisco.tfmodisco_workflow.workflow
reload(modisco.tfmodisco_workflow.workflow)
import modisco.aggregator
reload(modisco.aggregator)
import modisco.cluster
reload(modisco.cluster.core)
reload(modisco.cluster.phenograph.core)
reload(modisco.cluster.phenograph.cluster)
import modisco.core
reload(modisco.core)
import modisco.coordproducers
reload(modisco.coordproducers)
import modisco.metaclusterers
reload(modisco.metaclusterers)

tfmodisco_results = modisco.tfmodisco_workflow.workflow.TfModiscoWorkflow(
                    #Slight modifications from the default settings
                    sliding_window_size=15,
                    flank_size=5,
                    seqlets_to_patterns_factory=
                     modisco.tfmodisco_workflow.seqlets_to_patterns.TfModiscoSeqletsToPatternsFactory(
                        trim_to_window_size=15,
                        initial_flank_to_add=5,
                        kmer_len=5, num_gaps=1,
                        num_mismatches=0,
                        final_min_cluster_size=60,
                        verbose=False)
                )(
                task_names=["task0","task1","task2","task3","task4"],
                contrib_scores=task_to_scores,
                hypothetical_contribs=task_to_scores,
                one_hot=onehot_data)

print("**************** workflow done *********************")

import os
import h5py
import modisco.util
reload(modisco.util)
os.system('rm results.hdf5')
grp = h5py.File("results.hdf5")
tfmodisco_results.save_hdf5(grp)

end = datetime.datetime.now()
print ("start time " + str(start))
print ("end time " + str(end))
print("**************** result saved *********************")
