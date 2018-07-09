#!/usr/bin/env python
from __future__ import print_function, division

import logging
logging.basicConfig(
        format='%(asctime)s %(levelname)-5s %(message)s',
        level=logging.DEBUG,
        datefmt='%Y-%m-%d %H:%M:%S')

logging.debug("start run_tf_modisco")
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

# ### Class to calculate average deeplift scores among overlapping sequences
class AverageScores:

    def __init__(self, avg_hyp_scores_list, avg_scores_seq_list,
                 chrom=None, avg_st=0, avg_en=0, core_size=400):
        ''' 
        constructor
        self.avg_hyp_scores_list : newly calculated average hyp scores list
        self.avg_scores_seq_list : new sequence list (after merging/averaging overlaps)
        self.avg_st       : start of the new interval
        self.avg_en       : end   of the new interval
        self.avg_scores   : new scores array, will be add to avg_hyp_scores_list
        self.avg_counts   : counts array, has the same len as the self.avg_scores
        self.avg_seq      : new sequence, will be add to avg_scores_seq_list
        self.max_seq_size : max seq size this obj has seen
        '''
        self.avg_hyp_scores_list = avg_hyp_scores_list
        self.avg_scores_seq_list        = avg_scores_seq_list
        self.chrom      = chrom      # rest init at new interval
        self.avg_st     = avg_st
        self.avg_en     = avg_en
        self.avg_seq    = ""
        self.avg_scores = np.zeros((0, 4))
        self.avg_counts = np.zeros(0)
        self.max_seq_size = 0
        self.core_size    = core_size
        logging.debug("AverageScores: %s", 
                      "append to seq_list" if self.avg_scores_seq_list != None else "-" )

    def close_interval(self):

        if self.avg_st == self.avg_en: return # nothing to close

        ### generate a sequence between avg_st, avg_en ###

        self.avg_scores /= self.avg_counts[:, None] # calculate the average

        self.avg_hyp_scores_list.append(self.avg_scores)
        if self.avg_scores_seq_list != None:
            self.avg_scores_seq_list.append(self.avg_seq)

        if self.avg_en - self.avg_st > self.max_seq_size : # i'm keeping track of the max sequence size
            self.max_seq_size = self.avg_en - self.avg_st

        '''
        logging.debug("close interval: %s:%d-%d, len=%d max=%d %s", 
                      self.chrom, self.avg_st, self.avg_en, 
                      self.avg_en - self.avg_st, self.max_seq_size,
                      "append to seq_list" if self.avg_scores_seq_list != None else "-" )
        # write the interval to file
        if self.ofh :
            self.ofh.write(self.chrom + "\t" + str(self.avg_st) + "\t" + \
                           str(self.avg_en) + "\n")
        logging.debug("scores = ")
        logging.debug(self.avg_scores)
        logging.debug("counts = ")
        logging.debug(self.avg_counts)
        '''

        # reset for calculating the avg of the next set of overlapping intervals
        self.avg_st = self.avg_en
        self.avg_scores = np.zeros((0, 4))
        self.avg_counts = np.zeros(0)
        self.avg_seq    = ""

    def append_interval(self, chrom, st, en, scores_hyp, seq):
        '''
        process the next interval <st, en>

        Parameters:
        ----------
        st:     start of the next interval to be processed
        en:     end   of the interval
        scores_hyp: input hyp scores array of the next interval, offset st : en

        Returns:
        ----------
        self.avg_en: the new end

                self.avg_st     self.avg_en
                          v     v
        self.avg_scores:  SSSSSS
               O=Overlap    OOOOAA    A=Append
        scores_hyp:         HHHHHH
                            ^     ^
                            st    en
        '''
        
        assert(st >= self.avg_st) # intervals must be sorted
        if chrom != self.chrom or st > self.avg_en : # handle the non-overlapping case
            self.chrom  = chrom
            self.avg_st = st
            self.avg_en = st

        '''
        length of the overlap = self.avg_en - st
        overlap is self.avg_scores[st - self.avg_st : self.avg_en - self.avg_st]
        overlap is scores_hyp[0 : self.avg_en - st]
        '''
        if st > self.avg_st : # non-empty overlap
            self.avg_scores[st - self.avg_st : self.avg_en - self.avg_st] += \
                scores_hyp[0:self.avg_en - st]
            self.avg_counts[st - self.avg_st : self.avg_en - self.avg_st] += \
                np.ones(self.avg_en - st)
        
        '''
        part to append, len = en - self.avg_en
        scores_hyp[self.avg_en-st : en-st    ]
        '''
        if en > self.avg_en : # non-empty append
            self.avg_scores = np.concatenate((self.avg_scores,
                                              scores_hyp[self.avg_en - st : en - st]), 
                                              axis = 0)
            self.avg_counts = np.concatenate((self.avg_counts, np.ones(en-self.avg_en)), 
                                              axis = 0)
            if self.avg_scores_seq_list != None:
                self.avg_seq    += seq[self.avg_en - st : en - st]
            self.avg_en = en

        return self.avg_en


    def process_one_interval(self, chrom, st, en, scores_hyp, seq):
        
        seq_size = en - st
        left     = int((seq_size - self.core_size) / 2)
        right    = left + self.core_size

        if chrom != self.chrom or st > self.avg_en: # start a new interval
            self.close_interval()

        #self.append_interval(chrom, st, en, scores_hyp, seq)
        self.append_interval(chrom, st + left, st + right, scores_hyp[left:right], 
                             seq[left:right])

        #logging.debug("processed interval %s:%d-%d, reduce to %d-%d",
        #              chrom, st, en, st+left, st+right)

def average_scores(in_tsv, hyp_scores_all, avg_hyp_scores_list, seq_list, avg_scores_seq_list):

    logging.debug("in_tsv = " + in_tsv)
    with open(in_tsv,'r') as tsvin:
        avg = AverageScores(avg_hyp_scores_list, avg_scores_seq_list)
        for idx, line in enumerate(tsvin):
            row = line.split()
            chrom = row[0]
            st    = int(row[1])
            en    = int(row[2])

            scores_hyp = hyp_scores_all[idx]
            avg.process_one_interval(chrom, st, en, scores_hyp, seq_list[idx])


# ### testing code for AverageScores
def gen_tests(seq_size, seq_stride):
    letters = ['A','C', 'T', 'G']
    tsv_list = []
    seq_list = []
    hyp_scores_all = np.zeros((10, seq_size, 4))
    for i in range(5):
        chrom ='chr1'
        st = 100 + seq_stride * i
        en = st + seq_size
        tsv_list.append([chrom, st, en])
        seq_list.append(letters[i%4]*seq_size)
        hyp_scores_all[i].fill(0.1 * (i+1))

    for i in range(5,10):
        chrom ='chr1'
        st = 500 + seq_stride * i
        en = st + seq_size
        tsv_list.append([chrom, st, en])
        seq_list.append(letters[i%4]*seq_size)
        hyp_scores_all[i].fill(0.1 * (i+1))

    return tsv_list, hyp_scores_all, seq_list


def do_test() :
    seq_size   = 10
    seq_stride = 2
    tsv_list, hyp_scores_all, seq_list = gen_tests(seq_size, seq_stride)

    logging.debug(tsv_list)
    logging.debug(hyp_scores_all)
    avg_hyp_scores_list = []
    avg_scores_seq_list        = []
    logging.debug("gen_tests DONE")

    avg = AverageScores(avg_hyp_scores_list, avg_scores_seq_list=avg_scores_seq_list,
                        core_size=8)
    for i in range(len(tsv_list)):
        chrom, st, en = tsv_list[i]
        scores_hyp = hyp_scores_all[i]
        seq        = seq_list[i]
        avg.process_one_interval(chrom, st, en, scores_hyp, seq)
    avg.close_interval()

    #logging.debug(tsv_list)
    for i in range(len(avg_hyp_scores_list)) :
        ar = avg_hyp_scores_list[i]
        logging.debug("i=%d shape of array = %s", i, str(ar.shape))
        logging.debug(ar)
        seq = avg_scores_seq_list[i]
        logging.debug(seq)

    logging.debug("average score DONE")

#do_test()
#quit()

import sys
if len(sys.argv) != 5:
    print("Syntax: ", sys.argv[0] , " <score prefix> <sequence fa file> <sequence tsv> <number of tasks> ")
    quit()

score_prefix = sys.argv[1]                 # "./scores/hyp_scores_task_"
input_name   = sys.argv[2]                 # subset.fa, sequences 
input_tsv    = sys.argv[3]                 # subset.tsv
num_tasks    = int(sys.argv[4])            # 

logging.debug("method file prefix is %s, input seq file is %s, input tsv is %s, number of tasks is %d", 
              score_prefix, input_name, input_tsv, num_tasks)

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
for i in range(num_tasks):
    loc_i = score_prefix + str(i) + ".npy"
    scores_loc.append(loc_i)
    task_names.append("task" + str(i))

# scores & their one-hot encodings

avg_scores_seq_list        = []
avg_scores_onehot_list     = []
for t in range(num_tasks):
    avg_hyp_scores_list    = []
    avg_target_scores_list = []

    task = task_names[t]
    hyp_scores_all = np.load(scores_loc[t])
    average_scores(input_tsv, hyp_scores_all, avg_hyp_scores_list, fasta_sequences,
                   avg_scores_seq_list = avg_scores_seq_list if t==0 else None)

    for i in range(len(avg_hyp_scores_list)):
        onehot_seq = one_hot_encode_along_channel_axis(avg_scores_seq_list[i])
        target_scores = avg_hyp_scores_list[i] * onehot_seq
        avg_target_scores_list.append(target_scores)
        if t == 0:
            avg_scores_onehot_list.append(onehot_seq)

    task_to_hyp_scores[task] = avg_hyp_scores_list
    task_to_scores[task]     = avg_target_scores_list

    logging.debug("shape of hyp_score " + str(task_to_hyp_scores['task0'][0].shape))
    logging.debug("shape of score " + str(task_to_scores['task0'][0].shape))

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
                            one_hot=avg_scores_onehot_list)


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



