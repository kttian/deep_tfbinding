# extract seqlets from results.hdf5, write into seqlets.tsv

from __future__ import print_function, division

import logging
import sys
import os
#import matplotlib as mpl
#mpl.use('Agg')

import numpy as np
import modisco
import sys
import argparse
TFNET_ROOT = os.getenv('TFNET_ROOT', '/home/ktian/kundajelab/tfnet/')
sys.path.append(TFNET_ROOT + "/scripts")
    
logging.basicConfig(
        format='%(asctime)s %(levelname)-5s %(message)s',
        level=logging.INFO,
        datefmt='%Y-%m-%d %H:%M:%S')

logging.info(" ".join(sys.argv))
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

'''
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
    args = parser.parse_args(args)
    return args

args = parse_args()

score_prefix = args.scores
input_name   = args.fasta
input_tsv    = args.tsv
start_task   = args.start_task
end_task     = args.end_task
target_fdr   = args.fdr
'''


#convert the motifs to log-odds space
def log_odds_space(pwm, background,pseudocount):
    new_pwm = []
    for pos_values in np.transpose(pwm,(1,0)):
        if sum(pos_values)==0:
            new_pwm.append(pos_values)
        else:
            pos_values = pos_values+pseudocount/(1+pseudocount*4)
            new_pwm.append(np.log(pos_values) - np.log(background))
    return np.array(new_pwm).transpose(1,0)

score_prefix = "../scores/hyp_scores_task_"
input_name   = "../interpret.fa"
input_tsv    = "../interpret.tsv"

start_task = 0
end_task   = 1

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
merged_tsv_list        = []
num_tasks = end_task - start_task
for t in range(num_tasks):
    merged_hyp_scores_list     = []
    merged_contrib_scores_list = []
    task = task_names[t]
    hyp_scores_all = np.load(scores_loc[t])
    merge_overlaps(input_tsv, hyp_scores_all, merged_hyp_scores_list, fasta_sequences,
                   merged_seq_list = merged_seq_list if t==0 else None,
                   merged_tsv_list = merged_tsv_list if t==0 else None                   )

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


# In[2]:


#save the scores and tsv
'''
np.savez("merged_scores.npz", *merged_hyp_scores_list)
with open("merged.tsv", 'w') as fh:
    for tsv in merged_tsv_list:
        fields = [str(f) for f in tsv]
        fh.write("\t".join(fields) + "\n")
'''


# ## Load the saved hdf5 file

# Load the results object from the saved file

# In[3]:


try:
    reload  # Python 2.7
except NameError:
    try:
        from importlib import reload  # Python 3.4+
    except ImportError:
        from imp import reload  # Python 3.0 - 3.3

import h5py
import numpy as np
from collections import OrderedDict
import modisco.util
import modisco.core
reload(modisco.core)
import modisco.metaclusterers
reload(modisco.metaclusterers)
import modisco.coordproducers
reload(modisco.coordproducers)
import modisco.tfmodisco_workflow.seqlets_to_patterns
import modisco.tfmodisco_workflow
reload(modisco.tfmodisco_workflow.seqlets_to_patterns)
reload(modisco.tfmodisco_workflow)
from modisco.tfmodisco_workflow import workflow
reload(workflow)


# read snp.txt format: chr:pos

# load the motifs from the results
track_set = modisco.tfmodisco_workflow.workflow.prep_track_set(
                task_names=task_names,
                contrib_scores=task_to_scores,
                hypothetical_contribs=task_to_hyp_scores,
                one_hot=merged_onehot_list)

grp = h5py.File("results.hdf5","r")
loaded_tfmodisco_results =    workflow.TfModiscoResults.from_hdf5(grp, track_set=track_set)
grp.close()


# In[4]:


tf_res = loaded_tfmodisco_results


# In[5]:


from modisco import affinitymat
reload(affinitymat.core)
reload(affinitymat)
from modisco import hit_scoring
reload(hit_scoring.fast_hit_scoring)
reload(hit_scoring)
from collections import OrderedDict

task_names = loaded_tfmodisco_results.task_names

seqlet_size_to_score_with = 25

metacluster_idx_to_scorer = OrderedDict()

all_pattern_scorers = []
all_pattern_names = []

for metacluster_name in    sorted(loaded_tfmodisco_results
           .metacluster_idx_to_submetacluster_results.keys()):
    submetacluster_results =(
        loaded_tfmodisco_results
            .metacluster_idx_to_submetacluster_results[metacluster_name])
    activity_pattern = submetacluster_results.activity_pattern
    relevant_task_names = [task_name for (task_name,x) in
                           zip(task_names, activity_pattern) if np.abs(x) != 0]
    
    patterns_in_submetacluster =        submetacluster_results.seqlets_to_patterns_result.patterns
 

    for pattern_idx, pattern in        enumerate(submetacluster_results.
                   seqlets_to_patterns_result.patterns):
        metacluster_idx = int(metacluster_name.split("_")[1])
        all_pattern_names.append("metacluster_"+str(metacluster_idx)
                             +",pattern_"+str(pattern_idx))


print(all_pattern_names)


# In[6]:


all_patterns = [x for y in
                  sorted(loaded_tfmodisco_results
                  .metacluster_idx_to_submetacluster_results.keys())
                  for x in
                   loaded_tfmodisco_results
                   .metacluster_idx_to_submetacluster_results[y]
                   .seqlets_to_patterns_result.patterns]
seqlets_to_score = []
seqlets_to_score_true_labels = []
for i,pattern in enumerate(all_patterns):
    seqlets_to_score.extend(pattern.seqlets)
    seqlets_to_score_true_labels.extend(
        [i for x in pattern.seqlets])

MAX_LOC = 250000000L
    
def seq_to_bed(seq):
    st = seq.coor.start
    en = seq.coor.end
    local = int((st + en)/2) # fix later using center of gravity

    tsv = merged_tsv_list[seq.coor.example_idx]
    location = int(tsv[1]) + local
    if location > MAX_LOC:
        print("location=%d, tsv_start=%d, ex=%d, st=%d, en=%d" %(location, int(tsv[1]), seq.coor.example_idx, seq.coor.start, seq.coor.end))
        print(tsv)

    return tsv[0], int(tsv[1]) + st, int(tsv[1]) + en


with open("seqlets.tsv", "w") as fh:
    for i, seq in enumerate(seqlets_to_score):
        chrom, st, en = seq_to_bed(seq)
        fh.write(chrom + "\t" + str(st) + "\t" + str(en) + "\n")

