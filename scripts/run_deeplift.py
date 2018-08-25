# ## This notebook demonstrates how to get importance scores for use with TF-MoDISco using DeepLIFT
# 
# It relies on the TAL-GATA example model in the DeepLIFT repository

from __future__ import print_function, division
import logging
import sys
import datetime
import gzip
import numpy as np
import psutil
import os

from sys import getsizeof
from pympler.asizeof import asizeof

logging.basicConfig(
        format='%(asctime)s %(levelname)-5s %(message)s',
        level=logging.DEBUG,
        datefmt='%Y-%m-%d %H:%M:%S')

logging.info(" ".join(sys.argv))

start = datetime.datetime.now()
print ("start time " + str(start))


def log_mem_usage(step=100, msg=""):
    GB = 2**30

    logging.debug(msg)
    logging.debug(psutil.virtual_memory())
    pid = os.getpid()
    py = psutil.Process(pid)
    rss = py.memory_info()[0]/GB  # rss use in GB
    vms = py.memory_info()[1]/GB  # vms use in GB
    logging.debug('memory use: vms=%f G, rss=%f G', vms, rss)

    if step >= 1:
        global fasta_sequences
        global deeplift_model
        global keras_model
        logging.debug("size of fasta_sequences = %f (%f) G",
                      asizeof(fasta_sequences)/GB, getsizeof(fasta_sequences)/GB)
        logging.debug("size of deeplift_model  = %f (%f) G",
                      asizeof(deeplift_model)/GB, getsizeof(deeplift_model)/GB)
        logging.debug("size of keras_model     = %f (%f) G",
                      asizeof(keras_model)/GB, getsizeof(keras_model)/GB)
    
    if step >= 2:
        global hyp_scores_all
        global hyp_scores
        global contrib_scores
        logging.debug("size of hyp_scores_all = %f G", getsizeof(hyp_scores_all)/GB)
        logging.debug("size of hyp_scores     = %f G", getsizeof(hyp_scores)/GB)
        logging.debug("size of contrib_scores = %f G", getsizeof(contrib_scores)/GB)

log_mem_usage(0, "memory check location 1")

#this is set up for 1d convolutions where examples
#have dimensions (len, num_channels)
#the channel axis is the axis for one-hot encoding.
def one_hot_encode_along_channel_axis(sequence):
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
            
if len(sys.argv) < 4 or len(sys.argv) > 5:
    print("Syntax: ", sys.argv[0] , " <model name> <sequence file> {<start_task> <end_task(exclusive)> | <end task(exclusive)>}")
    quit()

keras_model_weights = sys.argv[1] + "Weights.h5"
keras_model_json    = sys.argv[1] + "Json.json"
input_file          = sys.argv[2]    # subset.txt, sequences without fasta header line ">"
if len(sys.argv) == 4:
    start_task      = 0
    end_task        = int(sys.argv[3])

else:
    start_task      = int(sys.argv[3])
    end_task        = int(sys.argv[4])

logging.info("loading models from " + keras_model_json + " " + keras_model_weights)
logging.info("input sequence file is " + input_file + ", range of tasks are " + str(start_task)+":"+str(end_task))

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
logging.debug("len of input sequences = %d", len(fasta_sequences))

log_mem_usage(0, "memory check location 2")

# ### Load the keras model

# In[4]:

import deeplift
from keras.models import model_from_json

#load the keras model
keras_model = model_from_json(open(keras_model_json).read())
keras_model.load_weights(keras_model_weights)



# ## Prepare the deeplift models
# 
# ### Model conversion
# 
# Convert the keras models to a deeplift model capable of computing importance
# scores using DeepLIFT.

# In[5]:

import deeplift.conversion.kerasapi_conversion as kc
from collections import OrderedDict

deeplift_model = kc.convert_model_from_saved_files(
                     h5_file=keras_model_weights,
                     json_file=keras_model_json)


'''
# ### Sanity checks
# To ensure that the conversion happend correctly, ensure that the models give
# identical predictions
# 
# If you are using a functional model, see this issue for how to adapt the
# code: https://github.com/kundajelab/deeplift/issues/54

# In[6]:

#make sure predictions are the same as the original model
from deeplift.util import compile_func
deeplift_prediction_func    = compile_func(
                                  [deeplift_model.get_layers()[0].get_activation_vars()],
                                  deeplift_model.get_layers()[-1].get_activation_vars())

original_model_predictions  = keras_model.predict(onehot_data, batch_size=200)

converted_model_predictions = deeplift.util.run_function_in_batches(
                                  input_data_list=[onehot_data],
                                  func=deeplift_prediction_func,
                                  batch_size=200,
                                  progress_update=None)
logging.debug("maximum difference in predictions: %e",
      np.max(np.array(converted_model_predictions)-np.array(original_model_predictions)))

assert np.max(np.array(converted_model_predictions)-np.array(original_model_predictions)) < 10**-5
predictions = converted_model_predictions
'''

# ## Compute importance scores
# 
# ### Compile the DeepLIFT contribution scoring function
# Using the deeplift model, we obtain the functions capable of computing the
# importance scores.

# In[7]:

from deeplift.util import get_shuffle_seq_ref_function
from deeplift.dinuc_shuffle import dinuc_shuffle #function to do a dinucleotide shuffle

contribs_func = deeplift_model.get_target_contribs_func(find_scores_layer_idx=0,
                                                        target_layer_idx=-2)

# ### Adapt the scoring function to work with multiple shuffled references
# 
# For each sequence, we generate a collection of reference sequences by
# shuffling the original sequence, and then average the scores over these
# multiple shuffled versions of the sequence.

# In[8]:

contribs_many_refs_func = get_shuffle_seq_ref_function(
    #score_computation_function: is the original function to compute scores
    #shuffle_func: is the function that shuffles the sequence. On real genomic
    #    data, a dinuc shuffle is advisable due to the strong bias against CG
    #    dinucleotides
    score_computation_function=contribs_func,
    shuffle_func=dinuc_shuffle,
    one_hot_func=lambda x: np.array([one_hot_encode_along_channel_axis(seq)
                                     for seq in x]))


# ### Compile the "hypothetical" contribution scoring function
# 
# Hypothetical contribution scores are estimates of the contributions that
# different bases *would* have if it were present in the sequence. They are
# found by looking at the DeepLIFT mutlipliers (which satisfy the equation
# `difference_from_reference*multiplier = contribution`) and substituting in
# what the different values for `difference_from_reference` would be if other
# bases were present in the sequence. The hypothetical contributions can act
# as "autocompletes" of the motifs by revealing the preference of the network
# for bases that are not present in the sequence. For a longer discussion, see
# https://github.com/kundajelab/tfmodisco/issues/5

# In[9]:

from deeplift.util import get_hypothetical_contribs_func_onehot
import os

multipliers_func = deeplift_model.get_target_multipliers_func(find_scores_layer_idx=0,
                                                              target_layer_idx=-2)
hypothetical_contribs_func = get_hypothetical_contribs_func_onehot(multipliers_func)

#Once again, we rely on multiple shuffled references
hypothetical_contribs_many_refs_func = get_shuffle_seq_ref_function(
    score_computation_function=hypothetical_contribs_func,
    shuffle_func=dinuc_shuffle,
    one_hot_func=lambda x: np.array([one_hot_encode_along_channel_axis(seq)
                                     for seq in x]))


# ### Obtain the scores

# In[10]:

num_refs_per_seq = 10
all_tasks = range(start_task, end_task)

log_mem_usage(1, "memory check location 3")

for task_idx in all_tasks:

    logging.debug("On task %d",task_idx)

    block_size = 10000
    hyp_scores_all = np.zeros((0,1000,4))
    num_block = int((len(fasta_sequences) + block_size - 1) / block_size )


    for block_id in range(num_block):
        st = block_id * block_size
        seq_block = fasta_sequences[st:st+block_size]
        onehot_block = np.array([one_hot_encode_along_channel_axis(seq)
                                for seq in seq_block])

        hyp_scores = hypothetical_contribs_many_refs_func(
                task_idx=task_idx,
                input_data_sequences=seq_block,
                num_refs_per_seq=num_refs_per_seq,
                batch_size=200,
                progress_update=10000,
            )

        if task_idx == -1: # no need to do this
            contrib_scores = np.sum(contribs_many_refs_func(
                    task_idx=task_idx,
                    input_data_sequences=seq_block,
                    num_refs_per_seq=num_refs_per_seq,
                    batch_size=200,
                    progress_update=10000,
                ),axis=2)[:,:,None] * onehot_block

            # ### Sanity check the hypothetical contributions
            # 
            # We make sure that the "hypothetical" contributions of the bases that are
            # actually present in the sequence are the same as the actual
            # contributions of those bases. This should be true by design, so
            # technically you don't have to compute `task_to_contrib_scores[task_idx]`
            # separately in the step above; you can just obtain it using
            # `task_to_contrib_scores[task_idx] =
            # task_to_hyp_contrib_scores[task_idx]*onehot_data`.

            max_diff = np.max(np.abs(hyp_scores * onehot_block - \
                                     contrib_scores))
            logging.debug("task %d max diff: %e",task_idx, max_diff)
            assert max_diff < 1e-6 #assert the difference is within numerical precision

        else:
            contrib_scores = hyp_scores * onehot_block

        logging.debug("task %d block %d hyp_total shape= %s",
                      task_idx, block_id, str(hyp_scores_all.shape))

        # now concatentate the scores
        hyp_scores_all = np.concatenate((hyp_scores_all, hyp_scores))

        log_mem_usage(msg="memory check location 4")


    #task_to_contrib_scores[task_idx] = contrib_scores
    #task_to_hyp_contrib_scores[task_idx] = hyp_scores

    # ### Mean-normalize the hypothetical contributions at each base
    # 
    # This is more of a personal preference. I (Av Shrikumar) find the
    # hypothetical contributions easier to digest visually when they sum up to
    # 0 at each position. I suspect the normalization improves TF-MoDISco
    # results, but I haven't rigorously evaluated this.

    #hyp_scores = (hyp_scores - np.mean(hyp_scores, axis=-1)[:,:,None])

    # save to files
    #filename = "contrib_scores_task_" + str(task_idx) + ".npy"
    #print("saving contrib_scores to " + filename)
    #np.save(filename, contrib_scores)

    os.system("mkdir -p scores")
    filename = "scores/hyp_scores_task_" + str(task_idx) + ".npy"
    logging.info("saving hyp_scores_all to " + filename + ", shape = " + str(hyp_scores_all.shape))
    np.save(filename, hyp_scores_all.astype(np.float32))

# ### Visualize the contributions and hypothetical contributions on a few sequences
# 
# It's good to do this as another sanity check. The sequences picked here match the sequences visualized in the corresponding DeepLIFT genomics notebook.

# In[13]:

# In[ ]:


end = datetime.datetime.now()
logging.debug ("start time " + str(start))
logging.debug ("end time " + str(end))
