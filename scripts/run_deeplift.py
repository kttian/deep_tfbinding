from __future__ import print_function

import datetime
start = datetime.datetime.now()
print ("start time " + str(start))

# ### Load the keras model
import deeplift
import deeplift.conversion.kerasapi_conversion as kc
from keras.models import model_from_json

#load the keras model
    
#keras_model_weights = "../results/model_files/record_2_model_9BnmY_modelWeights.h5"
#keras_model_json = "../results/model_files/record_2_model_9BnmY_modelJson.json"
#input_file = "test.fa"

import sys
if len(sys.argv) != 3:
    print("Syntax: ", sys.argv[0] , " <model name> <sequence file>")
    quit()

keras_model_weights = sys.argv[1] + "Weights.h5"
keras_model_json    = sys.argv[1] + "Json.json"
input_file          = sys.argv[2]

print("loading models from ", keras_model_json, " ", keras_model_weights)
keras_model = model_from_json(open(keras_model_json).read())
keras_model.load_weights(keras_model_weights)

print("input sequence file is ", input_file)

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

sequences = []
fasta = fasta_iter(input_file)
i = 0
for header, seq in fasta:   
    if i % 10 == 0
        sequences.append(seq)
    i = i+1
print(len(sequences))

import numpy as np

#this model was trained on data one-hot encoded as a 2d image, with the row-axis being the axis
#for one-hot encoding.
def one_hot_encode_along_row_axis(sequence):
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
    #zeros_array should be an array of dim 4xlen(sequence), filled with zeros.
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
            
onehot_data = np.array([one_hot_encode_along_row_axis(seq) for seq in sequences])
#print(onehot_data.shape())

# ## Prepare the deeplift models
# 
# ### Model conversion
# 
# Convert the keras models to deeplift models capable of computing importance scores using DeepLIFT (with 3 different variants: rescale on the conv layers and revealcancel on the fully-connected layers (the genomics default), rescale on all layers, and revealcancel on all layers), gradients and guided backprop

from deeplift.layers import NonlinearMxtsMode
reload(deeplift.layers)
reload(deeplift.conversion.kerasapi_conversion)
from collections import OrderedDict

method_to_model = OrderedDict()
for method_name, nonlinear_mxts_mode in [
    #The genomics default = rescale on conv layers, revealcance on fully-connected
    ('rescale_conv_revealcancel_fc', NonlinearMxtsMode.DeepLIFT_GenomicsDefault)]:

    method_to_model[method_name] = kc.convert_model_from_saved_files(
        h5_file=keras_model_weights,
        json_file=keras_model_json,
        nonlinear_mxts_mode=nonlinear_mxts_mode)


# ### Sanity checks
# To ensure that the conversion happend correctly, ensure that the models give identical predictions

#make sure predictions are the same as the original model
from deeplift.util import compile_func
model_to_test = method_to_model['rescale_conv_revealcancel_fc']
deeplift_prediction_func = compile_func([model_to_test.get_layers()[0].get_activation_vars()],
                                         model_to_test.get_layers()[-1].get_activation_vars())
original_model_predictions = keras_model.predict(onehot_data, batch_size=200)
converted_model_predictions = deeplift.util.run_function_in_batches(
                                input_data_list=[onehot_data],
                                func=deeplift_prediction_func,
                                batch_size=200,
                                progress_update=None)
print("maximum difference in predictions:",np.max(np.array(converted_model_predictions)-np.array(original_model_predictions)))
assert np.max(np.array(converted_model_predictions)-np.array(original_model_predictions)) < 10**-5
predictions = converted_model_predictions


# ## Compute importance scores
# 
# ### Compile various scoring functions
# Using the deeplift models, we obtain the functions capable of computing the importance scores.

print("Compiling scoring functions")
method_to_scoring_func = OrderedDict()
for method,model in method_to_model.items():
    print("Compiling scoring function for: "+method)
    method_to_scoring_func[method] = model.get_target_contribs_func(find_scores_layer_idx=0,
                                                                    target_layer_idx=-2)

# ### Call scoring functions on the data
# 
# In the cell below, a reference representing 40% GC content is used

background = OrderedDict([('A', 0.3), ('C', 0.2), ('G', 0.2), ('T', 0.3)])

from collections import OrderedDict
method_to_task_to_scores = OrderedDict()
'''
for method_name, score_func in method_to_scoring_func.items():
    print("on method",method_name)
    method_to_task_to_scores[method_name] = OrderedDict()
    for task_idx in [0,1,2,3,4]:
        scores = np.array(score_func(
                    task_idx=task_idx,
                    input_data_list=[onehot_data],
                    input_references_list=[
                     np.array([background['A'],
                               background['C'],
                               background['G'],
                               background['T']])[None,None,:]],
                    batch_size=200,
                    progress_update=None))
        assert scores.shape[2]==4
        scores = np.sum(scores, axis=2)
        method_to_task_to_scores[method_name][task_idx] = scores
'''


# ## Using multiple shuffled references
# 
# As an alternative to using a flat reference based on GC content (which can sometimes produce artefacts), we propose averaging the scores produced using mutliple references which are produced by shuffling the original sequence. We find in practice that this can give more robust results. Note that in general, the optimal choice of reference is an area of active research.


reload(deeplift.util)
from deeplift.util import get_shuffle_seq_ref_function
#from deeplift.util import randomly_shuffle_seq
from deeplift.dinuc_shuffle import dinuc_shuffle #function to do a dinucleotide shuffle

rescale_conv_revealcancel_fc_many_refs_func = get_shuffle_seq_ref_function(
    #score_computation_function is the original function to compute scores
    score_computation_function=method_to_scoring_func['rescale_conv_revealcancel_fc'],
    #shuffle_func is the function that shuffles the sequence
    #technically, given the background of this simulation, randomly_shuffle_seq
    #makes more sense. However, on real data, a dinuc shuffle is advisable due to
    #the strong bias against CG dinucleotides
    shuffle_func=dinuc_shuffle,
    one_hot_func=lambda x: np.array([one_hot_encode_along_row_axis(seq) for seq in x]))


num_refs_per_seq=10 #number of references to generate per sequence
method_to_task_to_scores['rescale_conv_revealcancel_fc_multiref_'+str(num_refs_per_seq)] = OrderedDict()
for task_idx in [0,1,2,3,4]:
#for task_idx in [0]:
    method_to_task_to_scores['rescale_conv_revealcancel_fc_multiref_'+str(num_refs_per_seq)][task_idx] =        np.squeeze(np.sum(rescale_conv_revealcancel_fc_many_refs_func(
            task_idx=task_idx,
            input_data_sequences=sequences,
            num_refs_per_seq=num_refs_per_seq,
            batch_size=200,
            progress_update=10000,
        ),axis=2))
    for method_name in [
                        #'rescale_conv_revealcancel_fc',
                        'rescale_conv_revealcancel_fc_multiref_10'
                        ]:
        scores = method_to_task_to_scores[method_name][task_idx]
        filename = method_name + "_task_" + str(task_idx) + ".npy"
        print("saving to " + filename)
        np.save(filename, scores)

end = datetime.datetime.now()
print ("start time " + str(start))
print ("end time " + str(end))

