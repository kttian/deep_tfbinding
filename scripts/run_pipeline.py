#!/usr/bin/env python
# to run the modisco pipeline

import logging
logging.basicConfig(
        format='%(asctime)s %(levelname)-5s %(message)s',
        level=logging.DEBUG,
        datefmt='%Y-%m-%d %H:%M:%S')

import sys
logging.info(" ".join(sys.argv))

import os

ROOT_DIR    = os.getenv('TFNET_ROOT', "../../") 
genomeDir   = ROOT_DIR + "/genome/"
resultDir   = ROOT_DIR + "/results/"
templateDir = resultDir + "/templates/"

num_task = 3

import sys
if len(sys.argv) > 2:
    print("Syntax: ", sys.argv[0] , " [start]")
    quit()

if len(sys.argv) == 2:
    start = int(sys.argv[1])
else:
    start = 1 # start with normal training

os.system("mkdir -p logs")

if start <= -1:
#0 prepare_data with union of positives (no background) and train a model
    os.system("cp -r " + templateDir + "/config .")
    os.system("ln -s " + templateDir + "/make_hdf5_yaml .")
    os.system("cp -f config/hyperparameter_configs_list.yaml.from_scratch config/hyperparameter_configs_list.yaml")
    quit()

if start <= 0:
    os.system("python prepare_data.py --no-bg > logs/pre_prepare.log 2>&1")
    os.system("momma_dragonn_train > logs/pre_train.log 2>&1")

    os.system("mkdir -p pretrain")
    os.chdir("model_files")
    os.system("ln -s record_1_*Json.json  record_1_Json.json")
    os.system("ln -s record_1_*Weights.h5 record_1_Weights.h5")
    os.chdir("..")
    os.system("mv model_files pretrain")

    os.system("gunzip -c splits/test.txt.gz | sed 's/:/\t/; s/-/\t/' | sort -k1,1 -k2,2n > subset_nobg.tsv")
    os.system("bedtools getfasta -fi " + genomeDir + "hg19.fa -bed subset_nobg.tsv -fo subset_nobg.fa")
    os.system("mv splits pretrain")

    os.system("mv runs_perf-metric-auROC.db pretrain/")
    os.system("mv label* pretrain/")
    os.system("mv _tmp_* pretrain/")

    os.system("cp -f config/hyperparameter_configs_list.yaml.pretrain config/hyperparameter_configs_list.yaml")
    print("step 0 pre_train done")

#1 prepare_data with background for the main training
if start <= 1:
    os.system("python prepare_data.py > logs/prepare.log 2>&1")
    print("step 1 prepare_data done")

#2 train to continue from pre-trained data
if start <= 2:
    os.system("momma_dragonn_train > logs/train.log 2>&1")
    os.chdir("model_files")
    os.system("ln -s record_1_*Json.json  record_1_Json.json")
    os.system("ln -s record_1_*Weights.h5 record_1_Weights.h5")
    os.chdir("..")
    print("step 2 training done")

#3 deeplift
if start <= 3:
    # use the test set as the subset for deeplift
    # os.system("gunzip -c splits/test.txt.gz | perl -lane 'if ($.%2==1) {print}' | sed 's/:/\t/; s/-/\t/' | sort -k1,1 -k2,2n > splits/subset.tsv") # select half of testset
    #os.system("gunzip -c splits/test.txt.gz | sed 's/:/\t/; s/-/\t/' | sort -k1,1 -k2,2n > splits/subset.tsv")
    #os.system("bedtools getfasta -fi " + genomeDir + "hg19.fa -bed splits/subset.tsv -fo subset.fa")

    os.system("python $TFNET_ROOT/scripts/run_deeplift.py model_files/record_1_ subset_nobg.fa 5 > logs/deeplift.log 2>&1")
    print("step 3 deeplift done")

#4 modisco
if start <= 4:
    os.system("python $TFNET_ROOT/scripts/run_tfmodisco.py scores/hyp_scores_task_ subset_nobg.fa subset_nobg.tsv 5 > logs/modisco.log 2>&1")

"""
'''
Input for prepare_data.py
-------
list of positive  peaks for each task
list of ambiguous peaks for each task
fasta of the whole genome
sizes of the chromosomes

Step 1: prepare_data.py
-----------------------------------------------------------------------
'''
# run prepare_data.py (executable), which calls label_regions script to create labeled bins

'''
Outputs of prepare_data.py, inputs for momma_dragonn_train
-------
inputs.fa           # the input sequences for momma_dragonn_train
splits/train.txt.gz # training set
splits/valid.txt.gz # validation set
splits/test.txt.gz  # testing set


Step 2: momma_dragonn_train
-----------------------------------------------------------------------
'''
# train momma dragonn model

# re-name momma dragonn model files

'''
Outputs of momma_dragonn_train (after rename)
------
model_files/record_1_Json.json  # model architecture
model_files/record_1_Weights.h5 # model weights

Additional input for deeplift
------
subset.fa      # fasta file for the sequences. either 10% of inputs.fa, or test set (chr1)

'''

'''
Step 3: run_deeplift.py
-----------------------------------------------------------------------
'''
# run deeplift

'''
Output from deeplift
------
importance scores
scores/hyp_scores_task_0.npy
...
scores/hyp_scores_task_N.npy
where N= number of tasks

Additional input for modisco:
------
subset.txt    # sequences corresponding to subset.fa, excluding the fasta header lines with '>'

Step 4: run_tfmodisco.py
-----------------------------------------------------------------------
'''
# run tf modisco

"""

