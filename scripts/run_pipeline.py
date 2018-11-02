#!/usr/bin/env python
# to run the modisco pipeline

import logging
import sys
import os
import argparse

logging.basicConfig(
        format='%(asctime)s %(levelname)-5s %(message)s',
        level=logging.DEBUG,
        datefmt='%Y-%m-%d %H:%M:%S')

logging.info(" ".join(sys.argv))

ROOT_DIR    = os.getenv('TFNET_ROOT', "../../") 
genomeDir   = ROOT_DIR + "/genome/"
resultDir   = ROOT_DIR + "/results/"
templateDir = resultDir + "/templates/"

num_task = 5
end = 100

def parse_args(args = None):
    parser = argparse.ArgumentParser('run_pipeline.py',
                                     description='run pipe line for TF binding training',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--tfs', type=str, help="List of transcription factors, separated by ','")
    parser.add_argument('--cells', type=str, default=None, help="List of cell-lines, separated by ','")
    parser.add_argument('--start', type=int, default=20, help="start stage")
    parser.add_argument('--end', type=int, default=100, help="end stage")
    parser.add_argument('--start-task', type=int, default=0, help="start tast")
    parser.add_argument('--end-task', type=int, default=5, help="end task")
    parser.add_argument('--fdr', type=float, default=0.01, help="target FDR")
    parser.add_argument('--data-dir', type=str, default=None, help="DataDir")
    parser.add_argument('--stride', type=int, default=10, help="stride")
    parser.add_argument('--expr', type=str, default=None, help="Experiment Id")
    parser.add_argument('--pf', type=bool, default=False, help="Use pfasta")
    parser.add_argument('--min-seqlets', type=int, default=-1, help="min seqlets")
    parser.add_argument('--test-chroms', type=str, default=None, help="heldout the comma separated list of chroms for test and validation")
    args = parser.parse_args(args)
    return args

args = parse_args()

tfs        = args.tfs
cell_lines = args.cells
start      = args.start
end        = args.end
start_task = args.start_task
end_task   = args.end_task
num_tasks  = end_task - start_task

'''
if len(sys.argv) > 6:
    print("Syntax: ", sys.argv[0] , " <TF names> <cell_lines> <number of tasks> [start]")
    quit()

tf = sys.argv[1]
cell_lines = sys.argv[2]
num_tasks = int(sys.argv[3])
if num_tasks == 0:
    print("Syntax: ", sys.argv[0] , " <TF names> <cell_lines> <number of tasks> [start]")
    quit()

if len(sys.argv) == 5:
    start = int(sys.argv[4])
elif len(sys.argv) == 6:
    start = int(sys.argv[4])
    end   = int(sys.argv[5])
else:
    start = 20 # start with pre-train
'''

os.system("mkdir -p logs")

logging.debug("tf=%s, num_tasks=%d start=%d end=%d" % (tfs, num_tasks, start, end))

if args.data_dir != None:
    data_dir_str = " --data-dir " + args.data_dir + " "
else:
    data_dir_str = ""

if args.expr != None:
    expr_str = " --expr " + args.expr + " "
else:
    expr_str = ""

if args.cells != None:
    cell_str = " --cells " + args.cells + " "
else:
    cell_str = ""

if args.pf :
    hdf5_str = ""
    hdf5_suffix = ""
else:
    hdf5_str = " --hdf5 True "
    hdf5_suffix = ".h5"

if args.min_seqlets < 0 :
    min_seqlets_str = ""
else:
    min_seqlets_str = " --min-seqlets " + str(args.min_seqlets) + " "

if args.test_chroms != None:
    test_chroms_str = " --test-chroms " + args.test_chroms + " "
else:
    test_chroms_str = ""

#-------------------------------
if start <= 10:
#0 prepare_data with union of positives (no background) and train a model
    os.system("cp -r " + templateDir + "/config" + hdf5_suffix + " .")
    if hdf5_suffix != "":
        os.system("mv config" + hdf5_suffix + " config")
    os.system("ln -s " + templateDir + "/make_hdf5_yaml .")
    #os.system("cp -f config/hyperparameter_configs_list_for_pretrain.yaml config/hyperparameter_configs_list.yaml")
    os.system("rm config/hyperparameter_configs_list.yaml")
    os.system("cat config/hyperparameter_configs_list_for_pretrain.yaml | sed -e 's/_MESSAGE_/" + tfs + "/g; s/_OUT_DIM_/" + str(num_tasks) + "/g' > config/hyperparameter_configs_list.yaml")
    #sys.exit()

#-------------------------------
if start <= 20 and end > 20:
    cmd = "python $TFNET_ROOT/scripts/prepare_data_pf.py --tfs " + tfs + cell_str + expr_str + data_dir_str + " --no-bg True --stride " + str(args.stride) + test_chroms_str
    if num_tasks >1 :
        cmd += hdf5_str # for single task, there is no need to pre-train, but we need the tsv for interpretation
    cmd += " > logs/pre_prepare.txt 2>&1"
    logging.info(cmd)
    os.system(cmd)

    #os.system("gunzip -c splits/train.tsv.gz splits/valid.tsv.gz | sort -k1,1 -k2,2n > interpret.tsv")

    print("step 20 prepare for pre_train done")

#-------------------------------
if start <= 30 and end > 30:
    if num_tasks > 1:
        rc = os.system("momma_dragonn_train > logs/pre_train.txt 2>&1")
        print("momma_dragonn_train returned ", rc)
        if rc != 0:
            sys.exit()

        os.chdir("model_files")
        os.system("cp record_1_*Json.json  record_1_Json.json")
        os.system("cp record_1_*Weights.h5 record_1_Weights.h5")
        os.chdir("..")
        os.system("mkdir -p pretrain")
        os.system("mv model_files pretrain")
        os.system("mv -f *.hdf5 pretrain")
        os.system("mv splits pretrain")
        os.system("mv runs_perf-metric-auROC.db pretrain/")
        os.system("mv label* pretrain/")
        os.system("mv _tmp_* pretrain/")

        os.system("cat config/hyperparameter_configs_list_for_finetune.yaml | sed -e 's/_MESSAGE_/" + tfs + "/g; s/_OUT_DIM_/" + str(num_tasks) + "/g' > config/hyperparameter_configs_list.yaml")
        print("step 30 pre_train done")

#1 prepare_data with background for the main training
#-------------------------------
if start <= 40 and end > 40:
    cmd = "python $TFNET_ROOT/scripts/prepare_data_pf.py --tfs " + tfs + cell_str + expr_str + data_dir_str + " --stride " + str(args.stride) + test_chroms_str
    cmd += hdf5_str + " > logs/prepare.txt 2>&1"
    os.system(cmd)
    print("step 40 prepare_data done")

#2 train to continue from pre-trained data
#-------------------------------
if start <= 50 and end > 50:
    rc = os.system("momma_dragonn_train > logs/train.txt 2>&1")
    print("momma_dragonn_train returned ", rc)
    if rc != 0:
        sys.exit()

if start <= 52 and end > 52:
    os.chdir("model_files")
    os.system("cp record_1_*Json.json  record_1_Json.json")
    os.system("cp record_1_*Weights.h5 record_1_Weights.h5")
    os.chdir("..")

    os.system("mkdir -p finetune")
    os.system("cp -rp model_files finetune")
    os.system("mv -f *.hdf5 finetune")
    os.system("mv splits finetune")
    os.system("mv runs_perf-metric-auROC.db finetune/")
    os.system("mv label* finetune/")
    os.system("mv _tmp_* finetune/")

    print("step 50 training done")

#-------------------------------
if start <= 55 and end > 55:
    if args.test_chroms != None:
        test_chroms = args.test_chroms
    else:
        test_chroms = "chr1,chr2"
    test_chroms_list = test_chroms.split(",")
    if len(test_chroms_list) > 2: # more than one test chrom
        test_chrom  = "'" + test_chroms_list[0] + "\t|" + test_chroms_list[1] + "\t'"
    else:
        test_chrom  = "'" + test_chroms_list[0] + "\t'"
    cmd = "python $TFNET_ROOT/scripts/pick_summit.py --tfs " + tfs + cell_str + expr_str + data_dir_str
    cmd += " | grep -v -P " + test_chrom + " | sort -k1,1 -k2,2n > interpret.tsv"
    os.system(cmd)

    os.system("bedtools getfasta -fi " + genomeDir + "hg19.fa -bed interpret.tsv -fo interpret.fa")

#-------------------------------
if start <= 60 and end > 60:
    cmd = "python $TFNET_ROOT/scripts/prepare_data_pf.py --tfs " + tfs + cell_str + expr_str + data_dir_str + test_chroms_str + " --test-only True --bg-stride=50 > logs/test.txt 2>&1"
    os.system(cmd)

if start <= 65 and end > 65:
    cmd = "python $TFNET_ROOT/../momma_dragonn/scripts/momma_dragonn_eval --valid_data_loader_config config/test_data_loader_config_pf.yaml >> logs/test.txt 2>&1"
    os.system(cmd)
    os.system("python $TFNET_ROOT/scripts/analyze_data.py > logs/analyze.txt")
    os.system("mv counts.tsv model_files")
    os.system("paste model_files/counts.tsv model_files/record_1_Tsv.tsv > TF.tsv")
    os.system("cp -rp model_files finetune")

    print("step 60 testing done")


#3 deeplift
#-------------------------------
if start <= 70 and end > 70:
    os.system("python $TFNET_ROOT/scripts/run_deeplift.py model_files/record_1_ interpret.fa " + str(num_tasks) + " > logs/deeplift.txt 2>&1")
    print("step 70 deeplift done")

#4 modisco
#-------------------------------
if start <= 80 and end > 80:
    os.system("python $TFNET_ROOT/scripts/run_tfmodisco.py --scores scores/hyp_scores_task_ --fasta interpret.fa --tsv interpret.tsv " + \
        " --start-task " + str(args.start_task) + " --end-task " + str(args.end_task) +  " --fdr " + str(args.fdr) +  min_seqlets_str + " > logs/modisco.txt 2>&1")

"""
'''
Input for prepare_data.py
-------
list of positive  peaks for each task
list of ambiguous peaks for each task
fasta of the whole genome
sizes of the chromosomes

Step 0: prepare_data.py and pretrain
Step 1: prepare_data.py and train
-----------------------------------------------------------------------
'''
# run prepare_data.py (executable), which calls label_regions script to create labeled bins

'''
Outputs of prepare_data.py, inputs for momma_dragonn_train
-------
inputs.fa           # the input sequences for momma_dragonn_train
train.hdf5 # training set
valid.hdf5 # validation set
test.hdf5  # testing set


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

