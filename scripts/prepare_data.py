#!/usr/bin/env python

import os


"""
TFNET_ROOT is set to the root of TFNET dir
 +-- scripts       such as label_region, prepare_data.py, run_deeplift.py, run_modisco.py, run_pipeline.py etc
 +-- genome        hg19.fa hg19.chrom.sizes
 +-- ENCODE_data   intervals
 +-- results       results of the pipeline
      +-- ZNF143
      +-- CTCF
      +-- SIX5
"""
ROOT_DIR   = os.getenv('TFNET_ROOT', "../../") 
scriptDir  = ROOT_DIR + "/scripts/"
dataDir    = ROOT_DIR + "/ENCODE_data/"
genomeDir  = ROOT_DIR + "/genome/"
resultsDir = "./"
logDir     = resultsDir + "log/"
tmpDir     = "./tmp/"

positives = []
ambiguous = []

import logging
logging.basicConfig(
        format='%(asctime)s %(levelname)-5s %(message)s',
        level=logging.DEBUG,
        datefmt='%Y-%m-%d %H:%M:%S')


# must run from resultsDir

# loop through the TFs
#for tf in ['ZNF143']:

def process_tf(tf):
    header = "id"
    tmp_file_list = []
    tmp_empty_file = tmpDir + "_tmp_empty_file"
    tmp_file_list.append(tmp_empty_file)
    os.system("touch " + tmp_empty_file)
    for cell, exp, rep in [['GM12878','ENCSR000DZL',3], 
                           ['GM12878','ENCSR936XTK',2], 
                           ['H1-hESC','ENCSR000EBW',3],
                           ['HeLa-S3','ENCSR000ECO',0],
                           ['K562'   ,'ENCSR000EGP',3]]:
        positive = dataDir + cell + "-" + tf + "-human-" + exp + "-optimal_idr.narrowPeak.gz"
        positives.append(positive)
        header += "\t" + exp
        merged = tmp_empty_file
        # merged = ""
        if rep == 1: # bit0 == 1, has rep1
            rep1 = cell + "-" + tf + "-human-" + exp + "-rep1.narrowPeak.gz"
            merged = dataDir+rep1
        if rep == 2: # bit1 == 1, has rep2
            rep2 = cell + "-" + tf + "-human-" + exp + "-rep2.narrowPeak.gz"
            merged = dataDir+rep2
        if rep == 3: # bit1,bit2 == 1, has both rep1 and rep2
            rep1 = cell + "-" + tf + "-human-" + exp + "-rep1.narrowPeak.gz"
            rep2 = cell + "-" + tf + "-human-" + exp + "-rep2.narrowPeak.gz"
            tmp_merged = tmpDir + "_tmp_" + cell + "-" + tf + "-human-" + exp + "-merged.narrowPeak.gz"
            tmp_file_list.append(tmp_merged)
            merged = tmp_merged
            os.system("pigz -d -c " + dataDir+rep1 + " " + dataDir+rep2 + \
                      " | cut -f 1-3 | bedtools sort | bedtools merge | pigz -c > " + merged)
        ambiguous.append(merged)

    if positives != []:
        positives_str = " --positives " + ','.join(positives)
    else:
        positives_str = ""

    if ambiguous != []:
        ambiguous_str = " --ambiguous " + ','.join(ambiguous)
    else:
        ambiguous_str = ""
    
    background_str = " --background " + genomeDir + "hg19.tsv "

    # call tfdragonn labelregions
    #
    # labels_multitask_gz = "tflabel.intervals_file.tsv.gz"
    # cmd = "tfdragonn labelregions " + positives_str + ambiguous_str + \
    #       " --genome hg19 --prefix tflabel" #       + " --background background "

    labels_multitask_gz = "label.intervals_file.tsv.gz"
    cmd = scriptDir + "label_regions " + positives_str + ambiguous_str + \
          " --genome hg19 --prefix label " + background_str + " --stride 20"
    logging.debug(cmd)
    os.system(cmd)

    labels_multitask    = labels_multitask_gz[:-3]

    os.system("pigz -d -c " + labels_multitask_gz +  " > " + labels_multitask)

    tmp_labels_wo_title = tmpDir + "_tmp_labels_without_title.txt"

    # use sed to change each line from
    # chr<\t>start<\t>end<\t>label1<\t>label2 ...
    # to 
    # chr:start-end<\t>label1<\t>label2 ...
    os.system("cat " + labels_multitask + " | sed 's/\t/:/; s/\t/-/' > " + tmp_labels_wo_title)
    tmp_file_list.append(tmp_labels_wo_title)

    os.system("bedtools getfasta -fi " + genomeDir + "hg19.fa -bed " + labels_multitask + " -fo inputs.fa")

    #make the final inputs labels files from the shuffled lines (tfdragonn shuffles already)
    os.system("echo " + header + " > labels.txt")
    os.system("cat " + tmp_labels_wo_title + " >> labels.txt")

    #remove the intermediate files
    for tmp_file in tmp_file_list:
        os.system("rm -f " + tmp_file)

    logging.info("split and make hdf5")
    os.system("mkdir -p splits")

    #make the splits
    valid_chrom = "chr2"
    test_chrom  = "chr1"

    os.system("cat labels.txt | grep " + valid_chrom + ": | pigz -c > splits/valid.txt.gz")
    os.system("cat labels.txt | grep " + test_chrom  + ": | pigz -c > splits/test.txt.gz")
    cmd = "cat labels.txt | grep -v \"" + test_chrom  + ":\|"+ valid_chrom + ":\|^id" + "\" | pigz -c > splits/train.txt.gz"
    os.system(cmd)

    os.system("pigz -f labels.txt")
    #os.system("pigz -f inputs.fa")

    os.system("make_hdf5 --yaml_configs make_hdf5_yaml/* --output_dir .")

    logging.info("prepare_data done")

# guarantee to clean up tmp dir
import contextlib
import tempfile
import shutil
@contextlib.contextmanager
def make_temp_directory():
    temp_dir = tempfile.mkdtemp(dir = ".", prefix = "_tmp_")
    try:
        yield temp_dir
    finally:
        shutil.rmtree(temp_dir)

if __name__ == '__main__':

    with make_temp_directory() as temp_dir:
        global tmpDir
        tmpDir = temp_dir + "/"
        process_tf('ZNF143')

