#!/usr/bin/env python


import logging
logging.basicConfig(
        format='%(asctime)s %(levelname)-5s %(message)s',
        level=logging.DEBUG,
        datefmt='%Y-%m-%d %H:%M:%S')


import sys
import os
logging.info(" ".join(sys.argv))



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
#tmpDir     = "./tmp/"

import sys
if len(sys.argv) < 2:
    print("Syntax: ", sys.argv[0] , " <TF name> <cell-lines>")
    quit()

tf = sys.argv[1]

cell_lines = sys.argv[2:]

# must run from resultsDir

# loop through the TFs
#for tf in ['ZNF143']:

import gzip
def process_files(in_names, out_name, bin_size):
    with open(out_name, 'w') as tsvout:
        print("out_name= " + out_name)
        for in_name in in_names:
            with gzip.open(in_name,'r') as tsvin:
    
                for cnt, line in enumerate(tsvin):
                
                    fields = line.split('\t')
                    if len(fields) < 10:
                        continue

                    chrom = fields[0]
                    start = int(fields[1])
                    end   = int(fields[2])
                    peak  = int(fields[9])
                    left  = start + peak - int(bin_size/2)

                    tsvout.write(chrom + "\t" + str(left) + "\t" + str(left + bin_size) + "\n")
    
        
def process_tf(tf, cells_set):

            #           -4   -3    -2          -1
    #neutrophil-CTCF-human-ENCSR785YRL-optimal_idr.narrowPeak.gz
    #neutrophil-CTCF-human-ENCSR785YRL-rep1.narrowPeak.gz
    #neutrophil-CTCF-human-ENCSR785YRL-rep2.narrowPeak.gz

    import glob
    tf_files = glob.glob(dataDir + "*-" + tf + "-human-*-optimal*")

    count = 0
    task_list = []
    for path_name in tf_files:
        fn = os.path.basename(path_name)
        fn_list = fn.split('-')
        exp  = fn_list[-2]
        tf   = fn_list[-4]
        cell = '-'.join(fn_list[:-4])
        if cells_set != None and len(cell_set) != 0: # select cells only in the specified set
            if not cell in cells_set:
                continue
        task_list.append([cell, tf, exp])
        print(path_name)
        count = count + 1

    #print(task_list)
    print(count)

    positives = []
    for cell, tf, exp in task_list:

        positive = dataDir + cell + "-" + tf + "-human-" + exp + "-optimal_idr.narrowPeak.gz"
        if not os.path.isfile(positive):
            print("ERR does not exist: ", positive)
        positives.append(positive)

    out_name = tf + "_summit.tsv"
    process_files(positives, out_name, 1000)


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

        cell_set = set(cell_lines)
        process_tf(tf, cell_set)

