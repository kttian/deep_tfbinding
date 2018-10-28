#!/usr/bin/env python


import sys
import logging
import os
import argparse

logging.basicConfig(
        format='%(asctime)s %(levelname)-5s %(message)s',
        level=logging.DEBUG,
        datefmt='%Y-%m-%d %H:%M:%S')


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

positives = []
ambiguous = []

use_hdf5=False

# must run from resultsDir

# loop through the TFs
#for tf in ['ZNF143']:

def process_tf(tfs, cell_set=None, expr=None, test_chroms="chr1,chr2"):
            #           -4   -3    -2          -1
    #neutrophil-CTCF-human-ENCSR785YRL-optimal_idr.narrowPeak.gz
    #neutrophil-CTCF-human-ENCSR785YRL-rep1.narrowPeak.gz
    #neutrophil-CTCF-human-ENCSR785YRL-rep2.narrowPeak.gz

    if expr == None:
        expr_str = ""
    else:
        expr_str = expr

    import glob
    tf_files = []
    for tf in tfs:
        tf_files.extend(glob.glob(dataDir + "*-" + tf + "-human-" + expr_str + "*-optimal*"))

    count = 0
    task_list = []
    for path_name in tf_files:
        fn = os.path.basename(path_name)
        fn_list = fn.split('-')
        exp  = fn_list[-2]
        tf   = fn_list[-4]
        cell = '-'.join(fn_list[:-4])
        if cell_set != None and len(cell_set) != 0: # select cells only in the specified set
            if not cell in cell_set:
                continue
        task_list.append([cell, tf, exp])
        print(path_name)
        count = count + 1

    #print(task_list)
    print(count)

    header = "id"
    tmp_file_list = []
    tmp_empty_file = tmpDir + "_tmp_empty_file"
    tmp_file_list.append(tmp_empty_file)
    os.system("touch " + tmp_empty_file)
    has_ambiguous = False
    for cell, tf, exp in task_list:
        positive = dataDir + cell + "-" + tf + "-human-" + exp + "-optimal_idr.narrowPeak.gz"
        if not os.path.isfile(positive):
            print("ERR does not exist: ", positive)
        positives.append(positive)
        header += "\t" + exp
        merged = tmp_empty_file
        # merged = ""
        rep1 = cell + "-" + tf + "-human-" + exp + "-rep1.narrowPeak.gz"
        rep2 = cell + "-" + tf + "-human-" + exp + "-rep2.narrowPeak.gz"

        has_rep1 = os.path.isfile(dataDir + rep1)
        has_rep2 = os.path.isfile(dataDir + rep2)

        if has_rep1:
            merged = dataDir+rep1
            has_ambiguous = True
        if has_rep2:
            merged = dataDir+rep2
            has_ambiguous = True
        if has_rep1 and has_rep2:
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

    if has_ambiguous:
        ambiguous_str = " --ambiguous " + ','.join(ambiguous)
    else:
        #ambiguous_str = ""
        ambiguous_str = " --ambiguous " + ','.join(ambiguous)
    
    background_str = " --background " + genomeDir + "hg19.tsv "

    # call tfdragonn labelregions
    #
    # labels_multitask_gz = "tflabel.intervals_file.tsv.gz"
    # cmd = "tfdragonn labelregions " + positives_str + ambiguous_str + \
    #       " --genome hg19 --prefix tflabel" #       + " --background background "

    labels_multitask_gz = "label.intervals_file.tsv.gz"
    cmd = scriptDir + "label_regions " + positives_str + ambiguous_str + \
          " --genome hg19 --prefix label " + " --stride " + str(args.stride) + " --test-chroms " + args.test_chroms

    if args.bg_stride >= 0:
        cmd += " --bg-stride " + str(args.bg_stride) + " "

    if args.no_bg == False:
        cmd += background_str

    if args.test_only:
        cmd += " --test-only True "

    logging.debug(cmd)
    os.system(cmd)

    labels_multitask    = labels_multitask_gz[:-3]

    os.system("pigz -d -c " + labels_multitask_gz +  " > " + labels_multitask)

    #make the splits
    os.system("mkdir -p splits")
    test_chroms_list = test_chroms.split(",")
    valid_chrom = "'" + test_chroms_list[-1] + "\t'"
    if len(test_chroms_list) > 2: # more than one test chrom
        test_chrom  = "'" + test_chroms_list[0] + "\t|" + test_chroms_list[1] + "\t'"
    else:
        test_chrom  = "'" + test_chroms_list[0] + "\t'"
    os.system("cat " + labels_multitask + " | grep -P " + valid_chrom + " | pigz -c > splits/valid.tsv.gz")
    os.system("cat " + labels_multitask + " | grep -P " + test_chrom  + " | pigz -c > splits/test.tsv.gz")
    cmd =     "cat " + labels_multitask + " | grep -v -P " + test_chrom + " | grep -v -P " + valid_chrom + " | grep -v -P '^id' | pigz -c > splits/train.tsv.gz"
    os.system(cmd)

    os.system("echo '" + header + "' > headers.txt")
    if args.hdf5 != True:
        logging.info("prepare_data done")
        return

    # following preparations are for hdf5
    tmp_labels_wo_title = tmpDir + "_tmp_labels_without_title.txt"
    # use sed to change each line from
    # chr<\t>start<\t>end<\t>label1<\t>label2 ...
    # to 
    # chr:start-end<\t>label1<\t>label2 ...
    os.system("cat " + labels_multitask + " | sed 's/\t/:/; s/\t/-/' > " + tmp_labels_wo_title)
    tmp_file_list.append(tmp_labels_wo_title)

    os.system("bedtools getfasta -fi " + genomeDir + "hg19.fa -bed " + labels_multitask + " -fo inputs.fa")

    #make the final inputs labels files from the shuffled lines (tfdragonn shuffles already)
    os.system("echo '" + header + "' > labels.txt")
    os.system("cat " + tmp_labels_wo_title + " >> labels.txt")

    #remove the intermediate files
    for tmp_file in tmp_file_list:
        os.system("rm -f " + tmp_file)

    logging.info("split and make hdf5")
    os.system("mkdir -p splits")

    #make the splits
    test_chroms_list = test_chroms.split(",")
    valid_chrom = "'" + test_chroms_list[-1] + ":'"
    if len(test_chroms_list) > 2: # more than one test chrom
        test_chrom  = "'" + test_chroms_list[0] + ":|" + test_chroms_list[1] + ":'"
    else:
        test_chrom  = "'" + test_chroms_list[0] + ":'"
    os.system("cat labels.txt | grep -P " + valid_chrom + " | pigz -c > splits/valid.txt.gz")
    os.system("cat labels.txt | grep -P " + test_chrom  + " | pigz -c > splits/test.txt.gz")
    cmd =     "cat labels.txt | grep -v -P " + test_chrom + " | grep -v -P " + valid_chrom + " | grep -v -P '^id' | pigz -c > splits/train.txt.gz"
    print(cmd)
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

def parse_args(args = None):
    parser = argparse.ArgumentParser('prepare_data_pf.py',
                                     description='prepare data for TF binding training',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--tfs', type=str, default=None, help="List of transcription factors, separated by ','")
    parser.add_argument('--cells', type=str, default=None, help="List of cell-lines, separated by ','")
    parser.add_argument('--no-bg', type=bool, default=False, help="No background regions")
    parser.add_argument('--bg-stride', type=int, default=-1, help="Stride in background regions")
    parser.add_argument('--test-only', type=bool, default=False, help="only prepare data for test and validation")
    parser.add_argument('--hdf5', type=bool, default=False, help="produce hdf5 as well")
    parser.add_argument('--data-dir', type=str, default=None, help="DataDir")
    parser.add_argument('--stride', type=int, default=10, help="stride in positive regions")
    parser.add_argument('--expr', type=str, default=None, help="Experiment Id")
    parser.add_argument('--test-chroms', type=str, default="chr1,chr2", help="heldout the comma separated list of chroms for test and validation")
    args = parser.parse_args(args)
    return args

if __name__ == '__main__':

    args = parse_args()

    if args.data_dir != None:
        dataDir = args.data_dir + "/"
        logging.info("dataDir="+dataDir)

    with make_temp_directory() as temp_dir:
        global tmpDir
        tmpDir = temp_dir + "/"

        tfs = args.tfs.split(',')
        if args.cells == '-' or args.cells == None:
            cell_set = None
        else:
            cell_lines = args.cells.split(',')
            cell_set = set(cell_lines)

        process_tf(tfs, cell_set, args.expr, args.test_chroms)

