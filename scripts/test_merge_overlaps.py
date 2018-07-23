import numpy as np
from merge_overlaps import MergeOverlaps

# ### testing code for MergeOverlaps
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

    print(tsv_list)
    print(hyp_scores_all)
    merged_hyp_scores_list = []
    merged_seq_list        = []
    print("gen_tests DONE")

    merged = MergeOverlaps(merged_hyp_scores_list, merged_seq_list=merged_seq_list,
                           core_size=8)
    for i in range(len(tsv_list)):
        chrom, st, en = tsv_list[i]
        hyp_scores = hyp_scores_all[i]
        seq        = seq_list[i]
        merged.process_one_interval(chrom, st, en, hyp_scores, seq)
    merged.close_interval()

    #print(tsv_list)
    for i in range(len(merged_hyp_scores_list)) :
        ar = merged_hyp_scores_list[i]
        print("i=%d shape of array = %s" % (i, str(ar.shape)))
        print(ar)
        seq = merged_seq_list[i]
        print(seq)

    print("average score DONE")

do_test()
