#sample the genome randomly to get candidate negatives

#get 200000 random samples from hg19
python /users/annashch/anna_utils/seq_utils/sample_genome_randomly.py --region_n 200000 --outf candidate_negatives.bed --main_chroms_only

#select candidate negatives so that they don't overlap the positives.
bedtools intersect -v -a candidate_negatives.bed -b HepG2-JUND-summit.bed > confirmed_negatives.bed

#assign positive and negative labels
python /users/annashch/anna_utils/seq_utils/add_fixed_labels_to_bed.py --bed HepG2-JUND-summit.bed --labels 1 --outf positives.bed
python /users/annashch/anna_utils/seq_utils/add_fixed_labels_to_bed.py --bed confirmed_negatives.bed --labels 0 --outf negatives.bed
#concatenative positives and negatives into a single file 
cat positives.bed negatives.bed > combined.bed

#run the dinucleotide matching code on "combined.bed" -- it has labels of "1" for positives and "0" for negatives.
#in the first pass, generate the matrix of dinucleotide frequencies. 
python gen_dinucleotide_freqs.py --bed_path combined.bed --ratio_neg_to_pos 2.45446958912  --outf combined.freqs --ref_fasta /users/jocelins/hg19.fa

#now that the frequency matrix has been generated,  subselect negatives 
python /srv/scratch/annashch/deeplearning/form_inputs/code/gc_dinuc_balanced/gen_dinucleotide_freqs.py --bed_path combined.bed --ratio_neg_to_pos 2.45446958912  --outf JUND_dinuc.bed --ref_fasta /users/jocelins/hg19.fa --dinuc_freqs combined.freqs

