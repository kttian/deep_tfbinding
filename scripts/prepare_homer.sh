ratio=1.5
homer_dir="homer/"
#SCRIPT_PATH=/srv/scratch/annashch/deeplearning/form_inputs/code/gc_dinuc_balanced/
SCRIPT_PATH="$TFNET_ROOT/scripts/"


orig_positives=label.intervals_file.tsv
positives=$homer_dir/positives.bed
combined=$homer_dir/combined.bed
dinuc0=$homer_dir/dinuc.bed.0
dinuc=$homer_dir/dinuc.bed.0.shuf

#sample the genome randomly to get candidate negatives

mkdir -p $homer_dir/homer_out/

$TFNET_ROOT/scripts/prepare_data_pf.py --tfs SPI1 --cells GM12878 --no-bg True --stride 20
python $TFNET_ROOT/scripts/pick_summit.py --tfs SPI1 --cells GM12878  | grep -v -P 'chr1\t' | bedtools sort > interpret.tsv
bedtools slop -i interpret.tsv -g $TFNET_ROOT/genome/hg19.chrom.sizes -b -300 | bedtools merge > $homer_dir/peaks.bed # 400 surround each peak then merge


#get 200000 random samples from hg19
#python /users/annashch/anna_utils/seq_utils/sample_genome_randomly.py --region_n 200000 --outf $homer_dir/candidate_negatives.bed --main_chroms_only
python $SCRIPT_PATH/sample_genome_randomly.py --region_n 500000 --outf $homer_dir/candidate_negatives.bed --main_chroms_only

#select candidate negatives so that they don't overlap the positives.
bedtools intersect -v -a $homer_dir/candidate_negatives.bed -b $orig_positives > $homer_dir/confirmed_negatives.bed

#assign positive and negative labels
python /users/annashch/anna_utils/seq_utils/add_fixed_labels_to_bed.py --bed $homer_dir/peaks.bed --labels 1 --outf $homer_dir/positives.bed
python /users/annashch/anna_utils/seq_utils/add_fixed_labels_to_bed.py --bed $homer_dir/confirmed_negatives.bed --labels 0 --outf $homer_dir/negatives.bed
#concatenative positives and negatives into a single file 

#cat $orig_positives | grep -v '\-1' > $positives
cat $positives $homer_dir/negatives.bed > $homer_dir/combined.bed

#run the dinucleotide matching code on "combined.bed" -- it has labels of "1" for positives and "0" for negatives.
#in the first pass, generate the matrix of dinucleotide frequencies. 

python $SCRIPT_PATH/gen_dinucleotide_freqs.py --bed_path $combined --ratio_neg_to_pos $ratio --outf $homer_dir/combined.freqs --ref_fasta /users/jocelins/hg19.fa

#now that the frequency matrix has been generated,  subselect negatives 
python $SCRIPT_PATH/gen_dinucleotide_freqs.py --bed_path $combined --ratio_neg_to_pos $ratio --outf $homer_dir/dinuc.bed --ref_fasta /users/jocelins/hg19.fa --dinuc_freqs $homer_dir/combined.freqs #--seed $(($i+42)) &

cat $dinuc0 | bedtools sort | uniq -u | shuf > $dinuc

cat $dinuc | grep '1$' | bedtools sort | bedtools merge > $homer_dir/pos.bed
cat $dinuc | grep '0$' | bedtools sort | bedtools merge > $homer_dir/neg.bed

bedtools getfasta -fi $TFNET_ROOT/genome/hg19.fa -bed $homer_dir/pos.bed -fo $homer_dir/pos.fa
bedtools getfasta -fi $TFNET_ROOT/genome/hg19.fa -bed $homer_dir/neg.bed -fo $homer_dir/neg.fa

findMotifs.pl $homer_dir/pos.fa human $homer_dir/homer_out -fasta $homer_dir/neg.fa -p 8 &> logs/homer.txt &

dinuc_shuf=$homer_dir/dinuc.bed.0.shuf
dinuc_tmp=$homer_dir/_tmp_shuf.bed

python /users/annashch/anna_utils/seq_utils/add_fixed_labels_to_bed.py --bed interpret.tsv --labels 1 --outf $homer_dir/_tmp_pos.bed
cat $dinuc_shuf | grep '0$'  > $homer_dir/_tmp_neg.bed
cat $homer_dir/_tmp_pos.bed $homer_dir/_tmp_neg.bed | bedtools sort | uniq -u | shuf > $dinuc_tmp

bedtools getfasta -fi $TFNET_ROOT/genome/hg19.fa -bed $dinuc_tmp -fo inputs.fa
cat headers.txt $dinuc_tmp > labels.tsv
cat labels.tsv | sed -e 's/\t/:/; s/\t/-/' > labels.txt
cat labels.txt | pigz -c > labels.txt.gz

valid_chrom="chr2:"
test_chrom="chr1:"
cat labels.txt | grep "chr2:" | pigz -c > splits/valid.txt.gz
cat labels.txt | grep "chr1:"  | pigz -c > splits/test.txt.gz
cat labels.txt | grep -v "chr1:\|chr2:\|^id" | pigz -c > splits/train.txt.gz

ln -s $TFNET_ROOT/results/templates/make_hdf5_yaml .
make_hdf5 --yaml_configs make_hdf5_yaml/* --output_dir .

rm $homer_dir/_tmp_pos.bed $homer_dir/_tmp_neg.bed $homer_dir/_tmp_shuf.bed

rm $homer_dir/peaks.bed

