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

cat $dinuc | grep -P 'chr1\t' | pigz -c > splits/test.tsv.gz
cat $dinuc | grep -P 'chr2\t' | pigz -c > splits/valid.tsv.gz
cat $dinuc | grep -v -P 'chr1\t' | grep -v -P 'chr2\t' | pigz -c > splits/train.tsv.gz


bedtools getfasta -fi $TFNET_ROOT/genome/hg19.fa -bed $dinuc -fo inputs.fa
cat headers.txt $dinuc > labels.txt
cat labels.txt | sed -e 's/\t/:/; s/\t/-/' | pigz -c > labels.txt.gz
pigz -c -d splits/test.tsv.gz | sed -e 's/\t/:/; s/\t/-/' | pigz -c > splits/test.txt.gz
pigz -c -d splits/valid.tsv.gz | sed -e 's/\t/:/; s/\t/-/' | pigz -c > splits/valid.txt.gz
pigz -c -d splits/train.tsv.gz | sed -e 's/\t/:/; s/\t/-/' | pigz -c > splits/train.txt.gz

make_hdf5 --yaml_configs make_hdf5_yaml/* --output_dir .


rm $homer_dir/peaks.bed

