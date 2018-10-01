
homer_dir="homer/"
SCRIPT_PATH="$TFNET_ROOT/scripts/"

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


