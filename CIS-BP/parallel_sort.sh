#parallel -k --pipepart -a FIMO-CIS-BP-Homo_sapiens.hg19.fimo.txt.annotated --block 1g cut -f 2-4 > out
#parallel -k --pipepart -a cut_out_n --block 1g python calc_loc.py > out

for i in `seq 1 22` X Y
do
    chrom=chr$i
    #echo $chrom
    sort -n -k 2 data/loc_$chrom.txt > sort_$chrom.txt&
done

