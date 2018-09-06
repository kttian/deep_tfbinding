#parallel -k --pipepart -a FIMO-CIS-BP-Homo_sapiens.hg19.fimo.txt.annotated --block 1g cut -f 2-4 > out
#parallel -k --pipepart -a cut_out_n --block 1g python calc_loc.py > out

for i in `seq 3 22` X Y x y
do
    chromstr=chr$i
    chrom=chr$i
    #echo $chrom
    grep -P $chromstr cut_out_n > out_$chrom.txt &
done

grep -P "chr1\t" cut_out_n >out_chr1.txt&
grep -P "chr2\t" cut_out_n >out_chr2.txt&
