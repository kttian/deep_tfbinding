#grep 'JUND$' FIMO-CIS-BP-Homo_sapiens.hg19.fimo.txt.annotated_n > JUND.txt &
#grep 'SPI1$' FIMO-CIS-BP-Homo_sapiens.hg19.fimo.txt.annotated_n > SPI1.txt &
#grep 'RELA$' FIMO-CIS-BP-Homo_sapiens.hg19.fimo.txt.annotated_n > RELA.txt &
#grep 'STAT1$' FIMO-CIS-BP-Homo_sapiens.hg19.fimo.txt.annotated_n > STAT1.txt &
#grep 'CTCF$' FIMO-CIS-BP-Homo_sapiens.hg19.fimo.txt.annotated_n > CTCF.txt &

for i in JUND SPI1 RELA STAT1 CTCF
do
#cat $i.txt | python ../calc_loc.py > ${i}_loc.txt&
cat ${i}_loc.txt | sort -n -k2 > ${i}_sort.txt&
done

