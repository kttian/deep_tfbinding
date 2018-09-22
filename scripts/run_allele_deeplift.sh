#cat $TFNET_ROOT/results/nandi/bQTL/analysis/bQTL_all_SNPs/SPI1_10k.txt | python $TFNET_ROOT/scripts/snp_to_bed.py 1000 > SPI1_10k.tsv
#ln -s SPI1_10k.tsv interpret.tsv
#bedtools getfasta -fi $TFNET_ROOT/genome/hg19.fa -bed interpret.tsv -fo interpret.fa
#ln -s /home/ktian/kundajelab/tfnet/results/nandi/SPI1/SPI1_GM12878_refine_18_09_04/model_files .
python $TFNET_ROOT/scripts/run_deeplift.py model_files/record_1_ scores/interpret.fa0 1 > logs/deeplift0.txt 2>&1
mv scores/hyp_scores_task_0.npy scores/hyp_scores_task_0.npy0
python $TFNET_ROOT/scripts/run_deeplift.py model_files/record_1_ scores/interpret.fa1 1 > logs/deeplift1.txt 2>&1
mv scores/hyp_scores_task_0.npy scores/hyp_scores_task_0.npy1
python $TFNET_ROOT/scripts/run_deeplift.py model_files/record_1_ scores/interpret.fa2 1 > logs/deeplift2.txt 2>&1
mv scores/hyp_scores_task_0.npy scores/hyp_scores_task_0.npy2
python $TFNET_ROOT/scripts/run_deeplift.py model_files/record_1_ scores/interpret.fa3 1 > logs/deeplift3.txt 2>&1
mv scores/hyp_scores_task_0.npy scores/hyp_scores_task_0.npy3


#python get_scores.py 
