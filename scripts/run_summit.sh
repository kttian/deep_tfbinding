#python /home/ktian/kundajelab/tfnet/scripts/pick_summit.py NANOG H1-hESC | bedtools sort > NANOG_summit.tsv
#grep -P 'chr1\t' NANOG_summit.tsv > NANOG_summit_test.tsv 
#grep -v -P 'chr1\t' NANOG_summit.tsv > NANOG_summit_train_valid.tsv 
#bedtools getfasta -fi /home/ktian/kundajelab/tfnet/genome/hg19.fa -bed NANOG_summit.tsv -fo NANOG_summit.fa
bedtools getfasta -fi /home/ktian/kundajelab/tfnet/genome/hg19.fa -bed NANOG_summit_train_valid.tsv -fo NANOG_summit_train_valid.fa
python /home/ktian/kundajelab/tfnet/scripts/run_deeplift.py model_files/record_1_ NANOG_summit_train_valid.fa 1 > logs/deeplift.log 2>&1
python run_tfmodisco.py scores/hyp_scores_task_ NANOG_summit_train_valid.fa NANOG_summit_train_valid.tsv 1234 1 > logs/modisco.log 2>&1
