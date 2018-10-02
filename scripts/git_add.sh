modisco_run="modisco.run1"
mkdir $modisco_run
mkdir $modisco_run/logs
mv logs/modisco.txt $modisco_run/logs
mv results.hdf5 $modisco_run
mv figures $modisco_run

git add $modisco_run/logs $modisco_run/figures tfmodisco* SNP*

git add config
git add run_all.sh prepare_homer.sh
git add *counts.tsv
git add headers.txt

git add finetune/model_files/record_1_model*
git add finetune/model_files/*.tsv
git add finetune/model_files/*.png
git add finetune/runs*

git add logs/pipe.txt
git add logs/pre_prepare.txt
git add logs/prepare.txt
git add logs/analyze.txt
git add logs/test.txt
git add logs/homer.txt
git add logs/prep_homer.txt
