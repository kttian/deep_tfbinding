#This is an interpretation run with the CTCF model that Georgi trained. 
#The model is at ~/oak/marinovg/various/2018-04-20-K562-TFs/tfdragonn-SequenceClassifier-K562-CTCF-human-ChIP-seq-HA

ln -s ~/oak/marinovg/various/2018-04-20-K562-TFs/tfdragonn-SequenceClassifier-K562-CTCF-human-ChIP-seq-HA/model.arch.yaml .
ln -s ~/oak/marinovg/various/2018-04-20-K562-TFs/tfdragonn-SequenceClassifier-K562-CTCF-human-ChIP-seq-HA/model.weights.h5 .
ln -s ~/oak/marinovg/various/2018-04-20-K562-TFs/labelregions-K562-CTCF-ChIP-seq-HA.intervals_file.tsv.gz original_labels.tsv.gz
