DataDir = "./"
header_file = DataDir + "headers.txt"
valid_file  = DataDir + "splits/valid.tsv.gz"
import gzip

with open(header_file, 'r') as f:
    head = f.readline()
    head = head.split()[1:]
#print(len(head))

import pandas as pd

#df = pd.read_csv(valid_file, sep='\t', index_col=0, header=None, compression='gzip')
df0 = pd.read_csv(valid_file, sep='\t', index_col=None, header=None, compression='gzip')
df = df0.iloc[:,3:]
df.columns = head

#print(df.head(2))


counts_list = []
for c in range(df.shape[1]):
    counts = df.iloc[:, c].value_counts().sort_index(ascending=False)
    counts_list.append(counts)

counts_df = pd.concat(counts_list, axis=1).sort_index(ascending=False)
#print(counts_df)
sums = counts_df.sum()
#print (sums)

#for r in range(counts_df.shape[0]):
#        percent = counts_df.iloc[r] / sums
# In[ ]:

percent_df = counts_df / sums

new_index = [ int(1), int(0), int(-1)]
counts_df  = counts_df.reindex(new_index)
percent_df = percent_df.reindex(new_index)

#print (percent_df)

all=pd.concat([counts_df, percent_df])


all_t = all.transpose()
all_t.columns=['pos', 'neg', 'amb', 'pos%', 'neg%', 'amb%']
#print(all_t)
all_t.to_csv("counts.tsv", sep='\t')

