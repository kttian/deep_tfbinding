
import gzip

#with open('labels.txt', 'r') as f:
with gzip.open('labels.txt.gz', 'r') as f:
    head = f.readline()
    head = head.split()[1:]
#print(len(head))

import pandas as pd

df = pd.read_csv("splits/valid.txt.gz", sep='\t', index_col=0, header=None, 
                 compression='gzip')
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

#print (percent_df)

all=pd.concat([counts_df, percent_df])


all_t = all.transpose()
all_t.columns=['pos', 'neg', 'amb', 'pos%', 'neg%', 'amb%']
#print(all_t)
all_t.to_csv("counts.tsv", sep='\t')

