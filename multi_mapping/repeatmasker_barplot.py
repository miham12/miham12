import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#unique part 
unique_dna=pd.read_csv('SRR3633289_unique_repeats.bed',sep='\t',header=None,usecols=[3,5],names=['read_id','repeat_id'])
unique_dna[['repeat_id', 'file']] = unique_dna['repeat_id'].str.split('\|file:', expand=True)

unique_repeatmasker=unique_dna[unique_dna['file']=='repeatmasker']
del unique_dna

unique_repeatmasker[['repeat_id', 'class','family']]=unique_repeatmasker['repeat_id'].str.split(';_;', expand=True)
del unique_repeatmasker['repeat_id']
del unique_repeatmasker['family']
unique_repeatmasker['times_mapped']='unique'


#multi part 
multi_dna=pd.read_csv('SRR3633289_multi_repeats.bed',sep='\t',header=None,usecols=[3,5],names=['read_id','repeat_id'])

#all_ids_unique=list(unique_dna['read_id'])
#all_ids_multi=list(multi_dna['read_id'])
#all_ids=list(set(all_ids_unique+all_ids_multi))

multi_dna[['repeat_id', 'file']] = multi_dna['repeat_id'].str.split('\|file:', expand=True)

multi_repeatmasker=multi_dna[multi_dna['file']=='repeatmasker']
del multi_dna

multi_repeatmasker[['repeat_id', 'class','family']]=multi_repeatmasker['repeat_id'].str.split(';_;', expand=True)
del multi_repeatmasker['repeat_id']
del multi_repeatmasker['family']

read_id_counts = multi_repeatmasker['read_id'].value_counts().reset_index()
read_id_counts.columns = ['read_id', 'count']
# Добавим новый столбец times_mapped
multi_repeatmasker = multi_repeatmasker.merge(read_id_counts, on='read_id', how='left')
multi_repeatmasker['times_mapped'] = pd.cut(multi_repeatmasker['count'], bins=[2, 9, float('inf')], labels=['2_9', '10_more'], right=False)
del multi_repeatmasker['count']

repeatmasker = pd.concat([unique_repeatmasker, multi_repeatmasker])
del unique_repeatmasker 
del multi_repeatmasker
plt.figure(dpi=1200)
plt.figure(figsize=(15, 9))
ax=sns.countplot(data=repeatmasker, x='class', hue='times_mapped',palette='Set2')

ax.axes.set_title("DNA reads annotation via repeatmasker",fontsize=20)
ax.set_xlabel("Repeatmasker class",fontsize=16)
ax.set_ylabel("Counts",fontsize=16)
ax.tick_params(labelsize=13,labelrotation=90)
plt.legend(fontsize = "18", loc ="upper right")
plt.yscale('log')

plt.tight_layout() 
plt.savefig('SRR3633289_intersected_repets.jpg', dpi=500)