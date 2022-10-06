import pandas as pd
import tarfile

genes_tar_file = 'genes_referenceset.tar'
file = tarfile.open(genes_tar_file, 'w')

df = pd.read_csv('ReferenceSet.csv')

for (gene), group in df.groupby(['GeneID']):
    group_values = list(group['referenceSet'].value_counts())
    if len(group_values) > 1:
        if all(value >=6 for value in group_values):
            group.to_csv('referencesets/'f'{gene}.csv', index=False)
            file.add('referencesets/'f'{gene}.csv')
file.close()
