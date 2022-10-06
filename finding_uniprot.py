import gzip
import pandas as pd
import csv

refset = pd.read_csv('ReferenceSet.csv')

gene_ids = refset['GeneID'].unique()
gene_ids = [str(gene) for gene in gene_ids]


fh = open('uniprot_ids.csv', 'w', newline='')
writer = csv.writer(fh)
writer.writerow(['GeneID', 'Uniprot Number'])

with gzip.open('HUMAN_9606_idmapping_selected.tab.gz', 'rb') as f:
    for line in f:
        line = line.decode('utf-8')
        splitted_line = line.split('\t')
        gene_id = splitted_line[2]
        if gene_id in gene_ids:
            uniprot_number = splitted_line[0]
            writer.writerow([gene_id, uniprot_number])


