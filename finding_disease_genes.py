import pandas as pd 
import re

df = pd.read_csv('variant_summary.txt', sep='\t+')
df = df[['#AlleleID', 'Type', 'Name', 'GeneID', 'GeneSymbol', 
         'ClinicalSignificance', 'RS# (dbSNP)', 'Assembly', 'ReviewStatus',]]

df = df[df['Type'] == 'single nucleotide variant']

df = df[df['Assembly'] == 'GRCh38']

df = df[df.Name.str.contains("\(p\.[a-zA-Z]{3}\d+[a-zA-Z]{3}\)$")]

df = df[(df['ReviewStatus'] == 'criteria provided, single submitter') |
       (df['ReviewStatus'] == 'criteria provided, multiple submitters, no conflicts') |
       (df['ReviewStatus'] == 'criteria provided, conflicting interpretations') |
       (df['ReviewStatus'] == 'reviewed by expert panel') |
       (df['ReviewStatus'] == 'practice guideline')]

df = df[(df['ClinicalSignificance'] == 'Pathogenic') |
       (df['ClinicalSignificance'] == 'Benign') |
       (df['ClinicalSignificance'] == 'Likely pathogenic') |
       (df['ClinicalSignificance'] == 'Likely benign') |
       (df['ClinicalSignificance'] == 'Benign/Likely benign') |
       (df['ClinicalSignificance'] == 'Pathogenic/Likely pathogenic')]

df_pathogenic = df[df['ClinicalSignificance'] == 'Pathogenic']

disease_genes_set = df_pathogenic.GeneSymbol.unique()

disease_genes = df[df['GeneSymbol'].isin(disease_genes_set)]

hgvsp = []
for name in disease_genes.Name:
    x = ''.join(re.findall("p\.[a-zA-Z]{3}\d+[a-zA-Z]{3}", name))
    hgvsp.append(x)

disease_genes['hgvsp'] = hgvsp

disease_genes.sort_values(by=['GeneSymbol'], inplace=True)


disease_genes = disease_genes[['GeneSymbol', 'GeneID', 'hgvsp', 'ClinicalSignificance']]

disease_genes = disease_genes[~disease_genes.GeneSymbol.str.contains(';')]
disease_genes = disease_genes.reset_index(drop=True)

disease_genes.to_csv('disease_genes.csv')

