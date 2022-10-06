import pandas as pd
import pickle

disease_genes_df = pd.read_csv('disease_genes.csv', index_col=0, low_memory=False)
disease_genes_dict = dict.fromkeys(tuple(zip(list(disease_genes_df['GeneSymbol']), list(disease_genes_df['GeneID']))))

for gene in disease_genes_dict:
    df = disease_genes_df[disease_genes_df['GeneSymbol'] == gene[0]]
    disease_genes_dict[gene] = tuple(zip(list(df['hgvsp']), list(df['ClinicalSignificance'])))

disease_genes_file = open('disease_genes_dictionary.pkl', 'wb')
pickle.dump(disease_genes_dict, disease_genes_file)
disease_genes_file.close()