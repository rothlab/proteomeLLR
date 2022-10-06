import pickle

fh = open('disease_genes_dictionary.pkl', 'rb')
disease_genes_dict = pickle.load(fh)

print(disease_genes_dict)