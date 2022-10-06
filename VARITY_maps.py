import pandas as pd
import gzip
import Bio.SeqUtils
import csv

amino_acids_dict = Bio.SeqUtils.IUPACData.protein_letters_1to3
df = pd.read_csv('uniprot_ids.csv')
uniprot_numbers = df['Uniprot Number']
uniprot_numbers = list(uniprot_numbers)

with gzip.open('varity_all_predictions.tar.gz', 'rb') as f:
    for line in f:
        line = line.decode('utf-8')
        splitted_line = line.split('\t')
        uniprot_number = splitted_line[4]
        if uniprot_number in uniprot_numbers:
            
            fh = open(f"varity_maps/{uniprot_number}_map.csv", 'a', newline='')
            writer = csv.writer(fh)

            hgvsp = f"p.{amino_acids_dict[splitted_line[6]]}{splitted_line[5]}{amino_acids_dict[splitted_line[7]]}"
            score = splitted_line[8]

            writer.writerow([hgvsp, score])
