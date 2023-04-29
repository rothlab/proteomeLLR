#!/usr/bin/env python3

"""
This script takes a csv with LLRs across 0.1 VARITY score intervals and a GenCC submissions file 
(from https://search.thegencc.org/download) and outputs the same csv with an additional column 
with inheritance data

NOTE: all columns but 'gene_symbol', 'disease_title', 'classification_title' and 'moi_title'
were removed from the GenCC file beforehand

See "LLR_intervals.csv" for sample input csv
See "gencc-submissions" for sample GenCC file

See "LLR_intervals_w_inheritance.csv" for sample output
"""

import argparse

def make_inheritance_dict(input_csv: str, gencc_file: str) -> dict:

	llr_genes = {}

	with open(input_csv, "r") as llr_intervals:
		llr_intervals.readline()
		for line in llr_intervals:
			gene = line.split(",")[0].strip()
			llr_genes[gene] = {}

	with open(gencc_file, "r") as gencc:

		gencc.readline()

		for line in gencc:
			gene_info = line.split("\t")
			if (gene_info[0].strip() in llr_genes) and (gene_info[3].strip() not in llr_genes[gene_info[0].strip()]): 
				llr_genes[gene_info[0].strip()][gene_info[3].strip()] = []
				if (gene_info[3].strip() in llr_genes[gene_info[0].strip()]) and (gene_info[2].strip() not in llr_genes[gene_info[0].strip()][gene_info[3].strip()]):
					if gene_info[2].strip() in ["Definitive", "Strong"]:
						llr_genes[gene_info[0].strip()][gene_info[3].strip()].append(gene_info[2].strip())	

	return llr_genes


def write_output(input_csv: str, inheritance_dict: dict, output_name: str) -> None:

	with open(output_name + ".csv", "w") as output:
		with open(input_csv, "r") as llr_intervals:

			output.write(llr_intervals.readline().strip() + ",Inheritance\n")

			for line in llr_intervals:
				gene = line.split(",")[0].strip()

				if list(inheritance_dict[gene].keys()) == ['Autosomal recessive']:
					inheritance = "AR"
				elif list(inheritance_dict[gene].keys()) == ['Autosomal dominant']:
					inheritance = "AD"
				elif list(inheritance_dict[gene].keys()).sort() == ['Autosomal dominant', 'Autosomal recessive']:
					if "Definitive" in inheritance_dict[gene]['Autosomal dominant'] and "Definitive" not in inheritance_dict[gene]['Autosomal recessive']:
						inheritance = "AD"
					elif "Definitive" in inheritance_dict[gene]['Autosomal recessive'] and "Definitive" not in inheritance_dict[gene]['Autosomal dominant']:
						inhertiance = "AR"
				else:
					inheritance = "other"

				output.write(line.strip() + f",{inheritance}\n")


if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Add inheritance column to LLR csv from GENCC data')
	parser.add_argument("-g", "--gencc_file", help="Path to GENCC  file", type=str)
	parser.add_argument("-i", "--input_csv", help="Path to LLR file", type=str)
	parser.add_argument("-o", "--output", help="Output filename", type=str)
	args = parser.parse_args()

	
	inheritance_dict = make_inheritance_dict(args.input_csv, args.gencc_file)

	write_output(args.input_csv, inheritance_dict, args.output)
			


"""
Types of inheritance: ['Autosomal recessive', 'Autosomal dominant', 'Semidominant', 'Unknown', 
'X-linked', 'X-linked recessive', 'X-linked dominant', 'Somatic mosaicism', 'Autosomal dominant inheritance with paternal imprinting', 'Mitochondrial']
Types of classification: Definitive, Disputed Evidence, Limited, Moderate, No Known Disease Relationship, Refuted Evidence, Strong, Supportive
"""
