#!/usr/bin/env python3

"""
This script takes a csv with LLRs across 0.1 VARITY score intervals and the ouptut from 
the DAVID Functional Annotation Tool and outputs the same csv with an additional column 
with 'functional cluster annotation number'

See "LLR_intervals.csv" for sample input csv
See "david_functional_annotation_clustering_GO_MF_all_terms_low_stringency.txt" for sample DAVID output

See "LLR_intervals_w_david_functional_annotation.csv" for sample output
"""

import argparse

def make_cluster_dict(annotation_file: str) -> dict:

	cluster_dict = {}
	cluster_count = 1

	with open(annotation_file, "r") as clusters:	
		for cluster in (clusters.read().split("\n\t")):
			if cluster != "":
				if cluster_count not in cluster_dict:
					cluster_dict[cluster_count] = []
				for category in cluster.split("\n"):
					if "Annotation Cluster" not in category and ("PValue" not in category):
						genes_in_category = category.split("\t")[5]
						for gene in genes_in_category.strip("\"").split(","):
							if gene.strip() != "" and gene.strip() not in cluster_dict[cluster_count]:
								cluster_dict[cluster_count].append(gene.strip())
				cluster_dict[cluster_count].sort()
				cluster_count +=1

	return cluster_dict

def write_output(input_csv: str, cluster_dict: dict, output_name: str) -> None:
	with open(output_name + ".csv", "w") as output:
		with open(input_csv, "r") as llr_list:
			output.write(llr_list.readline().strip() + ",functional_cluster\n")
			for line in llr_list:
				gene = line.split(",")[0].strip()
				for cluster in cluster_dict:
					if gene in cluster_dict[cluster]:
						output.write(line.strip() + f",{cluster}\n")


if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Add DAVID annotation cluster column to LLR csv')
	parser.add_argument("-a", "--annotation_file", help="Path to DAVID annotation cluster file", type=str)
	parser.add_argument("-i", "--input_csv", help="Path to LLR file", type=str)
	parser.add_argument("-o", "--output", help="Output filename", type=str)
	args = parser.parse_args()

	cluster_dict = make_cluster_dict(args.annotation_file)

	write_output(args.input_csv, cluster_dict, args.output)
