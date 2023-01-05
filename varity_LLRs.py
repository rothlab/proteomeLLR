#!/usr/bin/env python3

import re
from xopen import xopen
import os
import argparse
import logging
from TileSeqMut import help_functions


# def uniprot_to_geneid(id_mapping: str):
#     """
#     Create a dictionary from the id_mapping file 'id_mapping' where the keys are UniProt IDs and 
#     the values are corresponding GeneIDs.
#     """

#     mapping_dict = {}
#     with xopen(id_mapping, "rt") as mapping:
#         for line in mapping:
#             mapping_dict[line.strip().split("\t")[-1]] = line.strip().split("\t")[0]

#     return mapping_dict


# def get_varity_supported_ids(varity_supported_ids: str):

#     supported_ids = {}

#     with xopen(varity_supported_ids, "rt") as ids:
#         ids.readline()
#         for line in ids:
#             supported_ids[line.strip().split("\t")[0]] = ""

#     return supported_ids


def uniprot_to_genesymbol(varity_supported_ids: str):
    """
    Return dictionary where the keys are UniProt IDs and the values are corresponding 
    gene symbols.
    """

    mapping_dict = {}
    with xopen(varity_supported_ids, "rt") as mapping:
        for line in mapping:
            mapping_dict[line.strip().split("\t")[0]] = line.strip().split("\t")[1]

    return mapping_dict


def get_clinvar_entries(variant_dict: dict, clinvar_db: str):
    """
    Modify dictionary 'variant_dict' to include variants from the clinvar database file 'clinvar_db'.

    Each variant is associated with a list - [Clinical Significance, Positive/Negative, Source, maf, num_hom]
    """

    missense_pattern = re.compile("^p\\.\\w{3}\\d+\\w{3}$")

    with open(clinvar_db) as clinvar:

        for line in clinvar:
            info = line.strip().split("\t")
            if (info[1] == "single nucleotide variant") and (missense_pattern.match(info[2].split(" ")[-1].strip("()"))) and ("Ter" not in info[2]) and ("no assertion" not in info[24]) and (info[6] not in ["Uncertain significance", "Conflicting interpretations of pathogenicity"]):
                gene_symbol = info[4].split(";")[0] #figure out to take first or last
                hgvs = info[2].split(" ")[-1].strip("()")
                if "pathogenic" in info[6].lower():
                    ref_set = "Positive"
                elif "benign" in info[6].lower():
                    ref_set = "Negative"
                
                if not gene_symbol in variant_dict:
                    variant_dict[gene_symbol] = {}
                    variant_dict[gene_symbol]["ref_set"] = []
                variant_dict[gene_symbol][hgvs] = [info[6], ref_set, "Clinvar", "NA", "NA"]
                variant_dict[gene_symbol]["ref_set"] =  variant_dict[gene_symbol]["ref_set"] + [ref_set]

    # Name (2), ClinicalSignificance (6), ClinSigSimple (7), ReviewStatus (24)



def get_frequencies(variant_dict: dict, args, pop_freq: str):
    """
    Update the variant dictionary "variant_dict" with the 'AF_popmax' allele frequency for  
    each variant from the cleaned gnomAD file 'pop_freq'
    """

    missense_pattern = re.compile("^p\\.\\w{3}\\d+\\w{3}")

    with xopen(pop_freq, "rt") as all_freqs:

        for freq in all_freqs: 
            num_hom = None
            maf = None
            gene_symbol = None

            freq_info = freq.strip().split(";")
            for info in freq_info:
                if info.startswith("controls_nhomalt"):
                    num_hom = info.split("=")[-1]
                elif info.startswith("AF_popmax"):
                    maf = info.split("=")[-1]
                elif "|" in info:
                    vep = info.split("|")
                    if missense_pattern.match(vep[1]):
                        gene_symbol = vep[0]
                        hgvs = missense_pattern.match(vep[1]).group(0)

            #check hgvs formatting in gnomad file

            # take max of MAF and num_hom from variants that cause the same protein-level change
            if gene_symbol in variant_dict and hgvs in variant_dict[gene_symbol]:

                if maf is not None and variant_dict[gene_symbol][hgvs][3] is not None:
                    variant_dict[gene_symbol][hgvs][3] = str(max(float(maf), float(variant_dict[gene_symbol][hgvs][3])))
                elif variant_dict[gene_symbol][hgvs][3] is None:
                    variant_dict[gene_symbol][hgvs][3] = maf
                if num_hom is not None and variant_dict[gene_symbol][hgvs][4] is not None:
                    variant_dict[gene_symbol][hgvs][4] = str(max(int(num_hom), int(variant_dict[gene_symbol][hgvs][4])))
                elif variant_dict[gene_symbol][hgvs][4] is None:
                    variant_dict[gene_symbol][hgvs][4] = num_hom
                

            elif gene_symbol in variant_dict and ((maf is not None and float(maf) > args.maf_cutoff) or (num_hom is not None and int(num_hom) > 0)) and (hgvs != ""):
                # if gene_symbol not in variant_dict:
                #     variant_dict[gene_symbol] = {}
                #     variant_dict[gene_symbol]["ref_set"] = []
                variant_dict[gene_symbol][hgvs] = ["Predicted benign", "Negative", "GnomAD", maf, num_hom]
                variant_dict[gene_symbol]["ref_set"] =  variant_dict[gene_symbol]["ref_set"] + ["Negative"]


def get_varity_predictions(varity_predictions: str, prediction_dict: dict):
    """
    Update dictionary 'prediction_dict' with all varity predictions from file 
    'varity_predictions', where the keys are uniprot gene identifiers and the 
    values are sub-dictionaries, where the keys are variants within the gene  
    and the values are their corresponding VARITY score.
    """

    code_conversion = {'C': 'Cys', 'D': 'Asp', 'S': 'Ser', 'Q': 'Gln', 'K': 'Lys',
     'I': 'Ile', 'P': 'Pro', 'T': 'Thr', 'F': 'Phe', 'N': 'Asn', 
     'G': 'Gly', 'H': 'His', 'L': 'Leu', 'R': 'Arg', 'W': 'Trp', 
     'A': 'Ala', 'V':'Val', 'E': 'Glu', 'Y': 'Tyr', 'M': 'Met'}

    with open(varity_predictions, "r") as predictions: 
        predictions.readline()
        for line in predictions:
            prediction = line.strip().split("\t")
            uniprot_id = prediction[4]
            hgvs = f"p.{code_conversion[prediction[6]]}{prediction[5]}{code_conversion[prediction[7]]}"
            score = prediction[8]

            if uniprot_id not in prediction_dict:
                prediction_dict[uniprot_id] = {}
            if hgvs not in prediction_dict[uniprot_id]:
                prediction_dict[uniprot_id][hgvs] = []
            prediction_dict[uniprot_id][hgvs].append(float(score))


def filter_variants(variant_dict, id_mapping_dict, prediction_dict, min_ref, main_log):
    """
    Update dictionary 'variant_dict' to remove genes with insufficient reference variants
    and genes that do not have VARITY predictions
    """

    to_remove = []

    for gene in variant_dict:
        #check if there sufficient positive and negative reference variants
        if (variant_dict[gene]["ref_set"].count("Positive") < min_ref) or (variant_dict[gene]["ref_set"].count("Negative") < min_ref): 
            main_log.info(f"Insufficient reference variants for {gene}")
            to_remove.append(gene)
        #check if the gene is in the id_mapping file and if the id_mapping of the gene id is in the prediction files
        elif (gene not in id_mapping_dict) or (id_mapping_dict[gene] not in prediction_dict): # or (gene not in varity_ids):
            main_log.info(f"Varity predictions unavailable for {gene}")
            to_remove.append(gene)

    for gene in to_remove:
        variant_dict.pop(gene) 


def make_refsets(variant_dict, id_mapping_dict, prediction_dict, main_log):
    """
    Output a reference file for each gene in dictionary 'variant_dict'
    """

    for gene in variant_dict:
        with open(f"{gene}_reference.csv", "w") as output:
            main_log.info(f"Writing to {gene}_reference.csv")
            output.write("hgvsp,Clinical Significance,referenceSet,source,maf,hom\n")
            for variant in variant_dict[gene]:
                if variant != "ref_set" and variant != "GeneID": ######## FIX THIS ###########
                    output.write(",".join([variant] + variant_dict[gene][variant]) + "\n")


def make_maps(variant_dict, prediction_dict, id_mapping_dict, main_log):
    """
    Output a map file for each gene in dictionary 'variant_dict'
    """
    for gene in variant_dict:
        uniprot_id = id_mapping_dict[gene]

        with open(f"{gene}_map.csv", "w") as output:
            main_log.info(f"Writing to {gene}_map.csv")
            output.write("hgvs_pro,score\n")
            for hgvs in prediction_dict[uniprot_id]: 
                prediction_list = prediction_dict[uniprot_id][hgvs]
                output.write(hgvs + "," + str(sum(prediction_list)/len(prediction_list)) + "\n")


def main(args):
    if not args.output:
        if not os.path.isdir(os.path.join(os.getcwd(), "refsets")):
            output_directory = os.mkdir(os.path.join(os.getcwd(), "refsets"))
        else: 
            output_directory = os.path.join(os.getcwd(), "refsets")
    elif args.output and not os.path.isdir(args.output):
        raise FileNotFoundError(f"Output directory not found: {args.output}")
    else:
        output_directory = args.output

    clinvar_db = args.clinvar
    pop_freq = args.gnomad
    varity_supported_ids = args.varity_ids
    varity_predictions = args.varity_predictions

    main_log_name = os.path.join(output_directory, "main.log")
    main_log = help_functions.logginginit("info", main_log_name)

    variant_dict = {}
    prediction_dict = {}

    #main_log.info(f"Mapping Uniprot IDs to GeneIDs")
    #id_mapping = "/home/rothlab/skathirg/varity/idmapping.txt"
    #id_mapping_dict = uniprot_to_geneid(id_mapping)

    main_log.info(f"Mapping varity supported UniProt IDs to gene symbols")
    id_mapping_dict = uniprot_to_genesymbol(varity_supported_ids)


    #main_log.info(f"Getting varity supported ids")
    #varity_ids = get_varity_supported_ids(varity_supported_ids)
    #main_log.info(f"{len(varity_ids)} varity supported ids")

    main_log.info(f"Getting varity predictions")
    get_varity_predictions(varity_predictions, prediction_dict)

    main_log.info(f"Getting ClinVar variants from {clinvar_db}")
    get_clinvar_entries(variant_dict, clinvar_db)

    main_log.info(f"Getting variant frequencies from {pop_freq}")
    get_frequencies(variant_dict, args, pop_freq)

    filter_variants(variant_dict, id_mapping_dict, prediction_dict, args.min_ref, main_log)

    os.chdir(output_directory)
    
    make_refsets(variant_dict, id_mapping_dict, prediction_dict, main_log)

    make_maps(variant_dict, prediction_dict, id_mapping_dict, main_log)

    main_log.shutdown()


if __name__ == "__main__":

    uniprot_geneid_mapping = "/home/rothlab/skathirg/varity/id_mapping.tab.gz"

    parser = argparse.ArgumentParser(description='Generate reference sets for calcLLR function')
    parser.add_argument("--clinvar", help="Path to ClinVar variant summary file", type=str, default="variant_summary.txt")
    parser.add_argument("--gnomad", help="Path to GnomAD variant summary file", type=str, default="gnomad_frequencies.vcf.gz")
    parser.add_argument("--varity_ids", help="List of VARITY supported ids", type=str, default="varity_supported_ids.txt")
    parser.add_argument("--varity_predictions", help="VARITY predictions", type=str, default="varity_all_predictions.txt")
    parser.add_argument("--num_hom", help="GnomAD homozygote threshold", type=int, default=1)
    parser.add_argument("--min_ref", help="Minimum number of reference pathogenic and benign variants ", type=int, default=11)
    parser.add_argument("--maf_cutoff", help="MAF threshold to be considered benign", type=float, default=0.001)
    parser.add_argument("-o", "--output", help="Output directory where reference files should be written", type=str)
    args = parser.parse_args()

    main(args)

"""
,#AlleleID (0),Type (1),Name (2),GeneID (3),GeneSymbol (4),ClinicalSignificance (6),ClinSigSimple,Assembly (16),Chromosome (18),Start (19),Stop 
(20),ReferenceAllele (21),AlternateAllele (22)

"""
