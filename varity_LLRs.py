#!/usr/bin/env python3

import re
from xopen import xopen
import os
import argparse
import logging
from TileSeqMut import help_functions # would have to implement this on my own to remove this dependency


def uniprot_to_geneid(id_mapping: str):

    mapping_dict = {}
    with xopen(id_mapping, "rt") as mapping:
        for line in mapping:
            mapping_dict[line.strip().split("\t")[-1]] = line.strip().split("\t")[0]

    return mapping_dict


def get_clinvar_entries(variant_dict: dict, clinvar_db: str):

    missense_pattern = re.compile("^p\\.\\w{3}\\d+\\w{3}$")

    with open(clinvar_db) as clinvar:

        for line in clinvar:
            info = line.strip().split("\t")
            if (info[1] == "single nucleotide variant") and (missense_pattern.match(info[2].split(" ")[-1].strip("()"))) and ("Ter" not in info[2]) and ("no assertion" not in info[24]) and (info[6] not in ["Uncertain significance", "Conflicting interpretations of pathogenicity"]):
                gene_symbol = info[4].split(";")[-1]
                hgvsp = info[2].split(" ")[-1].strip("()")
                if "pathogenic" in info[6].lower():
                    ref_set = "Positive"
                elif "benign" in info[6].lower():
                    ref_set = "Negative"
                
                if not gene_symbol in variant_dict:
                    variant_dict[gene_symbol] = {}
                    variant_dict[gene_symbol]["ref_set"] = []
                variant_dict[gene_symbol][hgvsp] = [info[6], ref_set, "Clinvar", "NA", "NA"] # add info[24] for review status
                variant_dict[gene_symbol]["ref_set"] =  variant_dict[gene_symbol]["ref_set"] + [ref_set]
                variant_dict[gene_symbol]["GeneID"] = info[3]

    # Name (2), ClinicalSignificance (6), ClinSigSimple (7), ReviewStatus (24)



def get_frequencies(variant_dict: dict, args, pop_freq: str):
    """
    """

    with xopen(pop_freq, "rt") as all_freqs:
        all_freqs.readline()
        all_freqs.readline()
        all_freqs.readline()

        for freq in all_freqs: 
            num_hom = "NA"
            maf = "NA"

            freq_info = freq.strip().split(";")
            for info in freq_info:
                if info.startswith("nhomalt"):
                    num_hom = info.split("=")[-1]
                elif info.startswith("AF_popmax"):
                    maf = info.split("=")[-1]
                elif "|" in info:
                    vep = info.split("|")
                    gene_symbol = vep[0]
                    hgvsp = vep[1]


            #check hgvsp formatting in gnomad file

            #take max of MAF and num_hom from variants that cause the same protein-level change
            if gene_symbol in variant_dict and hgvsp in variant_dict[gene_symbol]:

                if maf != "NA" and variant_dict[gene_symbol][hgvsp][3] != "NA":
                    variant_dict[gene_symbol][hgvsp][3] = str(max(float(maf), float(variant_dict[gene_symbol][hgvsp][3])))
                elif variant_dict[gene_symbol][hgvsp][3] == "NA":
                    variant_dict[gene_symbol][hgvsp][3] = maf
                if num_hom != "NA" and variant_dict[gene_symbol][hgvsp][4] != "NA":
                    variant_dict[gene_symbol][hgvsp][4] = str(max(int(num_hom), int(variant_dict[gene_symbol][hgvsp][4])))
                elif variant_dict[gene_symbol][hgvsp][4] == "NA":
                    variant_dict[gene_symbol][hgvsp][4] = num_hom
                

            elif gene_symbol in variant_dict and ((maf != "NA" and float(maf) > args.maf_cutoff) or (num_hom != "NA" and int(args.num_hom) > 1)):
                # if gene_symbol not in variant_dict:
                #     variant_dict[gene_symbol] = {}
                #     variant_dict[gene_symbol]["ref_set"] = []
                variant_dict[gene_symbol][hgvsp] = ["Predicted benign", "Negative", "GnomAD", maf, num_hom]
                variant_dict[gene_symbol]["ref_set"] =  variant_dict[gene_symbol]["ref_set"] + ["Negative"]


def get_varity_supported_ids(varity_supported_ids: str):

    supported_ids = {}

    with xopen(varity_supported_ids, "rt") as ids:
        ids.readline()
        for line in ids:
            supported_ids[line.strip().split("\t")[0]] = ""

    return supported_ids


def get_varity_predictions(varity_predictions: str, prediction_dict: dict):

    code_conversion = {'C': 'Cys', 'D': 'Asp', 'S': 'Ser', 'Q': 'Gln', 'K': 'Lys',
     'I': 'Ile', 'P': 'Pro', 'T': 'Thr', 'F': 'Phe', 'N': 'Asn', 
     'G': 'Gly', 'H': 'His', 'L': 'Leu', 'R': 'Arg', 'W': 'Trp', 
     'A': 'Ala', 'V':'Val', 'E': 'Glu', 'Y': 'Tyr', 'M': 'Met'}

    with open(varity_predictions, "r") as predictions: 
        predictions.readline()
        for line in predictions:
            prediction = line.strip().split("\t")
            uniprot_id = prediction[4]
            hgvsp = f"p.{code_conversion[prediction[6]]}{prediction[5]}{code_conversion[prediction[7]]}"
            score = prediction[8]

            if uniprot_id not in prediction_dict:
                prediction_dict[uniprot_id] = []
            prediction_dict[uniprot_id] = prediction_dict[uniprot_id] + [f"{hgvsp}\t{score}"]


def filter_variants(variant_dict, id_mapping_dict, varity_ids, prediction_dict, min_ref, main_log):

    to_remove = []

    for gene in variant_dict:
        if (variant_dict[gene]["ref_set"].count("Positive") < min_ref) or (variant_dict[gene]["ref_set"].count("Negative") < min_ref):
            main_log.info(f"Insufficient reference variants for {gene}")
            to_remove.append(gene)
        elif (variant_dict[gene]["GeneID"] not in id_mapping_dict) or (gene not in varity_ids) or (id_mapping_dict[variant_dict[gene]["GeneID"]] not in prediction_dict):
            main_log.info(f"Varity predictions unavailable for {gene}")
            to_remove.append(gene)

    for gene in to_remove:
        variant_dict.pop(gene) 


def make_refsets(variant_dict, id_mapping_dict, varity_ids, prediction_dict, main_log):

    for gene in variant_dict:
        with open(f"{gene}_reference.csv", "w") as output:
            main_log.info(f"Writing to {gene}_reference.csv")
            output.write("hgvsp\tClinical Significance\treferenceSet\tsource\tmaf\thom\n")
            for variant in variant_dict[gene]:
                if variant != "ref_set" and variant != "GeneID":
                    output.write("\t".join([variant] + variant_dict[gene][variant]) + "\n")


def make_maps(variant_dict, prediction_dict, varity_ids, id_mapping_dict, main_log):
    for gene in variant_dict:
        uniprot_id = id_mapping_dict[variant_dict[gene]["GeneID"]]

        with open(f"{gene}_map.csv", "w") as output:
            main_log.info(f"Writing to {gene}_map.csv")
            output.write("hgvsp\tscore\n")
            for prediction in prediction_dict[uniprot_id]:
                output.write(prediction + "\n")


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

    main_log.info(f"Mapping Uniprot IDs to GeneIDs")
    id_mapping = "/home/rothlab/skathirg/varity/idmapping.txt"
    id_mapping_dict = uniprot_to_geneid(id_mapping)

    main_log.info(f"Getting varity supported ids")
    varity_ids = get_varity_supported_ids(varity_supported_ids)
    main_log.info(f"{len(varity_ids)} varity supported ids")

    main_log.info(f"Getting varity predictions")
    get_varity_predictions(varity_predictions, prediction_dict)

    main_log.info(f"Getting ClinVar variants from {clinvar_db}")
    get_clinvar_entries(variant_dict, clinvar_db)

    main_log.info(f"Getting variant frequencies from {pop_freq}")
    get_frequencies(variant_dict, args, pop_freq)

    filter_variants(variant_dict, id_mapping_dict, varity_ids, prediction_dict, args.min_ref, main_log)

    os.chdir(output_directory)
    
    make_refsets(variant_dict, id_mapping_dict, varity_ids, prediction_dict, main_log)

    make_maps(variant_dict, prediction_dict, varity_ids, id_mapping_dict, main_log)

    main_log.shutdown()


if __name__ == "__main__":

    uniprot_geneid_mapping = "/home/rothlab/skathirg/varity/id_mapping.tab.gz"

    parser = argparse.ArgumentParser(description='Generate reference sets for calcLLR function')
    parser.add_argument("--clinvar", help="Path to ClinVar variant summary file", type=str, default="variant_summary.txt")
    parser.add_argument("--gnomad", help="Path to GnomAD variant summary file", type=str, default="gnomad_frequencies.vcf.gz")
    parser.add_argument("--varity_ids", help="List of VARITY supported ids", type=str, default="varity_supported_ids.txt")
    parser.add_argument("--varity_predictions", help="VARITY predictions", type=str, default="varity_all_predictions.txt")
    parser.add_argument("--num_hom", help="GnomAD homozygote threshold", type=int, default=1)
    parser.add_argument("--min_ref", help="Minimum number of reference pathogenic and benign variants ", type=int, default=6)
    parser.add_argument("--maf_cutoff", help="MAF threshold to be considered benign", type=float, default=0.001)
    parser.add_argument("-o", "--output", help="Output directory where reference files should be written", type=str)
    args = parser.parse_args()

    main(args)


    #id_mapping_dict = uniprot_to_geneid(uniprot_geneid_mapping)

    #variant_dict = {}

    # with open(varity_supported_ids, "r") as varity_ids:
    #     varity_ids.readline()

    #     for line in varity_ids:
    #         if uniprot_id not in id_mapping_dict:
    #             gene_symbol = id_mapping_dict[uniprot_id]
                
    #         else:
    #             print(f"Uniprot ID {uniprot_id} not in mapping file .... skipping.")

    # get_clinvar_entries(variant_dict, clinvar_db)

    # get_frequencies(variant_dict, args, pop_freq)

    # os.chdir(output_directory)

    # for gene in variant_dict:
    #     if variant_dict[gene]["ref_set"].count("Positive") >= 6 and variant_dict[gene]["ref_set"].count("Negative") >= 6:
    #         with open(gene + "_reference.csv", "w") as output:
    #             print(f"Writing to {gene}_reference.csv")
    #             output.write("hgvsp\tClinical Significance\treferenceSet\tsource\tmaf\thom\n")
    #             for variant in variant_dict[gene]:
    #                 if variant != "ref_set":
    #                     output.write("\t".join([variant] + variant_dict[gene][variant]) + "\n")
    #     else:
    #         print(f"Insufficient reference variants for {gene}")




"""

[row], 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 
,#AlleleID (0),Type (1),Name (2),GeneID (3),GeneSymbol (4),HGNC_ID (5),ClinicalSignificance (6),ClinSigSimple (7),LastEvaluated (8),RS# (dbSNP) (9),nsv/esv (dbVar) (10),RCVaccession (11),PhenotypeIDS (12),PhenotypeList (13),Origin (14),OriginSimple (15),Assembly (16),
ChromosomeAccession (17),Chromosome (18),Start (19),Stop (20),ReferenceAllele (21),AlternateAllele (22),Cytogenetic (23),ReviewStatus (24),NumberSubmitters (25),Guidelines (26),TestedInGTR (27),OtherIDs (28),SubmitterCategories (29),VariationID (30),PositionVCF (31),
ReferenceAlleleVCF (32), AlternateAlleleVCF (33)


vep=G|missense_variant|MODERATE|OR4F5|ENSG00000186092|Transcript|ENST00000335137|protein_coding|1/1||ENST00000335137.3:c.44A>G|ENSP00000334393.3:p.Glu15Gly|44|44|15|E/G|gAa/gGa|rs781394307|1||1
Allele|Consequence    |IMPACT  |SYMBOL|Gene         |Feature_type|Feature       |BIOTYPE       |EXON|INTRON|HGVSc             |HGVSp|cDNA_position


,#AlleleID (0),Type (1),Name (2),GeneID (3),GeneSymbol (4),ClinicalSignificance (6),ClinSigSimple,Assembly (16),Chromosome (18),Start (19),Stop (20),ReferenceAllele (21),AlternateAllele (22)



"""
