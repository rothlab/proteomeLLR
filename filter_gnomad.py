#!/usr/bin/env python3

"""
This script produces a cleaned, condensed version of the gnomad frequencies file where:

controls_nhomalt = Count of homozygous individuals in samples in the controls subset

AF_popmax = Maximum allele frequency across populations (excluding samples of Ashkenazi, Finnish, and indeterminate ancestry)

vep = Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position
|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL
|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|GMAF|AFR_MAF|AMR_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|ExAC_MAF
|ExAC_Adj_MAF|ExAC_AFR_MAF|ExAC_AMR_MAF|ExAC_EAS_MAF|ExAC_FIN_MAF|ExAC_NFE_MAF|ExAC_OTH_MAF|ExAC_SAS_MAF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|
HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info
"""

from xopen import xopen

gnomad_file = "gnomad.exomes.r2.1.1.sites.vcf.gz"
filtered_gnomad_file = "gnomad_frequencies.vcf"

if __name__ == "__main__":
	
	with xopen(gnomad_file, "rt") as gnomad_db:
		with open(filtered_gnomad_file, "w") as output:
			for entry in gnomad_db:
				if "missense_variant" in entry: 
					variant_info = entry.strip().split("\t")[-1].split(";")
					for info in variant_info:
						if (info.startswith("controls_nhomalt=")) or (info.startswith("AF_popmax")):
							output.write(info + ";")
						elif info.startswith("vep="):
							for vep_entry in info.split(","):
								if "missense" in vep_entry and "YES" in vep_entry: # missense variant on canonical transcript
									vep_info = vep_entry.strip().split("|")
									output.write(vep_info[3] + "|" + vep_info[11].split(":")[-1])
					output.write("\n")

"""
Example VEP section (entries split by commas):

vep=T|missense_variant|MODERATE|CCDC18|ENSG00000122483|Transcript|ENST00000334652|protein_coding|17/30||ENST00000334652.5:c.77C>T
|ENSP00000334084.5:p.Ser26Leu|2714|77|26|S/L|tCa/tTa||1||1||SNV|1|HGNC|30370|||||ENSP00000334084||F5H8K2|UPI000049DCA2||deleterious(0.01)
|probably_damaging(0.94)|Coiled-coils_(Ncoils):Coil&hmmpanther:PTHR18875&hmmpanther:PTHR18875:SF4||||||||||||||||||||||||||||||,T|missense_variant
|MODERATE|CCDC18|ENSG00000122483|Transcript|ENST00000338949|protein_coding|16/27||ENST00000338949.4:c.1457C>T|ENSP00000344380.4:p.Ser486Leu|2582
|1457|486|S/L|tCa/tTa||1||1||SNV|1|HGNC|30370|||||ENSP00000344380||E9PDA0|UPI00015E0BFC||deleterious(0.01)|benign(0.227)
|Coiled-coils_(Ncoils):Coil&hmmpanther:PTHR18875&hmmpanther:PTHR18875:SF4||||||||||||||||||||||||||||||,T|missense_variant|MODERATE|CCDC18
|ENSG00000122483|Transcript|ENST00000370276|protein_coding|17/29||ENST00000370276.1:c.2350C>T|ENSP00000359299.1:p.Ser784Leu|2350|2351|784|S/L
|tCa/tTa||1||1|cds_start_NF|SNV|1|HGNC|30370|YES||||ENSP00000359299|||UPI0001F78148||deleterious(0.01)|possibly_damaging(0.615)
|Coiled-coils_(Ncoils):Coil&hmmpanther:PTHR18875&hmmpanther:PTHR18875:SF4||||||||||||||||||||||||||||||

"""
