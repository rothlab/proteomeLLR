import gzip
import re
import pickle
import csv
import pandas as pd
import operator

fh = open('disease_genes_dictionary.pkl', 'rb')
disease_genes_dict = pickle.load(fh)

fh = open('ReferenceSet.csv', 'w', newline='')
writer = csv.writer(fh)
writer.writerow(['GeneSymbol', 'GeneID', 'hgvsp', 'ClinicalSignificance', 'referenceSet', 'controls_AC', 'controls_AN', 'maf', 'hom', 'source'])

with gzip.open('gnomad.exomes.r2.1.1.sites.1.vcf.gz' , 'rb') as f:
    for line in f:
        line = line.decode("utf-8")
        if line.startswith('#'):
            pass
        else:
            line_list = line.split('\t')
            info = line_list[7]
            info_list = info.split(';')
            for item in info_list:
                if item.startswith('vep='):
                    veps_list = item.split(',')
                    for vep in veps_list:
                        if '|missense_variant|' in vep and '|YES|' in vep:
                            controls_ac = ''
                            controls_an = ''
                            controls_af = -1.0
                            controls_nhomalt = -1.0
                            source = ''
                            clinical_significance = -1
                            reference_set = ''
                            hgvsp = re.findall("p\.[a-zA-Z]{3}\d+[a-zA-Z]{3}", vep)[0]
                            gene_symbol = vep.split('|')[3]
                            for gene in disease_genes_dict.keys():
                                if gene_symbol == gene[0]:
                                    gene_id = gene[1]
                                    source = 'Clinvar'
                                    hgvsps = list(map(operator.itemgetter(0), disease_genes_dict[gene]))
                                    if hgvsp in hgvsps:
                                        clinical_significance = [item[1] for item in disease_genes_dict[gene] if item[0] == hgvsp]
                                        clinical_significance = clinical_significance[0]
                                    else:
                                        source = 'GnomAD'
                                        clinical_significance = 'NA'
                                    for inf in info_list:
                                        if inf.startswith('controls_AC='):
                                            controls_ac = inf[12:]
                                        elif inf.startswith('controls_AN='):
                                            controls_an = inf[12:]
                                        elif inf.startswith('controls_AF='):
                                            controls_af = float(inf[12:])
                                        elif inf.startswith('controls_nhomalt='):
                                            controls_nhomalt = float(inf[17:])
                                        
                                    if source == 'Clinvar':
                                        if clinical_significance == 'Pathogenic' or clinical_significance == 'Likely pathogenic' or clinical_significance == 'Pathogenic/Likely pathogenic':
                                            reference_set = 'Positive'
                                        else:
                                            reference_set = 'Negative'
                                        writer.writerow([gene_symbol, gene_id, hgvsp, clinical_significance, reference_set, controls_ac, controls_an, controls_af, controls_nhomalt, source])
                                        break
                                    elif source == 'GnomAD':
                                        if controls_af >= 0.001 or controls_nhomalt >= 1.0:
                                            reference_set = 'Negative'
                                            writer.writerow([gene_symbol, gene_id, hgvsp, clinical_significance, reference_set, controls_ac, controls_an, controls_af, controls_nhomalt, source])
                                            break