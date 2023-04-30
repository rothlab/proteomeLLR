This repository contains scripts to:
1) Generate **reference variant sets** for all genes with VARITY predictions
2) Add additional annotations to CSV files with log-likelihood ratios (LLRs) at 0.1 VARITY score intervals
3) Cluster genes based on inheritance and functional activity, and attempt to build a naive bayes model to classify genes based on their LLR curves

## Generating Reference Sets

To generate reference variant sets for genes with VARITY predictions, you will need: 
1. A **ClinVar variant summary file** (which can be found at https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/)
2. A **GnomAD variant summary file** (which can be found at https://gnomad.broadinstitute.org/downloads/)
3. A list of **VARITY supported IDs** (from http://varity.varianteffect.org/)
4. All **VARITY predictions** (from http://varity.varianteffect.org/)    

>**NOTE**: The GnomAD file must be filtered to speed up the time it takes to make reference sets. To do this:
>1. Change the name of [gnomad_file](https://github.com/rothlab/proteomeLLR/blob/e3d64daf983332e78d9e6207647d32af71e351f4/filter_gnomad.py#L19) to the path of your input GnomAD variant summary file and [filtered_gnomad_file](https://github.com/rothlab/proteomeLLR/blob/e3d64daf983332e78d9e6207647d32af71e351f4/filter_gnomad.py#L20) to your desired output name
>2. Run `python filter_gnomad.py`

You can run the script **varity_LLRs.py** to generate reference sets. 

    Required arguments:
	    --clinvar              Path to ClinVar variant summary file
	    --gnomad               Path to filtered GnomAD variant summary file
	    --varity_ids           Path to list of VARITY supported IDs
	    --varity_predictions   Path to VARITY predictions
	    -o, --output           Output directory where reference sets should be written
	
	Optional arguments: 
	    --num_hom              GnomAD homozygote threshold (default: 1)
	    --min_ref              Minimum number of reference pathogenic and benign variants (default: 11)    
	    --maf_cutoff           MAF threshold to be considered benign (default: 0.001)

Example usage: 

    python varity_LLRs.py --clinvar data/variant_summary.txt --gnomad data/gnomad_frequencies.vcf.gz --varity_ids data/varity_supported_ids.txt --varity_predictions data/varity_all_predictions.txt --output reference_sets --num_hom 2

Once the reference sets are produced, you can modify the [calcLLR.R](https://github.com/rothlab/tileseqMave/blob/master/inst/scripts/calcLLR.R) script within TileseqMave to output LLRs at 0.1 VARITY score intervals into a single output called **LLR_intervals.csv**. To do so, add the following code to the end of calcLLR.R and update TileseqMave locally. 

    interval_scores = sapply(seq(0,1,0.1), llrObj$llr)
    cat(paste(strsplit(outprefix, "_")[[1]][1], ","), file="LLR_intervals.csv", append=TRUE)
    cat(interval_scores, file="LLR_intervals.csv", sep=",", append=TRUE)
    cat("\n", file="LLR_intervals.csv", append=TRUE)

You can then write a quick script such as the one below to run this for all genes with reference sets: 

    REFDIR=path_to_directory_with_reference_sets
    REFS=($(ls $REFDIR | grep reference.csv))
    cd $REFDIR
    for REF in ${REFS[@]}
    do
	    GENE=$(echo $REF | sed 's/_reference.csv//')
	    MAP=$(echo $REF | sed 's/_reference.csv/_map.csv/')
	    tsm calcLLR $MAP $REF --kernel gaussian #can change parameters here as needed
	done

## Adding Additional Annotations

Additional information can be added to LLR_intervals.csv, such as the mode of inheritance of diseases associated with genes in the file or functional annotations for each gene. Two scripts have been written to add such annotations: **[add_functional_annotations.py](https://github.com/rothlab/proteomeLLR/blob/main/add_functional_annotations.py "add_functional_annotations.py")** and **[add_inheritance_annotations.py](https://github.com/rothlab/proteomeLLR/blob/main/add_inheritance_annotations.py "add_inheritance_annotations.py")**. 

Example add_inheritance_annotations.py usage:

    python add_inheritance_annotations.py -i path_to_input_csv -g path_to_gencc_file -o output_filename

 - *path_to_input_csv* refers to the path of LLR_intervals.csv
 - *path_to_gencc_file* refers to a GenCC TSV annotation file (found here https://search.thegencc.org/download)
	 >**NOTE**: Modify the GenCC file to only include the 'gene_symbol', 'disease_title', 'classification_title' and 'moi_title' columns
 - *output_filename* refers to the desired output filename

Example add_inheritance_annotations.py usage:

    python add_functional_annotations.py -i path_to_input_csv -g path_to_david_annotation_cluster_file -o output_filename

 - *path_to_input_csv* refers to the path of LLR_intervals.csv
 - *path_to_david_annotation_cluster_file* refers to the DAVID functional annotation output file (can be made here: https://david.ncifcrf.gov/tools.jsp)
 - *output_filename* refers to the desired output filename 

## Clustering

The script [**clusterLLR.R**](https://github.com/rothlab/proteomeLLR/blob/main/clusterLLR.R) takes LLR_intervals.csv and will attempt to cluster the genes based on their LLR curves and output PCA and UMAP plots

Additionally, if provided with an LLR_intervals.csv file with additional annotations, it attempts to build a classification model using naive bayes to see if VARITY performs better for certain subgroups of genes (ex. inheritance, different functional sub classes of genes). You can change the input LLR_intervals.csv file directly within the script. 