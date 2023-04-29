# proteomeLLR
This repository contains scripts to:
1) Generate **reference variant sets** for all genes with VARITY predictions
2) Add additional annotations to 
3) Perform clustering


## Generating Reference Sets

To generate reference variant sets for genes with VARITY predictions, you will need: 
1. A ClinVar variant summary file (which can be found at https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/)
2. A GnomAD variant summary file (which can be found at https://gnomad.broadinstitute.org/downloads/)
3. A list of VARITY supported IDs (from http://varity.varianteffect.org/)
4. All VARITY predictions (from http://varity.varianteffect.org/)    

>**NOTE**: The GnomAD file must be filtered to speed up the time it takes to generate reference sets:
>1. Change the name of [gnomad_file](https://github.com/rothlab/proteomeLLR/blob/e3d64daf983332e78d9e6207647d32af71e351f4/filter_gnomad.py#L19) to your input GnomAD variant summary file and [filtered_gnomad_file](https://github.com/rothlab/proteomeLLR/blob/e3d64daf983332e78d9e6207647d32af71e351f4/filter_gnomad.py#L20) to your desired output name
>2. Run `python filter_gnomad.py`

You can run the script **varity_LLRs.py** to generate reference sets. 

    Required arguments:
	    --clinvar              Path to ClinVar variant summary file
	    --gnomad               Path to filtered GnomAD variant summary file
	    --varity_ids           Path to list of VARITY supported IDs
	    --varity_predictions   Path to VARITY predictions
	    -o, --output           Output directory where reference sets should be written
	Optional arguments: 
		--num_hom			   GnomAD homozygote threshold (default: 1)
		--min_ref              Minimum number of reference pathogenic and benign variants (default: 11)    
		--maf_cutoff           MAF threshold to be considered benign (default: 0.001)

Once the reference sets are produced, you can modify the [calcLLR.R](https://github.com/rothlab/tileseqMave/blob/master/inst/scripts/calcLLR.R) script within TileseqMave to output the LLRs at 0.1 VARITY score intervals into a single output called LLR_intervals.csv. To do so, add the following code to the end of TileseqMave and update TileseqMave locally. 

    interval_scores = sapply(seq(0,1,0.1), llrObj$llr)
    cat(paste(strsplit(outprefix, "_")[[1]][1], ","), file="LLR_intervals.csv", append=TRUE)
    cat(interval_scores, file="LLR_intervals.csv", sep=",", append=TRUE)
    cat("\n", file="LLR_intervals.csv", append=TRUE)