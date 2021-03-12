#!/bin/bash

# Run the first R script (process a single sample from one patient in order to compute the loci that will be used for epiallele inference)

Rscript path_to_R_script/script1_loci_assembly.R \
	--source_dir path_to_R_source_code \
	--input_dir path_to_input_data \
	--sample_id exampleID \
	--output_dir path_to_output_directory 
