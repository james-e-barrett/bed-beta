#!/bin/bash

# Run the second R script (combine tumour and normal samples from one patient and infer epialleles)

Rscript path_to-script/script2_epiallele_inference.R \
	--source_dir path_to_source \
	--tumour_input_dir /path_to_R1/example_R1_Z.Rdata,/path_to_R2/example_R2_Z.Rdata \
	--normal_input_dir /paht_to_N/example_N_Z.Rdata \
	--sample_id example \
	--output_dir path_to_output