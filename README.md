NOTE: this is a beta version and not under active development.

Code to run the bayesian epiallele detection (BED) pipeline.

## Step 1: convert to sam file

Each sample needs to converted to a name sorted .sam file. 

``samtools view -h -o path_to_sam_files/LTX001/LTX001_N_top_n6dupsRemoved_sorted.sam /home/regmjeb@ad.ucl.ac.uk/cifs/CI_FeberTanic_TracerX/release_170808/LTX001/RRBS/Nugen/LTX001_N_top_n6dupsRemoved_sorted.bam``

## Step 2: run perl script

This will process each read and pull out the relevant methylation calls (taking into acount indels etc.). The script outputs one text file for each chr. These text files are the inputs for the next step. There is a version for both single and paired end data. This is definitely not the most efficient way to do this but it's what we hacked together at the time...

## Step 3: run the first R script

Use `bed_script1.sh` to run the first R script. This script will load the outputs from the previous script. This code will chop the reads into distinct loci. The output is a .Rdata list structure that contains the the reads (represented as 0's and 1's where 1 means methylated) at each locus.

## Step 4: run the second R script

Use `bed_script2.sh` to run the second R script. This script will take all of the tumour regions and the matched normal region, combine them, and infer what epialleles are present at each locus.

## Step 5: use script3_data_analysis.R

This is the final bit of the analysis where a phylogenetic tree was obtained and some heatmaps etc. This will depend on what's needed for the analysis.
