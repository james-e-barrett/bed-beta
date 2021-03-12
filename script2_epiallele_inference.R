    Bayesian epiallele detection
    Copyright (C) 2019 James E. Barrett (regmjeb@ucl.ac.uk)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

#----------------------------------------#
# Parse command line args
#----------------------------------------#

library(getopt)
library(parallel)

spec = matrix( c(
  "source_dir", "s" , 1 , "character" ,
  "tumour_input_dir", "t" , 1 , "character",
  "normal_input_dir", "n" , 1 , "character",
  "output_dir", "o" , 1 , "character" ,
  "sample_id", "i" , 1 , "character" ,
  "help", "h" , 0 , "logical" 
), ncol=4 , byrow=TRUE )

opt = getopt(spec , opt = commandArgs(TRUE))

if( !is.null(opt$help) ) { 
  cat(getopt(spec , usage=TRUE))
  q(status=1)
}

if( is.null(opt$source_dir) || is.null(opt$normal_input_dir) || is.null(opt$tumour_input_dir) || is.null(opt$output_dir) ) { 
  cat(getopt(spec , usage=TRUE))
  q() 
}

if ( is.null(opt$sample_id)){
  opt$sample_id <- '_'
}


# ----------------------------------------#
# Prepare for parallelisation
# ----------------------------------------#

NCORES <- 22
cl <- makeCluster(NCORES)


#----------------------------------------#
# Source R_files
#----------------------------------------#

setwd(opt$source_dir)

for (src in dir('R_files')){
  source(paste('R_files/',src,sep=''))
}

clusterExport(cl=cl, varlist=ls())

#----------------------------------------#
# Combine tumour regions and normal sample
#----------------------------------------#
# merged list structure should be indexed as Z[[chr]][[r]][[i]]

cat("\nBegin merge samples...\n")

CHR <- 22
total_samples <- length(strsplit(opt$tumour_input_dir,split=',')[[1]]) + 1
Z_merge <- vector('list', CHR)

# load & merge tumour samples
for (r in 1:(total_samples-1)){
  load(strsplit(opt$tumour_input_dir,split=',')[[1]][r])
  for (chr in 1:CHR){
    Z_merge[[chr]][[r]] <- Z[[chr]]
  }
}

# load & merge normal sample
load(opt$normal_input_dir)
for (chr in 1:CHR){
  Z_merge[[chr]][[total_samples]] <- Z[[chr]]
}

Z <- Z_merge
save(Z,file=paste(opt$output_dir,paste(opt$sample_id,'_merged','_Z.Rdata',sep=""),sep="/"))

cat("Merge samples complete.\n\n")

#----------------------------------------#
# Quality control plots
#----------------------------------------#

cat("Begin QC plots...\n")

epiallele_qc(Z, opt)

cat("QC plots complete.\n\n")


#----------------------------------------#
# Compute loci depths and indices
#----------------------------------------#

cat("Begin indexing epiallele loci...\n")

Edepth <- parLapply(cl, Z, fun = epiallele_depth)
setwd(opt$output_dir)
save('Edepth',file=paste(opt$output_dir,paste(opt$sample_id,'_Edepth.Rdata',sep=""),sep="/"))

cat("Indexing epiallele loci complete.\n\n")

#----------------------------------------#
# Infer epialleles
#----------------------------------------#

cat("Begin epiallele inference...\n")

input <- vector('list',length(Z))
for (chr in 1:length(Z)){
    input[[chr]] <- list(Z=Z[[chr]],e.depth=Edepth[[chr]])
}

Edata <- parLapply(cl, input, fun = epiallele)
#Edata <- lapply(input, FUN = epiallele)
save('Edata',file=paste(opt$output_dir,paste(opt$sample_id,'_Edata.Rdata',sep=""),sep="/"))

cat("Epiallele inference complete.\n\n")
