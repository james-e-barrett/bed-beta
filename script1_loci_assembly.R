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
   "input_dir", "i" , 1 , "character",
   "output_dir", "o" , 1 , "character" ,
   "sample_id", "d" , 1 , "character" ,
   "help", "h" , 0 , "logical" 
   ), ncol=4 , byrow=TRUE )

opt = getopt(spec , opt = commandArgs(TRUE))

if( !is.null(opt$help) ) { 
   cat(getopt(spec , usage=TRUE))
   q(status=1)
}

if( is.null(opt$source_dir) || is.null(opt$input_dir) || is.null(opt$output_dir) ) { 
   cat(getopt(spec , usage=TRUE))
   q() 
}

if ( is.null(opt$sample_id)){
   opt$sample_id <- '_'
}

#----------------------------------------#
# Prepare for parallelisation
#----------------------------------------#

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
# Load output from processed sam files
#----------------------------------------#

file.cols <- vector('list', 7)
file.cols[[1]] <- 'name'
file.cols[[2]] <- 'flag1'
file.cols[[3]] <- 'flag2'
file.cols[[4]] <- 'start'
file.cols[[5]] <- 'end'
file.cols[[6]] <- 'z'
file.cols[[7]] <- 'Z'


str <- strsplit(opt$input_dir,split='/')
sample_id <- str[[1]][length(str[[1]])]

# list of files in input directory
file_list <- dir(opt$input_dir)

# vector of missing chr files
missing_chr <- match(paste(sample_id,'_chr',seq(1:22),'.txt',sep=''),
                     dir(opt$input_dir))
  
if(sum(is.na(missing_chr)) >0){
  warning('At least one chr file missing from input directory')
}

# load the chr files
setwd(opt$input_dir)
cat("Begin read input chr files...\n")
sam_data <- vector('list',sum(!is.na(missing_chr)))
for (chr in which(!is.na(missing_chr))){
  sam_data[[chr]] <- scan(paste(sample_id,'_chr',chr,'.txt',sep=''),
                          what = file.cols)  
}
cat("Read input chr files complete.\n\n")


#----------------------------------------#
# Assemble loci from raw data
#----------------------------------------#

cat("Begin data assembly...\n")

setwd(opt$source_dir)
Y <- parLapply(cl, sam_data, fun = data_assembly)

save(Y,file=paste(opt$output_dir,paste(opt$sample_id,"Y.Rdata",sep="_"),sep="/"))

cat("Data assembly complete.\n\n")


#----------------------------------------#
# Filter and split loci
#----------------------------------------#

N.MIN <- 10
D.MIN <- 5

cat("Begin split reads...\n")

Z <- parLapply(cl, Y, fun = split_reads, N.MIN, D.MIN)

save(Z,file=paste(opt$output_dir,paste(opt$sample_id,"Z.Rdata",sep="_"),sep="/"))

cat("Split reads complete.\n\n")


#----------------------------------------#
# Logging
#----------------------------------------#

epiallele_logging(sam_data, Y, Z, N.MIN, D.MIN, opt)

