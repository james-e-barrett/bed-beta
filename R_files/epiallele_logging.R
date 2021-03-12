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


epiallele_logging <- function(sam_data, Y, Z, N.MIN, D.MIN, opt){
    
   filename <- paste(opt$output_dir,'/',
                     opt$sample_id,"_epiallele_log.txt",sep="")
   
    # Log total number of reads from sam file
    total_reads <- 0
    for (chr in 1:length(sam_data)){
        total_reads <- total_reads + length(sam_data[[chr]][[1]])
    }
    
    #==============================================#
    # Before filtering and splitting
    #==============================================#
    
    # Log number of loci inferred
    total_loci <- 0
    for (chr in 1:length(Y)){
        total_loci <- total_loci + length(Y[[chr]])
    }
    
    loci_depth <- numeric(total_loci)
    loci_width <- numeric(total_loci)
    counter <- 1
    for (chr in 1:length(Y)){
        for (i in 1:length(Y[[chr]])){
            loci_depth[counter] <- nrow(Y[[chr]][[i]])
            loci_width[counter] <- ncol(Y[[chr]][[i]])
            counter <- counter + 1
        }
    }
    
    #==============================================#
    # Create log file
    #==============================================#
    
    # Create log file
    str <- paste('=======================',
                 '\nEpiallele log file',
                 '\n=======================',
                 '\n\nGenerated on: ', Sys.time(),
                 '\nsource_dir: ', opt$source_dir,
                 '\ninput_dir: ', opt$input_dir,
                 '\noutput_dir: ', opt$output_dir,
                 '\n\n=======================',
                 '\nBefore filtering',
                 '\n=======================',
                 '\n\nTotal reads from sam file: ', total_reads,
                 '\nTotal loci: ', total_loci,sep='')
    
    write(str, file=filename, append = FALSE)
    
    write('\nSummary of locus depth:', file=filename, append = TRUE)
    capture.output(summary(loci_depth), file=filename, append = TRUE)
    
    write('\nSummary of locus width:', file=filename, append = TRUE)
    capture.output(summary(loci_width), file=filename, append = TRUE)
    
    str <- paste('\n=======================',
                 '\nAfter filtering',
                 '\n=======================',sep='')
    write(str, file=filename, append = TRUE)
    
    #==============================================#
    # After filtering and splitting
    #==============================================#
    
    # Log number of loci inferred
    total_loci <- 0
    for (chr in 1:length(Z)){
        total_loci <- total_loci + length(Z[[chr]])
    }
    
    loci_depth <- numeric(total_loci)
    loci_width <- numeric(total_loci)
    counter <- 1
    for (chr in 1:length(Z)){
        for (i in 1:length(Z[[chr]])){
            loci_depth[counter] <- nrow(Z[[chr]][[i]])
            loci_width[counter] <- ncol(Z[[chr]][[i]])
            counter <- counter + 1
        }
    }
    
    discard <- as.character(100*(1- sum(loci_depth)/total_reads))
    
    str <- paste('\nMinimum locus depth: ',N.MIN,
                 '\nMinimum locus width: ',D.MIN,
                 '\n\nTotal loci after filtering and splitting: ', total_loci,
                 '\nNumber of reads discarded: ',
                 total_reads-sum(loci_depth),
                 ' (',substr(discard,1,4),'%)',sep='')
    
    write(str, file=filename, append = TRUE)
    
    write('\nSummary of locus depth:', file=filename, append = TRUE)
    capture.output(summary(loci_depth), file=filename, append = TRUE)
    
    write('\nSummary of locus width:', file=filename, append = TRUE)
    capture.output(summary(loci_width), file=filename, append = TRUE)
    
    str <- paste('\n=======================',
                 '\nEnd log',
                 '\n=======================',sep='')
    write(str, file=filename, append = TRUE)
    
    
}