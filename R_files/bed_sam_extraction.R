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


bed_sam_extraction <- function(sam_data){
   
   # Total reads in sam_data
   TOTAL.READS <- length(sam_data[[1]])
   
   # data structure to be returned
   #     name = name of read from sam file
   #     flag = flag from sam file
   #     chr = chromosome number
   #     start = genome coords of first CpG on read
   #     end =genome coords of last CpG on read
   #     z = genome coords of unmethylated CpG sites
   #     Z = genome coords of methylated CpG sites
   
   rtn_data <- list(name=vector('list',TOTAL.READS),
                    flag=vector('list',TOTAL.READS),
                    chr=vector('list',TOTAL.READS),
                    start=vector('list',TOTAL.READS),
                    end=vector('list',TOTAL.READS), 
                    z=vector('list',TOTAL.READS), 
                    Z=vector('list',TOTAL.READS))
   
   # Loop though all reads
   for (i in 1:TOTAL.READS){
      
      name <- as.character(sam_data[[1]][i]) # read name
      CIGAR <- as.character(sam_data[[6]][i]) # CIGAR string
      READ.LENGTH <- nchar(sam_data[[14]][i])-5 
      POS <- as.numeric(sam_data[[4]][i]) # start position of read
      m <- substr(sam_data[[14]][i],6,READ.LENGTH+5) # meth call string
      
      # Get the absolute coordinates of each nucleotide
      ar <- bed_absolute_reference(CIGAR, POS, READ.LENGTH)
      
      # Pull out the indices where there is a 'z' or 'Z'
      z.index <- gregexpr('[z]',m)[[1]]
      Z.index <- gregexpr('[Z]',m)[[1]]
      
      # Pull out the correspondng genome coords of z and Z
      if (!any(z.index==-1)){z.loc <- ar[z.index]} else {z.loc <- NA}
      if (!any(Z.index==-1)){Z.loc <- ar[Z.index]} else {Z.loc <- NA}
      
      rtn_data$name[i] <- sam_data[[1]][i]
      rtn_data$flag[i] <- as.numeric(sam_data[[2]][i])
      rtn_data$chr[i] <- sam_data[[3]][i]
      rtn_data$start[i] <- min(z.loc,Z.loc)
      rtn_data$end[i] <- max(z.loc,Z.loc)
      rtn_data$z[[i]] <- z.loc
      rtn_data$Z[[i]] <- Z.loc
      
   } # End loop over reads
   
   return(rtn_data)
}