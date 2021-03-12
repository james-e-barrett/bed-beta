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

bed_absolute_reference <- function(CIGAR, POS, READ.LENGTH){
   
#    CIGAR = A compact string that (partially) summarizes the alignment of the raw 
#            sequence read to the reference genome.
#    POS   = Left-most position within the reference genome where the alignment occurs
#    READ.LENGTH  = Template Length.
   
   # Total number of elements in CIGAR
   N.CIGAR <- nchar(CIGAR)
   
   # An empty vector of length READ.LENGTH
   absolute.reference <- numeric(READ.LENGTH)
   
   # Which elements of CIGAR contain a letter
   letter.location <- gregexpr('[MID]',CIGAR)[[1]]
   
   # d = total number of letters in CIGAR
   d <- length(letter.location)
   
   # An index to keep track of our current position on the read
   read.position <- 1
   
   # An index to keep track of our current position on the reference genome
   reference.position <- 1
   
   # Begin loop over the d letters
   for (mu in seq(1,d)){
   
      # Here we want to extract the number corresponding to the mu-th letter
      # The start position of the number in the CIGAR string
      if (mu>1){start = letter.location[mu-1]+1} else {start = 1}
      
      # The end position of the number in the CIGAR string
      stop = letter.location[mu]-1
      
      # len = the number itself corresponding to the d-th letter
      len <- as.numeric(substr(CIGAR,start=start,stop=stop))
      
      # If mu-th letter is 'M' continue to add indices
      if(substr(CIGAR,letter.location[mu],letter.location[mu])=='M'){
         absolute.reference[read.position:(read.position+len-1)] <- seq(reference.position,reference.position+len-1)
         read.position <- read.position + len
         reference.position <- reference.position + len
      }
      # If mu-th letter is 'I' add NA and advance read.position
      if(substr(CIGAR,letter.location[mu],letter.location[mu])=='I'){
         absolute.reference[read.position:(read.position+len-1)] <- rep(NA,len)
         read.position <- read.position + len
      }
      # If mu-th letter is 'D' advance reference.position only
      if(substr(CIGAR,letter.location[mu],letter.location[mu])=='D'){
         reference.position <- reference.position + len
      }
   } # End loop over mu
   
   # Shift indices to the starting postion (POS)
   absolute.reference <- absolute.reference + (POS-1)
   
   return(absolute.reference)
}