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


data_assembly <- function(sam_data){
  
  #    Flags 99 and 147 mean that the genome coordinate of a CpG aligns with the C (forward reads).
  #    Flags 83 and 163 mean that the coordinate alingns with the G (reverse reads).
  #    Therefore reverse reads will have 1 bp subtracted from start, end, and coords.
  #    After this forward and reverse reads are pooled together.
  
  #    sam_data[[1]] <- 'name'
  #    sam_data[[2]] <- 'flag1'
  #    sam_data[[3]] <- 'flag2'
  #    sam_data[[4]] <- 'start'
  #    sam_data[[5]] <- 'end'
  #    sam_data[[6]] <- 'z'
  #    sam_data[[7]] <- 'Z'
  
  # source files (seems to be necessary for parallelisation)
  for (src in dir('R_files')){
    source(paste('R_files/',src,sep=''))
  }
  
  # ----------------------------------------------------------- #
  # Extract relevant strand
  #-------------------------------------------------------------#
  
  data <- sam_data
  data[[4]] <- as.numeric(data[[4]])
  data[[5]] <- as.numeric(data[[5]])
  
  # Remove NA starts and ends
  na <- which(is.na(data[[5]])|is.na(data[[5]]))
  for (mu in c(1:7)){data[[mu]] <- data[[mu]][-na]}
  names(data) <- c('name','flag1','flag2','start','end','z','Z')   
  
  # Subtract 1 from reverse reads
  data[[4]][which(data[[2]]=="83")] <- data[[4]][which(data[[2]]=="83")] - 1
  data[[5]][which(data[[2]]=="83")] <- data[[5]][which(data[[2]]=="83")] - 1
  
  
  # ----------------------------------------------------------- #
  # Generate partition maps
  #-------------------------------------------------------------#
  
  start <- data[[4]]
  end <- data[[5]]
  ord <- order(start)
  start <- start[ord]
  end <- end[ord]
  for (mu in c(1:7)){data[[mu]] <- data[[mu]][ord]}
  partition <- chr_partition(start, end)
  
  
  # ----------------------------------------------------------- #
  # Assemble datasets
  #-------------------------------------------------------------#
  
  # No of partitions
  P <- partition[length(partition)]
  
  data.Y <- vector('list',P)
  
  cat('Assembling datasets...\n')
  pB <- txtProgressBar(min=1,max=P, width =50L, style = 3)
  
  for (p in 1:P){
    
    setTxtProgressBar(pB, p)
    
    index <- which(partition==p)
    
    on.index <- NULL
    off.index <- NULL
    
    for(mu in 1:length(index)){
      split.on <- as.numeric(strsplit(data$Z[index[mu]], split=':')[[1]])
      split.off <- as.numeric(strsplit(data$z[index[mu]], split=':')[[1]])
      if (data$flag1[index[mu]] == "83"){
        split.on <- split.on - 1
        split.off <- split.off - 1
      }
      on.index <- union(on.index, split.on)
      off.index <- union(off.index, split.off)
    }
    on.index <- on.index[!is.na(on.index)]
    off.index <- off.index[!is.na(off.index)]
    
    CpG.loci <- union(on.index, off.index)
    CpG.loci <- CpG.loci[order(CpG.loci)]
    
    Y <- matrix(NA, length(index), length(CpG.loci))
    colnames(Y) <- CpG.loci
    
    for(mu in seq(1, length(index))){
      split.on <- as.numeric(strsplit(data$Z[index[mu]], split=':')[[1]])
      split.off <- as.numeric(strsplit(data$z[index[mu]], split=':')[[1]])
      if (data$flag1[index[mu]] == "83"){
        split.on <- split.on - 1
        split.off <- split.off - 1
      }
      Y[mu,match(split.on,CpG.loci)] <- 1
      Y[mu,match(split.off,CpG.loci)] <- 0 
    }
    
    data.Y[[p]] <- Y
  }
  
  close(pB)
  return(data.Y)
  
}