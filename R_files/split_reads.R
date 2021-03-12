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


split_reads <- function(Y, N.MIN, D.MIN){
   
  # source files (seems to be necessary for parallelisation)
  for (src in dir('R_files')){
    source(paste('R_files/',src,sep=''))
  }
  
   M <- length(Y)
   Z <- vector('list',5*M)
   counter <-1
   
   cat("Splitting and filtering reads...\n")
   
   pB <- txtProgressBar(min=1,max=M, width =50L, style = 3)
   
   for (m in 1:M){
      setTxtProgressBar(pB, m)
      
      if ((nrow(Y[[m]])>=N.MIN) & (ncol(Y[[m]])>=D.MIN)){
         
         # Downsample very deep regions to 10,000 reads
         if (nrow(Y[[m]]) > 10000){
            Y.reduced <- Y[[m]][sample(1:nrow(Y[[m]]),10000,replace=F),]
            results <- partition2(Y.reduced, N.MIN, D.MIN)
         } else{
            results <- partition2(Y[[m]], N.MIN, D.MIN)   
         }
         n.z <- length(results)
         if(!is.null(results)){
            for (i in 1:n.z){
               Z[[counter]] <- results[[i]]
               counter <- counter + 1
            }
         }
      }
   }
   
   close(pB)
   
   Z.temp <- vector('list',counter-1)
   for (i in 1:(counter-1)){Z.temp[[i]] <- Z[[i]]}
   Z <- Z.temp
   
   return(Z)
   
}