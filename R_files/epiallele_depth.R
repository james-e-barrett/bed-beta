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


epiallele_depth <- function(Z){
   
   R <- length(Z)
   
   # ------------------------------------------------------------------ #
   # Match loci and generate matrix of sequencing depth
   # ------------------------------------------------------------------ #
   
   # Generate list of start locations
   start <- vector('list',R)
   start.union <- NULL
   for (r in 1:R){
      start[[r]] <- numeric(length(Z[[r]]))
      for (mu in 1:length(Z[[r]])){start[[r]][mu] <- colnames(Z[[r]][[mu]])[1]}
      start.union <- union(start.union,start[[r]])
   }
   
   # Index the overlapping regions
   N <- length(start.union)
   e.N <- matrix(NA,N,R)
   e.ind <- matrix(NA,N,R)
   
   pB <- txtProgressBar(min=1,max=N, width =50L, style = 3)
   for (n in 1:N){
   
      setTxtProgressBar(pB, n)
      
      # Index which regions are matched
      ind <- rep(NA,R)
      for (r in 1:R){
         if(any(start[[r]] == start.union[n])){ind[r] <- which(start[[r]] == start.union[n])}
      }
      
      e.ind[n,] <- ind
      for (r in 1:R){
         if(!is.na(ind[r])){
            e.N[n,r] <- nrow(Z[[r]][[ind[r]]])
         }
         
      }
      
   } # END loop
   close(pB)
   
   
   return(list(N = e.N, ind = e.ind))
     
}