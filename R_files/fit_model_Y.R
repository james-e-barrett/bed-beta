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


fit_model_Y <- function(Y){
   
   # Wrapper code for fitting a model to DNAm sequence read data. Returns a fitted model 
   # (see MAPxw.R) and the number of subtypes.
   #
   # Inputs:
   #    Y:       N by 2*D character matrix where each row is a barcode pair and 
   #             each element must equal T,C,G,A. D=7 is the barcode length.
   #    counts:  numeric vector of length N containing the read count for each barcode pair
   #    
   #
   # Outputs:
   #    model:  list where element q is the output from the LL function
   #    rho:    no of latent components
   
   
   # Total number of reads
   N <- nrow(Y)
   # Preallocate a list for the output
   results <- list(model=vector('list',N), BIC=rep(NA,N))
   
   H <- hamming(Y)
   fit <- hclust(as.dist(H), method = 'ward.D')
   
   # ---------- Begin loop over N -------- #
   for (q in 1:N){
      
      # BREAK if there's only one barcode pair observed  
      if (N==1){
         results$rho <- 1
         break
      }
      
      # Fit model using LL function
      results$model[[q]] <- MAPxw(Y, q, fit)
      results$BIC[q] <- results$model[[q]]$BIC
      results$rho <- which.min(results$BIC) # Current best estimate for corrected read count
      
      # BREAK if there are no mismatches and q=1. This means the model has achieved a 
      # perfect fit so we can stop.
      # BREAK if there are no mismatches. 
      if (results$model[[q]]$mismatches==0){
         results$rho <- which.min(results$BIC[1:q])
         break
      }
      
      # Search at least until q=10, then check to see if we've hit the minimum BIC score.
      # The BIC is non monotonic so don't search in the last five values of q.
      if ((q>=10) & (which.min(results$BIC[1:q])<(q-5))) break
   }
   # ---------- End loop over N -------- #
   
   return(results)
}