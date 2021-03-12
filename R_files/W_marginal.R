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


W_marginal <- function(Y,X,W){
      
   Q <- nrow(X)
   N <- nrow(Y)
   
   N.plus <- numeric(Q)
   N.minus <- numeric(Q)
   
   for (q in seq(1,Q)){
      N.q <- sum(W==q)
      term1 <- colSums(Y[W==q,,drop=FALSE], na.rm = TRUE)
      term3 <- colSums(!is.na(Y[W==q,,drop=FALSE]), na.rm = TRUE)
      term2 <- term3 - term1
      N.plus[q] <- sum(term1[X[q,]==1]) + sum(term2[X[q,]==0])
      N.minus[q] <- sum(term3) - N.plus[q]
   }
   
   matches <- sum(N.plus)
   mismatches <- sum(N.minus)
   
   # change w
   w.margin <- matrix(NA,N,Q)

   for (i in 1:N){
      w <- numeric(Q)
      for (mu in 1:Q){

         # subtract MAP matches and add current matches
         w.matches <- matches - sum(Y[i,] == X[W[i],], na.rm = TRUE) + 
            sum(Y[i,] == X[mu,], na.rm = TRUE)
         w.mismatches <- mismatches - sum(Y[i,] != X[W[i],], na.rm = TRUE) + 
            sum(Y[i,] != X[mu,], na.rm = TRUE)
   
         #beta <- max(mismatches,1)/(matches+mismatches)
         beta <- 0.05

         w[mu] <- w.matches*log(1-beta) + w.mismatches*log(beta)
      }
      w.margin[i,] <- exp(w - mean(w))/sum(exp(w - mean(w)))
   }
   
   return(w.margin)
}
