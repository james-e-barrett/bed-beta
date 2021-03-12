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


MAPxw <- function(Y, Q, fit){
  
  MIN.N <- 25
  MIN.CpG <- 10
  
  d <- ncol(Y)
  N <- nrow(Y)
  #    if (N <= Q) return(NULL)
  #    if (d < MIN.CpG) return(NULL)
  #    if (N < MIN.N) return(NULL)
  
  #    H <- hamming(Y)
  #    fit <- hclust(as.dist(H), method = 'ward.D')
  W.new <- cutree(fit, k=Q)
  
  
  #---------------------------------------------------------------#
  # MAP X and W solutions
  #---------------------------------------------------------------#
  
  MAX <- 10
  counter <- 0
  
  while (counter < MAX){
    
    X <- matrix(0,Q,d)
    for (q in seq(1,Q)){
      X[q,colSums((Y[W.new==q,,drop=FALSE]), na.rm = TRUE)/colSums(!is.na(Y[W.new==q,,drop=FALSE]))>0.5] <- 1
    }
    
    W.old <- W.new
    for (i in seq(1,N)){
      W.new[i] <- which.max(colSums(t(X)==Y[i,],na.rm = TRUE))
    }
    
    if (identical(W.old,W.new)) break
    counter <- counter + 1
  }
  W <- W.new
  
  #---------------------------------------------------------------#
  # alpha and gamma (used to find MAP beta)
  #---------------------------------------------------------------#
  
  N.plus <- numeric(Q)
  N.minus <- numeric(Q)
  
  for (q in seq(1,Q)){
    N.q <- sum(W==q)
    term1 <- colSums(Y[W==q,,drop=FALSE], na.rm = TRUE)
    term3 <- colSums(!is.na(Y[W==q,,drop=FALSE]), na.rm = TRUE)
    term2 <- term3 - term1
    N.plus[q] <- sum(term1[X[q,]==1]) + sum(term2[X[q,]==0]) # counts matches between Y and X. sum is over mu
    N.minus[q] <- sum(term3) - N.plus[q]
  }
  
  matches <- sum(N.plus)
  mismatches <- sum(N.minus)
  
  beta <- max(mismatches,1)/(matches+mismatches) # Noise hyperparameter
  #beta <- 0.05
  
  # Log likelihood
  LL <- matches*log(1-beta) + mismatches*log(beta)
  # Bayes information criterion score
  BIC <- -2*LL + 2*Q*(d)
  
  # Return a list with any relevant quantities
  return(list(W=W,X=X,matches=matches,mismatches=mismatches,beta=beta,BIC=BIC))
  
  
}