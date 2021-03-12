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


epiallele_distance <- function(Z, ind, rseed){
   
   set.seed(314*rseed)
   
   R <- length(ind)   
   
   # ------------------------------------------------------------------ #
   # Downsample to median
   # ------------------------------------------------------------------ #
   
   N <- rep(NA,R)
   for (r in 1:R){if(!is.na(ind[r])){N[r] <- nrow(Z[[r]][[ind[r]]])}}
   med <- round(median(N,na.rm=TRUE))
   
   # ------------------------------------------------------------------ #
   # Skip regions of extreme coverage (more than 10,000 reads)
   # ------------------------------------------------------------------ #
   if (sum(N,na.rm=T)>10000){ return(NA) }
   
   
   # ------------------------------------------------------------------ #
   # Fit model
   # ------------------------------------------------------------------ #
   
   # See which loci overlap within the matched regions
   overlap.col <- colnames(Z[[which(!is.na(ind))[1]]][[ind[which(!is.na(ind))[1]]]])
   for (r in 2:R){
      if(!is.na(ind[r])){overlap.col <- intersect(overlap.col, colnames(Z[[r]][[ind[r]]]))}
   }
   
   # dataset of all reads from pooling the eight regions
   Ze <- NULL
   for (r in 1:R){
      if(!is.na(N[r])){
         if (N[r] > med){
            #Ze <- rbind(Ze, Z[[r]][[ind[r]]][sample(1:N[r],med,replace=T),!is.na(match(colnames(Z[[r]][[ind[r]]]),overlap.col)),drop=FALSE])
            Ze <- rbind(Ze, Z[[r]][[ind[r]]][,!is.na(match(colnames(Z[[r]][[ind[r]]]),overlap.col)),drop=FALSE])
         } else {
            Ze <- rbind(Ze, Z[[r]][[ind[r]]][,!is.na(match(colnames(Z[[r]][[ind[r]]]),overlap.col)),drop=FALSE])
         }
      }
   }
   
   # region indexes which tumour region each row of Ze comes from
   region <- NULL
   for (r in 1:R){
      if(!is.na(N[r])){
         if(N[r] > med){
            #region <- c(region, rep(r,med))
            region <- c(region, rep(r,N[r]))
         } else {
            region <- c(region, rep(r,N[r]))
         }
      }
   }
   
   model <- fit_model_Y(Ze)
   
   # ------------------------------------------------------------------ #
   # Discard low frequency epialleles
   # ------------------------------------------------------------------ # 
   
   MIN.FREQ <- 0.01
   
   n <- nrow(Ze)
   
   # epiallele frequencies
   x <- as.numeric(summary(as.factor(model$model[[model$rho]]$W)))/n
   
   # epialleles with sufficient frequency
   X <- model$model[[model$rho]]$X[x>MIN.FREQ,,drop=F]
   
   # If Q=1 break
   if(nrow(X)==1){
      result = list(X = X,epi.profile.margin = matrix(N,ncol=1))
      return(result)
   }
   
   # it may happen than no epialleles exceed min freq
   if(nrow(X)==0){return(NA)}
   
   # compute new W
   W <- numeric(n)
   for (i in seq(1,n)){
      W[i] <- which.max(colSums(t(X)==Ze[i,],na.rm = TRUE))
   }
   
   # ------------------------------------------------------------------ #
   # Compute distance
   # ------------------------------------------------------------------ # 
   
   w.margin <- W_marginal(Ze, X, W)
   epi.profile.margin <- matrix(NA, R, nrow(X))
   for (r in 1:R){
      if(ncol(w.margin)==1){
         epi.profile.margin[r,] <- sum(w.margin[region==r,])
      }else{
         epi.profile.margin[r,] <- colSums(w.margin[region==r,])
      }
   }
   
   result = list(X = X,epi.profile.margin = epi.profile.margin)
   return(result)
   
   
}

