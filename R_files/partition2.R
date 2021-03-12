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


partition2 <- function(Y, N.MIN, D.MIN){
   
   a <- matrix(0,nrow=nrow(Y),ncol=ncol(Y))
   a[!is.na(Y)] <- 1
   H <- hamming(a)
   fit <- hclust(as.dist(H), method = 'ward.D')
   
   # ------------------------------------------------------------------------ #
   # Determine appropriate number of clusters
   # ------------------------------------------------------------------------ #
   
   cumsum.old <- sum(a)/(ncol(a)*nrow(a))
   
   Q <- 2
   while(T){
      W <- cutree(fit, k=Q)
      cumsum <- 0
      for (q in 1:Q){
         b <- a[W==q,,drop=FALSE]
         if (nrow(b)>1){c <- colSums(b)} else {c <-b}
         
         cumsum <- cumsum + mean(c[which(c!=0)])/sum(W==q)
      }
      cumsum <- cumsum/Q
      diff <- (cumsum-cumsum.old)/cumsum.old
      
      # BREAK if there are no mismatches. 
      if (cumsum > 0.75){
         #Q <- Q - 1
         break
      } else {
         
         cumsum.old <- cumsum
         Q <- Q + 1
      }
   }
   
   result <- vector('list',Q)
   W <- cutree(fit, k=Q)
   for (q in 1:Q){
      b <- a[W==q,,drop=FALSE]
      if (nrow(b)>1){c <- colSums(b)} else {c <-b}
      result[[q]] <- Y[W==q,which(c!=0),drop=FALSE]
   }
   
   # ------------------------------------------------------------------------ #
   # filter for n and d
   # ------------------------------------------------------------------------ #
   
   n.z <- length(result)
   result.tmp <- vector('list',n.z)
   counter <- 1
   for (i in 1:n.z){
      z <- result[[i]][,which(colSums(is.na(result[[i]]))/nrow(result[[i]]) < 0.25),drop=FALSE]
      D.PASS <- (ncol(z) >= D.MIN)
      N.PASS <- (nrow(z) >= N.MIN)
      NA.PASS <- (sum(is.na(z))/(ncol(z)*nrow(z)) < 0.25)
      if (N.PASS & D.PASS & NA.PASS){
         result.tmp[[counter]] <- z
         counter <- counter + 1
      }
   }
   Q <- counter - 1
   if (Q == 0){
      return(NULL)
   }
   result <- vector('list',Q)
   for (i in 1:Q){result[[i]] <- result.tmp[[i]]}
   
   
   
   # ------------------------------------------------------------------------ #
   # stitch together any regions with same start site
   # ------------------------------------------------------------------------ #
   start <- numeric(Q)
   for (q in 1:Q){start[q] <- colnames(result[[q]])[1]}
   start.u <- unique(start)
   
   results.tmp <- vector('list',length(start.u))
   
   for (q in 1:length(start.u)){
      ind <- which(start==start.u[q])
      if (length(ind)==1){
         results.tmp[[q]] <- result[[ind]]
      }else{
         overlap.col <- colnames(result[[ind[1]]])
         for (mu in 2:length(ind)){overlap.col <- intersect(overlap.col,colnames(result[[ind[mu]]]))}
         
         # bind regions together
         #results.tmp[[q]] <- NULL
         for (mu in 1:length(ind)){
            results.tmp[[q]] <- rbind(results.tmp[[q]], result[[ind[mu]]][,!is.na(match(colnames(result[[ind[mu]]]), overlap.col)),drop=F])
         }
      }
      
   } # end loop over unique start sites
   result <- results.tmp
   # ------- END stitch together any Z's with same start site -------- #
   
   return(result)
   
}

