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


epiallele <- function(input){
   
   
   for (src in dir('R_files')){
      source(paste('R_files/',src,sep=''))
   }
   
   Z <- input$Z
   e.depth <- input$e.depth
   
   # ------------------------------------------------------------------ #
   # Fit model to regions with sufficient depth
   # ------------------------------------------------------------------ #
   
   regions <- which(apply(e.depth$N, 1, FUN='median', na.rm=T)>=100)
   
   Ed <- vector('list',length(regions))
   pB <- txtProgressBar(min=1,max=length(regions), width =50L, style = 3)
   counter <- 1
   for (n in regions){ 
      
      #setTxtProgressBar(pB, counter)
      cat(n,'\n')
      set.seed(32*n)
      
      # tryCatch({
      #    result <- epiallele_distance(Z, e.depth$ind[n,], n)
      #    Ed[[counter]]$dist <- result$epi.profile.margin
      #    Ed[[counter]]$epi <- result$X
      #    Ed[[counter]]$ind <- e.depth$ind[n,]
      #    Ed[[counter]]$N <- e.depth$N[n,]
      # }, error=function(e) NULL)
      
      result <- epiallele_distance(Z, e.depth$ind[n,], n)
      if(all(!is.na(result))){
         Ed[[counter]]$dist <- result$epi.profile.margin
         Ed[[counter]]$epi <- result$X
         Ed[[counter]]$ind <- e.depth$ind[n,]
         Ed[[counter]]$N <- e.depth$N[n,]
      } else {
         Ed[[counter]]$dist <- NA
         Ed[[counter]]$epi <- NA
         Ed[[counter]]$ind <- e.depth$ind[n,]
         Ed[[counter]]$N <- e.depth$N[n,]
      }
      
      counter <- counter + 1
   }
   close(pB)
   
   return(Ed)
   
}
