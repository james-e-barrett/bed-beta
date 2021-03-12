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


normal_decontamination <- function(Ed, purity.estimate){
   
   Ed.norm <- Ed
   
   # normalise Ed
   for (r in 1:nrow(Ed)){Ed.norm[r,] <- Ed[r,]/sum(Ed[r,])}
   
   N <- rowSums(Ed, na.rm = TRUE)
   Ed.norm[N==0,] <- rep(NA,ncol(Ed.norm))
   
   # Subtract epiallales that can be attributed to normal tissue
   
   Ed.decon <- matrix(0,8,ncol(Ed))
   for(mu in 1:ncol(Ed)){
      Ed.decon[1:7,mu] <- (Ed.norm[1:7,mu] - (1-purity.estimate)*Ed.norm[8,mu])/purity.estimate
   }
   
   Ed.decon[Ed.decon<0] <- 0
   Ed.decon[Ed.decon>1] <- 1
   
   Ed.decon <- Ed.decon/rowSums(Ed.decon)
   Ed.decon[8,] <- rep(NA,ncol(Ed.norm))
   
   return(Ed.decon)
   
}