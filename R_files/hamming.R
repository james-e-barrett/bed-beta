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


hamming <- function(Y){
   
   d <- ncol(Y)
#       N <- nrow(Y)
#       
#       H <- matrix(0,N,N)
#       
#       for (i in seq(1,N)){
#          for (j in seq(1,N)){
#             H[i,j] <- sum(Y[i,] != Y[j,], na.rm = TRUE)
#          }
#       }
   
   B <- Y
   B[is.na(Y)] <- 0
   B1 <- (B==1)
   
   C <- Y
   C[is.na(Y)] <- 1
   C0 <- (C==0)
   
   H <- (B1%*%t(C0) + C0%*%t(B1))/d
   
   return(H)
   
}
