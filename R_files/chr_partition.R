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


chr_partition <- function(start,end){
   
   M <- length(start)
   
   partition <- rep(NA,M)
   
   right <- end[1]  
   counter <- 1

#    l-------r
#    s--------e   
#       
#    l-------r
#       s--------e   
# 
#    l-----------r
#      s---e
#             s--------e
   
   for (m in 1:M){
      
      if (start[m] <= right){
         partition[m] <- counter
         if (right < end[m]) right <- end[m]
      } else {
         counter <- counter + 1
         partition[m] <- counter
         right <- end[m]
      }   
   }
   
   return(partition)
   
}