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


plot_meth <- function(cells){
   
   y <- cells
   D <- hamming(cells)
   cluster = hclust(as.dist(D), method = 'ward.D')
   y <- y[cluster$order,]
   
   #pdf('y1.pdf', width = 0.5*ncol(y), height=max(5,0.1*nrow(y)))
   
   x <- seq(1,ncol(y))
   z <- rep(1,ncol(y))
   plot(x,z,type='b',ylim=c(0,nrow(y)),col='white')
   
   
   for (i in seq(1,nrow(y))){
      
      x <- seq(1,ncol(y))
      z.on <- which(y[i,]==1)
      z.ind <- !is.na(y[i,])
      z <- rep(i,ncol(y))
      points(x[z.ind],z[z.ind],type='b',ylim=c(0,nrow(y)))
      points(z.on,rep(i,length(z.on)),type='p', pch=19)
   }
   
   #dev.off()
   
}