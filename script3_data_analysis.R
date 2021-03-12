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

#----------------------------------------#
# Source R_files
#----------------------------------------#

for (src in dir('R_files')){
    source(paste('R_files/',src,sep=''))
}

library(ggplot2)

#----------------------------------------#
# Clean out Null and NA fields
#----------------------------------------#

for (chr in 1:22){
  counter <- 1
  while(counter <= length(Edata[[chr]])){
    
    if(is.null(Edata[[chr]][[counter]])){
      
      Edata[[chr]][[counter]] <- NULL
      
    } else {counter <- counter+1}
  }
}


for (chr in 1:22){
    counter <- 1
    while(counter <= length(Edata[[chr]])){
        
        if (is.null(Edata[[chr]][[counter]])){
            counter = counter + 1
        } else {
            if(is.na(Edata[[chr]][[counter]]$dist)){
                
                Edata[[chr]][[counter]] <- NULL
                
            } else {counter <- counter+1}
        }
        
    }
}

#=============================================#
# Summary stats
#=============================================#

Q <- vector('list',22)
for (chr in 1:22){
  Q[[chr]] <- numeric(length(Edata[[chr]]))
  for (i in 1:length(Edata[[chr]])){
    Q[[chr]][i] <- ncol(Edata[[chr]][[i]]$dist)
  }
}
summary(as.factor(unlist(Q)))
plot(as.factor(unlist(Q)))


#=============================================#
# Plot epiallele distributions
#=============================================#
plot_epiallele_distribution(Edata[[1]][[3]])

#=============================================#
# Compute overall distance matrix
#=============================================#

R <- 8

# distance matrix
D <- matrix(0,R,R)
# counts no of non missing D entries
counter <- matrix(0,R,R)

for (chr in 1:22){
  for (i in 1:length(Edata[[chr]])){
    
    if (!is.null(Edata[[chr]][[i]])){
      
      NCOL <- (ncol(Edata[[chr]][[i]]$dist) > 1)
      R.NA <- (sum(is.na(Edata[[chr]][[i]]$N)) < 3)
      NORMAL.MISSING <- !is.na(Edata[[chr]][[i]]$N[R])
      
      if (NCOL & R.NA & NORMAL.MISSING){
        
        #Ed <- Edata[[chr]][[i]]$dist
        Ed <- normal_decontamination(Ed, purity.estimate)
        
        # normalisation
        for (r in 1:nrow(Ed)){Ed[r,] <- Ed[r,]/sum(Ed[r,])}
        
        d <- as.matrix(dist(Ed))
        counter[!is.na(d)] <- counter[!is.na(d)] +1
        d[is.na(d)] <- 0

        D <- D + d
      }
    }
  }
}


D <- (D/counter)
D <- D[1:R,1:R]
plot(fastme.bal(D))

#=============================================#
# Heatmap
#=============================================#

# IMPORTANT: recompute the matrix D above with normal_decontamination commented
# and for (r in 1:nrow(Ed)){Ed[r,] <- Ed[r,]/sum(Ed[r,])} uncommented

d <- 0
for (chr in 1:22){
   for (i in 1:length(Edata[[chr]])){
      
      NCOL <- (ncol(Edata[[chr]][[i]]$dist) > 1)
      R.NA <- (sum(is.na(Edata[[chr]][[i]]$N)) < 3)
      NORMAL.MISSING <- !is.na(Edata[[chr]][[i]]$N[R])
      
      if (NCOL & R.NA & NORMAL.MISSING){
         d <- d + ncol(Edata[[chr]][[i]]$dist)
      }
   }
}


data <- matrix(NA,R,d)
pB <- txtProgressBar(min=1,max=22, width =50L, style = 3)
counter <- 1
for (chr in 1:22){
   setTxtProgressBar(pB, chr)
   for (i in 1:length(Edata[[chr]])){
      
      NCOL <- (ncol(Edata[[chr]][[i]]$dist) > 2)
      R.NA <- (sum(is.na(Edata[[chr]][[i]]$N)) < 3)
      NORMAL.MISSING <- !is.na(Edata[[chr]][[i]]$N[R])
      
      if (NCOL & R.NA & NORMAL.MISSING){
         #Ed <- Edata[[chr]][[i]]$dist
         Ed <- normal_decontamination(Ed, purity.estimate)
         for (r in 1:nrow(Ed)){Ed[r,] <- Ed[r,]/sum(Ed[r,])}
         data[,counter:(counter+ncol(Edata[[chr]][[i]]$dist)-1)] <- Ed
         counter <- counter + ncol(Edata[[chr]][[i]]$dist)
      }
   }
}
close(pB)


# Heatmap ----------------- #

if (!require("gplots")) {
   install.packages("gplots", dependencies = TRUE)
   library(gplots)
}
if (!require("RColorBrewer")) {
   install.packages("RColorBrewer", dependencies = TRUE)
   library(RColorBrewer)
}

my_palette <- colorRampPalette(c("cornsilk1", "darkseagreen2", "darkslategray3", "dodgerblue4"))(n = 299)
#my_palette <- colorRampPalette(c("maroon","sienna1", "seashell2", "lightsteelblue2","royalblue3"))(n = 299)

#data <- data[-8,]

data.sd <- apply(data,2,'sd')
ord <- order(data.sd, decreasing=TRUE)
data <- data[,ord]
data <- data[,1:200]
rownames(data) <- c('R1','R2','R3','R4','R5','R6','R7','Normal')
#'data <- data[c(1,2,7,4,3,5,6,8),]


cluster = hclust(as.dist(D), method = "ward.D")
cluster2 = hclust(dist(t(data)), method = "ward.D")

#heatmap.2(data, Rowv = as.dendrogram(cluster),Colv = as.dendrogram(cluster2), dendrogram="both",trace="none",col=my_palette)

heatmap.2(data, Rowv = NULL,Colv = as.dendrogram(cluster2), dendrogram="none",trace="none",col=my_palette)



#=============================================#
# Purity estimate
#=============================================#

d <- 0
for (chr in 1:22){
   for (i in 1:length(Edata[[chr]])){
      
      NCOL <- (ncol(Edata[[chr]][[i]]$dist) > 1)
      R.NA <- (sum(is.na(Edata[[chr]][[i]]$N)) < 3)
      NORMAL.MISSING <- !is.na(Edata[[chr]][[i]]$N[R])
      
      if (NCOL & R.NA & NORMAL.MISSING){
         d <- d + 1
      }
   }
}

purity <- matrix(NA,R,d)
counter <- 1
for (chr in 1:22){
   for (i in 1:length(Edata[[chr]])){
      
      # NCOL <- (ncol(Edata[[chr]][[i]]$dist) > 1)
      # R.NA <- (sum(is.na(Edata[[chr]][[i]]$N)) < 3)
      # NORMAL.MISSING <- !is.na(Edata[[chr]][[i]]$N[R])
      # 
      # if (NCOL & R.NA & NORMAL.MISSING){
      #    Ed <- Edata[[chr]][[i]]$dist
      #    Ed.decon <- normal_decontamination(Ed)
      #    purity[,counter] <- rowSums(abs(Ed.decon))/2
      #    counter <- counter + 1
      # }
     NCOL <- (ncol(Edata[[chr]][[i]]$dist) > 1)
     R.NA <- (sum(is.na(Edata[[chr]][[i]]$N)) < 3)
     NORMAL.MISSING <- !is.na(Edata[[chr]][[i]]$N[R])
     
     if (NCOL & R.NA & NORMAL.MISSING){
       Ed <- Edata[[chr]][[i]]$dist
       Ed <- Ed/rowSums(Ed)
       purity[,counter] <- 0.5*rowSums(abs(t(t(Ed)-Ed[R,])))

       counter <- counter + 1
       
     }
   }
}
# proportion of reads that CANNOT be attributed to normal
rowMeans(purity,na.rm=T)

d <- density(purity[1,],na.rm=T,bw=0.02)
qplot(d$x,d$y, geom='line') + ylim(0,2) + xlab(expression(xi)) + ylab('Density') +xlim(-0.05,1)
ggsave(file='density_R1.pdf',width=10,height=5,units=c('cm'))

d <- density(purity[2,],na.rm=T,bw=0.02)
qplot(d$x,d$y, geom='line') + ylim(0,2) + xlab(expression(xi)) + ylab('Density')  +xlim(-0.05,1) 
ggsave(file='density_R2.pdf',width=10,height=5,units=c('cm'))

d <- density(purity[3,],na.rm=T,bw=0.02)
qplot(d$x,d$y, geom='line') + ylim(0,2) + xlab(expression(xi)) + ylab('Density')  +xlim(-0.05,1)
ggsave(file='density_R3.pdf',width=10,height=5,units=c('cm'))

d <- density(purity[4,],na.rm=T,bw=0.02)
qplot(d$x,d$y, geom='line') + ylim(0,2) + xlab(expression(xi)) + ylab('Density')  +xlim(-0.05,1)
ggsave(file='density_R4.pdf',width=10,height=5,units=c('cm'))

d <- density(purity[5,],na.rm=T,bw=0.02)
qplot(d$x,d$y, geom='line') + ylim(0,2) + xlab(expression(xi)) + ylab('Density')  +xlim(-0.05,1)
ggsave(file='density_R5.pdf',width=10,height=5,units=c('cm'))

d <- density(purity[6,],na.rm=T,bw=0.02)
qplot(d$x,d$y, geom='line') + ylim(0,2) + xlab(expression(xi)) + ylab('Density')  +xlim(-0.05,1)
ggsave(file='density_R6.pdf',width=10,height=5,units=c('cm'))

d <- density(purity[7,],na.rm=T,bw=0.02)
qplot(d$x,d$y, geom='line') + ylim(0,2) + xlab(expression(xi)) + ylab('Density')  +xlim(-0.05,1)
ggsave(file='density_R7.pdf',width=10,height=5,units=c('cm'))

#=============================================#
# Disorder
#=============================================#

d <- 0
for (chr in 1:22){
   for (i in 1:length(Edata[[chr]])){
      
      NCOL <- (ncol(Edata[[chr]][[i]]$dist) > 1)
      R.NA <- (sum(is.na(Edata[[chr]][[i]]$N)) < 3)
      NORMAL.MISSING <- !is.na(Edata[[chr]][[i]]$N[R])
      
      if (NCOL & R.NA & NORMAL.MISSING){
         d <- d + 1
      }
   }
}

entropy <- matrix(NA,R,d)
counter <- 1
for (chr in 1:22){
   for (i in 1:length(Edata[[chr]])){
      
      NCOL <- (ncol(Edata[[chr]][[i]]$dist) > 1)
      R.NA <- (sum(is.na(Edata[[chr]][[i]]$N)) < 3)
      NORMAL.MISSING <- !is.na(Edata[[chr]][[i]]$N[R])
      
      if (NCOL & R.NA & NORMAL.MISSING){
         Ed <- Edata[[chr]][[i]]$dist
         for (r in 1:nrow(Ed)){Ed[r,] <- Ed[r,]/sum(Ed[r,])}
         entropy[,counter] <- -rowSums(Ed*apply(Ed, 2, log2))
         counter <- counter + 1
      }
   }
}

bdat <- data.frame(out = c(entropy[1,], entropy[2,], entropy[3,], entropy[4,], entropy[5,]),
                   hl = c(rep('R1',d), rep('R2',d), rep('R3',d), rep('R4',d), rep('N',d)))
#boxplot(out~hl,bdat, main = 'pHER3', ylab = 'pHER3', xlab = 'HER1 and HER3')
p <- ggplot(bdat, aes(factor(hl), out))
p + geom_boxplot() + ylab('Entropy')
ggsave(file='disorder.pdf',width=20,height=8,units=c('cm'))


