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


plot_epiallele_distribution <- function(E){
  
  # number of epialleles   
  N <- nrow(E$epi)
  
  # number of sampling regions   
  R <- length(E$ind)
  
  # Normalise epiallele disttribution
  Ed <- E$dist
  for (r in 1:nrow(Ed)){Ed[r,] <- Ed[r,]/sum(Ed[r,])} 
  
  e.index <- 1:ncol(Ed)
  
  # prepare data frame for plotting
  
  epiallele <- NULL
  for (n in 1:N){
    epiallele <- c(epiallele,
                   rep(paste(as.character(E$epi[n,]),collapse=''),R))
  }
  
  data <- data.frame(epiallele=epiallele,
                     prop <- as.vector(Ed), 
                     region=rep(as.character(1:R),N))
  
  # plot distribution
  ggplot(data,aes(x=epiallele,y=prop,fill=factor(region)))+ geom_bar(stat="identity",position="dodge")+ xlab("")+ylab("Proportion of tissue sample") + theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank()) + scale_fill_brewer(palette="Set2") + coord_flip(ylim=c(0,1))
  
  #ggsave(file='fig2_R.pdf',width=10,height=11.5,units=c('cm'))
  
}