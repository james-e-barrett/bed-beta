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


epiallele_qc <- function(Z,opt){
    
    R <- length(Z[[1]])
    
    #====================================#
    # CpG distribution
    #====================================#
    
    cpg_x <- matrix(NA, nrow=512, ncol=R)
    cpg_y <- matrix(NA, nrow=512, ncol=R)
    loci_x <- matrix(NA, nrow=512, ncol=R)
    loci_y <- matrix(NA, nrow=512, ncol=R)
    
    for (r in 1:R){
        
        total_cols <- 0 
        total_loci <- 0
        
        for (chr in 1:length(Z)){
            total_loci <- total_loci + length(Z[[chr]][[r]])
            for (i in 1:length(Z[[chr]][[r]])){
                total_cols <- total_cols + ncol(Z[[chr]][[r]][[i]])
            }
        }
        
        mean_cpg <- rep(NA, total_cols)
        mean_loci <- rep(NA, total_loci)
        strt <- 1
        stp <- 1
        for (chr in 1:length(Z)){
            for (i in 1:length(Z[[chr]][[r]])){
                stp <- strt + ncol(Z[[chr]][[r]][[i]]) - 1
                mean_cpg[strt:stp] <- colMeans(Z[[chr]][[r]][[i]],na.rm=TRUE)
                mean_loci[i] <- mean(colMeans(Z[[chr]][[r]][[i]],na.rm=TRUE))
                strt <- stp + 1
            }    
        }
        
        d_cpg <- density(mean_cpg, bw=0.02, na.rm=TRUE)    
        d_loci <- density(mean_loci, bw=0.02, na.rm=TRUE)    
    
        cpg_x[,r] <- d_cpg$x
        cpg_y[,r] <- d_cpg$y
        loci_x[,r] <- d_loci$x
        loci_y[,r] <- d_loci$y
    }
    
    setwd(opt$output_dir)
    # png(filename = paste(opt$sample_id,'_QC_CpG_mean.png',sep=''),
    #     width = 800, height = 400, units = 'px',type="cairo")
    bitmap(file = paste(opt$sample_id,'_QC_CpG_mean.png',sep=''),
           type='png16m',
           width = 800,
           height = 400,
           units = 'px')
    
    matplot(cpg_x,cpg_y, type='l',
            xlim = c(-0.05,1.05),
            col = 1:R,
            xlab = 'Methylation', 
            ylab= 'Density',
            main=paste(opt$sample_id,
                       ' Mean methylation per CpG',sep=''))
    legend('top',
           col=1:R,
           lty=rep(1,R),
           legend = c(paste('R',1:(R-1),sep=''),'N'))
    grid()
    
    dev.off()
    
    # png(filename = paste(opt$sample_id,'_QC_loci_mean.png',sep=''),
    #     width = 800, height = 400, units = 'px',type="cairo")
    bitmap(file = paste(opt$sample_id,'_QC_loci_mean.png',sep=''),
           type='png16m',
           width = 800,
           height = 400,
           units = 'px')
    
    matplot(loci_x,loci_y, type='l',
            xlim = c(-0.05,1.05),
            col = 1:R,
            xlab = 'Methylation', 
            ylab= 'Density',
            main=paste(opt$sample_id,
                       ' Mean methylation per locus',sep=''))
    legend('top',
           col=1:R,
           lty=rep(1,R),
           legend = c(paste('R',1:(R-1),sep=''),'N'))
    grid()
    
    dev.off()
       
}
