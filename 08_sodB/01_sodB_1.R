############################################################################
####
#### R script for Fujita (2019)
####
#### Visualizing community assembly
#### 2019.12.09 Fujita
#### 2022.02.24 Toju
#### R 4.1.2
#### Set working directory of 'MTS' folder -- setwd('~/Desktop/Microbiome_TimeSeries/MTS/')
#### 
############################################################################

setwd("/Users/toju/Dropbox/YAJIMA_Daii/Statistics_R")

library(Rtsne)

set.seed(123)

## -- Loading Function and Library
source('functions/functions.R')
load.lib( c('vegan', 'ggplot2', 'tidyr', 'cowplot','RColorBrewer', 'scales', 'ggsnippets'))
library(parallel)
library(RcppThread)
library(RcppEigen)
library(doParallel)
library(foreach)
library(ggplot2)
library(ggsci)
library(imputeTS)
library(rEDM)
library(reshape2)

# -- Create directory to save
dir <- make.dir('./08_sodB/01_output')

# -- Load data table

data  <- read.table('Table/sodB_mat.txt', header=T, row.names=1)

info <- data[,1]
comm <- data[,-1]

Nonpathogenic <- colSums(subset(comm, info=="nonpathogenic"))
Pathogenic <- colSums(subset(comm, info=="pathogenic"))
commcl <- rbind(Nonpathogenic, Pathogenic)

rowSums(commcl)

l1 <- paste("S_1", formatC(1:128, width=3, flag="0"), sep="")
l2 <- paste("S_2", formatC(1:128, width=3, flag="0"), sep="")
l3 <- paste("S_3", formatC(1:128, width=3, flag="0"), sep="")
l4 <- paste("S_4", formatC(1:128, width=3, flag="0"), sep="")
l5 <- paste("S_5", formatC(1:128, width=3, flag="0"), sep="")
Sample.ID <- c(l1, l2, l3, l4, l5)

#d1 <- data.frame(Sample.ID, Day=rep(1:128, times=5))
d1 <- data.frame(Sample.ID)
d2 <- data.frame(Sample.ID=colnames(commcl), t(commcl))
d3 <- merge(d1, d2, by="Sample.ID", all.x=T)
d3[is.na(d3)] <- 0


############################################################################

col.list <- c("skyblue3", "lightpink3")

############################################################################


df <- melt(d3)
df$category <- factor(df$variable, levels=colnames(d3))

g <- list()

for (i in 1:5) {
	df2 <- df[grep(sprintf("S_%s", i), df$Sample.ID),]
	
g[[i]] <- ggplot(df2,
	aes (
		x = Sample.ID,
		y = as.numeric(value),
		fill = category
	)
)

g[[i]] <- g[[i]] + geom_bar(stat = "identity") + scale_fill_manual(values = col.list) +
	theme(axis.text.x = element_blank()) +
	theme(legend.position = "right") + 
	labs(x= "Day", subtitle=sprintf("Tank %s", i), y= "Number of sodB reads")
}

ggsave(plot=plot_grid(plotlist=g, nrow=5, byrow=FALSE),
       filename=sprintf("%s/Timeseries_sodB.pdf", dir$figdir), h=11.69, w=8.27)
    

