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

# -- Create directory to save
dir <- make.dir('./02_NMDS_tSNE/01_output')

# -- Load data table

info <- read.table('Table/Sample.List_1.txt', header=T)
data  <- read.table('Table/seqtab.16SrRNAcount.pro.txt', header=T)

merged.mat <- merge(info, data, by='Sample.ID', all = FALSE, sort = TRUE)
rownames(merged.mat) <- merged.mat$Sample.ID
sml  <- merged.mat[, 1:3]
comm <- merged.mat[, 4:ncol(merged.mat)]

#taxa <- read.table('Table/Taxa_list.pro_delspace.txt', row.name=1, header=T)
taxa <- readRDS('Table/Taxa_list.rds')

l='Family'
comm.new <- t(Taxa.mat(comm[rownames(sml),], taxa, l))

##### NMDS

sink(sprintf('%s/tSNE_Relative_Family_log.txt', dir$tabledir))
nmds <- metaMDS(comm.new/rowSums(comm.new), distance='bray')
#nmds <- metaMDS(comm.new, distance='bray')

df <- cbind(sml, nmds$points[rownames(sml), ])

g1 <- ggplot(df) + 
	geom_point(aes(x=MDS1, y=MDS2, color=as.character(Tank), fill=after_scale(alpha(color, 0.7)), size=as.numeric(Day)), alpha=0.7, shape=21)+
    scale_size_area(max_size=4, breaks=c(1, 25, 50, 75, 100, 125), guide=guide_legend(ncol=1, title.position = 'top'), name="Day") + labs(color="Tank", x="NMDS 1", y="NMDS 2")
                     
plot(g1)

ggsave(plot=g1, filename=sprintf('%s/NMDS_Relative_Family.pdf', dir$figdir), h=5, w=6)


#### tSNE


#tsne.r <- Rtsne(comm.new, ,perplexity=50, check_duplicates = FALSE, verbose=TRUE)
tsne.r <- Rtsne(comm.new/rowSums(comm.new), check_duplicates = FALSE, verbose=TRUE)

df2 <- cbind(sml, tsne.r$Y)

g2 <- ggplot(df2) + 
	geom_point(aes(x=tsne.r$Y[,1], y=tsne.r$Y[,2], color=as.character(Tank), fill=after_scale(alpha(color, 0.7)), size=as.numeric(Day)), alpha=0.7, shape=21)+
    scale_size_area(max_size=4, breaks=c(1, 25, 50, 75, 100, 125), guide=guide_legend(ncol=1, title.position = 'top'), name="Day") + labs(color="Tank", x="Axis 1", y="Axis 2")

plot(g2)
               
ggsave(plot=g2, filename=sprintf('%s/tSNE_Relative_Family.pdf', dir$figdir), h=5, w=6)

sink()
