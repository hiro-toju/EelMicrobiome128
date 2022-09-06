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
dir <- make.dir('./02_NMDS_tSNE/04_output')

# -- Load data table

info <- read.table('Table/Sample.List_1.txt', header=T)
data  <- read.table('Table/seqtab.16SrRNAcount.pro.txt', header=T)

param <- as.matrix(read.table("Table/Parameters_3.txt", header=T, na.strings = "NA"))
param2 <- melt(param)
colnames(param2) <- c("Day", "Parameter", "value")
pH <- param2[grep("pH", param2$Parameter), ]
DO <- param2[grep("DO", param2$Parameter), ]
AS <- param2[grep("AS", param2$Parameter), ]
colnames(pH)[3] <- 'pH.value'
colnames(DO)[3] <- 'DO.value'
colnames(AS)[3] <- 'AS.value'

l1 <- paste("S_1", subset(formatC(pH$Day, width=3, flag="0"), pH$Parameter=="T1_pH"), sep="")
l2 <- paste("S_2", subset(formatC(pH$Day, width=3, flag="0"), pH$Parameter=="T2_pH"), sep="")
l3 <- paste("S_3", subset(formatC(pH$Day, width=3, flag="0"), pH$Parameter=="T3_pH"), sep="")
l4 <- paste("S_4", subset(formatC(pH$Day, width=3, flag="0"), pH$Parameter=="T4_pH"), sep="")
l5 <- paste("S_5", subset(formatC(pH$Day, width=3, flag="0"), pH$Parameter=="T5_pH"), sep="")	
	
Sample.ID <- c(l1, l2, l3, l4, l5)

env1 <- data.frame(cbind(pH=pH$pH.value, DO=DO$DO.value, Eel.Activity =AS$AS.value))
rownames(env1) <- Sample.ID
env2 <- na.omit(env1)


merged.mat <- merge(info, data, by='Sample.ID', all = FALSE, sort = TRUE)
rownames(merged.mat) <- merged.mat$Sample.ID
sml  <- merged.mat[, 1:3]
comm <- merged.mat[, 4:ncol(merged.mat)]

#taxa <- read.table('Table/Taxa_list.pro_delspace.txt', row.name=1, header=T)
taxa <- readRDS('Table/Taxa_list.rds')

l='Genus'
comm.new <- t(Taxa.mat(comm[rownames(sml),], taxa, l))

common <- intersect(rownames(comm.new), rownames(env2))

comm2 <- comm.new[common,]
env3 <- env2[common,]

##### NMDS

nmds <- metaMDS(comm2/rowSums(comm2), k=2, distance='bray')
nmds

envNMDS <- envfit(nmds, env3, permutations=10000)
envNMDS

data.scores = as.data.frame(scores(nmds)$sites)
species.scores = as.data.frame(scores(nmds, "species"))
en_coord_cont = as.data.frame(scores(envNMDS, "vectors"))* ordiArrowMul(envNMDS)

genlist <- c("Flavobacterium", "Cetobacterium", "Edaphobaculum", )

genvec <- species.scores[genlist,]

#ordiplot(nmds, type="n")
#orditorp(nmds, display="sites", pch=16)
#plot(envNMDS, p.max=0.01)

#intersect(common, rownames(sml))

#df <- cbind(sml[common,], nmds$points)
df <- cbind(sml[common,], data.scores)


g1 <- ggplot(df) + 
	geom_point(aes(x=NMDS1, y=NMDS2, color=as.character(Tank), fill=after_scale(alpha(color, 0.7)), size=as.numeric(Day)), alpha=0.7, shape=21)+
    scale_size_area(max_size=4, breaks=c(1, 25, 50, 75, 100, 125), guide=guide_legend(ncol=1, title.position = 'top'), name="Day") + labs(color="Tank", x="NMDS 1", y="NMDS 2")+
    
geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
       data = genvec, size =0.7, alpha = 0.5, colour = "red")+
       
geom_text(data = genvec, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
       fontface = "bold", label = row.names(genvec)) 
                     
plot(g1)

ggsave(plot=g1, filename=sprintf('%s/NMDS_Relative_Vectors_Genus.pdf', dir$figdir), h=5, w=6)



