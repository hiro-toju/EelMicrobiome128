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
dir <- make.dir('./02_NMDS_tSNE/03_output')

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

env1 <- data.frame(cbind(pH=pH$pH.value, DO=DO$DO.value, AS=AS$AS.value))
rownames(env1) <- Sample.ID
env2 <- na.omit(env1)


merged.mat <- merge(info, data, by='Sample.ID', all = FALSE, sort = TRUE)
rownames(merged.mat) <- merged.mat$Sample.ID
sml  <- merged.mat[, 1:3]
comm <- merged.mat[, 4:ncol(merged.mat)]

#taxa <- read.table('Table/Taxa_list.pro_delspace.txt', row.name=1, header=T)
taxa <- readRDS('Table/Taxa_list.rds')

l='Family'
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
en_coord_cont = as.data.frame(scores(envNMDS, "vectors"))* ordiArrowMul(envNMDS)


#ordiplot(nmds, type="n")
#orditorp(nmds, display="sites", pch=16)
#plot(envNMDS, p.max=0.01)

#intersect(common, rownames(sml))

#df <- cbind(sml[common,], nmds$points)
df <- cbind(sml[common,], data.scores)

write.table(df, file=sprintf('%s/NMDS_Relative_Representative_Family.txt', dir$tabledir), sep='\t', quote=F, row.names=F)


############################################################################

Unidentified <- comm.new[, colnames(comm.new)=="Unidentified"]
comm.exU <- comm.new[, colnames(comm.new)!="Unidentified"]

comm.exU.ordered <- comm.exU[, order(colSums(comm.exU), decreasing=TRUE)]
Others <- rowSums(comm.exU.ordered[, 31:ncol(comm.exU.ordered)])
comm.top30.OU <- cbind(comm.exU.ordered[, 1:30], Others, Unidentified)


# Setting colors

col1 <- brewer.pal(7, "Set1")
col2 <- brewer.pal(7, "Set2")
col3 <- brewer.pal(8, "Set3")
col4 <- brewer.pal(8, "Pastel1")

col.list <- rev(c(col2, col3, col1, col4, "grey90", "grey80"))

############################################################################

sel <- c("S_1029", "S_4024", "S_2028", "S_4111", "S_4073", "S_5065", "S_5041")


x <- data.frame(comm.top30.OU[sel,])

data <- cbind(Sample = rownames(x), x)

df <- melt(data)
df$category <- factor(df$variable, levels=rev(colnames(x)))

g <- ggplot(df,
	aes (
		x = factor(Sample, levels=sel),
		y = value,
		fill = category
	)
)

g1 <- g + geom_bar(stat = "identity", position='fill') + scale_fill_manual(values = col.list) +
	theme(axis.text.x = element_text(angle=45, hjust=1)) +
	theme(legend.position = "none") + 
	theme(axis.title.x = element_blank(), axis.title.y = element_blank())

plot(g1)

ggsave(plot=g1, filename=sprintf('%s/NMDS_Relative_Representative_Family.pdf', dir$figdir), h=3, w=4)


g2 <- g + geom_bar(stat = "identity", position='fill') + scale_fill_manual(values = col.list) +
	theme(axis.text.x = element_text(angle=45, hjust=1))+ 
	theme(axis.title.x = element_blank(), axis.title.y = element_blank())

plot(g2)

ggsave(plot=g2, filename=sprintf('%s/NMDS_Relative_Representative_Family_legend.pdf', dir$figdir), h=4, w=7)
