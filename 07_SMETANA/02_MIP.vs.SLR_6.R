############################################################################
####
#### 
#### 2022.05.20 Toju
#### R 4.1.2
#### 
############################################################################

setwd("/Users/toju/Dropbox/YAJIMA_Daii/Statistics_R")

set.seed(123)

## -- Loading Function and Library
source('functions/functions.R')

load.lib( c('vegan', 'ggplot2', 'tidyr', 'cowplot','RColorBrewer', 'scales', 'ggsnippets'))
library(SpiecEasi)
library(igraph)
library(tidyverse)
library(graphlayouts)
library(ggforce)
library(scatterpie)
library(ggraph)
library(oaqc)
library(concaveman)
library(graphlayouts)
library(doParallel)
library(foreach)
library(ggplot2)
library(ggsci)
library(reshape2)
library(RColorBrewer)
library(cowplot)
library(imputeTS)

n.ex.rep <- 5 # number of tanks
sample.th <- 30
nrep.star <- 100
nlatent <- 10

load(sprintf("05_SpiecEasi/01_output/RData/SpiecEasi_th%s_latent%s.RData", sample.th, nlatent))

# -- Create directory to save
dir <- make.dir('07_SMETANA/02_output')

mip <- data.frame(readRDS("Table/mipMat.rds"))
mro <- data.frame(readRDS("Table/mroMat.rds"))

taxa <- read.table("Table/Taxa_list.pro_delspace.txt", header=T)

t_1 <- read.table("03_Surrogate/01_output/Table/Surrogate_ASVs_All_th30_Tank1_10000.txt", header=T)
t_2 <- read.table("03_Surrogate/01_output/Table/Surrogate_ASVs_All_th30_Tank2_10000.txt", header=T)
t_3 <- read.table("03_Surrogate/01_output/Table/Surrogate_ASVs_All_th30_Tank3_10000.txt", header=T)
t_4 <- read.table("03_Surrogate/01_output/Table/Surrogate_ASVs_All_th30_Tank4_10000.txt", header=T)
t_5 <- read.table("03_Surrogate/01_output/Table/Surrogate_ASVs_All_th30_Tank5_10000.txt", header=T)

colnames(t_1)[2] <- 'pcor.AS_1'
colnames(t_2)[2] <- 'pcor.AS_2'
colnames(t_3)[2] <- 'pcor.AS_3'
colnames(t_4)[2] <- 'pcor.AS_4'
colnames(t_5)[2] <- 'pcor.AS_5'

t_12 <- merge(t_1[,1:2], t_2[,1:2], by='ASV.ID', all=T)
t_123 <- merge(t_12, t_3[,1:2], by='ASV.ID', all=T)
t_1234 <- merge(t_123, t_4[,1:2], by='ASV.ID', all=T)
t_12345 <- merge(t_1234, t_5[,1:2], by='ASV.ID', all=T)

pcor.AS_mean <- rowMeans(t_12345[,2:ncol(t_12345)], na.rm=TRUE)
t_all <- cbind(t_12345, pcor.AS_mean)

asvs <- merge(t_all, taxa, by.x='ASV.ID', by.y='ID', all=F)


######################################################
# MB and SLR values

mb.ave <- list()
mb.ave.posi <- list()
mb.ave.nega <- list()
slr.posi <- list()
slr.nega <- list()
ig.mb.posi.w <- list()
ig.mb.nega.w <- list()
ig.slr.posi.w <- list()
ig.slr.nega.w <- list()
ig.mb.posi.uw <- list()
ig.mb.nega.uw <- list()
ig.slr.posi.uw <- list()
ig.slr.nega.uw <- list()

for(i in 1:n.ex.rep){
mb.ave[[i]] <- as.matrix(symBeta(mb.beta[[i]], mode = "ave"))

rownames(mb.ave[[i]]) <- colnames(dcomm[[i]])
colnames(mb.ave[[i]]) <- colnames(dcomm[[i]])
rownames(slr.cor[[i]]) <- colnames(dcomm[[i]])
colnames(slr.cor[[i]]) <- colnames(dcomm[[i]])

diag(mb.ave[[i]]) <- NA
mb.ave.posi[[i]] <- mb.ave[[i]]
mb.ave.posi[[i]][which(mb.ave.posi[[i]] <= 0)] <- 0
mb.ave.nega[[i]] <- mb.ave[[i]]
mb.ave.nega[[i]][which(mb.ave.nega[[i]] >= 0)] <- 0
ig.mb.posi.w[[i]] <- graph.adjacency(mb.ave.posi[[i]], mode='undirected', diag=FALSE, weighted= TRUE)
ig.mb.nega.w[[i]] <- graph.adjacency(-mb.ave.nega[[i]], mode='undirected', diag=FALSE, weighted= TRUE)

mb.binary.posi <- mb.ave.posi[[i]]
mb.binary.posi[which(mb.binary.posi > 0)] <- 1
mb.binary.nega <- mb.ave.nega[[i]]
mb.binary.nega[which(mb.binary.nega < 0)] <- -1
ig.mb.posi.uw[[i]] <- graph.adjacency(mb.binary.posi, mode='undirected', diag=FALSE, weighted= NULL)
ig.mb.nega.uw[[i]] <- graph.adjacency(-mb.binary.nega, mode='undirected', diag=FALSE, weighted= NULL)

diag(slr.cor[[i]]) <- NA
slr.posi[[i]] <- slr.cor[[i]]
slr.posi[[i]][which(slr.posi[[i]] <= 0)] <- 0
slr.nega[[i]] <- slr.cor[[i]]
slr.nega[[i]][which(slr.nega[[i]] >= 0)] <- 0
ig.slr.posi.w[[i]] <- graph.adjacency(slr.posi[[i]], mode='undirected', diag=FALSE, weighted= TRUE)
ig.slr.nega.w[[i]] <- graph.adjacency(-slr.nega[[i]], mode='undirected', diag=FALSE, weighted= TRUE)

slr.binary.posi <- slr.posi[[i]]
slr.binary.posi[which(slr.binary.posi > 0)] <- 1
slr.binary.nega <- slr.nega[[i]]
slr.binary.nega[which(slr.binary.nega < 0)] <- -1
ig.slr.posi.uw[[i]] <- graph.adjacency(slr.binary.posi, mode='undirected', diag=FALSE, weighted= NULL)
ig.slr.nega.uw[[i]] <- graph.adjacency(-slr.binary.nega, mode='undirected', diag=FALSE, weighted= NULL)
}

######################################################

list.mip <- rownames(mip)

asvlist <- unique(c(colnames(mb.ave[[1]]), colnames(mb.ave[[2]]), colnames(mb.ave[[3]]), colnames(mb.ave[[4]]), colnames(mb.ave[[5]])))

m <- list()
s <- list()

for(i in 1:n.ex.rep){
ID <- rownames(mb.ave[[i]])
m[[i]] <- mb.ave[[i]][ID %in% asvlist, ID %in% asvlist]

ID <- rownames(slr.cor[[i]])
s[[i]] <- slr.cor[[i]][ID %in% asvlist, ID %in% asvlist]
}


mb.sel <- list()
slr.sel <- list()

for(i in 1:n.ex.rep){
mb.sel[[i]] <- m[[i]][rownames(m[[i]]) %in% list.mip, rownames(m[[i]]) %in% list.mip]
slr.sel[[i]] <- s[[i]][rownames(s[[i]]) %in% list.mip, rownames(s[[i]]) %in% list.mip]
}

##########################################################
# ASVs in all tanks
# with X_0002

common12 <- intersect(rownames(mb.sel[[1]]), rownames(mb.sel[[2]]))
common123 <- intersect(common12, rownames(mb.sel[[3]]))
common1234 <- intersect(common123, rownames(mb.sel[[4]]))
common1235 <- intersect(common1234, rownames(mb.sel[[5]]))

mb_1 <- as.numeric(mb.sel[[1]][common1235,'X_0002'])
mb_2 <- as.numeric(mb.sel[[2]][common1235,'X_0002'])
mb_3 <- as.numeric(mb.sel[[3]][common1235,'X_0002'])
mb_4 <- as.numeric(mb.sel[[4]][common1235,'X_0002'])
mb_5 <- as.numeric(mb.sel[[5]][common1235,'X_0002'])
mb <- cbind(mb_1, mb_2, mb_3, mb_4, mb_5)
MB.mean_X_0002 <- rowMeans(mb, na.rm=TRUE)
ASV.ID <- common1235
MB_X_0002 <- cbind(ASV.ID, MB.mean_X_0002)

slr_1 <- as.numeric(slr.sel[[1]][common1235,'X_0002'])
slr_2 <- as.numeric(slr.sel[[2]][common1235,'X_0002'])
slr_3 <- as.numeric(slr.sel[[3]][common1235,'X_0002'])
slr_4 <- as.numeric(slr.sel[[4]][common1235,'X_0002'])
slr_5 <- as.numeric(slr.sel[[5]][common1235,'X_0002'])
slr <- cbind(slr_1, slr_2, slr_3, slr_4, slr_5)
SLR.mean_X_0002 <- rowMeans(slr, na.rm=TRUE)
SLR_X_0002 <- cbind(ASV.ID, SLR.mean_X_0002)

mip_X_0002 <- mip[,2]
mro_X_0002 <- mro[,2]

ASV.ID <- rownames(mip)
MIP_X_0002 <- cbind(ASV.ID, mip_X_0002)

ASV.ID <- rownames(mro)
MRO_X_0002 <- cbind(ASV.ID, mro_X_0002)

out <- merge(asvs, MB_X_0002, by="ASV.ID", all=FALSE)
out <- merge(out, SLR_X_0002, by="ASV.ID", all=FALSE)
out <- merge(out, MIP_X_0002, by="ASV.ID", all=FALSE)
out <- merge(out, MRO_X_0002, by="ASV.ID", all=FALSE)

write.table(out, file=sprintf('%s/Indices_X_0002_ASVs_AllTanks.txt', dir$tabledir), sep='\t', quote=F, row.names=F)


##########################################################
# ASVs appearing in at least one tank
# with X_0002

asvlist <- unique(rownames(mb.sel[[1]]), rownames(mb.sel[[2]]), rownames(mb.sel[[3]]), rownames(mb.sel[[4]]), rownames(mb.sel[[5]]))

mb_1 <- as.numeric(mb.sel[[1]][asvlist,'X_0002'])
mb_2 <- as.numeric(mb.sel[[2]][asvlist,'X_0002'])
mb_3 <- as.numeric(mb.sel[[3]][asvlist,'X_0002'])
mb_4 <- as.numeric(mb.sel[[4]][asvlist,'X_0002'])
mb_5 <- as.numeric(mb.sel[[5]][asvlist,'X_0002'])
mb <- cbind(mb_1, mb_2, mb_3, mb_4, mb_5)
MB.mean_X_0002 <- rowMeans(mb, na.rm=TRUE)
ASV.ID <- common1235
MB_X_0002 <- cbind(ASV.ID, MB.mean_X_0002)

slr_1 <- as.numeric(slr.sel[[1]][asvlist,'X_0002'])
slr_2 <- as.numeric(slr.sel[[2]][asvlist,'X_0002'])
slr_3 <- as.numeric(slr.sel[[3]][asvlist,'X_0002'])
slr_4 <- as.numeric(slr.sel[[4]][asvlist,'X_0002'])
slr_5 <- as.numeric(slr.sel[[5]][asvlist,'X_0002'])
slr <- cbind(slr_1, slr_2, slr_3, slr_4, slr_5)
SLR.mean_X_0002 <- rowMeans(slr, na.rm=TRUE)
SLR_X_0002 <- cbind(ASV.ID, SLR.mean_X_0002)

mip_X_0002 <- mip[,2]
mro_X_0002 <- mro[,2]

ASV.ID <- rownames(mip)
MIP_X_0002 <- cbind(ASV.ID, mip_X_0002)

ASV.ID <- rownames(mro)
MRO_X_0002 <- cbind(ASV.ID, mro_X_0002)

out <- merge(asvs, MB_X_0002, by="ASV.ID", all=FALSE)
out <- merge(out, SLR_X_0002, by="ASV.ID", all=FALSE)
out <- merge(out, MIP_X_0002, by="ASV.ID", all=FALSE)
out <- merge(out, MRO_X_0002, by="ASV.ID", all=FALSE)

write.table(out, file=sprintf('%s/Indices_X_0002_ASVs_AllASVs.txt', dir$tabledir), sep='\t', quote=F, row.names=F)

