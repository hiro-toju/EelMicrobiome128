############################################################################
####
#### R script for Fujita (2019)
####
#### SpiecEasi
#### 2022.02.25 Toju
#### 2022.04.12 Toju
#### R 4.1.2
#### 
############################################################################

setwd("/Users/toju/Dropbox/YAJIMA_Daii/Statistics_R")

set.seed(123)

n.ex.rep <- 5 # number of tanks
sample.th <- 30
nrep.star <- 100
nlatent <- 10

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

load(sprintf("05_SpiecEasi/01_output/RData/SpiecEasi_th%s_latent%s.RData", sample.th, nlatent))

# -- Create directory to save
dir <- make.dir('07_SMETANA/03_output')

mip <- data.frame(readRDS("Table/mip.rds"))
mro <- data.frame(readRDS("Table/mro.rds"))



######################################################

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
######################################################
## MIP and MRO

diag(mro) <- NA
diag(mip) <- NA


######################################################

mro.v <- c(unlist(mro))
mip.v <- c(unlist(mip))
mro_mip <- data.frame(cbind(mro.v, mip.v))

bdr <- c(0, 2.5, 5, 50, 95, 97.5, 100)
quantile(mro.v, bdr*0.01, na.rm=T)
quantile(mip.v, bdr*0.01, na.rm=T)

df <- data.frame(MRO=mro_mip$mro.v)
g <- ggplot(df, aes(x=MRO)) + geom_histogram()
plot(g)

ggsave(g, filename=sprintf("%s/Histogram_MRO.pdf", dir$figdir), h=0.8, w=4)


df <- data.frame(MIP=mro_mip$mip.v)
g <- ggplot(df, aes(x=MIP)) + geom_histogram(binwidth=1)
plot(g)

ggsave(g, filename=sprintf("%s/Histogram_MIP.pdf", dir$figdir), h=0.8, w=4)






