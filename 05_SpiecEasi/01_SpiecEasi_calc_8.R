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
source('functions/rEDM075_surrogate_codes.R')

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

cluster = makeCluster(parallel::detectCores(logical = FALSE), type = "FORK")
registerDoParallel(cluster)

# -- Create directory to save
dir <- make.dir('05_SpiecEasi/01_output')

# -- Load data table

info <- read.table('Table/Sample.List_1.txt', header=T)
data  <- read.table('Table/seqtab.16SrRNAcount.pro.txt', header=T)
taxa <- readRDS('Table/Taxa_list.pro.rds')

merged.mat <- merge(info, data, by='Sample.ID', all = FALSE, sort = TRUE)
rownames(merged.mat) <- merged.mat$Sample.ID
sml  <- merged.mat[, 1:3]
mat <- as.matrix(merged.mat[, 4:ncol(merged.mat)])
dim(mat)

mat.binary <- mat
mat.binary[which(mat.binary > 0)] <- 1
comm <- t(subset(t(mat), rowSums(t(mat.binary)) >= sample.th))
dim(comm)

Samples.selected <- rownames(comm)

comm.out <- cbind(Samples.selected, comm)

write.table(comm.out, file=sprintf("%s/SpiecEasi.Matrix_th%s.txt", dir$tabledir, sample.th), sep="\t", quote=F, row.names=F)


############################################################################
# Threshold for each sample

dinfo <- list()
dmat <- list()
dcomm <- list()
asv.info <- list()
r.slr <- list()
r.mb <- list()

for(i in 1:n.ex.rep){ 	
    dinfo[[i]] <- sml[sml$Tank==i, ]
    dmat[[i]] <- comm[sml$Tank==i, ]
    
    dmat.binary <- dmat[[i]]
	dmat.binary[which(dmat.binary > 0)] <- 1
	dcomm[[i]] <- t(subset(t(dmat[[i]]), rowSums(t(dmat.binary)) >= sample.th))
	dim(dcomm[[i]])
	asv.list <- as.matrix(colnames(dcomm[[i]]), ncol=1)
	colnames(asv.list) <- c("ASV.ID")
	asv.info[[i]] <- merge(asv.list, taxa, by.x='ASV.ID', by.y="ID", all = FALSE, sort = TRUE)
}


for (i in 1:n.ex.rep) {
	tank <- i
	Samples <- rownames(dcomm[[i]])
	out <- cbind(Samples, dcomm[[i]])
	
	write.table(out, file=sprintf("%s/CommData_th%s_Tank%s.txt", dir$tabledir, sample.th, tank), sep="\t", quote=FALSE, row.names=FALSE)
	
	write.table(asv.info[[i]], file=sprintf("%s/ASV.Info_th%s_Tank%s.txt", dir$tabledir, sample.th, tank), sep="\t", quote=FALSE, row.names=FALSE)
}

############################################################################
# Spiec-Easi

start <- proc.time()[3]

r.mb <- foreach(i=1:n.ex.rep) %dopar% spiec.easi(dcomm[[i]], method='mb',lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=nrep.star, ncores=6))

end <- proc.time()[3]
print(end-start)


start <- proc.time()[3]

r.slr <- foreach(i=1:n.ex.rep) %dopar% spiec.easi(dcomm[[i]], method='slr', r=nlatent, lambda.min.ratio=1e-2,　nlambda=20, pulsar.params=list(rep.num=nrep.star, ncores=6))

end <- proc.time()[3]
print(end-start)



############################################################################
# SAVING Spiec-Easi results

saveRDS(r.mb, sprintf("%s/SpiecEasi.MB_th%s.rds", dir$rdsdir, sample.th))

saveRDS(r.slr, sprintf("%s/SpiecEasi.SLR_th%s_latent%s.rds", dir$rdsdir, sample.th, nlatent))

############################################################################

mb.beta <- list()
mb.posi <- list()
mb.nega <- list()
ig.mb.posi <- list()
ig.mb.nega <- list()
slr.cor <- list()
slr.posi <- list()
slr.nega <- list()
ig.slr.posi <- list()
ig.slr.nega <- list()
vsize <- list()
coord <- list()

for(i in 1:n.ex.rep){ 
vname <- colnames(dcomm[[i]])

mb.beta[[i]] <- as.matrix(getOptBeta(r.mb[[i]]))
mb.refit <- as.matrix(getRefit(r.mb[[i]]))
mb.beta[[i]][which(mb.refit==FALSE)] <- 0
diag(mb.beta[[i]]) <- 0
mb.posi[[i]] <- mb.beta[[i]]
mb.posi[[i]][which(mb.posi[[i]] <= 0)] <- 0
mb.nega[[i]] <- mb.beta[[i]]
mb.nega[[i]][which(mb.nega[[i]] >= 0)] <- 0
ig.mb.posi[[i]] <- graph.adjacency(mb.posi[[i]], mode='undirected', diag=FALSE, weighted=TRUE)
ig.mb.nega[[i]] <- graph.adjacency(-mb.nega[[i]], mode='undirected', diag=FALSE, weighted=TRUE)

slr.icov <- solve(as.matrix(r.slr[[i]]$est$icov[[getOptInd(r.slr[[i]])]]))
slr.cor[[i]] <- cov2cor(slr.icov)
slr.refit <- as.matrix(getRefit(r.slr[[i]]))
slr.cor[[i]][which(slr.refit==FALSE)] <- 0
diag(slr.cor[[i]]) <- 0
slr.posi[[i]] <- slr.cor[[i]]
slr.posi[[i]][which(slr.posi[[i]] <= 0)] <- 0
slr.nega[[i]] <- slr.cor[[i]]
slr.nega[[i]][which(slr.nega[[i]] >= 0)] <- 0
ig.slr.posi[[i]] <- graph.adjacency(slr.posi[[i]], mode='undirected', diag=FALSE, weighted=TRUE)
ig.slr.nega[[i]] <- graph.adjacency(-slr.nega[[i]], mode='undirected', diag=FALSE, weighted=TRUE)
}

######################################################

# MB vs. SLR plot

for(i in 1:n.ex.rep){ 
pdf(sprintf("%s/Plot_MBvsSLR_th%s_Tank%s_latent%s.pdf", dir$figdir, sample.th, i, nlatent))
plot(as.vector(mb.beta[[i]]), as.vector(slr.cor[[i]]))
abline(a=0, b=1, lty=2)
dev.off()
}

######################################################
# modules

mbet.mb <- list()
mbet.slr <- list()

for(i in 1:n.ex.rep){ 
mbet.mb[[i]] <- cluster_edge_betweenness(ig.mb.posi[[i]], weights = E(ig.mb.posi[[i]])$weight)
mbet.slr[[i]] <- cluster_edge_betweenness(ig.slr.posi[[i]], weights = E(ig.slr.posi[[i]])$weight)
}

######################################################

save.image(sprintf("%s/SpiecEasi_th%s_latent%s.RData", dir$rdatadir, sample.th, nlatent))

