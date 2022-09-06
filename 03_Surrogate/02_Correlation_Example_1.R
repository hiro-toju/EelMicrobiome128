
#setwd("/Users/toju/Dropbox/YAJIMA_Daii/Statistics_R")
setwd("/Users/toju/Dropbox/YAJIMA_Daii/Statistics_R")


set.seed(123)

n.ex.rep <- 5 # number of tanks
sample.th <- 30
n.boot <- 2 # (100000)

source('functions/functions.R')
source('functions/rEDM075_surrogate_codes.R')

# -- Create directory to save
dir <- make.dir('03_Surrogate/02_output')

# SugiharaLabバージョンのrEDMではtwin-surrogateが使えないので
# Hao Ye さんが公開されている過去のバージョンをインストールします.
# remotes::install_github("ha0ye/rEDM")

library(RcppThread)
library(RcppEigen)
library(doParallel)
library(foreach)
library(ggplot2)
library(ggsci)
library(imputeTS)
library(rEDM)
library(reshape2)

############################################################################
## Threshold for all tanks

info <- read.table('Table/Sample.List_1.txt', header=T)
data  <- read.table('Table/seqtab.16SrRNAcount.pro.txt', header=T)
taxa <- readRDS('Table/Taxa_list.rds')

param <- as.matrix(read.table("Table/Parameters_3.txt", header=T, na.strings = "NA"))

param2 <- melt(param)
colnames(param2) <- c("Day", "Parameter", "value")

pH <- param2[grep("pH", param2$Parameter), ]
DO <- param2[grep("DO", param2$Parameter), ]
AS <- param2[grep("AS", param2$Parameter), ]
colnames(pH)[3] <- 'pH.value'
colnames(DO)[3] <- 'DO.value'
colnames(AS)[3] <- 'AS.value'

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
#write.table(comm.out, file=sprintf("./Timeseries_graph/Surrogate/Table/Matrix_th%s.txt", sample.th), sep="\t", quote=F, row.names=F)


############################################################################
## Threshold for each tank

day <- as.matrix(1:128, ncol=1)
colnames(day) <- c('Day')
dinfo <- list()
dmat <- list()
dcomm <- list()
asv.info <- list()
r.slr <- list()
r.mb <- list()
dcomm2 <- list()

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
	Day <- as.numeric(gsub(sprintf('S_%s', i), '', Samples))	
	comm2 <- cbind(Samples, Day, dcomm[[i]])
	dcomm2[[i]] <- merge(day, comm2, by="Day", all.x=TRUE)
	
	write.table(dcomm2[[i]], file=sprintf("%s/CommData_th%s_Tank%s.txt", dir$tabledir, sample.th, tank), sep="\t", quote=FALSE, row.names=FALSE)
	
	write.table(asv.info[[i]], file=sprintf("%s/ASV.Info_th%s_Tank%s.txt", dir$tabledir, sample.th, tank), sep="\t", quote=FALSE, row.names=FALSE)
}

############################################################################
## Twin-surrogate functions

mcor <- function(tss) {
    cor <- cor(tss, method = 'spearman')
    ncor <- length(cor[1,])
    (mean(cor) - 1/ncor) * (ncor/(ncor-1))
}

shuffle <- function(tss) apply(tss, 2, function(x) {make_surrogate_twin(x, phase_lock = FALSE)})

synctest <- function (tss, nboot=1999) {
  mc <- mcor(tss)
  boosts <- foreach(i=1:nboot, .combine='c') %dopar% mcor(shuffle(tss))
  mcp <- (1/(nboot+1))*(1+sum(boosts >= mc))
  result <- c(mc, mcp)
  name <- c("Correlation", "p-value")
  names(result) <- name
  result
}


parcor.test <- function (x, y, z, nboot=1999) {
  r_xy <- cor(x, y, method='spearman')
  r_xz <- cor(x, z, method='spearman')
  r_yz <- cor(y, z, method='spearman')
  r_xy_z <- (r_xy - r_xz*r_yz) / (((1-r_xz^2)^0.5)*((1-r_yz^2)^0.5))
  
  boot_r_xy <- foreach(i=1:nboot, .combine='c') %dopar% cor(shuffle(cbind(x, y)), method='spearman')[1,2]
  boot_r_xz <- foreach(i=1:nboot, .combine='c') %dopar% cor(shuffle(cbind(x, z)), method='spearman')[1,2]
  boot_r_yz <- foreach(i=1:nboot, .combine='c') %dopar% cor(shuffle(cbind(y, z)), method='spearman')[1,2]
  boot_r_xy_z <- foreach(i=1:nboot, .combine='c') %dopar% (boot_r_xy[i] - boot_r_xz[i]*boot_r_yz[i]) / (((1-boot_r_xz[i]^2)^0.5)*((1-boot_r_yz[i]^2)^0.5))
  
  Pvalue <- (1/(nboot+1))*(1+sum(boot_r_xy_z >= r_xy_z))
  result <- c(r_xy_z, Pvalue)
  name <- c("Partial.Correlation", "p-value")
  names(result) <- name
  result
}

############################################################################
## Surrogate permutation

dcomm3 <- list()
for (j in 1:n.ex.rep) {
	dcomm3[[j]] <- dcomm2[[j]][, -c(1,2)]
	}
	


## pH 

pH.int <- pH
pH.int[,3] <- na_interpolation(pH.int[,3], option = "linear")

j <- 1

v <- subset(pH.int, pH.int$Parameter==sprintf('T%s_pH', j))
data.p <- foreach(i=1:ncol(dcomm3[[j]])) %dopar% cbind(v, na_interpolation(as.numeric(dcomm3[[j]][, i]), option = "linear"))

plot(data.p[[4]][,3], data.p[[4]][,4])

ex <- data.p[[4]]
colnames(ex)[4] <- "Abundance"

g1 <- ggplot(ex, aes(x= pH.value, y=Abundance, color="deepskyblue4")) + geom_point() + labs(title='X_0004', x='pH', y='Absolute abundance (16S rRNA copies/uL)') +theme(legend.position = "none")

ggsave(g1, filename=sprintf("%s/Correlations_Example.pdf", dir$figdir), h=3, w=3)
