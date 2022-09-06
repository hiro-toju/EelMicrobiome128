############################################################################
####
#### R script for Fujita (2019)
####
#### SpiecEasi
#### 2022.02.16 Kenta Suzuki
#### 2022.04.14 Hirokazu Toju 
#### R 4.1.2
#### 
############################################################################

# Repeat the process for lag = 1 to 5

setwd("/Users/toju/Dropbox/YAJIMA_Daii/Statistics_R")

set.seed(123)

n.ex.rep <- 5 # number of tanks
sample.th <- 30
n.boot <- 2
lag <- 1

source('functions/functions.R')
source('functions/rEDM075_surrogate_codes.R')

# -- Create directory to save
dir <- make.dir('04_Surrogate_Timelag/01_output')

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
taxa <- readRDS('Table/Taxa_list.pro.rds')

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

#	write.table(dcomm2[[i]], file=sprintf("%s/CommData_th%s_Tank%s.txt", dir$tabledir, sample.th, tank), sep="\t", quote=FALSE, row.names=FALSE)
	
#	write.table(asv.info[[i]], file=sprintf("%s/ASV.Info_th%s_Tank%s.txt", dir$tabledir, sample.th, tank), sep="\t", quote=FALSE, row.names=FALSE)
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
## Premaring timelag data

dcomm3 <- list()

for (j in 1:n.ex.rep) {
	dcomm3[[j]] <- dcomm2[[j]][, -c(1,2)]
	}

pH.int <- pH
pH.int[,3] <- na_interpolation(pH.int[,3], option = "linear")
colnames(pH.int)[3] <- 'pH.value'

DO.int <- DO
DO.int[,3] <- na_interpolation(DO.int[,3], option = "linear")
colnames(DO.int)[3] <- 'DO.value'

AS.int <- AS
AS.int[,3] <- na_interpolation(AS.int[,3], option = "linear")
colnames(AS.int)[3] <- 'AS.value'

	
############################################################################
## Surrogate permutation
	
cluster = makeCluster(parallel::detectCores(logical = FALSE), type = "FORK")
registerDoParallel(cluster)

############################################################################
## Environments + i day: dealy
## Community    + 0 day

## pH 

pH.result <- list()

start <- proc.time()[3]
for (j in 1:n.ex.rep) {
v <- subset(pH.int, pH.int$Parameter==sprintf('T%s_pH', j))[-c(1:lag),]

data.p <- foreach(i=1:ncol(dcomm3[[j]])) %dopar% cbind(v, na_interpolation(as.numeric(dcomm3[[j]][1:(nrow(dcomm3[[j]])-lag), i]), option = "linear"))

pH.r <- c()
pH.r <- foreach(i=1:ncol(dcomm3[[j]]), .combine='rbind', .errorhandling="pass") %dopar% synctest(as.matrix(data.p[[i]][, c(3,4)]), nboot=n.boot)

corr <- foreach(i=1:ncol(dcomm3[[j]]), .combine='rbind', .errorhandling="pass") %dopar% synctest(as.matrix(data.p[[i]][, c(3,4)]), nboot=1)

pH.result[[j]] <- data.frame(cbind(colnames(dcomm3[[j]]), corr[,1], pH.r[,2]))
colnames(pH.result[[j]]) <- c("Variable", "Correlation", "P")
pH.result[[j]]$P[which(pH.result[[j]]$P=="mcor(shuffle(tss))")] <- NA
P.opposite <-  1 - as.numeric(pH.result[[j]]$P)
pH.result[[j]] <- cbind(pH.result[[j]], P.opposite)
colnames(pH.result[[j]]) <- c("Variable", "r_pH", "P_pH", "P.opp_pH")
gc(reset = TRUE)
gc(reset = TRUE)
}
end <- proc.time()[3]
print(end-start)


## DO 

DO.result <- list()

start <- proc.time()[3]
for (j in 1:n.ex.rep) {
v <- subset(DO.int, DO.int$Parameter==sprintf('T%s_DO', j))[-c(1:lag),]

data.p <- foreach(i=1:ncol(dcomm3[[j]])) %dopar% cbind(v, na_interpolation(as.numeric(dcomm3[[j]][1:(nrow(dcomm3[[j]])-lag), i]), option = "linear"))

DO.r <- c()
DO.r <- foreach(i=1:ncol(dcomm3[[j]]), .combine='rbind', .errorhandling="pass") %dopar% synctest(as.matrix(data.p[[i]][, c(3,4)]), nboot=n.boot)

corr <- foreach(i=1:ncol(dcomm3[[j]]), .combine='rbind', .errorhandling="pass") %dopar% synctest(as.matrix(data.p[[i]][, c(3,4)]), nboot=1)

DO.result[[j]] <- data.frame(cbind(colnames(dcomm3[[j]]), corr[,1], DO.r[,2]))
colnames(DO.result[[j]]) <- c("Variable", "Correlation", "P")
DO.result[[j]]$P[which(DO.result[[j]]$P=="mcor(shuffle(tss))")] <- NA
P.opposite <-  1 - as.numeric(DO.result[[j]]$P)
DO.result[[j]] <- cbind(DO.result[[j]], P.opposite)
colnames(DO.result[[j]]) <- c("Variable", "r_DO", "P_DO", "P.opp_DO")
gc(reset = TRUE)
gc(reset = TRUE)
}
end <- proc.time()[3]
print(end-start)


## AS 

AS.result <- list()

start <- proc.time()[3]
for (j in 1:n.ex.rep) {
v <- subset(AS.int, AS.int$Parameter==sprintf('T%s_AS', j))[-c(1:lag),]

data.p <- foreach(i=1:ncol(dcomm3[[j]])) %dopar% cbind(v, na_interpolation(as.numeric(dcomm3[[j]][1:(nrow(dcomm3[[j]])-lag), i]), option = "linear"))

AS.r <- c()
AS.r <- foreach(i=1:ncol(dcomm3[[j]]), .combine='rbind', .errorhandling="pass") %dopar% synctest(as.matrix(data.p[[i]][, c(3,4)]), nboot=n.boot)

corr <- foreach(i=1:ncol(dcomm3[[j]]), .combine='rbind', .errorhandling="pass") %dopar% synctest(as.matrix(data.p[[i]][, c(3,4)]), nboot=1)

AS.result[[j]] <- data.frame(cbind(colnames(dcomm3[[j]]), corr[,1], AS.r[,2]))
colnames(AS.result[[j]]) <- c("Variable", "Correlation", "P")
AS.result[[j]]$P[which(AS.result[[j]]$P=="mcor(shuffle(tss))")] <- NA
P.opposite <-  1 - as.numeric(AS.result[[j]]$P)
AS.result[[j]] <- cbind(AS.result[[j]], P.opposite)
colnames(AS.result[[j]]) <- c("Variable", "r_AS", "P_AS", "P.opp_AS")
gc(reset = TRUE)
gc(reset = TRUE)
}
end <- proc.time()[3]
print(end-start)


############################################################################
## Partial AS : controlling pH

AS.result_Contr.pH <- list()

start <- proc.time()[3]
for (j in 1:n.ex.rep) {
v.AS <- subset(AS.int, AS.int$Parameter==sprintf('T%s_AS', j))
v.pH <- subset(pH.int, pH.int$Parameter==sprintf('T%s_pH', j))
AS.pH <- cbind(v.AS, v.pH$pH.value)[-c(1:lag),]

data.p <- foreach(i=1:ncol(dcomm3[[j]])) %dopar% cbind(AS.pH, na_interpolation(as.numeric(dcomm3[[j]][1:(nrow(dcomm3[[j]])-lag), i]), option = "linear"))

AS.r <- c()

AS.r <- foreach(i=1:ncol(dcomm3[[j]]), .combine='rbind', .errorhandling="pass") %dopar% parcor.test(data.p[[i]][, 5], data.p[[i]][, 3], data.p[[i]][, 4], nboot=n.boot)

corr <- foreach(i=1:ncol(dcomm3[[j]]), .combine='rbind', .errorhandling="pass") %dopar% parcor.test(as.matrix(data.p[[i]][, 5]), as.matrix(data.p[[i]][, 3]), as.matrix(data.p[[i]][, 4]), nboot=1)

AS.result_Contr.pH[[j]] <- data.frame(cbind(colnames(dcomm3[[j]]), corr[,1], AS.r[,2]))
colnames(AS.result_Contr.pH[[j]]) <- c("Variable", "Partial.Correlation", "P")

AS.result_Contr.pH[[j]]$P[grep("spearman", AS.result_Contr.pH[[j]]$P)] <- NA

P.opposite <-  1 - as.numeric(AS.result_Contr.pH[[j]]$P)
AS.result_Contr.pH[[j]] <- cbind(AS.result_Contr.pH[[j]], P.opposite)
colnames(AS.result_Contr.pH[[j]]) <- c("Variable", "Partial.r_AS.contr.pH", "P_AS.contr.pH", "P.opp_AS.contr.pH")
gc(reset = TRUE)
gc(reset = TRUE)
}
end <- proc.time()[3]
print(end-start)


############################################################################
## Partial AS : controlling DO

AS.result_Contr.DO <- list()

start <- proc.time()[3]
for (j in 1:n.ex.rep) {
v.AS <- subset(AS.int, AS.int$Parameter==sprintf('T%s_AS', j))
v.DO <- subset(DO.int, DO.int$Parameter==sprintf('T%s_DO', j))
AS.DO <- cbind(v.AS, v.DO$DO.value)[-c(1:lag),]

data.p <- foreach(i=1:ncol(dcomm3[[j]])) %dopar% cbind(AS.DO, na_interpolation(as.numeric(dcomm3[[j]][1:(nrow(dcomm3[[j]])-lag), i]), option = "linear"))

AS.r <- c()

AS.r <- foreach(i=1:ncol(dcomm3[[j]]), .combine='rbind', .errorhandling="pass") %dopar% parcor.test(data.p[[i]][, 5], data.p[[i]][, 3], data.p[[i]][, 4], nboot=n.boot)

corr <- foreach(i=1:ncol(dcomm3[[j]]), .combine='rbind', .errorhandling="pass") %dopar% parcor.test(as.matrix(data.p[[i]][, 5]), as.matrix(data.p[[i]][, 3]), as.matrix(data.p[[i]][, 4]), nboot=1)

AS.result_Contr.DO[[j]] <- data.frame(cbind(colnames(dcomm3[[j]]), corr[,1], AS.r[,2]))
colnames(AS.result_Contr.DO[[j]]) <- c("Variable", "Partial.Correlation", "P")

AS.result_Contr.DO[[j]]$P[grep("spearman", AS.result_Contr.DO[[j]]$P)] <- NA

P.opposite <-  1 - as.numeric(AS.result_Contr.DO[[j]]$P)
AS.result_Contr.DO[[j]] <- cbind(AS.result_Contr.DO[[j]], P.opposite)
colnames(AS.result_Contr.DO[[j]]) <- c("Variable", "Partial.r_AS.contr.DO", "P_AS.contr.DO", "P.opp_AS.contr.DO")
gc(reset = TRUE)
gc(reset = TRUE)
}
end <- proc.time()[3]
print(end-start)


############################################################################
## Summarizing the resutls


result.all <- list()

for(i in 1:n.ex.rep) {
	result.all[[i]] <- cbind(AS.result_Contr.pH[[i]][, 1:2], AS.result_Contr.DO[[i]][,2],  AS.result[[i]][,2], pH.result[[i]][,2], DO.result[[i]][,2])
	
	colnames(result.all[[i]]) <- c('ASV.ID', 'AS.result_Contr.pH', 'AS.result_Contr.DO',  'AS.result', 'pH.result', 'DO.result')
		
	write.table(result.all[[i]], file=sprintf("%s/Surrogate_ASVs_EnvDealy%s_All_th%s_Tank%s_%s.txt", dir$tabledir, lag, sample.th, i, n.boot), sep="\t", quote=FALSE, row.names=FALSE)
}


stopCluster(cluster)

saveRDS(result.all, sprintf("%s/Surrogate_ASVs_EnvDealy%s_All_th%s_Tank%s_%s.rds", dir$rdsdir, lag, sample.th, i, n.boot))


############################################################################

save.image(sprintf("%s/Surrogate_ASVs_EnvDealy%s_th%s_%s.RData", dir$rdatadir, lag, sample.th, n.boot))


