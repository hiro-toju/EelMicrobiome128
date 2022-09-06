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

source('functions/functions.R')


sample.th <- 30

# -- Create directory to save
dir <- make.dir('06_Picrust/01_output')

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

data <- t(comm)

out <- matrix(0, nrow=nrow(data), ncol=nrow(data))
diag(out) <- 1
name <- rownames(data)
colnames(out) <- name
out <- cbind(name, out)

write.table(out, file=sprintf("%s/square.th%s.txt", dir$tabledir, sample.th), quote=F, sep="\t", row.names=F)







