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
source('functions/functions.R')


set.seed(123)

library(vegan)
library(ggplot2)

# -- Create directory to save
dir <- make.dir('06_Picrust/03_output')

d <- read.table("Table/pred_metagenome_unstrat_descrip.tsv", header=T, row.name=1)
ko.name <- d$description

mat <- t(d[, -1])
mat2 <- subset(mat, rowSums(mat) > 0)

jac <- vegdist(mat2, method="jaccard", binary=TRUE)
pcoa <- cmdscale(jac, eig=TRUE)
ordiplot(pcoa, display='sites')

pcoa.cord <- data.frame(pcoa$points)
colnames(pcoa.cord) <- c('PCoA_1', 'PCoA_2')
X_0002 <- subset(pcoa.cord, rownames(pcoa.cord)=='X_0002')
X_0014 <- subset(pcoa.cord, rownames(pcoa.cord)=='X_0014')
X_0020 <- subset(pcoa.cord, rownames(pcoa.cord)=='X_0020')
X_0027 <- subset(pcoa.cord, rownames(pcoa.cord)=='X_0027')
X_0028 <- subset(pcoa.cord, rownames(pcoa.cord)=='X_0028')
X_0029 <- subset(pcoa.cord, rownames(pcoa.cord)=='X_0029')
X_0041 <- subset(pcoa.cord, rownames(pcoa.cord)=='X_0041')
X_0064 <- subset(pcoa.cord, rownames(pcoa.cord)=='X_0064')
X_0134 <- subset(pcoa.cord, rownames(pcoa.cord)=='X_0134')

g <- ggplot(NULL)
g <- g + geom_point(data=pcoa.cord, aes(x=PCoA_1, y=PCoA_2), na.rm=TRUE, col="grey40", cex=1) + theme(text = element_text(size = 11)) +labs(x= "PCoA 1", y="PCoA 2")
g <- g + geom_point(data= X_0002, aes(x=PCoA_1, y=PCoA_2), na.rm=TRUE, col="red", pch=18, cex=5)
g <- g + geom_point(data= X_0014, aes(x=PCoA_1, y=PCoA_2), na.rm=TRUE, col="magenta3", pch= 15, cex=4)
g <- g + geom_point(data= X_0020, aes(x=PCoA_1, y=PCoA_2), na.rm=TRUE, col="magenta3", pch=17, cex=5)
g <- g + geom_point(data= X_0027, aes(x=PCoA_1, y=PCoA_2), na.rm=TRUE, col="magenta3", pch=7, cex=5)
g <- g + geom_point(data= X_0028, aes(x=PCoA_1, y=PCoA_2), na.rm=TRUE, col="magenta3", pch=9, cex=5)
g <- g + geom_point(data= X_0029, aes(x=PCoA_1, y=PCoA_2), na.rm=TRUE, col="magenta3", pch=10, cex=7)
g <- g + geom_point(data= X_0041, aes(x=PCoA_1, y=PCoA_2), na.rm=TRUE, col="magenta3", pch=11, cex=5)
g <- g + geom_point(data= X_0064, aes(x=PCoA_1, y=PCoA_2), na.rm=TRUE, col="darkgoldenrod2", pch=12, cex=5)
g <- g + geom_point(data= X_0134, aes(x=PCoA_1, y=PCoA_2), na.rm=TRUE, col="blue", pch=13, cex=6)
plot(g)

ggsave(g, filename=sprintf("%s/Picrust_KO_PCA.pdf", dir$figdir), h=4, w=4)






