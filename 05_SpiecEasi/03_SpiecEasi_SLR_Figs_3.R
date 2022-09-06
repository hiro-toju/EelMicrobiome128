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
library(imputeTS)

load(sprintf("05_SpiecEasi/01_output/RData/SpiecEasi_th%s_latent%s.RData", sample.th, nlatent))

load("03_Surrogate/01_output/RData/Surrogate_ASVs_th30_1e+05.RData")


# -- Create directory to save
dir <- make.dir('05_SpiecEasi/03_output')

cluster = makeCluster(parallel::detectCores(logical = FALSE), type = "FORK")
registerDoParallel(cluster)

######################################################
## Highlighting ASVs

size.vec <- foreach(i=1:n.ex.rep) %dopar% colnames(dcomm[[i]])

for (i in 1:n.ex.rep) {
	size.vec[[i]][which(size.vec[[i]]=="X_0002")] <- 7 # Best Cetobacterium
	size.vec[[i]][which(size.vec[[i]]=="X_0014")] <- 5 # Paraclostridium
	size.vec[[i]][which(size.vec[[i]]=="X_0020")] <- 5 # Plesiomonas
	size.vec[[i]][which(size.vec[[i]]=="X_0027")] <- 5 # Edwardsiella
	size.vec[[i]][which(size.vec[[i]]=="X_0028")] <- 5 # Romboutsia
	size.vec[[i]][which(size.vec[[i]]=="X_0029")] <- 5 # Clostridium sensu stricto 1
	size.vec[[i]][which(size.vec[[i]]=="X_0041")] <- 5 # Turicibacter	
	size.vec[[i]][which(size.vec[[i]]=="X_0064")] <- 5 # Barnesiellaceae	
	size.vec[[i]][which(size.vec[[i]]=="X_0134")] <- 5 # Worst cvE6	
	size.vec[[i]][grep('X_', size.vec[[i]])] <- 3
	mode(size.vec[[i]]) <- 'numeric'
}

shape.vec <- foreach(i=1:n.ex.rep) %dopar% colnames(dcomm[[i]])

for (i in 1:n.ex.rep) {
	shape.vec[[i]][which(shape.vec[[i]]=="X_0002")] <- 18 # Best Cetobacterium
	shape.vec[[i]][which(shape.vec[[i]]=="X_0014")] <- 15 # Paraclostridium
	shape.vec[[i]][which(shape.vec[[i]]=="X_0020")] <- 17 # Plesiomonas
	shape.vec[[i]][which(shape.vec[[i]]=="X_0027")] <- 7 # Edwardsiella
	shape.vec[[i]][which(shape.vec[[i]]=="X_0028")] <- 9 # Romboutsia
	shape.vec[[i]][which(shape.vec[[i]]=="X_0029")] <- 10 # Clostridium sensu stricto 1
	shape.vec[[i]][which(shape.vec[[i]]=="X_0041")] <- 11 # Turicibacter	
	shape.vec[[i]][which(shape.vec[[i]]=="X_0064")] <- 12 # Barnesiellaceae	
	shape.vec[[i]][which(shape.vec[[i]]=="X_0134")] <- 13 # Worst cvE6	
	shape.vec[[i]][grep('X_', shape.vec[[i]])] <- 16
	mode(shape.vec[[i]]) <- 'numeric'
}

######################################################
## SLR


## partial AS controlling pH ##

pr.AS_pH <- list()
g <- list()

for (i in 1:n.ex.rep) {
colnames(AS.result_Contr.pH[[i]])[1:2] <- c("ASV.ID", "Correlation")
mode(AS.result_Contr.pH[[i]]$ASV.ID) <- "character"
d <- merge(AS.result_Contr.pH[[i]], taxa, by.x="ASV.ID", by.y="ID")
FDR <- p.adjust(apply(matrix(as.numeric(unlist(d[,3:4])), ncol=2), 1, FUN=min), method="fdr")
Cor <- as.numeric(unlist(d$Correlation))
x <- na.omit(data.frame(cbind(Cor, FDR)))
pr.AS_pH[[i]] <- cbind(d, FDR)
}


g <- list()

for(i in 1:n.ex.rep){ 
Correlation <- as.numeric(pr.AS_pH[[i]]$Correlation)
g[[i]] <- ggraph(ig.slr.posi[[i]],layout = "stress")+ # OR layout = "nicely"
  geom_edge_link0(width=0.2,colour="darkslategray4")+
  geom_node_point(aes(col= Correlation), shape=shape.vec[[i]], size=size.vec[[i]], show.legend=FALSE)+
  scale_colour_gradient2(low = "blue", mid = "grey90", high = "red", na.value = "white", guide = "colourbar")+
  guides()+
  theme_graph()+
  theme(legend.position = "right")+
  theme(text = element_text("Helvetica"))

ggsave(g[[i]], filename=sprintf("%s/SpiecEasi_pAScontrPH_th%s_Tank%s.pdf", dir$figdir, sample.th, i), h=5.5, w=16)
}

g <- list()

for(i in 1:n.ex.rep){ 
Correlation <- as.numeric(pr.AS_pH[[i]]$Correlation)
g[[i]] <- ggraph(ig.slr.posi[[i]],layout = "stress")+ # OR layout = "nicely"
  geom_edge_link0(width=0.2,colour="darkslategray4")+
  geom_node_point(aes(col= Correlation), shape=shape.vec[[i]], size=size.vec[[i]], show.legend=TRUE)+
  scale_colour_gradient2(low = "blue", mid = "grey90", high = "red", na.value = "white", guide = "colourbar")+
  guides()+
  theme_graph()+
  theme(legend.position = "right")+
  theme(text = element_text("Helvetica"))

ggsave(g[[i]], filename=sprintf("%s/SLR_Partial.Corr_AS_cotlrPH_th%s_Tank%s.pdf", dir$figdir, sample.th, i), h=6, w=24)
}


## partial AS controlling DO ##

pr.AS_DO <- list()
g <- list()

for (i in 1:n.ex.rep) {
colnames(AS.result_Contr.DO[[i]])[1:2] <- c("ASV.ID", "Correlation")
mode(AS.result_Contr.DO[[i]]$ASV.ID) <- "character"
d <- merge(AS.result_Contr.DO[[i]], taxa, by.x="ASV.ID", by.y="ID")
FDR <- p.adjust(apply(matrix(as.numeric(unlist(d[,3:4])), ncol=2), 1, FUN=min), method="fdr")
Cor <- as.numeric(unlist(d$Correlation))
x <- na.omit(data.frame(cbind(Cor, FDR)))
pr.AS_DO[[i]] <- cbind(d, FDR)
}


g <- list()

for(i in 1:n.ex.rep){ 
Correlation <- as.numeric(pr.AS_DO[[i]]$Correlation)
g[[i]] <- ggraph(ig.slr.posi[[i]],layout = "stress")+ # OR layout = "nicely"
  geom_edge_link0(width=0.2,colour="darkslategray4")+
  geom_node_point(aes(col= Correlation), shape=shape.vec[[i]], size=size.vec[[i]], show.legend=FALSE)+
  scale_colour_gradient2(low = "blue", mid = "grey90", high = "red", na.value = "white", guide = "colourbar")+
  guides()+
  theme_graph()+
  theme(legend.position = "right")+
  theme(text = element_text("Helvetica"))

ggsave(g[[i]], filename=sprintf("%s/SpiecEasi_pAScontrDO_th%s_Tank%s.pdf", dir$figdir, sample.th, i), h=5.5, w=16)
}

g <- list()

for(i in 1:n.ex.rep){ 
Correlation <- as.numeric(pr.AS_DO[[i]]$Correlation)
g[[i]] <- ggraph(ig.slr.posi[[i]],layout = "stress")+ # OR layout = "nicely"
  geom_edge_link0(width=0.2,colour="darkslategray4")+
  geom_node_point(aes(col= Correlation), shape=shape.vec[[i]], size=size.vec[[i]], show.legend=TRUE)+
  scale_colour_gradient2(low = "blue", mid = "grey90", high = "red", na.value = "white", guide = "colourbar")+
  guides()+
  theme_graph()+
  theme(legend.position = "right")+
  theme(text = element_text("Helvetica"))

ggsave(g[[i]], filename=sprintf("%s/SLR_Partial.Corr_AS_cotlrDO_th%s_Tank%s_caption.pdf", dir$figdir, sample.th, i), h=5.5, w=16)
}

stopCluster(cluster)
