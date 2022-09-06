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

#set.seed(123)
set.seed(12)

n.ex.rep <- 5 # number of tanks
sample.th <- 30
nrep.star <- 100
nlatent <- 10
nperm.mantel <- 1000

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

dir <- make.dir('05_SpiecEasi/01_output')

cluster = makeCluster(parallel::detectCores(logical = FALSE), type = "FORK")
registerDoParallel(cluster)

load(sprintf("%s/SpiecEasi_th%s_latent%s.RData", dir$rdatadir, sample.th, nlatent))

######################################################

# MB vs. SLR plot

for(i in 1:n.ex.rep){ 
pdf(sprintf("%s/Plot_MBvsSLR_th%s_Tank%s_latent%s.pdf", dir$figdir, sample.th, i, nlatent))
plot(as.vector(mb.beta[[i]]), as.vector(slr.cor[[i]]))
abline(a=0, b=1, lty=2)
dev.off()
}

sink(sprintf('%s/MB_vs_SLR_Mantel.txt', dir$tabledir))

mantel(mb.beta[[1]], slr.cor[[1]], permutations=nperm.mantel)
mantel(mb.beta[[2]], slr.cor[[2]], permutations=nperm.mantel)
mantel(mb.beta[[3]], slr.cor[[3]], permutations=nperm.mantel)
mantel(mb.beta[[4]], slr.cor[[4]], permutations=nperm.mantel)
mantel(mb.beta[[5]], slr.cor[[5]], permutations=nperm.mantel)

sink()
	
	
g <- list()

for (i in 1:n.ex.rep) {

Coefficient_MB <- as.vector(mb.beta[[i]])
Coefficient_SLR <- as.vector(slr.cor[[i]])
d <- data.frame(cbind(Coefficient_MB, Coefficient_SLR))

g[[i]] <- ggplot(d, aes(x= Coefficient_MB, y= Coefficient_SLR), color="lightpink3") + geom_point() + labs(title=sprintf('Tank %s', i), x='Coefficients of microbe-to-microbe interactions (MB)', y='Coefficients of microbe-to-microbe interactions (SLR)')

}

g2 <- gridExtra::grid.arrange(g[[1]], g[[2]], g[[3]], g[[4]], g[[5]], nrow = 2)

ggsave(g2, filename=sprintf("%s/MB_vs_SLR_th%s.pdf", dir$figdir, sample.th), h=6, w=10)




######################################################


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
# Visualizing igraphs Edge betweenness-based modules

colvec <- c(brewer.pal(5, "Set1"), brewer.pal(5, "Set2"), rep("grey70", times=30))

for(i in 1:n.ex.rep){ 
group <- membership(mbet.mb[[i]])
extop10 <- as.vector(rownames(sort(table(group), decreasing=F)[1:(length(unique(group))-10)]))
g.extop10 <- group
g.extop10[g.extop10 %in% extop10] <- 999

g <- list()

g[[i]] <- ggraph(ig.mb.posi[[i]],layout = "stress")+ # OR layout = "nicely"
  geom_edge_link0(width=0.2,colour="darkslategray4")+
  geom_node_point(aes(col=factor(g.extop10)), shape=shape.vec[[i]], size=size.vec[[i]], show.legend=FALSE)+
  scale_color_manual(values = colvec)+
  guides()+
  theme_graph()+
  theme(legend.position = "right")+
  theme(text = element_text("Helvetica"))

ggsave(g[[i]], filename=sprintf("%s/SpiecEasi_MB_betModule_th%s_Tank%s.pdf", dir$figdir, sample.th, i), h=5.5, w=16)
}

######################################################
# Setting colors Phylum

list.phylum <- unique(taxa$Phylum)
category <- sort(table(taxa$Phylum), decreasing=T)
category.name <- names(category)
category.exUn <- category.name[category.name!="Unidentified"]
category.others <- category.exUn[21:length(category.exUn)]
category.list <- c(category.exUn, "Unidentified")

col1 <- brewer.pal(8, "Set1")
col2 <- brewer.pal(4, "Set2")
col3 <- brewer.pal(7, "Pastel1")
#col4 <- brewer.pal(7, "Pastel1")
#col5 <- brewer.pal(3, "Pastel2")

col.vec <- c("skyblue", col1, col2, col3, rep("grey80", times = length(category.others)), "grey40")

col.list <- data.frame(cbind(category.list, col.vec))
colnames(col.list) <- c("Phylum", "Color")

######################################################
# Visualizing igraphs Phylum

for(i in 1:n.ex.rep){ 
node.merge <- merge(asv.info[[i]], col.list, by='Phylum', all = FALSE, sort = TRUE)
node <- node.merge[order(node.merge$ASV.ID),]
Phylum <- node$Phylum
selected <- data.frame(unique(Phylum))
colnames(selected) <- c("Phylum")
selected.list <- merge(selected, col.list, by='Phylum', all = FALSE, sort = TRUE)

selected.list2 <- col.list[col.list$Phylum %in% unique(selected.list$Phylum), ]
Phylum <- factor(Phylum, levels=selected.list2$Phylum)

g <- list()

g[[i]] <- ggraph(ig.mb.posi[[i]],layout = "stress")+ # OR layout = "nicely"
  geom_edge_link0(width=0.2,colour="darkslategray4")+
  geom_node_point(aes(col=Phylum), shape=shape.vec[[i]], size=size.vec[[i]], show.legend=FALSE)+
  scale_color_manual(values = selected.list2$Color)+
  guides()+
  theme_graph()+
  theme(legend.position = "right")+
  theme(text = element_text("Helvetica"))

ggsave(g[[i]], filename=sprintf("%s/SpiecEasi_MB_Phylum_th%s_Tank%s.pdf", dir$figdir, sample.th, i), h=5.5, w=16)
}


for(i in 1:n.ex.rep){ 
node.merge <- merge(asv.info[[i]], col.list, by='Phylum', all = FALSE, sort = TRUE)
node <- node.merge[order(node.merge$ASV.ID),]
Phylum <- node$Phylum
selected <- data.frame(unique(Phylum))
colnames(selected) <- c("Phylum")
selected.list <- merge(selected, col.list, by='Phylum', all = FALSE, sort = TRUE)

selected.list2 <- col.list[col.list$Phylum %in% unique(selected.list$Phylum), ]
Phylum <- factor(Phylum, levels=selected.list2$Phylum)

g <- list()

g[[i]] <- ggraph(ig.mb.posi[[i]],layout = "stress")+ # OR layout = "nicely"
  geom_edge_link0(width=0.2,colour="darkslategray4")+
  geom_node_point(aes(col=Phylum), shape=shape.vec[[i]], size=size.vec[[i]], show.legend=TRUE)+
  scale_color_manual(values = selected.list2$Color)+
  guides()+
  theme_graph()+
  theme(legend.position = "right")+
  theme(text = element_text("Helvetica"))
  
ggsave(g[[i]], filename=sprintf("%s/SpiecEasi_MB_Phylum_th%s_Tank%s_Caption.pdf", dir$figdir, sample.th, i), h=5.5, w=16)
}

######################################################
# Visualizing igraphs Family
# Setting colors

list.family <- unique(taxa$Family)
category <- sort(table(taxa$Family), decreasing=T)
category.name <- names(category)
category.exUn <- category.name[category.name!="Unidentified"]
category.others <- category.exUn[39:length(category.exUn)]
category.list <- c(category.exUn, "Unidentified")

col1 <- brewer.pal(8, "Set1")
col2 <- brewer.pal(7, "Set2")
col3 <- brewer.pal(8, "Set3")
col4 <- brewer.pal(8, "Pastel1")
col5 <- brewer.pal(6, "Pastel2")

col.vec <- c("skyblue", col1, col2, col3, col4, col5, rep("grey80", times = length(category.others)), "grey40")

col.list <- data.frame(cbind(category.list, col.vec))
colnames(col.list) <- c("Family", "Color")


######################################################
# Visualizing igraphs Family

for(i in 1:n.ex.rep){ 
node.merge <- merge(asv.info[[i]], col.list, by='Family', all = FALSE, sort = TRUE)
node <- node.merge[order(node.merge$ASV.ID),]
Family <- node$Family
selected <- data.frame(unique(Family))
colnames(selected) <- c("Family")
selected.list <- merge(selected, col.list, by='Family', all = FALSE, sort = TRUE)

selected.list2 <- col.list[col.list$Family %in% unique(selected.list$Family), ]
Family <- factor(Family, levels=selected.list2$Family)

g <- list()

g[[i]] <- ggraph(ig.mb.posi[[i]],layout = "stress")+ # OR layout = "nicely"
  geom_edge_link0(width=0.2,colour="darkslategray4")+
  geom_node_point(aes(col=Family), shape=shape.vec[[i]], size=size.vec[[i]], show.legend=FALSE)+
  scale_color_manual(values = selected.list2$Color)+
  guides()+
  theme_graph()+
  theme(legend.position = "right")+
  theme(text = element_text("Helvetica"))

ggsave(g[[i]], filename=sprintf("%s/SpiecEasi_MB_th%s_Tank%s_Family.pdf", dir$figdir, sample.th, i), h=5.5, w=18)

}


for(i in 1:n.ex.rep){ 
node.merge <- merge(asv.info[[i]], col.list, by='Family', all = FALSE, sort = TRUE)
node <- node.merge[order(node.merge$ASV.ID),]
Family <- node$Family
selected <- data.frame(unique(Family))
colnames(selected) <- c("Family")
selected.list <- merge(selected, col.list, by='Family', all = FALSE, sort = TRUE)

selected.list2 <- col.list[col.list$Family %in% unique(selected.list$Family), ]
Family <- factor(Family, levels=selected.list2$Family)

g <- list()

g[[i]] <- ggraph(ig.mb.posi[[i]],layout = "stress")+ # OR layout = "nicely"
  geom_edge_link0(width=0.2,colour="darkslategray4")+
  geom_node_point(aes(col=Family), shape=shape.vec[[i]], size=size.vec[[i]], show.legend=TRUE)+
  scale_color_manual(values = selected.list2$Color)+
  guides()+
  theme_graph()+
  theme(legend.position = "right")+
  theme(text = element_text("Helvetica"))
  
ggsave(g[[i]], filename=sprintf("%s/SpiecEasi_MB_th%s_Tank%s_Caption_Family.pdf", dir$figdir, sample.th, i), h=5.5, w=18)

}



######################################################
### Node-level correlations on networks: direct correlations with pH, DO, and AS

load("03_Surrogate/01_output/RData/Surrogate_ASVs_th30_1e+05.RData")

# -- Create directory to save
dir <- make.dir('05_SpiecEasi/01_output')

## pH ##

sr.pH <- list()
g <- list()

for (i in 1:n.ex.rep) {
colnames(pH.result[[i]])[1:2] <- c("ASV.ID", "Correlation")
mode(pH.result[[i]]$ASV.ID) <- "character"
d <- merge(pH.result[[i]], taxa, by.x="ASV.ID", by.y="ID")
FDR <- p.adjust(apply(matrix(as.numeric(unlist(d[,3:4])), ncol=2), 1, FUN=min), method="fdr")
Cor <- as.numeric(unlist(d$Correlation))
x <- na.omit(data.frame(cbind(Cor, FDR)))
g[[i]] <- ggplot(x, aes(x=Cor, y=FDR)) + geom_point() + labs(title=sprintf('Tank %s', i), x='Correlation with pH', y='False discovery rate (FDR)')
sr.pH[[i]] <- cbind(d, FDR)
}

g2 <- gridExtra::grid.arrange(g[[1]], g[[2]], g[[3]], g[[4]], g[[5]], nrow = 1)

ggsave(g2, filename=sprintf("%s/FDR_pH_th%s.pdf", dir$figdir, sample.th), h=4, w=16)


g <- list()

for(i in 1:n.ex.rep){ 
Correlation <- as.numeric(sr.pH[[i]]$Correlation)
g[[i]] <- ggraph(ig.mb.posi[[i]],layout = "stress")+ # OR layout = "nicely"
  geom_edge_link0(width=0.2,colour="darkslategray4")+
  geom_node_point(aes(col= Correlation), shape=shape.vec[[i]], size=size.vec[[i]], show.legend=FALSE)+
  scale_colour_gradient2(low = "blue", mid = "grey90", high = "red", na.value = "white", guide = "colourbar")+
  guides()+
  theme_graph()+
  theme(legend.position = "right")+
  theme(text = element_text("Helvetica"))

ggsave(g[[i]], filename=sprintf("%s/SpiecEasi_Corr_pH_th%s_Tank%s.pdf", dir$figdir, sample.th, i), h=5.5, w=16)
}


for(i in 1:n.ex.rep){ 
Correlation <- as.numeric(sr.pH[[i]]$Correlation)
g[[i]] <- ggraph(ig.mb.posi[[i]],layout = "stress")+ # OR layout = "nicely"
  geom_edge_link0(width=0.2,colour="darkslategray4")+
  geom_node_point(aes(col= Correlation), shape=shape.vec[[i]], size=size.vec[[i]], show.legend=TRUE)+
  scale_colour_gradient2(low = "blue", mid = "grey90", high = "red", na.value = "white", guide = "colourbar")+
  guides()+
  theme_graph()+
  theme(legend.position = "right")+
  theme(text = element_text("Helvetica"))

ggsave(g[[i]], filename=sprintf("%s/SpiecEasi_Corr_pH_th%s_Tank%s_caption.pdf", dir$figdir, sample.th, i), h=5.5, w=16)
}


## DO ##

sr.DO <- list()
g <- list()

for (i in 1:n.ex.rep) {
colnames(DO.result[[i]])[1:2] <- c("ASV.ID", "Correlation")
mode(DO.result[[i]]$ASV.ID) <- "character"
d <- merge(DO.result[[i]], taxa, by.x="ASV.ID", by.y="ID")
FDR <- p.adjust(apply(matrix(as.numeric(unlist(d[,3:4])), ncol=2), 1, FUN=min), method="fdr")
Cor <- as.numeric(unlist(d$Correlation))
x <- na.omit(data.frame(cbind(Cor, FDR)))
g[[i]] <- ggplot(x, aes(x=Cor, y=FDR)) + geom_point() + labs(title=sprintf('Tank %s', i), x='Correlation with DO', y='False discovery rate (FDR)')
sr.DO[[i]] <- cbind(d, FDR)
}

g2 <- gridExtra::grid.arrange(g[[1]], g[[2]], g[[3]], g[[4]], g[[5]], nrow = 1)

ggsave(g2, filename=sprintf("%s/FDR_DO_th%s.pdf", dir$figdir, sample.th), h=4, w=16)


g <- list()

for(i in 1:n.ex.rep){ 
Correlation <- as.numeric(sr.DO[[i]]$Correlation)
g[[i]] <- ggraph(ig.mb.posi[[i]],layout = "stress")+ # OR layout = "nicely"
  geom_edge_link0(width=0.2,colour="darkslategray4")+
  geom_node_point(aes(col= Correlation), shape=shape.vec[[i]], size=size.vec[[i]], show.legend=FALSE)+
  scale_colour_gradient2(low = "blue", mid = "grey90", high = "red", na.value = "white", guide = "colourbar")+
  guides()+
  theme_graph()+
  theme(legend.position = "right")+
  theme(text = element_text("Helvetica"))

ggsave(g[[i]], filename=sprintf("%s/SpiecEasi_Corr_DO_th%s_Tank%s.pdf", dir$figdir, sample.th, i), h=5.5, w=16)
}


for(i in 1:n.ex.rep){ 
Correlation <- as.numeric(sr.DO[[i]]$Correlation)
g[[i]] <- ggraph(ig.mb.posi[[i]],layout = "stress")+ # OR layout = "nicely"
  geom_edge_link0(width=0.2,colour="darkslategray4")+
  geom_node_point(aes(col= Correlation), shape=shape.vec[[i]], size=size.vec[[i]], show.legend=TRUE)+
  scale_colour_gradient2(low = "blue", mid = "grey90", high = "red", na.value = "white", guide = "colourbar")+
  guides()+
  theme_graph()+
  theme(legend.position = "right")+
  theme(text = element_text("Helvetica"))

ggsave(g[[i]], filename=sprintf("%s/SpiecEasi_Corr_DO_th%s_Tank%s_caption.pdf", dir$figdir, sample.th, i), h=5.5, w=16)
}


## AS ##

sr.AS <- list()
g <- list()

for (i in 1:n.ex.rep) {
colnames(AS.result[[i]])[1:2] <- c("ASV.ID", "Correlation")
mode(AS.result[[i]]$ASV.ID) <- "character"
d <- merge(AS.result[[i]], taxa, by.x="ASV.ID", by.y="ID")
FDR <- p.adjust(apply(matrix(as.numeric(unlist(d[,3:4])), ncol=2), 1, FUN=min), method="fdr")
Cor <- as.numeric(unlist(d$Correlation))
x <- na.omit(data.frame(cbind(Cor, FDR)))
g[[i]] <- ggplot(x, aes(x=Cor, y=FDR)) + geom_point() + labs(title=sprintf('Tank %s', i), x="Correlation with eels' activity level", y='False discovery rate (FDR)')
sr.AS[[i]] <- cbind(d, FDR)
}

g2 <- gridExtra::grid.arrange(g[[1]], g[[2]], g[[3]], g[[4]], g[[5]], nrow = 1)

ggsave(g2, filename=sprintf("%s/FDR_AS_th%s.pdf", dir$figdir, sample.th), h=4, w=16)


g <- list()

for(i in 1:n.ex.rep){ 
Correlation <- as.numeric(sr.AS[[i]]$Correlation)
g[[i]] <- ggraph(ig.mb.posi[[i]],layout = "stress")+ # OR layout = "nicely"
  geom_edge_link0(width=0.2,colour="darkslategray4")+
  geom_node_point(aes(col= Correlation), shape=shape.vec[[i]], size=size.vec[[i]], show.legend=FALSE)+
  scale_colour_gradient2(low = "blue", mid = "grey90", high = "red", na.value = "white", guide = "colourbar")+
  guides()+
  theme_graph()+
  theme(legend.position = "right")+
  theme(text = element_text("Helvetica"))

ggsave(g[[i]], filename=sprintf("%s/SpiecEasi_Corr_AS_th%s_Tank%s.pdf", dir$figdir, sample.th, i), h=5.5, w=16)
}

g <- list()

for(i in 1:n.ex.rep){ 
Correlation <- as.numeric(sr.AS[[i]]$Correlation)
g[[i]] <- ggraph(ig.mb.posi[[i]],layout = "stress")+ # OR layout = "nicely"
  geom_edge_link0(width=0.2,colour="darkslategray4")+
  geom_node_point(aes(col= Correlation), shape=shape.vec[[i]], size=size.vec[[i]], show.legend=TRUE)+
  scale_colour_gradient2(low = "blue", mid = "grey90", high = "red", na.value = "white", guide = "colourbar")+
  guides()+
  theme_graph()+
  theme(legend.position = "right")+
  theme(text = element_text("Helvetica"))

ggsave(g[[i]], filename=sprintf("%s/SpiecEasi_Corr_AS_th%s_Tank%s_caption.pdf", dir$figdir, sample.th, i), h=5.5, w=16)
}


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
g[[i]] <- ggplot(x, aes(x=Cor, y=FDR)) + geom_point() + labs(title=sprintf('Tank %s', i), x="Correlation with eels' activity level", y='False discovery rate (FDR)')
pr.AS_pH[[i]] <- cbind(d, FDR)
}

g2 <- gridExtra::grid.arrange(g[[1]], g[[2]], g[[3]], g[[4]], g[[5]], nrow = 1)

ggsave(g2, filename=sprintf("%s/FDR_pAScontrPH_th%s.pdf", dir$figdir, sample.th), h=4, w=16)


g <- list()

for(i in 1:n.ex.rep){ 
Correlation <- as.numeric(pr.AS_pH[[i]]$Correlation)
g[[i]] <- ggraph(ig.mb.posi[[i]],layout = "stress")+ # OR layout = "nicely"
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
g[[i]] <- ggraph(ig.mb.posi[[i]],layout = "stress")+ # OR layout = "nicely"
  geom_edge_link0(width=0.2,colour="darkslategray4")+
  geom_node_point(aes(col= Correlation), shape=shape.vec[[i]], size=size.vec[[i]], show.legend=TRUE)+
  scale_colour_gradient2(low = "blue", mid = "grey90", high = "red", na.value = "white", guide = "colourbar")+
  guides()+
  theme_graph()+
  theme(legend.position = "right")+
  theme(text = element_text("Helvetica"))

ggsave(g[[i]], filename=sprintf("%s/SpiecEasi_pAScontrPH_th%s_Tank%s_caption.pdf", dir$figdir, sample.th, i), h=5.5, w=16)
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
g[[i]] <- ggplot(x, aes(x=Cor, y=FDR)) + geom_point() + labs(title=sprintf('Tank %s', i), x="Correlation with eels' activity level", y='False discovery rate (FDR)')
pr.AS_DO[[i]] <- cbind(d, FDR)
}

g2 <- gridExtra::grid.arrange(g[[1]], g[[2]], g[[3]], g[[4]], g[[5]], nrow = 1)

ggsave(g2, filename=sprintf("%s/FDR_pAScontrDO_th%s.pdf", dir$figdir, sample.th), h=4, w=16)


g <- list()

for(i in 1:n.ex.rep){ 
Correlation <- as.numeric(pr.AS_DO[[i]]$Correlation)
g[[i]] <- ggraph(ig.mb.posi[[i]],layout = "stress")+ # OR layout = "nicely"
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
g[[i]] <- ggraph(ig.mb.posi[[i]],layout = "stress")+ # OR layout = "nicely"
  geom_edge_link0(width=0.2,colour="darkslategray4")+
  geom_node_point(aes(col= Correlation), shape=shape.vec[[i]], size=size.vec[[i]], show.legend=TRUE)+
  scale_colour_gradient2(low = "blue", mid = "grey90", high = "red", na.value = "white", guide = "colourbar")+
  guides()+
  theme_graph()+
  theme(legend.position = "right")+
  theme(text = element_text("Helvetica"))

ggsave(g[[i]], filename=sprintf("%s/SpiecEasi_pAScontrDO_th%s_Tank%s_caption.pdf", dir$figdir, sample.th, i), h=5.5, w=16)
}

######################################################
## Node-level data outputs

ASV.table <- list()

for(i in 1:n.ex.rep){ 
tank <- i
mode(pH.result[[i]]$ASV.ID) <- "character"
mode(DO.result[[i]]$ASV.ID) <- "character"
mode(AS.result[[i]]$ASV.ID) <- "character"
pH <- merge(pH.result[[i]], taxa, by.x="ASV.ID", by.y="ID")
DO <- merge(DO.result[[i]], taxa, by.x="ASV.ID", by.y="ID")
AS <- merge(AS.result[[i]], taxa, by.x="ASV.ID", by.y="ID")
Correlation.pH <- pH$Correlation
Correlation.DO <- DO$Correlation
Correlation.AS <- AS$Correlation
Betweenness <- betweenness(ig.mb.posi[[i]], normalized=TRUE)
Eigenvector_centrality <- centr_eigen(ig.mb.posi[[i]], normalized=TRUE)$vector
Degree <- degree(ig.mb.posi[[i]])
ASV.ID <- pH$ASV.ID
Module <- membership(mbet.mb[[i]])
table <- data.frame(cbind(ASV.ID, Module, Correlation.pH, Correlation.DO, Correlation.AS, Betweenness, Eigenvector_centrality, Degree))
mode(table$ASV.ID) <- "character"
table2 <- merge(table, taxa, by.x="ASV.ID", by.y="ID")
ASV.table[[i]] <- apply(table2,2,as.character)
write.table(ASV.table[[i]], file=sprintf("%s/ASV.Table_th%s_Tank%s.txt", dir$tabledir, sample.th, tank), sep="\t", quote=FALSE, row.names=FALSE)
}

saveRDS(sr.pH, sprintf("%s/ASV.Table.SpiecEasi.MB_th%s.rds", dir$rdsdir, sample.th))

tbl <- list()

for(i in 1:n.ex.rep){ 
	Tank <- rep(sprintf('Tank_%s', i), times=nrow(ASV.table[[i]]))
	tbl[[i]] <- cbind(Tank, ASV.table[[i]])
}

ASV.All <- data.frame(rbind(tbl[[1]], tbl[[2]], tbl[[3]], tbl[[4]], tbl[[5]]))

ASV.All[, 4:9] <- apply(ASV.All[, 4:9],2,as.numeric)


write.table(ASV.All, file=sprintf("%s/ASV.Table_All_th%s.txt", dir$tabledir, sample.th), sep="\t", quote=FALSE, row.names=FALSE)


######################################################
## Node-level data: Correlations to pH, DO, or AS

table.n <- ASV.All

g1 <- ggplot(table.n, aes(x=Degree, y=Correlation.pH, color=Tank, shape=Tank)) + geom_point() + labs(title='', x='Degree centrality', y='Correlation with pH') + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) + theme(text = element_text(size = 11))

g2 <- ggplot(table.n, aes(x=Betweenness, y=Correlation.pH, color=Tank, shape=Tank)) + geom_point() + labs(title='', x='Betweenness centrality', y='Correlation with pH') + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) + theme(text = element_text(size = 11))

g3 <- ggplot(table.n, aes(x=Eigenvector_centrality, y=Correlation.pH, color=Tank, shape=Tank)) + geom_point() + labs(title='', x='Eigenvector centrality', y='Correlation with pH') + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) + theme(text = element_text(size = 11))

g4 <- ggplot(table.n, aes(x=Degree, y=Correlation.DO, color=Tank, shape=Tank)) + geom_point() + labs(title='', x='Degree centrality', y='Correlation with DO') + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) + theme(text = element_text(size = 11))

g5 <- ggplot(table.n, aes(x=Betweenness, y=Correlation.DO, color=Tank, shape=Tank)) + geom_point() + labs(title='', x='Betweenness centrality', y='Correlation with DO') + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) + theme(text = element_text(size = 11))

g6 <- ggplot(table.n, aes(x=Eigenvector_centrality, y=Correlation.DO, color=Tank, shape=Tank)) + geom_point() + labs(title='', x='Eigenvector centrality', y='Correlation with DO') + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) + theme(text = element_text(size = 11))

g7 <- ggplot(table.n, aes(x=Degree, y=Correlation.AS, color=Tank, shape=Tank)) + geom_point() + labs(title='', x='Degree centrality', y="Correlation with eels' activity level") + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) + theme(text = element_text(size = 11))

g8 <- ggplot(table.n, aes(x=Betweenness, y=Correlation.AS, color=Tank, shape=Tank)) + geom_point() + labs(title='', x='Betweenness centrality', y="Correlation with eels' activity level") + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) + theme(text = element_text(size = 11))

g9 <- ggplot(table.n, aes(x=Eigenvector_centrality, y=Correlation.AS, color=Tank, shape=Tank)) + geom_point() + labs(title='', x='Eigenvector centrality', y="Correlation with eels' activity level") + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) + theme(text = element_text(size = 11))

g10 <- gridExtra::grid.arrange(g1, g2, g3, g4, g5, g6, g7, g8, g9, nrow = 3)

ggsave(g10, filename=sprintf("%s/Centrality_Correlation_th%s.pdf", dir$figdir, sample.th), h=16, w=16)

######################################################
## Comparison correlations with factors

g1 <- ggplot(table.n, aes(x=Correlation.pH, y=Correlation.DO, color=Tank, shape=Tank)) + geom_point() + labs(title='', x='Correlation with pH', y='Correlation with DO') + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) + theme(text = element_text(size = 11))

g2 <- ggplot(table.n, aes(x=Correlation.pH, y=Correlation.AS, color=Tank, shape=Tank)) + geom_point() + labs(title='', x='Correlation with pH', y="Correlation with eels' activity level") + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) + theme(text = element_text(size = 11))

g3 <- ggplot(table.n, aes(x=Correlation.DO, y=Correlation.AS, color=Tank, shape=Tank)) + geom_point() + labs(title='', x='Correlation with DO', y="Correlation with eels' activity level") + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) + theme(text = element_text(size = 11))


g10 <- gridExtra::grid.arrange(g1, g2, g3, nrow = 1)

ggsave(g10, filename=sprintf("%s/Comparison_Correlations_th%s.pdf", dir$figdir, sample.th), h=5, w=16)


######################################################
## Correlatinos between base env factors

param <- as.matrix(read.table("Table/Parameters_3.txt", row.name=1, header=T, na.strings = "NA"))

param2 <- melt(param)
colnames(param2) <- c("Day", "Parameter", "value")

d.pH <- param2[grep("pH", param2$Parameter), ]
DO <- param2[grep("DO", param2$Parameter), ][,3]
AS <- param2[grep("AS", param2$Parameter), ][,3]

env <- cbind(d.pH, DO, AS)
colnames(env) <- c('Day', 'Tank', 'pH', 'DO', 'AS')


g <- list()
for(i in 1:n.ex.rep){ 
	env.sub <- subset(env, env$Tank==sprintf('T%s_pH', i))
g[[i]] <- ggplot(na.omit(env.sub), aes(x=pH, y=AS)) + geom_point() + labs(title=sprintf('Tank %s', i), x='pH', y='Eel appetite score') + theme(text = element_text(size = 11))
}

g10 <- gridExtra::grid.arrange(g[[1]], g[[2]], g[[3]], g[[4]], g[[5]], nrow = 1)

ggsave(g10, filename=sprintf("%s/pH_vs_AS_th%s.pdf", dir$figdir, sample.th), h=4, w=16)


g <- list()
for(i in 1:n.ex.rep){ 
	env.sub <- subset(env, env$Tank==sprintf('T%s_pH', i))
g[[i]] <- ggplot(na.omit(env.sub), aes(x=pH, y=DO)) + geom_point() + labs(title=sprintf('Tank %s', i), x='pH', y='DO') + theme(text = element_text(size = 11))
}

g10 <- gridExtra::grid.arrange(g[[1]], g[[2]], g[[3]], g[[4]], g[[5]], nrow = 1)

ggsave(g10, filename=sprintf("%s/pH_vs_DO_th%s.pdf", dir$figdir, sample.th), h=4, w=16)


g <- list()
for(i in 1:n.ex.rep){ 
	env.sub <- subset(env, env$Tank==sprintf('T%s_pH', i))
g[[i]] <- ggplot(na.omit(env.sub), aes(x=DO, y=AS)) + geom_point() + labs(title=sprintf('Tank %s', i), x='DO', y='Eel appetite score') + theme(text = element_text(size = 11))
}

g10 <- gridExtra::grid.arrange(g[[1]], g[[2]], g[[3]], g[[4]], g[[5]], nrow = 1)

ggsave(g10, filename=sprintf("%s/DO_vs_AS_th%s.pdf", dir$figdir, sample.th), h=4, w=16)

######################################################
### Node-level correlations on networks: partial correlations with pH, DO, and AS


xx <- list()
for (i in 1:n.ex.rep) {
	xx[[i]] <- data.frame(cbind(Tank=rep(sprintf("Tank_%s", i), times=nrow(AS.result_Contr.pH[[i]])), ASV.ID=AS.result_Contr.pH[[i]]$ASV.ID, Partial.r.AS_ctlrPH=AS.result_Contr.pH[[i]]$Correlation))
}

pAS1 <- rbind(xx[[1]], xx[[2]], xx[[3]], xx[[4]], xx[[5]])

yy <- list()
for (i in 1:n.ex.rep) {
	colnames(AS.result_Contr.DO[[i]])[1:2] <- c("ASV.ID", "Correlation")
mode(AS.result_Contr.DO[[i]]$ASV.ID) <- "character"
	yy[[i]] <- data.frame(cbind(Tank=rep(sprintf("Tank_%s", i), times=nrow(AS.result_Contr.DO[[i]])), ASV.ID=AS.result_Contr.DO[[i]]$ASV.ID, Partial.r.AS_ctlrDO=AS.result_Contr.DO[[i]]$Correlation))
}

pAS2 <- rbind(yy[[1]], yy[[2]], yy[[3]], yy[[4]], yy[[5]])

table.n2_0 <- data.frame(cbind(table.n, pAS1, pAS2))
cbind(table.n2_0$ASV.ID, table.n2_0$ASV.ID.1, table.n2_0$ASV.ID.2)

table.n2 <- table.n2_0[,-c(17,18,20,21)]

table.n2[,c(4:9,17,18)] <- apply(table.n2[, c(4:9,17,18)],2,as.numeric)


#write.table(table.n2, file=sprintf("%s/ASV.Table_All_th%s_2.txt", dir$tabledir, sample.th), sep="\t", quote=FALSE, row.names=FALSE)

######################################################

######################################################
## Partial AS controlling pH and DO
#library(maptools)


g1 <- ggplot(table.n2, aes(x=Partial.r.AS_ctlrPH, y=Partial.r.AS_ctlrDO, color=Tank, shape=Tank)) + geom_point() + labs(title='', x='Partial correlation with eel appetite score (pH-effects controlled)', y='Partial correlation with eel appetite score (DO-effects controlled)') + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) + theme(text = element_text(size = 11))

g2 <- ggplot(table.n2, aes(x= Betweenness, y=Partial.r.AS_ctlrPH, color=Tank, shape=Tank)) + geom_point() + labs(title='', x='Betweenness centrality', y='Partial correlation with eel appetite score (pH-effects controlled)') + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) + theme(text = element_text(size = 11))

g3 <- ggplot(table.n2, aes(x= Eigenvector_centrality, y=Partial.r.AS_ctlrPH, color=Tank, shape=Tank)) + geom_point() + labs(title='', x='Eigenvector centrality', y='Partial correlation with eel appetite score (pH-effects controlled)') + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) + theme(text = element_text(size = 11))

g4 <- ggplot(table.n2, aes(x= Degree, y=Partial.r.AS_ctlrPH, color=Tank, shape=Tank)) + geom_point() + labs(title='', x='Degree centrality', y='Partial correlation with eel appetite score (pH-effects controlled)') + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) + theme(text = element_text(size = 11))

g5 <- ggplot(table.n2, aes(x= Correlation.pH, y=Partial.r.AS_ctlrPH, color=Tank, shape=Tank)) + geom_point() + labs(title='', x='Correlation with pH', y='Partial correlation with eel appetite score (pH-effects controlled)') + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) + theme(text = element_text(size = 11))

g6 <- ggplot(table.n2, aes(x= Correlation.DO, y=Partial.r.AS_ctlrPH, color=Tank, shape=Tank)) + geom_point() + labs(title='', x='Correlation with DO', y='Partial correlation with eel appetite score (pH-effects controlled)') + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) + theme(text = element_text(size = 11))

ggsave(g1, filename=sprintf("%s/PartialCorr.AS_pHvsDO_th%s.pdf", dir$figdir, sample.th), h=5, w=5)

g10 <- gridExtra::grid.arrange(g1, g2, g3, g4, g5, g6, nrow = 2)

ggsave(g10, filename=sprintf("%s/PartialCorr.AS_AllPlots_th%s.pdf", dir$figdir, sample.th), h=10, w=16)


######################################################

write.table(table.n2, file=sprintf('%s/Partial.Corr.AS.txt', dir$tabledir), row.name=F, quote=F, sep='\t')

write.table(table.n2, file=sprintf('%s/ASV.Table_All_th%s_partial.txt', dir$tabledir, sample.th), sep='\t', quote=F, row.names=F)


######################################################
## Visualizing controlled appetite score: Module level
## Using Partial.r.AS_ctlrPH as controlled correlation

r.module <- c()

for (i in 1:5) {
	d <- subset(table.n2, table.n2$Tank==sprintf('Tank_%s', i))
	for (j in 1:length(unique(d$Module))) {
		d2 <- subset(d, d$Module==j)
		Tank <- i
		Module <- j
		N.nodes <- nrow(d2)
		Mean.Contl.Corr.AS <- mean(d2$Partial.r.AS_ctlrPH)
		SD.Contl.Corr.AS <- sd(d2$Partial.r.AS_ctlrPH)
		r <- c(Tank, Module, N.nodes, Mean.Contl.Corr.AS, SD.Contl.Corr.AS)
		r.module <- rbind(r.module, r)
	}
}

r.module <- data.frame(r.module)
colnames(r.module) <- c('Tank', 'Module', 'N.nodes', 'Mean.Contl.Corr.AS', 'SD.Contl.Corr.AS')
SE.Contl.Corr.AS <- r.module$SD.Contl.Corr.AS / (r.module$N.nodes^0.5)
r.module <- cbind(r.module, SE.Contl.Corr.AS)

write.table(r.module, file=sprintf('%s/Controlled.AS.Module.txt', dir$tabledir), row.name=F, quote=F, sep='\t')

r.module.m4 <- subset(r.module, r.module$N.nodes >= 4)

g1 <- ggplot(r.module.m4, aes(x=N.nodes, y=Mean.Contl.Corr.AS, color=as.character(Tank), shape=as.character(Tank))) + geom_point() + labs(title='', x='Number of ASVs within the module', y='Mean partial correlation with eel appetite score') + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) + theme(text = element_text(size = 11)) 

g2 <- g1 + geom_errorbar(aes(ymin=Mean.Contl.Corr.AS-SE.Contl.Corr.AS, ymax=Mean.Contl.Corr.AS+SE.Contl.Corr.AS)) + geom_hline(yintercept=0, color='red', lty=2)

ggsave(g2, filename=sprintf("%s/PartialCorr.AS_Module_th%s.pdf", dir$figdir, sample.th), h=4, w=5)

######################################################
## Visualizing controlled appetite score: Node level: Consistency among tanks

freq <- data.frame(table(table.n2$ASV.ID))
colnames(freq) <- c('ASV.ID', 'Frequency')
table.n4 <- merge(table.n2, freq, by='ASV.ID')

ASV.5 <- subset(table.n4, table.n4$Frequency >= 5)

g1 <- ggplot(ASV.5, aes(x=reorder(x=ASV.ID, X=-Partial.r.AS_ctlrPH, FUN=mean), y=Partial.r.AS_ctlrPH, color=Tank, shape=Tank)) + geom_point() + labs(title='', x='ASV', y="Partial correlation with eels' activity level") + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) + theme(text = element_text(size = 11)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_hline(yintercept=0, color='red', lty=2)

g2 <- ggplot(ASV.5) + geom_boxplot(aes(x=reorder(x=ASV.ID, X=-Partial.r.AS_ctlrPH, FUN=mean), y=Partial.r.AS_ctlrPH)) + labs(title='', x='ASV', y="Partial correlation with eels' activity level") + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) + theme(text = element_text(size = 11)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_hline(yintercept=0, color='red', lty=2) 

g10 <- gridExtra::grid.arrange(g1, g2, nrow = 2)

ggsave(g10, filename=sprintf("%s/PartialCorr.AS_ASV_5tanks_th%s.pdf", dir$figdir, sample.th), h=10, w=16)


table.genus <- subset(table.n2, table.n2$Genus!="Unidentified")

g1 <- ggplot(table.genus, aes(x=reorder(x=Genus, X=-Partial.r.AS_ctlrPH, FUN=mean), y=Partial.r.AS_ctlrPH, color=Tank, shape=Tank)) + geom_point() + labs(title='', x='Genus', y="Partial correlation with eels' activity level") + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) + theme(text = element_text(size = 11)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_hline(yintercept=0, color='red', lty=2)

g2 <- ggplot(table.genus) + geom_boxplot(aes(x=reorder(x=Genus, X=-Partial.r.AS_ctlrPH, FUN=mean), y=Partial.r.AS_ctlrPH)) + labs(title='', x='Genus', y="Partial correlation with eels' activity level") + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) + theme(text = element_text(size = 11)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_hline(yintercept=0, color='red', lty=2)

g10 <- gridExtra::grid.arrange(g1, g2, nrow = 2)

ggsave(g10, filename=sprintf("%s/PartialCorr.AS_Genus_th%s.pdf", dir$figdir, sample.th), h=10, w=16)



table.family <- subset(table.n2, table.n2$Family!="Unidentified")

g1 <- ggplot(table.family, aes(x=reorder(x=Family, X=-Partial.r.AS_ctlrPH, FUN=mean), y=Partial.r.AS_ctlrPH, color=Tank, shape=Tank)) + geom_point() + labs(title='', x='Family', y="Partial correlation with eels' activity level") + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) + theme(text = element_text(size = 11)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_hline(yintercept=0, color='red', lty=2)

g2 <- ggplot(table.family) + geom_boxplot(aes(x=reorder(x=Family, X=-Partial.r.AS_ctlrPH, FUN=mean), y=Partial.r.AS_ctlrPH)) + labs(title='', x='Family', y="Partial correlation with eels' activity level") + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) + theme(text = element_text(size = 11)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_hline(yintercept=0, color='red', lty=2)

g10 <- gridExtra::grid.arrange(g1, g2, nrow = 2)

ggsave(g10, filename=sprintf("%s/PartialCorr.AS_Family_th%s.pdf", dir$figdir, sample.th), h=10, w=16)


table.order <- subset(table.n2, table.n2$Order!="Unidentified")

g1 <- ggplot(table.order, aes(x=reorder(x=Order, X=-Partial.r.AS_ctlrPH, FUN=mean), y=Partial.r.AS_ctlrPH, color=Tank, shape=Tank)) + geom_point() + labs(title='', x='Order', y="Partial correlation with eels' activity level") + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) + theme(text = element_text(size = 11)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_hline(yintercept=0, color='red', lty=2)

g2 <- ggplot(table.order) + geom_boxplot(aes(x=reorder(x=Order, X=-Partial.r.AS_ctlrPH, FUN=mean), y=Partial.r.AS_ctlrPH)) + labs(title='', x='Order', y="Partial correlation with eels' activity level") + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) + theme(text = element_text(size = 11)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_hline(yintercept=0, color='red', lty=2)

g10 <- gridExtra::grid.arrange(g1, g2, nrow = 2)

ggsave(g10, filename=sprintf("%s/PartialCorr.AS_Order_th%s.pdf", dir$figdir, sample.th), h=10, w=16)



table.phylum <- subset(table.n2, table.n2$Phylum!="Unidentified")

g1 <- ggplot(table.phylum, aes(x=reorder(x=Phylum, X=-Partial.r.AS_ctlrPH, FUN=mean), y=Partial.r.AS_ctlrPH, color=Tank, shape=Tank)) + geom_point() + labs(title='', x='Phylum', y="Partial correlation with eels' activity level") + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) + theme(text = element_text(size = 11)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_hline(yintercept=0, color='red', lty=2)

g2 <- ggplot(table.phylum) + geom_boxplot(aes(x=reorder(x=Phylum, X=-Partial.r.AS_ctlrPH, FUN=mean), y=Partial.r.AS_ctlrPH)) + labs(title='', x='Phylum', y="Partial correlation with eels' activity level") + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) + theme(text = element_text(size = 11)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_hline(yintercept=0, color='red', lty=2)

g10 <- gridExtra::grid.arrange(g1, g2, nrow = 2)

ggsave(g10, filename=sprintf("%s/PartialCorr.AS_Phylum_th%s.pdf", dir$figdir, sample.th), h=10, w=8)


######################################################


g1 <- ggplot(table.n2, aes(x= Correlation.pH, y=Correlation.AS, color=Tank, shape=Tank)) + geom_point() + labs(title='', x='Correlation with pH', y='Correlation with eel activity level') + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) 

g2 <- ggplot(table.n2, aes(x= Correlation.DO, y=Correlation.AS, color=Tank, shape=Tank)) + geom_point() + labs(title='', x='Correlation with DO', y='Correlation with eel eel activity level') + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) 

g3 <- ggplot(table.n2, aes(x=Partial.r.AS_ctlrPH, y=Partial.r.AS_ctlrDO, color=Tank, shape=Tank)) + geom_point() + labs(title='', x='Partial correlation with eel activity level (pH-effects controlled)', y='Partial correlation with eel activity level (DO-effects controlled)') + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) 

g10 <- gridExtra::grid.arrange(g1, g2, g3, nrow = 1)

ggsave(g10, filename=sprintf("%s/Correlations_MainGraphs_th%s.pdf", dir$figdir, sample.th), h=4, w=16)


######################################################
stopCluster(cluster)

save.image(sprintf("%s/SpiecEasi_th%s_latent%s.RData", dir$rdatadir, sample.th, nlatent))

######################################################
