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


load(sprintf("05_SpiecEasi/01_output/RData/SpiecEasi_th%s_latent%s.RData", sample.th, nlatent))

# -- Create directory to save
dir <- make.dir('05_SpiecEasi/02_output')

######################################################
## Visualizing controlled appetite score: Node level: Consistency among tanks

freq <- data.frame(table(table.n2$ASV.ID))
colnames(freq) <- c('ASV.ID', 'Frequency')
table.n4 <- merge(table.n2, freq, by='ASV.ID')

ASV.5 <- subset(table.n4, table.n4$Frequency >= 5)

g1 <- ggplot(ASV.5, aes(x=reorder(x=ASV.ID, X=-Partial.r.AS_ctlrPH, FUN=mean), y=Partial.r.AS_ctlrPH, color=Tank, shape=Tank)) + geom_point() + labs(title='', x='ASV', y='Partial correlation with eel appetite score') + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) + theme(text = element_text(size = 14)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + geom_hline(yintercept=0, color='red', lty=2)

g2 <- ggplot(ASV.5) + geom_boxplot(aes(x=reorder(x=ASV.ID, X=-Partial.r.AS_ctlrPH, FUN=mean), y=Partial.r.AS_ctlrPH), color="deepskyblue4") + labs(title='', x='ASV', y='Partial correlation with eel appetite score') +  theme(text = element_text(size = 11)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + geom_hline(yintercept=0, color='red', lty=2) 

g10 <- gridExtra::grid.arrange(g1, g2, nrow = 2)

ggsave(g10, filename=sprintf("%s/PartialCorr.AS_ASV_5tanks_th%s.pdf", dir$figdir, sample.th), h=8, w=16)



table.genus_0 <- subset(table.n4, table.n4$Genus!="Unidentified")
table.genus <- subset(table.genus_0, table.genus_0$Frequency >= 5)

g1 <- ggplot(table.genus, aes(x=reorder(x=Genus, X=-Partial.r.AS_ctlrPH, FUN=mean), y=Partial.r.AS_ctlrPH, color=Tank, shape=Tank)) + geom_point() + labs(title='', x='Genus', y='Partial correlation with eel appetite score') + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) + theme(text = element_text(size = 14)) + theme(axis.text.x = element_text(face='italic', angle = 90, hjust = 1)) + geom_hline(yintercept=0, color='red', lty=2)

g2 <- ggplot(table.genus) + geom_boxplot(aes(x=reorder(x=Genus, X=-Partial.r.AS_ctlrPH, FUN=mean), y=Partial.r.AS_ctlrPH), color="deepskyblue4") + labs(title='', x='Genus', y='Partial correlation with eel appetite score') + theme(text = element_text(size = 14)) + theme(axis.text.x = element_text(face='italic', angle = 45, hjust = 1)) + geom_hline(yintercept=0, color='red', lty=2)

g10 <- gridExtra::grid.arrange(g1, g2, nrow = 2)

ggsave(g10, filename=sprintf("%s/PartialCorr.AS_Genus_5tanks_th%s.pdf", dir$figdi, sample.th), h=10, w=10)



