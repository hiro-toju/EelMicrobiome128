############################################################################
####
#### R script for Fujita (2019)
####
#### Visualizing community assembly
#### 2019.12.09 Fujita
#### 2022.02.24 Toju
#### R 4.1.2
#### Set working directory of 'MTS' folder -- setwd('~/Desktop/Microbiome_TimeSeries/MTS/')
#### 
############################################################################

setwd("/Users/toju/Dropbox/YAJIMA_Daii/Statistics_R")

source('functions/functions.R')
library(ggplot2)
library(ggsci)
library(reshape2)
library(RColorBrewer)
library(cowplot)

# -- Create directory to save
dir <- make.dir('./01_Timeseries_Graph/02_output')

col.list <- brewer.pal(5, "Set2")

param <- as.matrix(read.table("Table/Parameters_3.txt", row.name=1, header=T, na.strings = "NA"))

param2 <- melt(param)
colnames(param2) <- c("Day", "Parameter", "value")

pH <- param2[grep("pH", param2$Parameter), ]
DO <- param2[grep("DO", param2$Parameter), ]
AS <- param2[grep("AS", param2$Parameter), ]


##########################################################################

g <- list()

g[[1]] <- ggplot(pH, aes(x=Day, y=value, color=Parameter))
g[[1]] <- g[[1]] + geom_line() + labs(x= "Day", y= "pH")
g[[1]] <- g[[1]] + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) + theme(text = element_text(size = 11))
plot(g[[1]])

g[[2]] <- ggplot(DO, aes(x=Day, y=value, color=Parameter))
g[[2]] <- g[[2]] + geom_line() + labs(x= "Day", y= "Dissolved oxygen (mg/L)")
g[[2]] <- g[[2]] + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) + theme(text = element_text(size = 11))

plot(g[[2]])

g[[3]] <- ggplot(AS, aes(x=Day, y=value, color=Parameter))
g[[3]] <- g[[3]] + geom_line() + labs(x= "Day", y= "Eel appetite score")
g[[3]] <- g[[3]] + scale_color_hue(name="Tank", labels=c("1", "2", "3", "4", "5")) + theme(text = element_text(size = 11))
plot(g[[3]])

ggsave(plot=plot_grid(plotlist=g, nrow=3, byrow=FALSE), filename=sprintf('%s/TimeSeries_Graph_Parameters.pdf', dir$figdir),h=6, w=8.27)
       
