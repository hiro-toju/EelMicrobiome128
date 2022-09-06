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

library(ggplot2)
library(ggsci)
library(reshape2)
library(RColorBrewer)
library(cowplot)
source('functions/functions.R')
load.lib( c('vegan', 'ggplot2', 'tidyr', 'cowplot','RColorBrewer', 'scales', 'ggsnippets'))

# -- Create directory to save
dir <- make.dir('./01_Timeseries_Graph/03_output')

info <- read.table('Table/Sample.List_1.txt', header=T)
data  <- read.table('Table/seqtab.16SrRNAcount.pro.txt', header=T)
param <- as.matrix(read.table("Table/Parameters_3.txt", row.name=1, header=T, na.strings = "NA"))


ASV_Aeromonas <- c("X_0107", "X_2013", "X_3052", "X_3402", "X_3437", "X_3905", "X_5377", "X_6414", "X_8079")
ASV_Edwardsiella <- c("X_0027", "X_5825")
ASV_Flavobacterium <- c("X_0276", "X_6694", "X_7916")
ASV_Mycobacterium <- c("X_1631", "X_3216", "X_4368")
ASV_Pseudomonas <- c("X_5280")
ASV_Renibacterium <- c("X_6511")

mat_Aeromonas <- matrix(NA, nrow=nrow(data), ncol=length(ASV_Aeromonas))
mat_Edwardsiella <- matrix(NA, nrow=nrow(data), ncol=length(ASV_Edwardsiella))
mat_Flavobacterium <- matrix(NA, nrow=nrow(data), ncol=length(ASV_Flavobacterium))
mat_Mycobacterium <- matrix(NA, nrow=nrow(data), ncol=length(ASV_Mycobacterium))
mat_Pseudomonas <- matrix(NA, nrow=nrow(data), ncol=length(ASV_Pseudomonas))
mat_Renibacterium <- matrix(NA, nrow=nrow(data), ncol=length(ASV_Renibacterium))

for (i in 1:length(ASV_Aeromonas)){mat_Aeromonas[, i] <- data[, colnames(data)==ASV_Aeromonas[i]]}
for (i in 1:length(ASV_Edwardsiella)){mat_Edwardsiella[, i] <- data[, colnames(data)==ASV_Edwardsiella[i]]}
for (i in 1:length(ASV_Flavobacterium)){mat_Flavobacterium[, i] <- data[, colnames(data)==ASV_Flavobacterium[i]]}
for (i in 1:length(ASV_Mycobacterium)){mat_Mycobacterium[, i] <- data[, colnames(data)==ASV_Mycobacterium[i]]}
for (i in 1:length(ASV_Pseudomonas)){mat_Pseudomonas[, i] <- data[, colnames(data)==ASV_Pseudomonas[i]]}
for (i in 1:length(ASV_Renibacterium)){mat_Renibacterium[, i] <- data[, colnames(data)==ASV_Renibacterium[i]]}

Aeromonas <- rowSums(mat_Aeromonas)
Edwardsiella <- rowSums(mat_Edwardsiella)
Flavobacterium <- rowSums(mat_Flavobacterium)
Mycobacterium <- rowSums(mat_Mycobacterium)
Pseudomonas <- rowSums(mat_Pseudomonas)
Renibacterium <- rowSums(mat_Renibacterium)

Disease <- cbind(Aeromonas, Edwardsiella, Flavobacterium, Mycobacterium, Pseudomonas, Renibacterium)
Sample.ID <- data$Sample.ID
Disease2 <- data.frame(cbind(Sample.ID, Disease))
All <- rowSums(Disease)

Disease3 <- merge(info, Disease2,  by='Sample.ID', all = FALSE, sort = TRUE)
rownames(Disease3) <- Disease3$Sample.ID


col.list <- rev(brewer.pal(6, "Pastel1"))


############################################################################

glist1 <- c()

for(i in unique(Disease3$Tank)){ 
    
    ## ============================================================ ##
    ## -- Extracting one treatment matrix
    
    label <- Disease3[Disease3$Tank==i, 1:3]
    comm <- Disease3[Disease3$Tank==i, 4:ncol(Disease3)]
    
    ## ======================================= ##
    ## -- Visualizing community assembly

        lf <- gather(cbind(label, comm), 
                     key, value, -c(1:(ncol(label))))
        
        lf$key <- factor(lf$key, levels=rev(colnames(comm)) )
         
        g1 <- ggplot(lf)+
            geom_area(aes(x=as.numeric(Day), y=as.numeric(value), fill=key), color='grey30', stat='identity', show.legend=FALSE, size=0.1)+
          #  facet_wrap(~replicate.id,scales='free_y',ncol=1, strip.position='left')+
            scale_fill_manual( values=col.list ) + 
            theme_bw(base_size=15)+
            labs(x= "Day", subtitle=sprintf("Tank %s", i), y= "Absolute abundance")+
            guides(fill=guide_legend(title='',ncol=3, reverse=TRUE))+
            theme(strip.placement = 'outside')+
            scale_x_continuous(expand=c(0,0))+
            scale_y_continuous(expand=c(0,0), position = "left") + theme(text = element_text(size = 11))
        glist1[[i]] <- g1
       
}    

ggsave(plot=plot_grid(plotlist=glist1, nrow=5, byrow=FALSE),
       filename=sprintf("%s/Disease_TimeSeries.pdf", dir$figdir), h=11.69, w=8.27)
       
############################################################################

# Color Caption

ggsave(plot=ggplot(lf)+
            geom_area(aes(x=as.numeric(Day), y=as.numeric(value), fill=key), color='grey30', stat='identity', show.legend=TRUE, size=0.1)+
            scale_fill_manual( values=col.list ) + theme(text = element_text(size = 14)),
            filename=sprintf('%s/Disease_TimeSeries_Caption.pdf', dir$figdir), h=5, w=8.27)
       
############################################################################

