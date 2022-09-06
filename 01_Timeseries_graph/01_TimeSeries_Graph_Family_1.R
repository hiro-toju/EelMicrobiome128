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

## -- Loading Function and Library
source('functions/functions.R')
load.lib( c('vegan', 'ggplot2', 'tidyr', 'cowplot','RColorBrewer', 'scales', 'ggsnippets'))

# -- Create directory to save
dir <- make.dir('./01_Timeseries_Graph/01_output')

# -- Load data table

info <- read.table('Table/Sample.List_1.txt', header=T)
data  <- read.table('Table/seqtab.16SrRNAcount.pro.txt', header=T)

merged.mat <- merge(info, data, by='Sample.ID', all = FALSE, sort = TRUE)
rownames(merged.mat) <- merged.mat$Sample.ID
sml  <- merged.mat[, 1:3]
comm <- merged.mat[, 4:ncol(merged.mat)]

#taxa <- read.table('Table/Taxa_list.pro_delspace.txt', row.name=1, header=T)
taxa <- readRDS('Table/Taxa_list.rds')

l='Family'
comm.new <- t(Taxa.mat(comm[rownames(sml),], taxa, l))
Unidentified <- comm.new[, colnames(comm.new)=="Unidentified"]
comm.exU <- comm.new[, colnames(comm.new)!="Unidentified"]

comm.exU.ordered <- comm.exU[, order(colSums(comm.exU), decreasing=TRUE)]
Others <- rowSums(comm.exU.ordered[, 31:ncol(comm.exU.ordered)])
comm.top30.OU <- cbind(comm.exU.ordered[, 1:30], Others, Unidentified)


# Setting colors

col1 <- brewer.pal(7, "Set1")
col2 <- brewer.pal(7, "Set2")
col3 <- brewer.pal(8, "Set3")
col4 <- brewer.pal(8, "Pastel1")

col.list <- rev(c(col2, col3, col1, col4, "grey90", "grey80"))


############################################################################

glist1 <- c()

for(i in unique(sml$Tank)){ #i=names(dlist)[2]
    
    ## ============================================================ ##
    ## -- Extracting one treatment matrix
    
    smlsub <- sml[sml$Tank==i, ]
    ts <- comm.top30.OU[sml$Tank==i, ]
    rel.ts <- ts/rowSums(ts)
    
    ## ======================================= ##
    ## -- Visualizing community assembly

         
        
        lf <- gather(cbind(smlsub, ts), 
                     key, value, -c(1:(ncol(smlsub))))
        
        lf$key <- factor(lf$key, levels=rev(colnames(comm.top30.OU)) )
         
        g1 <- ggplot(lf)+
            geom_area(aes(x=as.numeric(Day), y=value, fill=key), color='grey30', stat='identity', show.legend=FALSE, size=0.1)+
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
       filename=sprintf('%s/TimeSeriesGraph_Absolute_Family.pdf',
       dir$figdir), h=11.69, w=8.27)
       

############################################################################


############################################################################

glist1 <- c()

for(i in unique(sml$Tank)){ #i=names(dlist)[2]
    
    ## ============================================================ ##
    ## -- Extracting one treatment matrix
    
    smlsub <- sml[sml$Tank==i, ]
    ts <- comm.top30.OU[sml$Tank==i, ]
    rel.ts <- ts/rowSums(ts)
    
    ## ======================================= ##
    ## -- Visualizing community assembly

         
        
        lf <- gather(cbind(smlsub, ts), 
                     key, value, -c(1:(ncol(smlsub))))
        
        lf$key <- factor(lf$key, levels=rev(colnames(comm.top30.OU)) )

        g1 <- ggplot(lf)+
            geom_area(aes(x=as.numeric(Day), y=value, fill=key), position="fill", color='grey30', stat='identity', show.legend=FALSE, size=0.1)+
          #  facet_wrap(~replicate.id,scales='free_y',ncol=1, strip.position='left')+
            scale_fill_manual( values=col.list ) + 
            theme_bw(base_size=15)+
            labs(x= "Day", subtitle=sprintf("Tank %s", i), y= "Relative abundance")+
            guides(fill=guide_legend(title='',ncol=3, reverse=TRUE))+
            theme(strip.placement = 'outside')+
            scale_x_continuous(expand=c(0,0))+
            scale_y_continuous(expand=c(0,0), position = "left") + theme(text = element_text(size = 11))
        glist1[[i]] <- g1
       
}    

ggsave(plot=plot_grid(plotlist=glist1, nrow=5, byrow=FALSE),
       filename=sprintf('%s/TimeSeriesGraph_Relative_Family.pdf',
       dir$figdir), h=11.69, w=8.27)
       
ggsave(plot=ggplot(lf)+
            geom_area(aes(x=as.numeric(Day), y=value, fill=key), position="fill", color='grey30', stat='identity', show.legend=TRUE, size=0.1)+
            scale_fill_manual( values=col.list ) + theme(text = element_text(size = 14)),
            filename=sprintf('%s/TimeSeriesGraph_Caption_Family.pdf',
       dir$figdir), h=5, w=8.27)
       
############################################################################

