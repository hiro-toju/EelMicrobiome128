#!/bin/bash
####################################################################
## 											
## ---------- 	  	 Filtering Check 			    ------------- ##
##
## 											2022. 02. 28. by Fujita
####################################################################

## Input directory containing not/filtered fastq
beforedir=02_Cutadaptor_fastaq
afterdir=03_FilterTrimming_fastaFiles
seqtabPATH=04_Denoising/seqtab_rmChimera.rds

thread=32
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##

piplineID="06_CheckFiltering"

####################################################################

## ============= Remove and Make directories ==================== ##

## --Remove older directory
if ls ${piplineID}_QCreport* >/dev/null 2>&1 ; then
	rm -r ${piplineID}*
fi

## ------------------------------------------------------------- ##
## -- Making directory to save results
mkdir -p ${piplineID}
mkdir -p ${piplineID}

####################################################################

## -- Read count
seqkit stat $beforedir/* -j $thread -a -b > $piplineID/StaticsBefore.tsv
seqkit stat $afterdir/* -j $thread -a -b > $piplineID/StaticsAfter.tsv


####################################################################
## -- Visualize 

cat <<RRR >  log_and_script/$piplineID.R
####################################################################
## 											
## ----------      Filtering Check Visualization    ------------- ##
##
## 										     2022. 02. 28. by Fujita
####################################################################
beforedir="02_Cutadaptor_fastaq"
afterdir="03_FilterTrimming_fastaFiles"
seqtabPATH="${seqtabPATH}"
diectory <- c(beforedir, afterdir)

path <-"${piplineID}"

####################################################################
RRR

cat <<'RRR' >> log_and_script/$piplineID.R
####################################################################

## -- Load packages and functions
library(AnalysisHelper)
load.lib(c("ShortRead", "ggplot2", "cowplot", "gridExtra", "tidyr","RreadQC"))

## -- Load data table
beforeread <- read.table(sprintf("%s/StaticsBefore.tsv", path), header=1, row.names=1)
rownames(beforeread) <- takeSamplenames1(x=rownames(beforeread), sep="__")
afterread <-  read.table(sprintf("%s/StaticsAfter.tsv", path), header=1, row.names=1)
rownames(afterread) <- takeSamplenames1(x=rownames(beforeread), sep="_filt")

## -- Load data table
fastqB <- readFastq(beforedir)
fastqA <- readFastq(afterdir)

####################################################################
## -- Visualize total sequence read

## -- Take common sample name
samplename <- intersect(rownames(beforeread), rownames(afterread) )

beforetmp <- beforeread[samplename, ]
aftertmp <- afterread[samplename, ]

## -- For ggplot data.frame
df <- data.frame(sample=samplename, 
                 read.in=as.numeric( gsub(",", "", beforetmp$num_seqs)),
                 read.out=as.numeric( gsub(",", "", aftertmp$num_seqs)) )
dfsort <- df
dfsort$sample <- factor(dfsort$sample, levels=dfsort$sample[order(dfsort$read.out)] )

## -- Make & save figure
nonsort <- ggreadFun(df)
sort <- ggreadFun(dfsort)

ggsave(plot=plot_grid(nonsort, sort, ncol=1),
       filename=sprintf("%s/Sequence_read_yield.pdf", path),
       w=15, h=10)

####################################################################

## -- Claculate quality indices
QstatB <- staticsQ(dir=beforedir)
QstatA <- staticsQ(dir=afterdir)

Qlist <- lapply( list(QstatB, QstatA), function(x){ #x=QstatB
  
  ## |||||||||||||||||||||||||||||||||||||||||||| ##
  ## -- Visualize total Quality score per 96 samples
  
  x$facet<-1
  samples <- unique(x$sample)
  sep=96
  ## |||||||||||||||||||||||||||||||||||||||||||| ##
  if( length(samples) > sep ){
    
    split <- round( length(samples)/sep )
    block <- cbind( c(1, (1:(split-1)*sep) +1),  (1:split*sep))
    
    for(i in 1:nrow(block)){ 
      sep <- samples[block[i,1]:block[i,2]]
      x[x$sample%in%sep, "facet"] <- i
    }
  }
  
  return(x)
  ## |||||||||||||||||||||||||||||||||||||||||||| ##

})

Qmerge <- rbind( cbind(filt="Before filtering", Qlist[[1]]),
                 cbind(filt="After filtering", Qlist[[2]]))
Qmerge$filt <- factor(Qmerge$filt, levels=c("Before filtering", "After filtering"))
## =============================================================== ##


glist <- c()
for(i in unique(Qmerge$facet) ){
  
  Qsub <- Qmerge[Qmerge$facet == i, ]
  
  ggQmean <- ggQFun(Qsub, xaxis="sample", yaxis="Mean")+
             scale_y_sqrt()+
             facet_wrap(~filt, ncol=1)+
             labs(title="Quality score")
  ggmaxEE <- ggQFun(Qsub, xaxis="sample", yaxis="maxEE")+
             scale_y_log10(label=scales::comma)+
             facet_wrap(~filt, ncol=1)+
             labs(title="Maximum expected error")
  ggLen <- ggQFun(Qsub, xaxis="sample", yaxis="length")+
           facet_wrap(~filt, ncol=1)+
           labs(title="Sequence read length")
  
  glist[[i]] <- plot_grid(ggQmean, ggmaxEE, ggLen, nrow=1)
}

ggsave(plot=marrangeGrob(glist, nrow=1, ncol=1),
       filename=sprintf("%s/Quality_check.pdf", path),
       w=15, h=10)
## ============================================================= ##
## -- Quality score distribution

gdist <- plot_grid(ggQdistFun(Qmerge, "Mean")+facet_wrap(~filt)+
                     scale_x_sqrt(expand=c(0,0,0.05,0), labels = scales::comma),
                   ggQdistFun(Qmerge, "maxEE")+facet_wrap(~filt)+
                     scale_x_log10(expand=c(0,0,0.05,0), labels = scales::comma),
                   ggQdistFun(Qmerge, "length")+facet_wrap(~filt)+
                     scale_x_log10(expand=c(0,0,0.05,0), labels = scales::comma), ncol=1)
ggsave(plot=gdist,
       filename=sprintf("%s/Sequence_read_distribution.pdf", path),
       w=8, h=10)

## ============================================================= ##

sampleMeanB <- aggregate(QstatB[,4:6], by=list(sample=QstatB$sample), mean)
rownames(sampleMeanB) <- sampleMeanB[,1]
sampleMeanA <- aggregate(QstatA[,4:6], by=list(sample=QstatA$sample), mean)
rownames(sampleMeanA) <- sampleMeanA[,1]

staticsResult <- cbind(df, before=sampleMeanB[df[,1],-1], after=sampleMeanA[df[,1],-1])

colnames(staticsResult) <- c("sample", 
                             "Total read before filtering",
                             "Total read after filtering",
                             paste(c("maxEE", "Mean Quality score","Sequence length"), 
                                   c("before filtering")),
                             paste(c("maxEE", "Mean Quality score","Sequence length"), 
                                   c("after filtering")))
total <- c("Total/Mean", colSums(staticsResult[,2:3]), colMeans(staticsResult[,4:ncol(staticsResult)]))

if(!is.null( seqtabPATH)){
  
  library(vegan)
  seqtab <- readRDS(seqtabPATH)
  slope=rareslope(seqtab, sample=c(500, 1000, 2000, 4000))
  for(i in c(500, 1000, 2000, 4000)){ 
    row <- which(rowSums(subtab)<i)  
    if(length(row)>0){
      slope[names(row), sprintf("N%s",i)] <- NA
    }
  }
  
  lf <- gather( data.frame(rownames(slope), slope), key, value, -1)
  lf$key <- factor(lf$key, levels=paste("N",c(500, 1000, 2000, 4000), sep=""))
  g <- ggplot(lf, aes(x=value))+
        geom_histogram(fill="royalblue4", color="honeydew1")+
        facet_wrap(~key)+
        theme_bw()+
        scale_y_continuous(expand=c(0,0,0.01,0))+
        scale_x_continuous(expand=c(0,0,0.01,0))
    
  ggsave(plot=g,
         filename=sprintf("%s/Rareslope_distribution.pdf", savedir),
         w=6, h=4)
  subtab <- seqtab[rowSums(seqtab)>100,]
  
  pdf(sprintf("%s/RaresCurve.pdf", savedir))
  quickRareCurve(subtab, sample=100, label = FALSE,
                 col=rainbow(ncol(subtab)))
  dev.off()
}


write.csv(rbind(total, staticsResult), sprintf("%s/summary_quality.csv", path),
          row.names=FALSE)
write.csv(slope, sprintf("%s/RareCurve.csv", path))
         
####################################################################          
RRR

Rscript log_and_script/$piplineID.R 2>&1 | tee log_and_script/log${piplineID}.txt

