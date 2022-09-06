############################################################################
#### 2022. 2. 21. Fujita
#### 2022. 2. 23. Toju
####
#### R script for eel gut microbiome
#### STD DNA check
####
#### R 4.1.2
#### 
############################################################################

#setwd("/Users/toju/Dropbox/YAJIMA_Daii/Statistics_R")
setwd("/Users/toju/Documents/Eel_Yajima/Eel_culture_TimeSeries")

ran.seed <- 1234
set.seed(ran.seed)

## -- Loading Function and Library
library(AnalysisHelper) 
load.lib(c("ggplot2", "cowplot", "tidyr", "reshape2", "stringr"))

# -- Create directory to save
dir = make.dir('./00_DADA2/02_output')

# -- Load workspace and functions
taxa.print <- readRDS("Table/taxaPrint.rds")
rownames(taxa.print) <- NULL
seqtab <- readRDS("Table/seqtab_rmchimera.rds")

############################################################################
## -- Parameters

## -- nM of each STD DNA
stdMix <- c( 0.1, 0.05, 0.02, 0.01,0.005)

## -- Dilution rate of STD DNA
STD.dilutionRate <- 20000

## -- DNA extraction condition
lysis.Buffer.Volume <- 500
sample.Volume <- 250

## -- PCR condition
total.Volume <- 8
template.DNA.volume <- 2
std.DNA.volume <- 0.32

############################################################################
## -- Calculate DNA copy concentration in sample (copy/micro litter)

## -- STD DNA concentration in PCR master mix
dilutedSTD.inPCR <- STD.dilutionRate/(std.DNA.volume/total.Volume)

## -- Convert to DNA copy number from nM
std.copy.n_in_mastermix <- (stdMix/dilutedSTD.inPCR)*(6.02*10^13/1000000)
std.copy.n <- std.copy.n_in_mastermix* (lysis.Buffer.Volume+sample.Volume)/sample.Volume

## -- Calibration curve function
sampleRatio = total.Volume/template.DNA.volume
lm.coef.fun <- function(x) summary(lm(as.numeric(x) ~ std.copy.n + 0))$coefficients[1] * sampleRatio

############################################################################
####
#### F1. Collection of helper functions for DNA extration study
#### 2017.12.1 Ushio
#### Now published in https://github.com/ong8181/micDNA-beads
###

# ggplot function 1
PlotStyle <-  function(ggobject){
    return(ggobject + theme_bw() + theme(axis.text.x = element_text(angle=0),
                                         panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank(),
                                         axis.text = element_text(size=12),
                                         axis.title = element_text(size=12),
                                         panel.background=element_rect(colour="black", fill=NA, size=0.8)))
}

# ggplot function 2
PlotStyle2 <- function(ggobject){
    return(ggobject + theme(axis.text.x = element_text(angle = 90, hjust = 1),
                            axis.title.x = element_blank(),
                            legend.position = "none") +
               geom_jitter(shape = 16, size = 2, alpha = 0.8, width = 0.1, height = 0))
}

# Merge standard DNA sequences
MergeSTD <- function(std.i, std.data = std.table){
    index.std <- which(match(colnames(std.table), std.i) == 1)
    if(length(index.std) > 1){
        std.tmp <- rowSums(std.table[,index.std])
    }else{
        std.tmp <- std.table[,index.std]
    }
    return(std.tmp)
}

# To check correlation
adj.r.fun <- function(x) summary(lm(as.numeric(x) ~ std.copy.n + 0))$adj.r.squared

############################################################################

# -- Extract standard sequeces
detected.std.name <- unique(taxa.print[which(substr(taxa.print[,"Phylum"], 1, 7) == "STD_pro"), "Phylum"])

n.std.seq <- which(substr(taxa.print[,"Phylum"], 1, 7) == "STD_pro")
std.table <- seqtab[,n.std.seq]
std.taxa <- taxa.print[n.std.seq, "Phylum"]

# --  STD reads - copy number relationship
# --  Rename colnames
colnames(std.table) <- std.taxa
# --  Merge the same colnames
new.std.table <- data.frame(std_rank1 = MergeSTD(detected.std.name[1], std.data = std.table),
                            std_rank2 = MergeSTD(detected.std.name[2], std.data = std.table),
                            std_rank3 = MergeSTD(detected.std.name[3], std.data = std.table),
                            std_rank4 = MergeSTD(detected.std.name[4], std.data = std.table),
                            std_rank5 = MergeSTD(detected.std.name[5], std.data = std.table))

# -- Interpolate 0 into missing value
new.std.table[is.na(new.std.table)] <- 0

# --  Linear regression
r2.summary <- apply(new.std.table, 1, adj.r.fun)
coef.summary <- apply(new.std.table, 1, lm.coef.fun)
new.seqtab <- as.data.frame(seqtab[,-n.std.seq]) # Make seq table without standard DNA

# --  Visualize regression results
# --  1. R2 value distribution
g1 <- ggplot(data.frame(values = r2.summary), aes(x = values))
g1 <- g1 + geom_histogram() + geom_hline(yintercept = 0, linetype = 2) #+ xlim(0.9,1) + ylim(0,13)
g1 <- g1 + xlab(expression(paste("R"^{2}, " values"))) + ylab("Count") #+ scale_x_continuous(breaks=seq(0.8,1,by=0.05))
g1 <- PlotStyle(g1)

# --  2. Slope distribution
g2 <- ggplot(data.frame(values = coef.summary), aes(x = values))
g2 <- g2 + geom_histogram() + geom_hline(yintercept = 0, linetype = 2)
g2 <- g2 + xlab("Regression slope") + ylab("Count")
g2 <- PlotStyle(g2)

# --  3. Regression examples
max.slope <- as.numeric(c(new.std.table[which.max(coef.summary),], 0))
med.slope <- as.numeric(c(new.std.table[which.min(abs(coef.summary - median(coef.summary))),], 0))
min.slope <- as.numeric(c(new.std.table[which.min(coef.summary),], 0))
slope.summary <- melt(data.frame(copy =c(std.copy.n, 0),
                                 max_slope = max.slope,
                                 med_slope = med.slope,
                                 min_slope = min.slope), id.vars = "copy")
g3 <- ggplot(slope.summary, aes(x = copy, y = value, group = variable, colour = variable))
g3 <- g3 + geom_point(size = 2) + scale_color_manual(name = "Regression slope", values = c("red3", "darkred", "black"))
g3 <- g3 + geom_smooth(method = "lm", size = 0.5, se = F)+theme_bw(base_size=26)
	  
g3 <- g3 + xlab(expression(paste("Standard DNA copies (", mu, l^{-1}, ")"))) + ylab("Standard DNA reads")
g3 <- PlotStyle(g3) + theme(legend.position = c(0.25, 0.75), legend.text = element_text(size = 7), legend.title = element_text(size = 8))

# --  4. Read visualization
new.seqtab$sample <- rownames(new.seqtab)

seqtab.plot <- gather(new.seqtab, key, value, -ncol(new.seqtab))
seqtab.plot <- seqtab.plot[seqtab.plot$value > 0,]
seqtab.plot$value2 <- seqtab.plot$value + 0.5
std.table.plot <- data.frame(sample = rownames(new.seqtab),
                             max = apply(new.std.table, 1, max),
                             min = apply(new.std.table, 1, min) + 0.5) # data.frame for overwrite standard DNA reads

seqtab.plot$group <- sapply(strsplit(seqtab.plot$sample, "_"), "[", 3)
# --  Box plot + jitter plot version
g4 <- ggplot(seqtab.plot, aes(x = sample, y = value2))
g4 <- g4 + geom_boxplot(shape = 16, alpha = 0.5)
g4 <- g4 + geom_jitter(shape = 16, size = 1.5, alpha = 0.5, position = position_jitter(0.2))
g4 <- g4 + scale_y_log10() + ylab("Sequence reads") + xlab("Sample ID")
g4 <- PlotStyle(g4) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
# Overwrite standard DNA reads
g4 <- g4 + geom_segment(data = std.table.plot, aes(x = sample, xend = sample, y = min, yend = max),
                        colour = "red3", size = 1.5, alpha = 0.5)
g5 <- g4+facet_wrap(~group)

# -- Summarize visualization
top.row <- plot_grid(g1, g2, g3, ncol = 3, align = "h", labels = c("A", "B", "C"), rel_widths = c(1,1,1.3))
Fig.std <- plot_grid(top.row, g4, ncol = 1, labels = c("", "D"))

############################################################################

ggsave(plot= Fig.std,
		filename=sprintf("%s/STD_Check.pdf", dir$figdir),heigh=10,width=15)


############################################################################

