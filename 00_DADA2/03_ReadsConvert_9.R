############################################################################
####
#### R script for Fujita (2019)
####
#### reads conversion
#### 2018.4.18 Ushio created
#### 2019.10.31 Fujita modified
#### 2022. 2. 23. Toju
#### R 4.1.2
####
############################################################################

#setwd("/Users/toju/Dropbox/YAJIMA_Daii/Statistics_R")
setwd("/Users/toju/Documents/Eel_Yajima/Eel_culture_TimeSeries")

ran.seed <- 1234
set.seed(ran.seed)


## -- Loading Function and Library
library(AnalysisHelper) 
load.lib( c('seqinr', 'ggplot2'))

# -- Create directory to save
dir = make.dir('./00_DADA2/03_output')

# -- Load data table
# -- Load workspace and functions
taxa.print <- readRDS("Table/taxaPrint.rds")
seqtab <- readRDS("Table/seqtab_rmchimera.rds")



############################################################################
## -- Parameters

## -- nM of each STD DNA
stdMix <- c( 0.1, 0.05, 0.02, 0.01,0.005) # nM

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

# STD copies per uL of master mix
std.copy.n_in_mastermix <- (stdMix/dilutedSTD.inPCR)*(6.02*10^14/1000000) # Avogadro constant conersion from L to uL: 10^(23-9)

lm.coef.fun <- function(x) summary(lm(as.numeric(x) ~ std.copy.n_in_mastermix + 0))$coefficients[1]

############################################################################
####
#### F1. Collection of helper functions for DNA extration study
#### 2017.12.1 Ushio
#### Now published in https://github.com/ong8181/micDNA-beads
###

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
coef.summary <- apply(new.std.table, 1, lm.coef.fun)
new.seqtab <- as.data.frame(seqtab[,-n.std.seq])

############################################################################

# -- Convert sequence reads to copy numbers

coef <- apply(new.std.table, 1, lm.coef.fun)

# -- Coefficient check
all(coef > 0)

# -- Conversion
seqtab.conv <- (new.seqtab/coef) * ((lysis.Buffer.Volume+sample.Volume)/sample.Volume) * (total.Volume/template.DNA.volume)

taxa.wo.std <- taxa.print[-n.std.seq,]
dim(seqtab.conv)
dim(taxa.wo.std)

# -- Conversion of sample reads to calculated copy numbers
amp.factor <- 1/min(seqtab.conv[seqtab.conv != 0])
seqtab.conv.correct <- round(seqtab.conv*amp.factor) # Converted to integers

# -- Filt low quality sample (low quality means no-correlation with std.copy)

std.check <- apply(new.std.table, 1, function(x){  cor( x, as.matrix(std.copy.n_in_mastermix))  })
detectSTD <- apply(new.std.table, 1, function(x) any(x==0))
low.q.sample <- which(std.check < 0.8 | rowSums(new.seqtab) < 2000 | detectSTD) 
plot(rowSums(new.seqtab), std.check); abline(v=2000)
if(length(low.q.sample) > 0) {
	seqtab.conv.correct.filt <- seqtab.conv.correct[-low.q.sample, ]
} else {
	seqtab.conv.correct.filt <- seqtab.conv.correct
}


############################################################################

# --Excluding contamination sequence
seqtab.c <- seqtab.conv.correct.filt

# -- Excluding Chloroplast and Mitochondria
taxaBac <- taxa.print[taxa.print[,'Kingdom']%in%c('Bacteria','Archaea'),]
mito = which(taxaBac == "Mitochondria", arr.ind=T)
chlo = which(taxaBac == "Chloroplast", arr.ind=T)

if(length(c(mito, chlo)) > 0) taxaBac <- taxaBac[-c(mito[,1], chlo[,1]),]
taxaBac[is.na(taxaBac)] <- 'Unidentified'

# --Excluding standard DNA ASV
row.std <- grep('STD_', taxaBac[,'Phylum'])
taxa.print.no.std <- taxaBac[-row.std, ]
taxa.print.named <- as.data.frame(cbind(ID=sprintf("X_%04d", 1:nrow(taxa.print.no.std)), taxa.print.no.std))

# --Output ASV table
asv.ex.contami <- taxa.print.named
saveRDS(cbind(asv.ex.contami,rownames(asv.ex.contami)), "Table/Taxa_list_with_sequence.rds")

ASV.seq <- rownames(asv.ex.contami)
ASV.out <- cbind(asv.ex.contami, ASV.seq)
write.table(ASV.out, file=sprintf("%s/Taxa_list.pro.txt", dir$tabledir), sep="\t", quote=F, row.name=F)
saveRDS(ASV.out, "Table/Taxa_list.pro.rds")

rownames(asv.ex.contami) <- NULL
rownames(asv.ex.contami) <- asv.ex.contami[,1]
saveRDS(asv.ex.contami, "Table/Taxa_list.rds")

# -- Output Fasta file
taxa.print.rename3 <- as.data.frame(asv.ex.contami)
fasta.table  <- cbind(as.character(taxa.print.rename3 $ID), rownames(taxa.print.named ))
	
write.fasta(as.list(fasta.table[,2]), fasta.table[,1],'Table/16S.fasta')

## ------------------------------------------------------------------ ##
# -- Out the sequence reads matrix
seqtab.bacteria <- seqtab.c[, rownames(taxa.print.named)]
if(all(colnames(seqtab.bacteria)==rownames(taxa.print.named))) {
	    colnames(seqtab.bacteria)= taxa.print.named[,1]
	}else{
	    stop(cat('colnames(seqtab.bacteria ) and rownames(asv.ex.contami) is different'))
}

rownames(seqtab.bacteria) <- gsub(".fastq.gz", "", rownames(seqtab.bacteria))
saveRDS(seqtab.bacteria, "Table/seqtabBac.rds")

write.table(cbind(Sample.ID=rownames(seqtab.bacteria), seqtab.bacteria), file=sprintf("%s/seqtab.16SrRNAcount.pro.txt", dir$tabledir), sep="\t", quote=F, row.names=F)


## ------------------------------------------------------------------ ##

# Original read matrix

new.seqtab.out <- new.seqtab[, rownames(taxa.print.named)]
if(all(colnames(new.seqtab.out)==rownames(taxa.print.named))) {
	    colnames(new.seqtab.out)= taxa.print.named[,1]
	}else{
	    stop(cat('colnames(new.seqtab.out ) and rownames(asv.ex.contami) is different'))
}

rownames(new.seqtab.out) <- gsub(".fastq.gz", "", rownames(new.seqtab.out))

saveRDS(new.seqtab.out, "Table/seqtabBac_OriginalReads.rds")

write.table(cbind(Sample.ID=rownames(new.seqtab.out), new.seqtab.out), file=sprintf("%s/seqtab.txt", dir$tabledir), sep="\t", quote=F, row.names=F)



# -- Out the matrix to check std reads proportion
seqtab.check <- seqtab
colnames(seqtab.check) <- taxa.print[,'Phylum']
	
sum.std <- rowSums(seqtab.check[,grep('STD_',colnames(seqtab.check))])
sum.na <-  rowSums(seqtab.check[,which(is.na(colnames(seqtab.check)))])
sum.pro <- rowSums(seqtab.check[,-c(grep('STD_',colnames(seqtab.check)),which(is.na(colnames(seqtab.check))))])
sum.reads <- rowSums(seqtab.check)
seqtab.check2 <- cbind(rownames(seqtab.check),sum.reads,sum.std, sum.pro, sum.na, seqtab.check)
	
lf <- as.data.frame(rbind( cbind("Sum of STD DNA reads",sum.std), cbind("Sum of ASV reads",sum.pro),
		           cbind("Sum of NA sequence reads",sum.na),cbind("Sum of each sample reads",sum.reads)))
colnames(lf) = c("key","value")
gg <- ggplot(data=lf,aes(x=as.numeric(as.vector(value)))) +
		  geom_histogram(fill='white', color='indianred4',size=1)+
		  geom_hline(yintercept=0)+
		  facet_wrap(~key)+
		  theme_minimal(base_size=20)
		  
pdf(sprintf("%s/Reads.count.pdf",dir$figdir))
plot(gg)	
dev.off()

## ------------------------------------------------------------------ ##


############################################################################
