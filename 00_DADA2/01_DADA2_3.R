############################################################################
#### 2022. 2. 21. Fujita
#### 2022. 2. 23. Toju
####
#### R script for eel gut microbiome
#### Filtering to Clustering ASVs
####
#### R 4.1.2
#### 
############################################################################

#setwd("/Users/toju/Dropbox/YAJIMA_Daii/Statistics_R")
setwd("/Users/toju/Documents/Eel_Yajima/Eel_culture_TimeSeries")

ran.seed <- 1234
set.seed(ran.seed)

## -- Loading Function and Library
library(AnalysisHelper) # devtools::install_github("hiroakif93/R-functions/AnalysisHelper")
load.lib(c('dada2', 'seqinr', 'stringr', 'ShortRead', 'Biostrings', 'ggplot2'))

# -- Create directory to save
dir.create("Table")
dir = make.dir('./00_DADA2/01_output')

###########################################################################

# -- Defined the name of forword and reverse sequence
pattern.f <- "fastq.gz"

# -- Difined the directory containing the fastq files
rootpath <- 'fastq' 
rootfiles <- list.files(rootpath)

reference <- 'reference/silva_nr99_v138_train_set_append_stdDNA.fa'
############################################################################

for(f in rootfiles) {
   
   # -- Create directory to store filtered fastq files
   runpath <- sprintf('%s/%s', dir$rdsdir, f)

   filtFs <- file.path(sprintf('%s/filtered', runpath))
   fastqFs <- sort(list.files(paste(rootpath, f, sep='/'), pattern = pattern.f, full.names = TRUE))
   
   ## ============================================= ##
   ## -- Filtered fastq directory
   
   out <- filterAndTrim(fwd=fastqFs, filt=filtFs,
                        truncLen=c(240), maxEE=2, truncQ=11, maxN=0, rm.phix=TRUE,
                        compress=TRUE, verbose=TRUE, multithread=TRUE)
   
   more0 <- out[which(out[,2]>0),]
   
   filtFs.ex0 <- fastqFs[basename(fastqFs)%in%rownames(more0)]
   filtFs.ex0 <- paste(filtFs, basename(filtFs.ex0), sep='/')
   
   sample.names <- basename(filtFs.ex0)
   names(filtFs.ex0) <- sample.names
   
   ## ============================================= ##
   # --  Learn error rates
   errF <- learnErrors(filtFs.ex0, nbases=1e10, multithread=TRUE, verbose=TRUE)
   
   pdf(sprintf('%s/plotErrors_%s.pdf', dir$figdir, f))
   print(plotErrors(errF, nominalQ=TRUE))
   dev.off()
   
   ## ============================================= ##
   # --  Sample inference and merger of paired-end reads
   
   derepF <- derepFastq(filtFs.ex0)
   ddF <- dada(derepF, err=errF, multithread=TRUE)
   
   st.all <- makeSequenceTable(ddF)
   
   saveRDS(st.all, sprintf("%s/stall_%s.rds", dir$rdsdir, f))
   save.image(sprintf("%s/Rdata_%s.RData", dir$rdsdir, f))
}

############################################################################


# -- Merge multiple runs 
rds.file <- list.files(dir$rdsdir, full.names = TRUE)
rds.file <- rds.file[grep(".rds",rds.file)]
rdslist <- lapply(rds.file, readRDS)

# stmerge <- mergeSequenceTables(tables=rdslist, repeats = "sum")

# -- Remove chimeras
#seqtab <- removeBimeraDenovo(stmerge, method="consensus", multithread=TRUE)
seqtab <- removeBimeraDenovo(rdslist, method="consensus", multithread=TRUE)
saveRDS(seqtab, "Table/seqtab_rmchimera.rds")


# -- Assign taxonomy
taxa <- assignTaxonomy(seqtab, reference, multithread=FALSE)
taxa[is.na(taxa)] <- 'Unidentified'
saveRDS(taxa, "Table/taxaprint.rds")

write.table(cbind(Sequence=rownames(taxa), taxa), file=sprintf("%s/taxa.txt", dir$tabledir), sep='\t', quote=F, row.name=F)

############################################################################
