
########################################################################
## -- DADA2
## -- Running on R

## == Setting parameter
inputdir <- "03_FilterTrimming_fastaFiles" 
outputdir <- "04_Denoising"
minident=1
vsearchpath="/home/toju/miniconda3/bin/vsearch"

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
piplineID="04_Denoising"

start <- Sys.time()
print(start); cat("\n")
########################################################################


library(dada2)
library(seqinr)

dir.create(outputdir)


takeSamplenames <- function(x, sep=''){#x=list.files(inputdir)
  
  split <- do.call(rbind, strsplit(x, sep) )
  check.unique <- sapply(apply(split, 2, table), length)
  uniquename <- which(check.unique>1)
  if(length(uniquename)==1  ) {
    samplename <- split[,uniquename]
  }else{
    samplename <- apply(split[,uniquename], 1, paste, collapse="_")
  }
  return(samplename)
}
########################################################################

# Learn forward error rates
cat( sprintf("Start learnerror from %s...\n", Sys.time()))
errF <- invisible( learnErrors(inputdir, nbases=1e8, multithread=TRUE) )

pdf(sprintf('%s/plotErrors.pdf', outputdir))
print(plotErrors(errF, nominalQ=TRUE))
dev.off()

# Infer sequence variants
derepFs <- derepFastq(inputdir, verbose = TRUE)

# Name the derep-class objects by the sample names
samplenames <- takeSamplenames(list.files(inputdir), sep='__')


names(derepFs) <- samplenames
dadaFs <- invisible( dada(derepFs, err = errF, multithread = TRUE) )

# Construct sequence table and write to disk
st.all <- makeSequenceTable(dadaFs)
saveRDS(st.all, sprintf('%s/stall_no_rmchimera.rds', outputdir))

## =============================================================== ##
# Remove chimeras
seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)

seq.mat <- cbind(colnames(seqtab),sprintf('X_%s', formatC(1:ncol(seqtab), width = nchar(ncol(seqtab)), flag = "0"))) 
write.fasta(as.list(seq.mat[,1]), seq.mat[,2], sprintf("%s/nonchim_seq.fasta", outputdir) )

seqtab2 <- seqtab
colnames(seqtab2) <- seq.mat[,2]
saveRDS(seqtab2, sprintf('%s/seqtab_rmChimera.rds', outputdir))
write.csv(cbind(sample=rownames(seqtab2), seqtab2), sprintf('%s/seqtab_rmChimera.csv', outputdir), row.names=FALSE)

##################################################################################################

system2(command = vsearchpath, 
        args = c("--cluster_fast $input", sprintf("%s/nonchim_seq.fasta", outputdir),
                 "--id", minident,
                 "--mothur_shared_out", sprintf("%s/ASV_OTU_corestab_%s.txt", outputdir, minident),
                 "--centroids", sprintf("%s/OTUseq_%s.fasta", outputdir, minident),
                 "--msaout",  sprintf("%s/seqAlign_%s.txt", outputdir, minident) ) )

otu <- read.table(sprintf("%s/ASV_OTU_corestab_%s.txt", outputdir, minident),
                            header=TRUE, row.names=2)[,-c(1:2)]

## ========================================== ##

otutab <- matrix(0, ncol=ncol(otu), nrow=nrow(seqtab2),
                 dimnames=list(rownames(seqtab2), colnames(otu)))

for(i in 1:ncol(otu)){
  
  if( sum(otu[,i])>1 ){
    memberSeq <- rownames(otu)[which(otu[,i]>0)]
    otutab[,i] <- rowSums(seqtab2[,which(colnames(seqtab2) %in% memberSeq) ])
  }else{
    centroidSeq <- colnames(otu)[i]
    otutab[,i] <- seqtab2[, which(colnames(seqtab2) == centroidSeq) ]
  }
  
}
  
saveRDS(otutab, sprintf("%s/seqOTUtab.rds", outputdir))
write.csv(cbind(sample=rownames(otutab), otutab), sprintf('%s/seqOTUtab.csv', outputdir), row.names=FALSE)

##################################################################################################

finish <- Sys.time()
print(finish-start); cat("\n")
##################################################################################################
sink("log_and_script/Version.txt", append=TRUE)
sessionInfo()
sink()

