
########################################################################
## -- DADA2
## -- Running on R

## == Setting parameter
minLen=200
minQ=10
maxEE=2
truncR=0

inputdir="02_Cutadaptor_fastaq"

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
piplineID="03_FilterTrimming"

########################################################################


# -- Load library
library(dada2) 

if( length(grep('fastq.gz$', list.files(inputdir))) > 0 ){
  path = inputdir
}else{
  stop("The path directory is missing\n")
}

########################################################################

## -- Set random seeds (for reproduction)
ran.seed <- 1234
set.seed(ran.seed)

## -- Infer Sequence Variants
fnFs <- sort(list.files(path, pattern="fastq.gz", full.names = TRUE))
sample.names <- apply( sapply(strsplit(basename(fnFs), "_"), `[`, c(1,2)), 2, paste, collapse="_")

## -- Path to output directory 
filtFs <- gsub(inputdir, sprintf("%s_fastaFiles",piplineID),fnFs)

## -- Filtering process
filterAndTrim(fnFs, filtFs, maxN = 0, 
                       maxEE = maxEE, minLen = minLen, truncLen=truncR, minQ=minQ, truncQ=11,
                       rm.phix = TRUE, compress = TRUE, multithread = TRUE) 


sink("log_and_script/Version.txt", append=TRUE)
sessionInfo()
sink()

