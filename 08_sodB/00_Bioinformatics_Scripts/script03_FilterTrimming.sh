#!/bin/bash

####################################################################
## 											
## ---------- 	   Quality control with DADA2		------------- ##
##
## 											2022. 02. 28. by Fujita
####################################################################
## Input directory
fastadir=02_Cutadaptor_fastaq

## -- Options
minLen=200
minQ=10
PercLQ=0.2
trancWindow=5
maxEE=2
truncR=0
avgqual=15

options=""
thread=32

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
piplineID="03_FilterTrimming"

####################################################################

## ============= Remove and Make directories ==================== ##

## --Remove older directory
if ls ${piplineID}_QCreport* >/dev/null 2>&1 ; then
	rm -r ${piplineID}_QCreport*
fi

if [ -d ${piplineID}_fastaFiles ]; then
	rm -r ${piplineID}_fastaFiles
fi

## ------------------------------------------------------------- ##
## -- Making directory to save results
mkdir -p ${piplineID}_QCreport
mkdir -p ${piplineID}_fastaFiles
mkdir -p workfile

## ------------------------------------------------------------- ##

## -- Version check
echo "## ++++++ ${piplineID} +++++ ##" >> log_and_script/Version.txt
multiqc --version >> log_and_script/Version.txt
fastqc --version >> log_and_script/Version.txt

####################################################################
echo "Start at $timestamp"

## -- Quality check before FiltTrim
## --Running in parallel proccess
echo > 0101.txt

for f in `ls ${fastadir} | grep fastq.gz ` ; do
	echo "fastqc -o workfile ${fastadir}/${f} " >> 0101.txt
done

parallel -a 0101.txt
rm 0101.txt

## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
## -- Merge quality check results
cd workfile
multiqc --filename before --outdir BeforeFiltering ./ 
cd ../

cp -r workfile/BeforeFiltering ${piplineID}_QCreport
rm -r workfile/*

## ============================================================ ##

cat <<RRR > log_and_script/script${piplineID}.R

########################################################################
## -- DADA2
## -- Running on R

## == Setting parameter
minLen=$minLen
minQ=$minQ
maxEE=$maxEE
truncR=$truncR

inputdir="$fastadir"

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
piplineID="$piplineID"

########################################################################

RRR

cat <<'RRR' >> log_and_script/script${piplineID}.R

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

RRR

Rscript log_and_script/script${piplineID}.R

find ${piplineID}_fastaFiles -size -500c -delete
####################################################################

## ============================================================ ##
## Quality check after filtering

echo > 0101.txt
for f in `ls ${piplineID}_fastaFiles| grep fastq.gz`; do
	echo "fastqc -o workfile ${piplineID}_fastaFiles/${f} " >> 0101.txt
done

parallel -a 0101.txt
rm 0101.txt

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
cd workfile
multiqc --filename after --outdir AfterFiltering --tag _filt.fastq.gz ./ 
cd ../

cp -r workfile/AfterFiltering ${piplineID}_QCreport
rm -r workfile

echo "Finish at $timestamp"
####################################################################

