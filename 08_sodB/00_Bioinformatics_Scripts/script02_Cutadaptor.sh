#!/bin/bash

####################################################################
## 											
## ---------- 	  Check and cut sequence adaptor     ------------ ##
## 
## 					2022. 02. 28. by Fujita
####################################################################

inputdir=01_Demultiplexed_fastaq
outputdir=02_Cutadaptor_fastaq
fseq=NNNNNNTGTCRTTCGAATTACCTGC
rseq=NNNNNNTCGATGTARTARGCGTGTTCCCA
thread=32
cutadaptpath=/home/toju/miniconda3/bin/cutadapt
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
piplineID=02_Cutadaptor

####################################################################

## ============= Remove and Make directories ==================== ##

## -- Remove older directory
if [ -d $outputdir ]; then
  rm -r $outputdir
fi

## -- Version check
echo "## ++++++ ${piplineID} +++++ ##" >> log_and_script/Version.txt
cv=`cutadapt --version`
echo "cutadapt ${cv}" >> log_and_script/Version.txt

########################################################################
## -- Cutadaptor by R

cat <<RRR > log_and_script/script${piplineID}.R

inputdir="${inputdir}"
outputdir="${piplineID}_fastaq"
primerf="${fseq}"
primerr="${rseq}"
thread=${thread}
cutadaptpath="${cutadaptpath}"

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
piplineID="${piplineID}"

start <- Sys.time()
print(start); cat("\n")
########################################################################
RRR

cat <<'RRR' >> log_and_script/script02_Cutadaptor.R
## ===================== Definition of function =====================  ##
primerHits <- function(primer, fn) {
   # Counts number of reads in which the primer is found
   nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
  }

allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
}

Check.primer = function(fwd= fwd, rev =rvs,
                        path=path, input.fastq=fnFs){

  FWD.orients <- allOrients(fwd)
  REV.orients <- allOrients(rev)

  tmp1 <- sapply(FWD.orients, primerHits, fn = list.files(path, full.names = TRUE))
  tmp2 = sapply(REV.orients, primerHits, fn = list.files(path, full.names = TRUE))
  rbind(tmp1,  tmp2)
  
}

########################################################################

######  ====================== Main part ======================  #######

########################################################################

# -- Create directory to save
dir.create(outputdir, showWarnings = F) 

# -- Load library and function
library(seqinr) ; library(stringr);library(ShortRead)
library(Biostrings) ; library(dada2); library(doParallel)

# Difined the directory containing the fastq files after unzipping
path.dir <- inputdir                                   

if( length( list.files(path.dir) ) > 0 ){
  path = path.dir
}else{
  stop("The path directory is missing\n")
}

################################################################################################

# Infer Sequence Variants
fnFs <- sort(list.files(path, full.names = TRUE))
#sample.names <- apply( sapply(strsplit(basename(fnFs), "_"), , c(1,2)), 2, paste, collapse="_")

FWD.ForwardReads = Check.primer(fwd= primerf, rev=primerr, 
                                input.fastq=fnFs, path=path )


#if(sum(FWD.ForwardReads)!=0){ 
  
  cat(sprintf('The primer found in  %s samples.\n', sum(FWD.ForwardReads)))
  print(FWD.ForwardReads)
  
  fnFs.cut <- file.path(outputdir, basename(fnFs))
  
  RVS.RC <- dada2:::rc(primerr)
  
  # Trim FWD and the reverse-complement of REV off of R1 (forward reads)
  R1.flags <- paste("-g", primerf,  "-a", RVS.RC) 
  
  cl <- makeCluster(detectCores(logical=FALSE))
  registerDoParallel(cl)
  
  # Run Cutadapt
  tmp <- foreach(f=1:length(fnFs)) %dopar% {
    invisible(system2(cutadaptpath, args = c(R1.flags, "-n", 2,  "-j", 2, "--max-n",0, # -n 2 required to remove FWD and REV from reads
                                         "-o", fnFs.cut[f], "-m", 10, # output files
                                         fnFs[f])) ) # input files
  }

#}

finish <- Sys.time()
print(sprintf("Finish. %s", finish-start)); cat("\n")

RRR
Rscript log_and_script/script02_Cutadaptor.R 2>&1 | tee  log_and_script/log02_Cutadaptor.txt

find 02_Cutadaptor_fastaq -size -500c -delete

