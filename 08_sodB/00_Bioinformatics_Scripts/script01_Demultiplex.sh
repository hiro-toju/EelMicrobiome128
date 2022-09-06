#!/bin/bash

####################################################################
## 											
## ---------- 	    Demultiplexing by Clident       ------------- ##
## 
## 											2022. 02. 28. by Fujita
####################################################################

## -- Options
runname=Sod_B # Experience name
index1=index1.txt # Fasta file containing R index tag sequence(P7-Hamady)
index2=index2.txt # Fasta file containing F index tag sequence(P5-Hamady)

fprimer=Fprimer.txt # Fasta file containing forward primer file
rprimer=Rprimer.txt # Fasta file containing reverse primer file

demultidir=demultiplex

thread=32

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
piplineID=01_Demultiplex

####################################################################

## ============= Remove and Make directories ==================== ##

## -- Remove older directory
if [ -d ${demultidir} ]; then
  rm -r ${demultidir}
fi

## -- Version check
echo "## ++++++ ${piplineID} +++++ ##" > log_and_script/Version.txt
clsplitseq 2> log_and_script/tmp.txt
head -n 1 log_and_script/tmp.txt | cat >> log_and_script/Version.txt; rm log_and_script/tmp.txt

####################################################################
## -- Demultiplexed 
# Removing undeterminated sequence and cut primer by using Claident

echo "Start at $timestamp"
clsplitseq --runname=$runname --truncateN=enable --minqualtag=30 \
		   --index1file=$index1 --index2file=$index2 --primerfile=$fprimer --reverseprimerfile=$rprimer \
		   --numthreads=$thread --outputjump=DISABLE \
		   lane1_NoIndex_L001_R1_001.fastq.gz lane1_NoIndex_L001_R2_001.fastq.gz \
		   		lane1_NoIndex_L001_R3_001.fastq.gz lane1_NoIndex_L001_R4_001.fastq.gz \
		   $demultidir 2>&1 | tee log_and_script/log${piplineID}.txt
echo "Finish at $timestamp"

#####################################################################
