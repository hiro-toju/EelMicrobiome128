#!/bin/bash

####################################################################
## 											
## ------ Taxonomy Annotation with  QC auto method  ------------- ##
##
## 									2022. 02. 28. by Fujita
####################################################################

inputdir=04_Denoising
outputdir=05_TaxonomyAnnotation

referencepath=/home/toju/Desktop/Scripts/fujitaScripts_version20220228/referenceDB/NULL
thread=32
maxpopposer=0.05
minsoratio=19
NNC=5,90%

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
piplineID="05_TaxonomyAnnotation"

####################################################################

## ============= Remove and Make directories ==================== ##

## --Remove older directory
if ls ${piplineID}_QCreport* >/dev/null 2>&1 ; then
	rm -r ${piplineID}
fi

if ls log_and_script/${piplineID}_log* >/dev/null 2>&1 ; then
	rm log_and_script/log${piplineID}
fi

## ------------------------------------------------------------- ##
## -- Making directory to save results
mkdir -p ${piplineID}

## ------------------------------------------------------------- ##
## -- Version check
echo "## ++++++ ${piplineID} +++++ ##" >> log_and_script/Version.txt

####################################################################

## -- Make cache database
clmakecachedb \
--blastdb=overall_genus \
--numthreads=$thread \
${inputdir}/nonchim_seq.fasta \
${outputdir}/01_cachedb_species

## ============================================================= ##
## --Assign taxonomy based on QCauto method
clidentseq \
--method=QC \
--blastdb=${outputdir}/01_cachedb_species \
--numthreads=$thread \
${inputdir}/nonchim_seq.fasta \
${outputdir}/02_neighborhoods_QC.txt

## -- Strict threshold
classigntax \
--taxdb=overall_genus \
${outputdir}/02_neighborhoods_QC.txt \
${outputdir}/03_taxalist_QC_strict.tsv

## -- Relaxed threshold
classigntax \
--taxdb=overall_genus --maxpopposer=${maxpopposer} --minsoratio=${minsoratio} \
${outputdir}/02_neighborhoods_QC.txt \
${outputdir}/03_taxalist_QC_relax.tsv

## ============================================================= ##
## --Assign taxonomy based on (95%-)5-NN method
clidentseq \
--method=${NNC} \
--blastdb=${outputdir}/01_cachedb_species  \
--numthreads=32 \
${inputdir}/nonchim_seq.fasta \
${outputdir}/02_neighborhoods_NNC.txt

classigntax \
--taxdb=overall_genus --minnsupporter=1 \
${outputdir}/02_neighborhoods_NNC.txt \
${outputdir}/03_taxalist_NNC.tsv


## ============================================================= ##

clmergeassign \
--preferlower \
--priority=descend \
${outputdir}/03_taxalist_QC_strict.tsv \
${outputdir}/03_taxalist_QC_relax.tsv \
${outputdir}/03_taxalist_NNC.tsv \
${outputdir}/taxonomy_merged.tsv

## ============================================================= ##
# Fill blank cells of taxonomic assignment
clfillassign \
${outputdir}/taxonomy_merged.tsv \
${outputdir}/taxonomy_merged_filled.tsv

cat <<RRR > log_and_script/script05_claident_reformat.R

## ============================================================= ##

outputdir="${outputdir}"

## ============================================================= ##
RRR
cat <<RRR >> log_and_script/script05_claident_reformat.R

df <- read.table(sprintf("%s/taxonomy_merged_filled.tsv",outputdir), header=T, row.names=1, sep="\t")
df2 <- df[, intersect( colnames(df), c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species") )]

for( i in 1:ncol(df2) ) {
	first.letter <- substr(colnames(df2)[i], 1, 1)
	colnames(df2)[i] <- paste( toupper(first.letter), substr(colnames(df2)[i], 2, nchar(colnames(df2)[i])), sep="")
}

mat <- as.matrix(df2)
mat[which(mat=="")] <- "Unidentified"
mat <- gsub(" ", "_", mat)
mat <- gsub("unidentified", "Unidentified", mat)

df3 <- as.data.frame(mat)

for(i in 1:(ncol(df3)-1)){

	unident.row <- which(df3[,i]!="Unidentified")
	if( length(unident.row ) > 0){
		df3[unident.row,"identified"] <- as.character(df3[unident.row, i]   )	
	}

}

write.table(cbind(ID=rownames(df3), df3), sprintf("%s/taxonomy_list.txt",outputdir), row.names=F, quote=F, sep="\t")
saveRDS(cbind(ID=rownames(df3), df3), sprintf("%s/taxonomy_list.rds",outputdir))

RRR
Rscript log_and_script/script05_claident_reformat.R

