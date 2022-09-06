
#!/bin/bash
## ============================================= ##
## -- Making scripts
sh /home/toju/Desktop/Scripts/fujitaScripts_version20220422/Filtering2Annotation_support/02_cutadapt_R.sh NNNNNNTGTCRTTCGAATTACCTGC NNNNNNTCGATGTARTARGCGTGTTCCCA 32 /home/toju/miniconda3/bin/cutadapt

sh /home/toju/Desktop/Scripts/fujitaScripts_version20220422/Filtering2Annotation_support/03_readQC_dada2.sh 200 10 0.2 5 2 0 15 32
sh /home/toju/Desktop/Scripts/fujitaScripts_version20220422/Filtering2Annotation_support/04_DenoisingRemoveCHIMERA_dada2.sh 32 1 /home/toju/miniconda3/bin/vsearch
sh /home/toju/Desktop/Scripts/fujitaScripts_version20220422/Filtering2Annotation_support/05_TaxonomyAnnotation_claident.sh /home/toju/Desktop/Scripts/fujitaScripts_version20220228/referenceDB/NULL 32 0.05 19 5,90%
sh /home/toju/Desktop/Scripts/fujitaScripts_version20220422/Filtering2Annotation_support/DRA_helper_script01_rename_fastq.sh
sh /home/toju/Desktop/Scripts/fujitaScripts_version20220422/Filtering2Annotation_support/06_QualityCheck.sh 32

## ============================================= ##

