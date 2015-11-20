#!/bin/bash

########################################
# retrieve the allele frequency of 
# all the heterozygous germline variants 
########################################

## define parameters
source parameters.sh
echo "polym distribution "

while IFS=$'\t' read -r -a line
do

id=${line[0]}

echo $id

SNPfile=$cmDetectDir/data/$id/Blood/$id.haplotypeCaller.hetero.filtered.vcf

$java -Xmx2g -jar $snpEffDir/SnpSift.jar extractFields -e "." $SNPfile CHROM POS REF ALT GEN[0].AD[0] GEN[0].AD[1] > $SNPfile.tab

awk -v id=$id '{FS=OFS="\t"} NR==1 {print} NR>1 {print $0"\t"$6/($5+$6)"\t"id}' $SNPfile.tab > $SNPfile.tab.tmp

checkError "polym distribution" $id

mv $SNPfile.tab.tmp $SNPfile.tab 

done < $SampleList
