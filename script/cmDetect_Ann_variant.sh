#! /usr/bin/bash

#################################
# annotate the vcf file with:
# dbsnp
# COSMIC
# functional
# dbnsfp
################################


## setup parameters
source parameters.sh

SampleId=$1
vcffile=$2
outputDir=$3


## dbsnp
$java -Xmx2g -jar $snpEffDir/SnpSift.jar annotate -v $refDbsnp $vcffile > $outputDir/$SampleId.dbsnp.vcf
checkError "Annotation dbsnp" $SampleId

## COSMIC
$java -Xmx2g -jar $snpEffDir/SnpSift.jar annotate -v $refCosmic $outputDir/$SampleId.dbsnp.vcf > $outputDir/$SampleId.dbsnp.cosmic.vcf
checkError "Annotation COSMIC" $SampleId

## Fonctional 
$java -Xmx2g -jar $snpEffDir/snpEff.jar eff -c $snpEffDir/snpEff.config -v -hgvs -t -noLog -noStats -noMotif -noNextProt $genome $outputDir/$SampleId.dbsnp.cosmic.vcf > $outputDir/$SampleId.dbsnp.cosmic.hg19.vcf
checkError "Annotation Fonctional" $SampleId

## dbNSFP
$java -Xmx2g -jar $snpEffDir/SnpSift.jar dbnsfp -v -db $dbNsfp $outputDir/$SampleId.dbsnp.cosmic.hg19.vcf > $outputDir/$SampleId.dbsnp.cosmic.hg19.eff.vcf
checkError "Annotation dbNSFP"  $SampleId

## Extract vcf to txt
$java -Xmx2g -jar $snpEffDir/SnpSift.jar extractFields -e "."  $outputDir/$SampleId.dbsnp.cosmic.hg19.eff.vcf $header >> $outputDir/$SampleId.ann.tab
checkError "Annotation Extract vcf to txt" $SampleId

## clean up
mv $outputDir/$SampleId.dbsnp.cosmic.hg19.eff.vcf $outputDir/$SampleId.ann.vcf
rm $outputDir/$SampleId.dbsnp*vcf
