#! /usr/bin/bash

#######################################################################
# Please define all the necessary parameters for the cmDetect workflow
# give complete path
#######################################################################

## cmDetect directory 
cmDetectDir=

## id of samples in the analysis per line 
SampleList=$cmDetectDir/script/SampleList.txt

## memory use in GATK, exemple: 4g 
javaMem=

## files
# genome reference used for the alignment
refFasta=

# database in the reference genome version 
genome=
refDbsnp=
refCosmic=
dbNsfp=


## tools
gatk=GenomeAnalysisTK.jar
vcfsort=vcftools_0.1.12b/bin/vcf-sort
samtools=samtools-0.1.19/samtools
R=R-3.2.1/bin/R
snpEffDir=snpEff-4.1c
java=jre1.7.0_25/bin/java

## output file format
header="CHROM POS REF ALT ID ANN[0].GENEID GEN[0].AD[0] GEN[0].AD[1] GEN[0].GT GEN[0].GQ dbNSFP_1000Gp1_AF dbNSFP_ESP6500_AA_AF dbNSFP_ESP6500_EA_AF ANN[0].EFFECT ANN[0].IMPACT ANN[0].FEATUREID  ANN[0].HGVS_C ANN[0].HGVS_P dbNSFP_Uniprot_acc dbNSFP_Polyphen2_HVAR_pred dbNSFP_SIFT_pred"

## variant filter
depth=10
Taltdepth=5
Tfreq=0.1

###############
## Check error
###############
# checkError step id

checkError() {
        if [[ ! $? -eq 0 ]]; then
            echo "ERROR in step $1, $2" | tee -a $cmDetectDir/log.txt
            #exit 1
        fi
}
