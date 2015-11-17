#! /usr/bin/bash

###################################################
# ctDNA Mutation Detection pipeline (cmDetect)
# detecte Somatic mutations with supporting 
# reads in the whole blood sample due to
# ctDNA contamination
###################################################

###################################################
## Setup parameters
###################################################
source parameters.sh

mkdir -p $cmDetectDir/data
mkdir -p $cmDetectDir/result

###################################################
## run GATK on all the tumor/normal sample 
## the following step may be done parallely, please 
## verify all the samples are succesfully analysed 
## before the next step
###################################################

while IFS=$'\t' read -r -a line
do
id=${line[0]}
mkdir -p $cmDetectDir/data/$id
Tbam=${line[1]}
mkdir -p $cmDetectDir/data/$id/Tumor
bash $cmDetectDir/script/cmDetect_call_variant.sh $id $Tbam $cmDetectDir/data/$id/Tumor
Nbam=${line[2]}
mkdir -p $cmDetectDir/data/$id/Blood
bash $cmDetectDir/script/cmDetect_call_variant.sh $id $Nbam $cmDetectDir/data/$id/Blood
done < $SampleList

###################################################
## Annotate variants called in tumor using SnpEff 
## the following step may be done parallely, please
## verify all the samples are succesfully analysed
## before the next step
###################################################

while IFS=$'\t' read -r -a line
do
id=${line[0]}
bash $cmDetectDir/script/cmDetect_Ann_variant.sh $id $cmDetectDir/data/$id/Tumor/$id.haplotypeCaller.hetero.filtered.vcf $cmDetectDir/data/$id/Tumor
done < $SampleList

###################################################
## retrieve supporting reads of called variant in 
## tumor and normal sample
## the following step may be done parallely, please
## verify all the samples are succesfully analysed
## before the next step
###################################################

while IFS=$'\t' read -r -a line
do
id=${line[0]}
Tbam=${line[1]}
Nbam=${line[2]}
variantfile=$cmDetectDir/data/$id/Tumor/$id.ann.tab
awk -v Tfreq="$Tfreq" -v depth="$depth" -v Taltdepth="$Talrdepth" '{FS=OFS="\t"} NR==1{print} NR>1 && $8/($7+$8) >= Tfreq  && $7+$8 >= depth && $8 >= Taltdeth {print}' $variantfile > $variantfile.tmp
mv $variantfile.tmp $variantfile
bash $cmDetectDir/script/cmDetect_CR.sh $id $variantfile $Tbam $Nbam $cmDetectDir/data/$id/Tumor variant
done < $SampleList

###################################################
## do fisher's exact test
###################################################

while IFS=$'\t' read -r -a line
do
id=${line[0]}
$R --vanilla --args $id $cmDetectDir/data/$id/Tumor/$id.variant.Num < $cmDetectDir/script/cmDetect_FisherTest.R
done < $SampleList

###################################################
## do FDR
###################################################

$R --vanilla --args $SampleList $cmDetectDir/data $cmDetectDir/result < $cmDetectDir/script/cmDetect_FDR.R

###################################################
## retrieve supporting reads in tumor and normal 
## samples on Variants passed fisher FDR
## the following step may be done parallely, please
## verify all the samples are succesfully analysed
## before the next step
###################################################

while IFS=$'\t' read -r -a line
do
id=${line[0]}
Tbam=${line[1]}
Nbam=${line[2]}
variantfile=$cmDetectDir/result/all.variant.FDR001.txt
bash $cmDetectDir/script/cmDetect_CR.sh $id $variantfile $Tbam $Nbam $cmDetectDir/data/$id/Tumor/ SeqBias
done < $SampleList

###################################################
## generate frequency distribution of sample 
## specific polymorphism
###################################################

bash $cmDetectDir/script/cmDetect_Dist_Polym.sh 

###################################################
## Call ctDNA variant
###################################################

$R --vanilla --args $SampleList $cmDetectDir/result/all.variant.FDR001.txt $cmDetectDir/data $cmDetectDir/result <  $cmDetectDir/script/cmDetect_ctDNA_variant_filter.R

