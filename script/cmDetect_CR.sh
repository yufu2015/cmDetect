#!/bin/bash

######################################
# retrieve and count supporting reads 
# from bam files using samtools 
########################################

## define parameters
source parameters.sh

SampleId=$1
variantfile=$2
Tbamfile=$3
Nbamfile=$4
OutputDir=$5
filename=$6

## make bed file
awk '{FS=OFS="\t"} {print $1"\t"($2-1)"\t"$2}' $variantfile | sort -k1,1 | uniq >  $OutputDir/$SampleId.$filename.bed
awk '{FS=OFS="\t"} {print $1"_"$2"\t"$0 }' $variantfile | sort -k1,1 | uniq > $OutputDir/$SampleId.$filename
checkError "CR make bed file" $SampleId

## get supporting read using samtools
$samtools mpileup -q 20 -Q 20 -l $OutputDir/$SampleId.$filename.bed -f $refFasta $Tbamfile $Nbamfile  > $OutputDir/$SampleId.$filename.read
awk '{FS=OFS="\t"} $4>0 && $7>0 {$NF=""; $(NF-3)=""; print}' $OutputDir/$SampleId.$filename.read | awk -v OFS="\t" '$1=$1' > $OutputDir/$SampleId.$filename.read.tmp

checkError "CR get supporting read using samtools" $SampleId

mv $OutputDir/$SampleId.$filename.read.tmp $OutputDir/$SampleId.$filename.read

## add samtools results in variant file
awk '{FS=OFS="\t"}{ print $1"_"$2"\t"$4"\t"$5"\t"$6"\t"$7 }' $OutputDir/$SampleId.$filename.read |sort -k1,1 | uniq > $OutputDir/$SampleId.$filename.read.sort
join $OutputDir/$SampleId.$filename $OutputDir/$SampleId.$filename.read.sort -t $'\t' | cut -d $'\t' -f 2- > $OutputDir/$SampleId.$filename.merge

checkError "CR add samtools results in variant file" $SampleId

rm $OutputDir/$SampleId.$filename.read.sort

## remove indel for point mutation before count variant supporting reads
awk 'BEGIN{FS=OFS="\t"} length($3)==1 && length($4)==1{
    if($(NF-2)~/[-+]/){
        num=split($(NF-2),a,"\\-|\\+");
        for(i=2;i<=num;i++){
            numm=split(a[i],aa,"")
            if(aa[2]~/[atcgATCG]/){
                a[i] = substr(a[i],2+aa[1],numm-(1+aa[1]))
                }else{
                    n=aa[1]*10+aa[2]
                    a[i] = substr(a[i],3+n,numm-(2+n))}}
        result = a[1]
        for (i = 2; i <= num; i++){result = result "" a[i]}
        $(NF-2)=result}
    if($NF~/[-+]/){
        num=split($NF,a,"\\-|\\+")
        for(i=2;i<=num;i++){
            numm=split(a[i],aa,"")
            if(aa[2]~/[atcgATCG]/){
                 a[i] = substr(a[i],2+aa[1],numm-(1+aa[1]))
                 }else{
                     n=aa[1]*10+aa[2]
                     a[i] = substr(a[i],3+n,numm-(2+n))}}
        result = a[1]
        for (i = 2; i <= num; i++){result = result "" a[i]}
        $NF=result}
     print $0}' $OutputDir/$SampleId.$filename.merge > $OutputDir/$SampleId.$filename.merge.tmp

checkError "CR remove indels for SNP"  $SampleId

awk 'BEGIN{FS=OFS="\t"} !(length($3)==1 && length($4)==1) {print $0}' $OutputDir/$SampleId.$filename.merge >> $OutputDir/$SampleId.$filename.merge.tmp
mv $OutputDir/$SampleId.$filename.merge.tmp $OutputDir/$SampleId.$filename.merge

# count variant supporting reads
head -1  $variantfile | awk '{FS=OFS="\t"} {gsub("#","");print $0"\tTdepth\tTaltdepth\tNdepth\tNaltdepth\tid" }' > $OutputDir/$SampleId.$filename.Num 
awk -v id=$SampleId '{FS=OFS="\t"}{motive=$4;
                                   if(length($3)>1){
                                     motive="-" length($3)-1 substr($3,2,length($3))
                                   }
                                   if(length($4)>1){
                                     motive="+" length($4)-1 substr($4,2,length($4))
                                   }
                                   Count1=split(toupper($(NF-2)),a,motive)-1;
                                   Count2=split(toupper($NF),a,motive)-1;
                                   $(NF-2)=Count1; $NF=Count2; 
                                   print $0"\t"id}' $OutputDir/$SampleId.$filename.merge >> $OutputDir/$SampleId.$filename.Num 

checkError "CR count variant supporting reads" $SampleId
