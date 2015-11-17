#!/bin/Rscript

##################################
# for each variant, do fisher's 
# exact test using variant evidence
# in tumor and blood sample
##################################

args <- commandArgs(TRUE)
id=args[1]
inputfile=args[2]
print(id)

#read mpileup out file
SNP = read.csv(file =inputfile,h = T, as.is = T ,sep="\t")
                       
SNP$TFreq = SNP$Taltdepth/SNP$Tdepth
SNP$NFreq = SNP$Naltdepth/SNP$Ndepth

SNP$pvfisher = sapply(1:dim(SNP)[1], function(i) {fisher.test(matrix(c(SNP$Taltdepth[i],
                                                                       SNP$Tdepth[i] - SNP$Taltdepth[i],
                                                                       SNP$Naltdepth[i],
                                                                       SNP$Ndepth[i] - SNP$Naltdepth[i]),
                                                              nr=2),alternative="greater")$p.value})

write.table(SNP, file=paste(inputfile,".Fish",sep=""),sep="\t",row.names=F,quote=F)

