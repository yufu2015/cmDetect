#!/bin/Rscript

##########################
# fisher's exact test FDR
##########################

args <- commandArgs(TRUE)
SampleList=args[1]
dataDir=args[2]
OutputDir=args[3]

id=read.table(file=SampleList, sep="\t", h=F, as.is=T)
id=id[,1]
filenames=sapply(unlist(id), function(i) list.files(path=dataDir, pattern=paste(i, ".variant.Num.Fish",sep=""), full.names=T, recursive = T))
data <- do.call(rbind,lapply(filenames, function(file) read.csv(file=file, h = T,as.is = T ,sep="\t")))
data=data[which(data$Naltdepth>0),]
data$fisher.FDR = p.adjust(as.numeric(data$pvfisher), method="fdr")
write.table(data, sep="\t", file=paste(OutputDir,"/all.variant.FDR.txt",sep=""), quote=F, row.names=F)
write.table(data[which(data$fisher.FDR<0.01),], sep="\t", file=paste(OutputDir,"/all.variant.FDR001.txt",sep=""), quote=F, row.names=F)
