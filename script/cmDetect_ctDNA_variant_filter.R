#!/bin/Rscript

###################################################
# for each variant, calculate 
# empirical p value of being a polymorphisme
# MAF in the population
# p value of being an error due to sequencing bias
# keep only the variants that pass the filter
###################################################

args <- commandArgs(TRUE)
SampleList=args[1]
inputfile=args[2]
dataDir=args[3]
OutputDir=args[4]

SNP=read.table(file=inputfile, h=T, sep="\t", as.is=T, quote="")
id=read.table(file=SampleList, sep="\t", h=F, as.is=T)
id=id[,1]

#polym distribution
filenames=sapply(unlist(id), function(i) list.files(path=dataDir, pattern=paste(i, ".haplotypeCaller.hetero.filtered.vcf.tab",sep=""), 
                                                    full.names=T, recursive = T))
polym <- do.call(rbind,lapply(filenames, function(file) read.table(file=file, h = F,as.is = T ,sep="\t")))

#Seqbias
filenames=sapply(unlist(id), function(i) list.files(path=dataDir, pattern=paste(i, ".SeqBias.Num",sep=""),
                                                    full.names=T, recursive = T))
Seqbias <- do.call(rbind,lapply(filenames, function(file) read.table(file=file, h = T,as.is = T)))


# variant MAF in samples
SNP$MAF=sapply(1:dim(SNP)[1], function(i){t=Seqbias[intersect(intersect(which(Seqbias$CHROM==SNP$CHROM[i]), which(Seqbias$POS==SNP$POS[i])),
                                                                 which(Seqbias$ALT==SNP$ALT[i])),c("Tdepth.1", "Taltdepth.1", "Ndepth.1", "Naltdepth.1", "id.1")]
                                          if(length(which(duplicated(t)==T))!=0){t=t[-which(duplicated(t)==T),]}
                                          t$Tfreq.1=t$Taltdepth.1/t$Tdepth.1
                                          t$Nfreq.1=t$Naltdepth.1/t$Ndepth.1
                                          length(intersect(which(t$Tfreq.1>=0.1), which(t$Nfreq>=0.1)))/length(id)})

# pv polym
SNP$pvpolym=sapply(1:dim(SNP)[1],function(i) sum(polym[which(polym[,8]==SNP$id[i]),7] < SNP$NFreq[i])/sum(which(polym[,8]==SNP$id[i])))
SNP$pvpolym.FDR=p.adjust(as.numeric(SNP$pvpolym), method="fdr")

# Sequencing Bias
SNP$Seqbias=sapply(1:dim(SNP)[1],function(i){t=Seqbias[intersect(intersect(which(Seqbias$CHROM==SNP$CHROM[i]), which(Seqbias$POS==SNP$POS[i])),
                                                                 which(Seqbias$ALT==SNP$ALT[i])),c("Tdepth.1", "Taltdepth.1", "Ndepth.1", "Naltdepth.1", "id.1")]
                                             if(length(which(duplicated(t)==T))!=0){t=t[-which(duplicated(t)==T),]}
                                             j=which(t$id.1==SNP$id[i])
                                             pbinom(t$Naltdepth.1[j]-1, t$Ndepth.1[j], sum(t$Naltdepth.1[-j])/sum(t$Ndepth.1[-j]),lower.tail =F)})
SNP$Seqbias.FDR=p.adjust(as.numeric(SNP$Seqbias), method="fdr")

write.table(SNP, paste(OutputDir, "/ctDNA.variant.txt",sep=""), sep="\t", row.names=F, quote=F)

# filter SNP by database: keep variants not in dbsnp/in COSMIC with MAF(1000G) < 0.001
SNP=SNP[unique(c(which(regexpr("COSM",SNP$ID)!= -1), which(SNP$ID=="."))),]
f1000G=sapply(1:dim(SNP)[1], function(i) max(c(SNP$dbNSFP_1000Gp1_AF[i], SNP$dbNSFP_ESP6500_AA_AF[i], SNP$dbNSFP_ESP6500_EA_AF[i]))>0.001) 
SNP=SNP[which(f1000G==F),]
# variant with MAF<0.1
SNP=SNP[which(SNP$MAF<0.1),]

# pv polym < 0.01
SNP=SNP[which(SNP$pvpolym.FDR<0.01),]

# Seqbias filter
SNP$Seqbias.FDR=p.adjust(as.numeric(SNP$Seqbias), method="fdr")
SNP=SNP[which(SNP$Seqbias.FDR<0.05),]

write.table(SNP, paste(OutputDir, "/ctDNA.variant.final.txt",sep=""), sep="\t", row.names=F, quote=F)