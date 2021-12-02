library(ape)

gtf_full <- read.gff("GRCh37_latest_genomic.gtf", na.strings = c(".", "?"), GFF3 = FALSE)
gtf <- gtf_full[gtf_full[,"feature"]=="gene",]

geneNames <- unname(sapply(gtf[,"attributes"],function(i)strsplit(i,";")[[1]][[1]]))
geneNames <- unname(sapply(geneNames,function(i)strsplit(i,"gene_id ")[[1]][[2]]))
geneNames <- unname(sapply(geneNames,function(i)strsplit(i,"\"")[[1]][[2]]))

myBed <- cbind(as.character(gtf[,"seqname"]),gtf[,"start"],gtf[,"end"],unname(geneNames),rep(".",nrow(gtf)),as.character(gtf[,"strand"]))

write.table(x=myBed,file="GRCh37_latest_genomic.bed6",append=FALSE,quote=FALSE,sep="\t", row.names=FALSE,col.names=FALSE)