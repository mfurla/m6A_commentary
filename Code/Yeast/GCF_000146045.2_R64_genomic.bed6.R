library(ape)

gtf_full <- read.gff("GCF_000146045.2_R64_genomic.gff", na.strings = c(".", "?"), GFF3 = TRUE)
gtf <- gtf_full[gtf_full[,"type"]=="gene",]

geneNames <- unname(sapply(gtf[,"attributes"],function(i)strsplit(i,";")[[1]][[1]]))
geneNames <- unname(sapply(geneNames,function(i)strsplit(i,"ID=gene-")[[1]][[2]]))

myBed <- cbind(as.character(gtf[,"seqid"]),gtf[,"start"],gtf[,"end"],unname(geneNames),rep(".",nrow(gtf)),as.character(gtf[,"strand"]))

write.table(x=myBed,file="GCF_000146045.2_R64_genomic.bed6",append=FALSE,quote=FALSE,sep="\t", row.names=FALSE,col.names=FALSE)
