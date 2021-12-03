##### This code is aimed at sub-sampling the WT sample #####

### Libraries
library(parallel)

### Custom function
readFastaFile <- function(file)
{
	fastaTmp <- read.table(file,header=FALSE)
	fastaTmp <- apply(fastaTmp,2,as.character)
	fastaTmp <- fastaTmp[,1]
	
	fastaTmpNames <- fastaTmp[seq_along(fastaTmp)%%2==1]
	fastaTmp <- fastaTmp[seq_along(fastaTmp)%%2==0]
	names(fastaTmp) <- fastaTmpNames
	fastaTmp
}

### Same amount of k-mers
## Number of reads WT * (Number of k-mers KO / Number of k-mers WT)
newReads_WT_equal_KO <- round(1037785 * (11592309/15829173))

### Smaller amount of k-mers
## Number of reads WT * ((Number of k-mers WT * (Number of k-mers KO / Number of k-mers WT))/Number of k-mers WT)
newSites_WT_smaller_KO <- round(11592309*(11592309/15829173))
newReads_WT_smaller_KO <- round(1037785 * (newSites_WT_smaller_KO/15829173))

### Reads selection
## All WT reads
readsFasta <- readFastaFile("/WT/nanom6A_step1/result.feature.fa") # File not released on GitHub due to its size.

## WT = KO
set.seed(1)
newReads_WT_equal_KO <- sample(names(readsFasta),size=newReads_WT_equal_KO)
readsFasta_WT_equal_KO <- readsFasta[newReads_WT_equal_KO]

write.table(readsFasta_WT_equal_KO,"WT_equal_KO/nanom6A_step1/result.feature.fa",quote=FALSE,col.names=FALSE,row.names=TRUE,sep="\n") # File not released on GitHub due to its size.

## WT < KO
set.seed(2)
newReads_WT_smaller_KO <- sample(names(readsFasta),size=newReads_WT_smaller_KO)
readsFasta_WT_smaller_KO <- readsFasta[newReads_WT_smaller_KO]

write.table(readsFasta_WT_smaller_KO,"/WT_smaller_KO/nanom6A_step1/result.feature.fa",quote=FALSE,col.names=FALSE,row.names=TRUE,sep="\n") # File not released on GitHub due to its size.

### Indexes
## WT = KO
system("samtools faidx /WT_equal_KO/nanom6A_step1/result.feature.fa")

## WT < KO
system("samtools faidx /WT_smaller_KO/nanom6A_step1/result.feature.fa")

### Features selection (only the features from the final reads)
readsFasta_WT_equal_KO <- readFastaFile("WT_equal_KO/nanom6A_step1/result.feature.fa") # File not released on GitHub due to its size.
readsFasta_WT_equal_KO <- names(readsFasta_WT_equal_KO)

readsFasta_WT_smaller_KO <- readFastaFile("WT_smaller_KO/nanom6A_step1/result.feature.fa") # File not released on GitHub due to its size.
readsFasta_WT_smaller_KO <- names(readsFasta_WT_smaller_KO)

## WT features
result.feature <- readLines("WT/nanom6A_step1/result.feature.tsv") # File not released on GitHub due to its size.

index.list = split(seq_along(result.feature), ceiling(seq_along(result.feature)/1000000))

namesTmp <- list()

for(i in seq_along(index.list))
{
	print(i)
	j <- index.list[[i]]
	namesTmp <- append(namesTmp,mclapply(result.feature[j],function(q)strsplit(q,"[|]")[[1]][[1]],mc.cores=1))
}

names(result.feature) <- paste0(">",unlist(namesTmp))

## WT = KO features
result.feature_WT_equal_KO <- result.feature[names(result.feature)%in%readsFasta_WT_equal_KO]
write.table(as.data.frame(result.feature_WT_equal_KO),"WT_equal_KO/nanom6A_step1/result.feature.tsv",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE) # File not released on GitHub due to its size.

## WT < KO features
result.feature_WT_smaller_KO <- result.feature[names(result.feature)%in%readsFasta_WT_smaller_KO]
write.table(as.data.frame(result.feature_WT_smaller_KO),"WT_smaller_KO/nanom6A_step1/result.feature.tsv",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE) # File not released on GitHub due to its size.

