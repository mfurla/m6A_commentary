##### This code is aimed at creating sets of shuffled modification probabilities #####

#### For all the available 5-mers (WT > KO)
### Mod and unmod tables
## KO
filesKO <- list.files("KO/nanom6A_step2/")
modFilesKO <- filesKO[grep("[.]mod",filesKO)]
unmodFilesKO <- filesKO[grep("[.]unmod",filesKO)]
filesKO <- c(modFilesKO, unmodFilesKO)

## WT
filesWT <- list.files("WT/nanom6A_step2/")
modFilesWT <- filesWT[grep("[.]mod",filesWT)]
unmodFilesWT <- filesWT[grep("[.]unmod",filesWT)]
filesWT <- c(modFilesWT, unmodFilesWT)

### All data
## KO
tablesKO <- lapply(filesKO,function(file)
{
	print(file)
	read.table(paste0("KO/nanom6A_step2/",file),sep="\t",header=FALSE)
})
elementsKO <- sapply(tablesKO,nrow)
names(elementsKO) <- filesKO
elementsKO <- sapply(names(elementsKO),function(k)
{
	rep(k,elementsKO[[k]])
})
elementsKO <- unlist(elementsKO)
elementsKO <- sapply(strsplit(elementsKO,"[.]"),"[[",1)

tablesKO <- do.call(rbind,tablesKO)

## WT
tablesWT <- lapply(filesWT,function(file)
{
	print(file)
	read.table(paste0("WT/nanom6A_step2/",file),sep="\t",header=FALSE)
})
elementsWT <- sapply(tablesWT,nrow)
names(elementsWT) <- filesWT
elementsWT <- sapply(names(elementsWT),function(k)
{
	rep(k,elementsWT[[k]])
})
elementsWT <- unlist(elementsWT)
elementsWT <- sapply(strsplit(elementsWT,"[.]"),"[[",1)

tablesWT <- do.call(rbind,tablesWT)

### All probabilities shuffling
shuffledProbabilities <- c(tablesKO[,"V1"],tablesWT[,"V1"])
set.seed(1)
shuffledProbabilities <- shuffledProbabilities[sample(seq_along(shuffledProbabilities),size=length(shuffledProbabilities))]

### Subdivision
## KO vs WT
shuffledProbabilitiesKO <- shuffledProbabilities[1:nrow(tablesKO)]
tablesKO[,"V1"] <- shuffledProbabilitiesKO
tablesKO[,"V2"] <- 1 - tablesKO[,"V1"]

shuffledProbabilitiesWT <- shuffledProbabilities[nrow(tablesKO)+(1:nrow(tablesWT))]
tablesWT[,"V1"] <- shuffledProbabilitiesWT
tablesWT[,"V2"] <- 1 - tablesWT[,"V1"]

## By k-mer
# KO
splittedTablesKO <- split(1:nrow(tablesKO),elementsKO)
splittedTablesKO <- lapply(splittedTablesKO,function(i)tablesKO[i,])

# WT
splittedTablesWT <- split(1:nrow(tablesWT),elementsWT)
splittedTablesWT <- lapply(splittedTablesWT,function(i)tablesWT[i,])

### Split and save data
## KO
for(i in seq_along(splittedTablesKO))
{
	print(i)

	splitTmp <- splittedTablesKO[[i]][,"V2"]>0.5
	modTmp <- splittedTablesKO[[i]][splitTmp,]
	unmodTmp <- splittedTablesKO[[i]][!splitTmp,]

	write.table(x=modTmp
			  , file=paste0("KO_shuffled/nanom6A_step2/",names(splittedTablesKO)[[i]],".shuffledmod")
			  , quote=FALSE
			  , sep="\t"
			  , row.names=FALSE
			  , col.names=FALSE)
	
	write.table(x=unmodTmp
			  , file=paste0("KO_shuffled/nanom6A_step2/",names(splittedTablesKO)[[i]],".shuffledunmod")
			  , quote=FALSE
			  , sep="\t"
			  , row.names=FALSE
			  , col.names=FALSE)
}

## WT
for(i in seq_along(splittedTablesWT))
{
	print(i)
	
	splitTmp <- splittedTablesWT[[i]][,"V2"]>0.5
	modTmp <- splittedTablesWT[[i]][splitTmp,]
	unmodTmp <- splittedTablesWT[[i]][!splitTmp,]

	write.table(x=modTmp
			  , file=paste0("WT_shuffled/nanom6A_step2/",names(splittedTablesWT)[[i]],".shuffledmod")
			  , quote=FALSE
			  , sep="\t"
			  , row.names=FALSE
			  , col.names=FALSE)
	
	write.table(x=unmodTmp
			  , file=paste0("WT_shuffled/nanom6A_step2/",names(splittedTablesWT)[[i]],".shuffledunmod")
			  , quote=FALSE
			  , sep="\t"
			  , row.names=FALSE
			  , col.names=FALSE)
}

#### For a subset of the available 5-mers (WT = KO)
### Mod and unmod tables
## KO
filesKO <- list.files("/KO/nanom6A_step2/")
modFilesKO <- filesKO[grep("[.]mod",filesKO)]
unmodFilesKO <- filesKO[grep("[.]unmod",filesKO)]
filesKO <- c(modFilesKO, unmodFilesKO)

## WT
filesWT <- list.files("/WT_equal_KO/nanom6A_step2/")
modFilesWT <- filesWT[grep("[.]mod",filesWT)]
unmodFilesWT <- filesWT[grep("[.]unmod",filesWT)]
filesWT <- c(modFilesWT, unmodFilesWT)

### All data
## KO
tablesKO <- lapply(filesKO,function(file)
{
	print(file)
	read.table(paste0("/KO/nanom6A_step2/",file),sep="\t",header=FALSE)
})
elementsKO <- sapply(tablesKO,nrow)
names(elementsKO) <- filesKO
elementsKO <- sapply(names(elementsKO),function(k)
{
	rep(k,elementsKO[[k]])
})
elementsKO <- unlist(elementsKO)
elementsKO <- sapply(strsplit(elementsKO,"[.]"),"[[",1)

tablesKO <- do.call(rbind,tablesKO)

## WT
tablesWT <- lapply(filesWT,function(file)
{
	print(file)
	read.table(paste0("/WT_equal_KO/nanom6A_step2/",file),sep="\t",header=FALSE)
})
elementsWT <- sapply(tablesWT,nrow)
names(elementsWT) <- filesWT
elementsWT <- sapply(names(elementsWT),function(k)
{
	rep(k,elementsWT[[k]])
})
elementsWT <- unlist(elementsWT)
elementsWT <- sapply(strsplit(elementsWT,"[.]"),"[[",1)

tablesWT <- do.call(rbind,tablesWT)

### All probabilities shuffling
shuffledProbabilities <- c(tablesKO[,"V1"],tablesWT[,"V1"])
set.seed(1)
shuffledProbabilities <- shuffledProbabilities[sample(seq_along(shuffledProbabilities),size=length(shuffledProbabilities))]

### Subdivision
## KO vs WT
shuffledProbabilitiesKO <- shuffledProbabilities[1:nrow(tablesKO)]
tablesKO[,"V1"] <- shuffledProbabilitiesKO
tablesKO[,"V2"] <- 1 - tablesKO[,"V1"]

shuffledProbabilitiesWT <- shuffledProbabilities[nrow(tablesKO)+(1:nrow(tablesWT))]
tablesWT[,"V1"] <- shuffledProbabilitiesWT
tablesWT[,"V2"] <- 1 - tablesWT[,"V1"]

## By k-mer
# KO
splittedTablesKO <- split(1:nrow(tablesKO),elementsKO)
splittedTablesKO <- lapply(splittedTablesKO,function(i)tablesKO[i,])

# WT
splittedTablesWT <- split(1:nrow(tablesWT),elementsWT)
splittedTablesWT <- lapply(splittedTablesWT,function(i)tablesWT[i,])

### Substitution and saving
## KO
for(i in seq_along(splittedTablesKO))
{
	print(i)

	splitTmp <- splittedTablesKO[[i]][,"V2"]>0.5
	modTmp <- splittedTablesKO[[i]][splitTmp,]
	unmodTmp <- splittedTablesKO[[i]][!splitTmp,]

	write.table(x=modTmp
			  , file=paste0("/KO_equal_shuffled/nanom6A_step2/",names(splittedTablesKO)[[i]],".shuffledequalmod")
			  , quote=FALSE
			  , sep="\t"
			  , row.names=FALSE
			  , col.names=FALSE)
	
	write.table(x=unmodTmp
			  , file=paste0("/KO_equal_shuffled/nanom6A_step2/",names(splittedTablesKO)[[i]],".shuffledequalunmod")
			  , quote=FALSE
			  , sep="\t"
			  , row.names=FALSE
			  , col.names=FALSE)
}

## WT
for(i in seq_along(splittedTablesWT))
{
	print(i)
	
	splitTmp <- splittedTablesWT[[i]][,"V2"]>0.5
	modTmp <- splittedTablesWT[[i]][splitTmp,]
	unmodTmp <- splittedTablesWT[[i]][!splitTmp,]

	write.table(x=modTmp
			  , file=paste0("/WT_equal_shuffled/nanom6A_step2/",names(splittedTablesWT)[[i]],".shuffledequalmod")
			  , quote=FALSE
			  , sep="\t"
			  , row.names=FALSE
			  , col.names=FALSE)
	
	write.table(x=unmodTmp
			  , file=paste0("/WT_equal_shuffled/nanom6A_step2/",names(splittedTablesWT)[[i]],".shuffledequalunmod")
			  , quote=FALSE
			  , sep="\t"
			  , row.names=FALSE
			  , col.names=FALSE)
}