### Libraries
library(seqinr)

### Custom functions
barplotText <- function(x,...)
{
	xTmp <- barplot(x,...)
	yTmp <- apply(x,2,function(c)c(c[1]/2,c[1]+c[2]/2,c[1]+c[2]+c[3]/2))
	
	for(i in seq_along(xTmp))
	{
		text(xTmp[[i]],yTmp[,i],round(x[,i],2),col="white")	
	}	
}

### Nanom6A ratios
## ratio_WT
ratio_WT = read.table("WT/nanom6A_step2/ratio.0.5.tsv",header=FALSE,sep="\t",fill=TRUE,col.names=1:100)

## ratio_KO
ratio_KO = read.table("KO/nanom6A_step2/ratio.0.5.tsv",header=FALSE,sep="\t",fill=TRUE,col.names=1:100)

## ratio_WT_shuffled
ratio_WT_shuffled = read.table("WT_shuffled/nanom6A_step2/ratio.0.5.tsv",header=FALSE,sep="\t",fill=TRUE,col.names=1:100)

## ratio_KO_shuffled
ratio_KO_shuffled = read.table("KO_shuffled/nanom6A_step2/ratio.0.5.tsv",header=FALSE,sep="\t",fill=TRUE,col.names=1:100)

## ratio_WT_equal
ratio_WT_equal = read.table("WT_equal_KO/nanom6A_step2/ratio.0.5.tsv",header=FALSE,sep="\t",fill=TRUE,col.names=1:100)

## ratio_WT_equal_shuffled
ratio_WT_equal_shuffled = read.table("WT_equal_shuffled/nanom6A_step2/ratio.0.5.tsv",header=FALSE,sep="\t",fill=TRUE,col.names=1:100)

## ratio_KO_equal_shuffled
ratio_KO_equal_shuffled = read.table("KO_equal_shuffled/nanom6A_step2/ratio.0.5.tsv",header=FALSE,sep="\t",fill=TRUE,col.names=1:100)

## ratio_WT_smaller
ratio_WT_smaller = read.table("WT_smaller_KO/nanom6A_step2/ratio.0.5.tsv",header=FALSE,sep="\t",fill=TRUE,col.names=1:100)

### m6A sites
## WT
m6Asites_WT = unique(unlist(apply(ratio_WT,1,function(r)paste0(r[[1]],sapply(sapply(r[-1][grepl("[|]",r[-1])],function(i)strsplit(i,"[|]")),"[[",1)))))

## KO
m6Asites_KO = unique(unlist(apply(ratio_KO,1,function(r)paste0(r[[1]],sapply(sapply(r[-1][grepl("[|]",r[-1])],function(i)strsplit(i,"[|]")),"[[",1)))))

## WT_shuffled
m6Asites_WT_shuffled = unique(unlist(apply(ratio_WT_shuffled,1,function(r)paste0(r[[1]],sapply(sapply(r[-1][grepl("[|]",r[-1])],function(i)strsplit(i,"[|]")),"[[",1)))))

## KO_shuffled
m6Asites_KO_shuffled = unique(unlist(apply(ratio_KO_shuffled,1,function(r)paste0(r[[1]],sapply(sapply(r[-1][grepl("[|]",r[-1])],function(i)strsplit(i,"[|]")),"[[",1)))))

## WT_equal
m6Asites_WT_equal = unique(unlist(apply(ratio_WT_equal,1,function(r)paste0(r[[1]],sapply(sapply(r[-1][grepl("[|]",r[-1])],function(i)strsplit(i,"[|]")),"[[",1)))))

## WT_equal_shuffled
m6Asites_WT_equal_shuffled = unique(unlist(apply(ratio_WT_equal_shuffled,1,function(r)paste0(r[[1]],sapply(sapply(r[-1][grepl("[|]",r[-1])],function(i)strsplit(i,"[|]")),"[[",1)))))

## KO_equal_shuffled
m6Asites_KO_equal_shuffled = unique(unlist(apply(ratio_KO_equal_shuffled,1,function(r)paste0(r[[1]],sapply(sapply(r[-1][grepl("[|]",r[-1])],function(i)strsplit(i,"[|]")),"[[",1)))))

## WT_smaller
m6Asites_WT_smaller = unique(unlist(apply(ratio_WT_smaller,1,function(r)paste0(r[[1]],sapply(sapply(r[-1][grepl("[|]",r[-1])],function(i)strsplit(i,"[|]")),"[[",1)))))

## WT vs KO
m6Asites_WT_vs_KO = c(table(union(m6Asites_WT,m6Asites_KO)%in%m6Asites_WT & union(m6Asites_WT,m6Asites_KO)%in%m6Asites_KO)[["TRUE"]]
				  ,table(union(m6Asites_WT,m6Asites_KO)%in%m6Asites_WT & !union(m6Asites_WT,m6Asites_KO)%in%m6Asites_KO)[["TRUE"]]
				  ,table(!union(m6Asites_WT,m6Asites_KO)%in%m6Asites_WT & union(m6Asites_WT,m6Asites_KO)%in%m6Asites_KO)[["TRUE"]])

## WT vs KO shuffled
m6Asites_WT_vs_KO_shuffled = c(table(union(m6Asites_WT_shuffled,m6Asites_KO_shuffled)%in%m6Asites_WT_shuffled & union(m6Asites_WT_shuffled,m6Asites_KO_shuffled)%in%m6Asites_KO_shuffled)[["TRUE"]]
							  ,table(union(m6Asites_WT_shuffled,m6Asites_KO_shuffled)%in%m6Asites_WT_shuffled & !union(m6Asites_WT_shuffled,m6Asites_KO_shuffled)%in%m6Asites_KO_shuffled)[["TRUE"]]
							  ,table(!union(m6Asites_WT_shuffled,m6Asites_KO_shuffled)%in%m6Asites_WT_shuffled & union(m6Asites_WT_shuffled,m6Asites_KO_shuffled)%in%m6Asites_KO_shuffled)[["TRUE"]])

## WT_equal vs KO
m6Asites_WT_equal_vs_KO = c(table(union(m6Asites_WT_equal,m6Asites_KO)%in%m6Asites_WT_equal & union(m6Asites_WT_equal,m6Asites_KO)%in%m6Asites_KO)[["TRUE"]]
				  ,table(union(m6Asites_WT_equal,m6Asites_KO)%in%m6Asites_WT_equal & !union(m6Asites_WT_equal,m6Asites_KO)%in%m6Asites_KO)[["TRUE"]]
				  ,table(!union(m6Asites_WT_equal,m6Asites_KO)%in%m6Asites_WT_equal & union(m6Asites_WT_equal,m6Asites_KO)%in%m6Asites_KO)[["TRUE"]])

## WT vs KO equal shuffled
m6Asites_WT_vs_KO_equal_shuffled = c(table(union(m6Asites_WT_equal_shuffled,m6Asites_KO_equal_shuffled)%in%m6Asites_WT_equal_shuffled & union(m6Asites_WT_equal_shuffled,m6Asites_KO_equal_shuffled)%in%m6Asites_KO_equal_shuffled)[["TRUE"]]
							  ,table(union(m6Asites_WT_equal_shuffled,m6Asites_KO_equal_shuffled)%in%m6Asites_WT_equal_shuffled & !union(m6Asites_WT_equal_shuffled,m6Asites_KO_equal_shuffled)%in%m6Asites_KO_equal_shuffled)[["TRUE"]]
							  ,table(!union(m6Asites_WT_equal_shuffled,m6Asites_KO_equal_shuffled)%in%m6Asites_WT_equal_shuffled & union(m6Asites_WT_equal_shuffled,m6Asites_KO_equal_shuffled)%in%m6Asites_KO_equal_shuffled)[["TRUE"]])

## WT_smaller vs KO
m6Asites_WT_smaller_vs_KO = c(table(union(m6Asites_WT_smaller,m6Asites_KO)%in%m6Asites_WT_smaller & union(m6Asites_WT_smaller,m6Asites_KO)%in%m6Asites_KO)[["TRUE"]]
				  ,table(union(m6Asites_WT_smaller,m6Asites_KO)%in%m6Asites_WT_smaller & !union(m6Asites_WT_smaller,m6Asites_KO)%in%m6Asites_KO)[["TRUE"]]
				  ,table(!union(m6Asites_WT_smaller,m6Asites_KO)%in%m6Asites_WT_smaller & union(m6Asites_WT_smaller,m6Asites_KO)%in%m6Asites_KO)[["TRUE"]])

## Plot
x = cbind("WT_vs_KO"=m6Asites_WT_vs_KO,"m6Asites_WT_vs_KO_shuffled"=m6Asites_WT_vs_KO_shuffled,"WT_equal_vs_KO"=m6Asites_WT_equal_vs_KO,"m6Asites_WT_vs_KO_equal_shuffled"=m6Asites_WT_vs_KO_equal_shuffled,"WT_smaller_vs_KO"=m6Asites_WT_smaller_vs_KO)
x = cbind(x,rep(NaN,3))
rownames(x) = c("Common","WT_only","KD_only")

pdf("comparativePlot.pdf",width=15,height=15)
par(mfrow=c(2,1))
barplotText(apply(x,2,function(i)i/sum(i)),ylab="Fraction of sites",col=c("black","darkgrey","lightgrey"),las=2)
legend("topright",pch=16,col=c("black","darkgrey","lightgrey"),legend=c("Common","WT only","KO only"))
dev.off()