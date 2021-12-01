### Libraries
library(pROC)

### Classification data
oligo1_GAACT.mod <- read.table("oligo1/nanom6A_step2/GAACT.mod",sep="\t",header=FALSE)
oligo1_GAACT.unmod <- read.table("oligo1/nanom6A_step2/GAACT.unmod",sep="\t",header=FALSE)
oligo1_AGACA.mod <- read.table("oligo1/nanom6A_step2/AGACA.mod",sep="\t",header=FALSE)
oligo1_AGACA.unmod <- read.table("oligo1/nanom6A_step2/AGACA.unmod",sep="\t",header=FALSE)
oligo1_GAACC.mod <- read.table("oligo1/nanom6A_step2/GAACC.mod",sep="\t",header=FALSE)
oligo1_GAACC.unmod <- read.table("oligo1/nanom6A_step2/GAACC.unmod",sep="\t",header=FALSE)
oligo1_GGACT.mod <- read.table("oligo1/nanom6A_step2/GGACT.mod",sep="\t",header=FALSE)
oligo1_GGACT.unmod <- read.table("oligo1/nanom6A_step2/GGACT.unmod",sep="\t",header=FALSE)

control_GAACT.mod <- read.table("control/nanom6A_step2/GAACT.mod",sep="\t",header=FALSE)
control_GAACT.unmod <- read.table("control/nanom6A_step2/GAACT.unmod",sep="\t",header=FALSE)
control_AGACA.mod <- read.table("control/nanom6A_step2/AGACA.mod",sep="\t",header=FALSE)
control_AGACA.unmod <- read.table("control/nanom6A_step2/AGACA.unmod",sep="\t",header=FALSE)
control_GAACC.mod <- read.table("control/nanom6A_step2/GAACC.mod",sep="\t",header=FALSE)
control_GAACC.unmod <- read.table("control/nanom6A_step2/GAACC.unmod",sep="\t",header=FALSE)
control_GGACT.mod <- read.table("control/nanom6A_step2/GGACT.mod",sep="\t",header=FALSE)
control_GGACT.unmod <- read.table("control/nanom6A_step2/GGACT.unmod",sep="\t",header=FALSE)

### ROCs
pdf("rocCurve.pdf",height=6,width=6)

## TRUEs (m6A+ RRACH motifs)
TRUEs <- c(oligo1_GGACT.mod[,2]
		  ,oligo1_GGACT.unmod[,2])

## FALSEs (m6A- RRACH motifs)
FALSEs <- c(oligo1_GAACT.mod[,2]
		   ,oligo1_GAACT.unmod[,2]
		   ,oligo1_AGACA.mod[,2]
		   ,oligo1_AGACA.unmod[,2]
		   ,oligo1_GAACC.mod[,2]
		   ,oligo1_GAACC.unmod[,2]
		   ,control_GAACT.mod[,2]
		   ,control_GAACT.unmod[,2]
		   ,control_AGACA.mod[,2]
		   ,control_AGACA.unmod[,2]
		   ,control_GAACC.mod[,2]
		   ,control_GAACC.unmod[,2]
		   ,control_GGACT.mod[,2]
		   ,control_GGACT.unmod[,2])

pmods_all <- c(TRUEs,FALSEs)
responses_all <- c(rep(1,length(TRUEs)),rep(0,length(FALSEs)))
roc_all <- roc(responses_all,pmods_all)
plot(roc_all,col=1,lwd=3)

legend("bottomright"
	  ,col=1
	  ,lwd=3
	  ,legend=c(paste0("AUC Oligo1+Control: ",round(roc_all$auc[[1]],2)," (2136 vs 28906)")))

dev.off()
