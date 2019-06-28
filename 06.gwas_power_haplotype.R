### comparison of single variant-based and haplotype-based GWAS : Figure 3 ###

## libraries ##
library(rrBLUP)
library(pROC)
source("inhouse.functions.R")

## data ##
load("data/geno.train.rda")
load("data/haplo.rda")
load("data/K.train.rda")
genetic.map = read.csv("data/genetic.map.csv")

## execution ##
Ncores = 2
cond.qtl = c(10, 25)
cond.h2 = c(0.3, 0.6)
n.rep = 25

n.hap = apply(haplo[,-(1:4)], 1, max, na.rm=TRUE)
bi.loci = names(n.hap[n.hap == 2])
bi.loci = haplo[is.element(rownames(haplo), bi.loci),]

H = t(bi.loci[,-(1:4)])
info = bi.loci[,1:4]

set.seed(1)
for (n.qtl in cond.qtl) {
	for (h2 in cond.h2) {
		for (i in 1:n.rep) {
			filename = paste("simQTLbi",n.qtl,h2*10,i,".rda",sep="_")
			simqtl = sim.maQTL(H, info, K.train, h2, n.qtl)
			save(simqtl, file=filename)
		}
	}
}

# n.qtl=10; h2=0.3; ct=1
# n.qtl=10; h2=0.6; ct=2
# n.qtl=25; h2=0.3; ct=3
# n.qtl=25; h2=0.6; ct=4

ct = 1
list.auc = vector(mode="list", length=4)
for (n.qtl in cond.qtl) {
	for (h2 in cond.h2) {
		exp.id = paste(n.qtl, h2*10, sep="_")
		names(list.auc)[[ct]] = exp.id
		mat = c()
		for (i in 1:n.rep) {
			filename = paste("simQTLbi",n.qtl,h2*10,i,".rda",sep="_")
			load(filename)
			res = GWAS(data.frame(gid=rownames(simqtl$val), y=simqtl$val$y), geno.train[,-(4:14)], fixed=NULL, K=K.train, n.PC=3, min.MAF=0.05, n.core=Ncores, P3D=TRUE, plot=FALSE)
			logP = res$y
			tru = c()
			qtl = simqtl$qtl
			for (j in 1:nrow(qtl)) {
				q = as.character(qtl$MARKER[j])
				ix = strsplit(q, "_")[[1]][1]
				if (ix=="Prom") {
					tru = c(tru, as.character(geno.train$MARKER[geno.train$PROM==q]))
				} else if (ix=="Down") {
					tru = c(tru, as.character(geno.train$MARKER[geno.train$DOWN==q]))
				} else {
					tru = c(tru, as.character(geno.train$MARKER[geno.train$GENE==q]))		
				}
			}
			Causative = res$MARKER %in% tru
			sv.roc = roc(Causative, logP, algorithm=2, quiet=FALSE, levels=c("FALSE", "TRUE"), direction="<")
			sv.auc = as.numeric(sv.roc$auc)
			hre = hapLMM(data.frame(gid=rownames(simqtl$val), y=simqtl$val$y), "y", haplo, K=K.train, fixed=NULL, n.PC=3)
			mht = hre$mht
			mht = mht[!is.na(mht$Score),]
			logP = -log10(mht$Score)
			Causative = mht$MARKER %in% simqtl$qtl$MARKER
			hap.roc = roc(Causative, logP, algorithm=2, quiet=FALSE, levels=c("FALSE", "TRUE"), direction="<")
			hap.auc = as.numeric(hap.roc$auc)
			mat = rbind(mat, c(sv.auc, hap.auc))
		}
		list.auc[[ct]] = mat
		ct = ct + 1	
	}
}

L = list.auc
M <- c()
for (i in 1:4) {
	x <- L[[i]]
	msv <- data.frame(
		cond=rep(names(L)[i], n.rep),
		GWAS=rep("SV", n.rep),
		AUC=x[,1]
	)
	mha <- data.frame(
		cond=rep(names(L)[i], n.rep),
		GWAS=rep("Hap", n.rep),
		AUC=x[,2]
	)
	M <- rbind(M,msv,mha)
}
gp <- ggplot(M, aes(y=AUC,x=GWAS, fill=GWAS)) + 
	geom_boxplot() + 
	labs(x="GWAS method",y="AUC") + 
	theme(axis.title=element_text(size=9)) + 
	theme(axis.text=element_text(size=9)) + 
	theme(strip.text=element_text(size=9,face="bold")) + 
	scale_fill_discrete("") + 
	theme(legend.position = 'none') +
	facet_grid(~cond, labeller=as_labeller(
		c("10_3"="QTL10 : 0.3",
			"10_6"="QTL10 : 0.6",
			"25_3"="QTL25 : 0.3",
			"25_6"="QTL25 : 0.6")))
ggsave(file="outputs/Fig3B.png", plot=gp, dpi=400, width = 4, height = 1.5)



n.hap = apply(haplo[,-(1:4)], 1, max, na.rm=TRUE)
ma.loci = names(n.hap[n.hap > 2])
ma.loci = haplo[is.element(rownames(haplo), ma.loci),]

H = t(ma.loci[,-(1:4)])
info = ma.loci[,1:4]

for (n.qtl in cond.qtl) {
	for (h2 in cond.h2) {
		for (i in 1:n.rep) {
			filename = paste("simQTLma",n.qtl,h2*10,i,".rda",sep="_")
			simqtl = sim.maQTL(H, info, K.train, h2, n.qtl)
			save(simqtl, file=filename)
		}
	}
}

# n.qtl=10; h2=0.3; ct=1
# n.qtl=10; h2=0.6; ct=2
# n.qtl=25; h2=0.3; ct=3
# n.qtl=25; h2=0.6; ct=4

ct = 1
list.auc = vector(mode="list", length=4)
for (n.qtl in cond.qtl) {
	for (h2 in cond.h2) {
		exp.id = paste(n.qtl, h2*10, sep="_")
		names(list.auc)[[ct]] = exp.id
		mat = c()
		for (i in 1:n.rep) {
			filename = paste("simQTLma",n.qtl,h2*10,i,".rda",sep="_")
			load(filename)
			res = GWAS(data.frame(gid=rownames(simqtl$val), y=simqtl$val$y), geno.train[,-(4:14)], fixed=NULL, K=K.train, n.PC=3, min.MAF=0.05, n.core=Ncores, P3D=TRUE, plot=FALSE)
			logP = res$y
			tru = c()
			qtl = simqtl$qtl
			for (j in 1:nrow(qtl)) {
				q = as.character(qtl$MARKER[j])
				ix = strsplit(q, "_")[[1]][1]
				if (ix=="Prom") {
					tru = c(tru, as.character(geno.train$MARKER[geno.train$PROM==q]))
				} else if (ix=="Down") {
					tru = c(tru, as.character(geno.train$MARKER[geno.train$DOWN==q]))
				} else {
					tru = c(tru, as.character(geno.train$MARKER[geno.train$GENE==q]))		
				}
			}
			Causative = res$MARKER %in% tru
			sv.roc = roc(Causative, logP, algorithm=2, quiet=FALSE, levels=c("FALSE", "TRUE"), direction="<")
			sv.auc = as.numeric(sv.roc$auc)
			hre = hapLMM(data.frame(gid=rownames(simqtl$val), y=simqtl$val$y), "y", haplo, K=K.train, fixed=NULL, n.PC=3)
			mht = hre$mht
			mht = mht[!is.na(mht$Score),]
			logP = -log10(mht$Score)
			Causative = mht$MARKER %in% simqtl$qtl$MARKER
			hap.roc = roc(Causative, logP, algorithm=2, quiet=FALSE, levels=c("FALSE", "TRUE"), direction="<")
			hap.auc = as.numeric(hap.roc$auc)
			mat = rbind(mat, c(sv.auc, hap.auc))
		}
		list.auc[[ct]] = mat
		ct = ct + 1	
	}
}

L = list.auc
M <- c()
for (i in 1:4) {
	x <- L[[i]]
	msv <- data.frame(
		cond=rep(names(L)[i], n.rep),
		GWAS=rep("SV", n.rep),
		AUC=x[,1]
	)
	mha <- data.frame(
		cond=rep(names(L)[i], n.rep),
		GWAS=rep("Hap", n.rep),
		AUC=x[,2]
	)
	M <- rbind(M,msv,mha)
}
gp <- ggplot(M, aes(y=AUC,x=GWAS, fill=GWAS)) + 
	geom_boxplot() + 
	labs(x="GWAS method",y="AUC") + 
	theme(axis.title=element_text(size=9)) + 
	theme(axis.text=element_text(size=9)) + 
	theme(strip.text=element_text(size=9,face="bold")) + 
	scale_fill_discrete("") + 
	theme(legend.position = 'none') +
	facet_grid(~cond, labeller=as_labeller(
		c("10_3"="QTL10 : 0.3",
			"10_6"="QTL10 : 0.6",
			"25_3"="QTL25 : 0.3",
			"25_6"="QTL25 : 0.6")))
ggsave(file="outputs/Fig3C.png", plot=gp, dpi=400, width = 4, height = 1.5)
