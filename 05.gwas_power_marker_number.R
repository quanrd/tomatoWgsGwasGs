### number of markers and GWAS power : Figure 3A ###

## libraries ##
library(rrBLUP)
library(pROC)
library(ggplot2)
source("inhouse.functions.R")

## data ##
load("data/geno.train.rda")
load("data/K.train.rda")
genetic.map = read.csv("data/genetic.map.csv")

## execution ##
geno.train = geno.train[geno.train$CHROM != "SL3.0ch00",]
geno.train = loessMapPos(genetic.map, geno.train)

G = t(geno.train[,-(1:14)])
info = geno.train[,1:14]

cond.qtl = c(10, 25)
cond.h2 = c(0.3, 0.6)
n.rep = 25
Ncore = 10
n.p = c(5000, 50000, 100000, 500000)

set.seed(1)
for (n.qtl in cond.qtl) {
	for (h2 in cond.h2) {
		for (i in 1:n.rep) {
			filename = paste("simQTLnmrks",n.qtl,h2*10,i,".rda",sep="_")
			simqtl = sim.QTL(G, info, K.train, h2, n.qtl)
			save(simqtl, file=filename)
		}
	}
}

Nmrks_list = vector(mode="list", length=4)
ct = 1
set.seed(1)
for (n.qtl in cond.qtl) {
	for (h2 in cond.h2) {
		AUC = c()
		Nmrks = c()
		Sim = c()
		names(Nmrks_list)[[ct]] = paste("Nqtl.",n.qtl,"_h2.",h2,sep="")
		for (i in 1:n.rep) {
			filename = paste("simQTLnmrks",n.qtl,h2*10,i,".rda",sep="_")
			load(filename)
			phenoIN = data.frame(gid=rownames(simqtl$val), y=simqtl$val$y)
			geno.base = geno.train[!is.element(geno.train$MARKER, simqtl$qtl$MARKER),]
			for (p in n.p) {
				gi = geno.base[sort(sample(1:nrow(geno.base), p)),]
				Ki = A.mat(t(gi[,-(1:14)]), n.core=Ncore)
				res = GWAS(phenoIN, gi[,-(4:14)], fixed=NULL, K=Ki, min.MAF=0.05, n.core=Ncore, P3D=TRUE, plot=FALSE)
				ans = reconstruct.res(res, simqtl$qtl)
				logP = res$y
				Causative = res$MARKER %in% ans[,1]
				not.detected = n.qtl - length(Causative[Causative])
				Causative = c(Causative, rep(TRUE, not.detected))
				logP = c(logP, rep(0, not.detected))
				val = roc(Causative, logP, algorithm=2, quiet=FALSE, levels=c("FALSE", "TRUE"), direction="<")$auc
				AUC = c(AUC, as.numeric(val))
				Nmrks = c(Nmrks, p)
				Sim = c(Sim, i)
			}
		}
		M = data.frame(AUC=AUC, N.markers=Nmrks, Sim=Sim)
		M$Sim = as.factor(M$Sim)
		M$N.markers = as.numeric(as.factor(M$N.markers))
		Nmrks_list[[ct]] = M
		ct = ct + 1
	}
}


M = c()
for (i in 1:length(Nmrks_list)) {
	cond.names = names(Nmrks_list)[i]
	Cond = rep(cond.names, 100)
	M = rbind(M, cbind(Cond, Nmrks_list[[i]]))
}
M$N.markers = as.numeric(as.factor(M$N.markers))
M$Sim = as.character(as.numeric(M$Sim))
xx = c(5,50,100,500)
for (i in 1:4) M$N.markers[M$N.markers==i] = xx[i]
M$N.markers = as.factor(M$N.markers)
g = ggplot(M, aes(x=N.markers, y=AUC, color=N.markers)) +
		geom_boxplot() +
		facet_grid( ~ Cond) +
		theme(legend.position='none')
ggsave(file = "outputs/Fig3A.png", plot = g, dpi = 400, width = 7, height = 2.5)

