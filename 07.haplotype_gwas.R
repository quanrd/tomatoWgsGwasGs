### haplotype-based GWAS : Figure 4 ###

## libraries ##
library(rrBLUP)
source("inhouse.functions.R")

## data ##
load("data/haplo.rda")
load("data/K.train.rda")
load("data/pheno.train.rda")
load("data/geno.train.rda")


### execution ###
thr.level = 5.393751
phenos = colnames(pheno.train)[-c(1,2)]

for (phenoi in phenos) {
	res = hapLMM(pheno.train, phenoi, haplo, K=K.train, fixed="year", n.PC=3)
	outfile = paste("outputs/fig4_", phenoi, "_res.rda", sep="")
	save(res, file=outfile)
	load(outfile)
	x = res$mht
	x = x[!is.na(x$Score),]
	x$Score = -log10(x$Score)
	x = x[x$Score > thr.level,]
	chroms = sort(unique(x$CHROM))
	haplo.list = vector(mode="list",length=length(chroms))
	names(haplo.list) = chroms
	for (i in 1:length(chroms)) {
		chrom = chroms[i]
		xi = x[x$CHROM==chrom,]
		if (nrow(xi)==1) haplo.list[[i]] = c(xi$MARKER, xi$Sig.Hap)
		if (nrow(xi)==1) print(c(xi$MARKER, xi$Sig.Hap))
		if (nrow(xi)==1) next
		haploi = haplo[is.element(haplo$MARKER,xi$MARKER),,drop=FALSE]
		hap.cov = list(); ct=1			
		fixq = xi[xi$Score==max(xi$Score),][1,]
		fixq = c(fixq$MARKER, fixq$Sig.Hap)
		print(fixq)
		hap.cov[[ct]] = fixq
		resi = hapLMM(pheno.train, phenoi, haploi, K.train, "year", 3, hap.cov=hap.cov)
		yi = resi$mht
		yi$Score = -log10(yi$Score)
		maxP = max(yi$Score)
		while (maxP > thr.level) {
			fixq = yi[yi$Score==maxP,][1,]
			fixq = c(fixq$MARKER, fixq$Sig.Hap)
			if (fixq[1]==hap.cov[[ct]][1]) break
			ct = ct + 1
			print(fixq)
			hap.cov[[ct]] = fixq
			resa = hapLMM(pheno.train, phenoi, haploi, K.train, "year", 3, hap.cov=hap.cov)
			yi = resa$mht
			yi$Score = -log10(yi$Score)
			maxP = max(yi$Score)
		}
		haplo.list[[i]] = hap.cov
	}
	outfile = paste("outputs/fig4_", phenoi, "_haplo.list.rda", sep="")
	save(haplo.list, file=outfile)
	sel.mrks = c()
	for (j in 1:length(haplo.list)) {
		g = haplo.list[[j]][[1]][1]
		vec = numeric(nrow(geno.train))
		vec[is.element(as.character(geno.train$PROM), g)] = 1
		vec[is.element(as.character(geno.train$GENE), g)] = 1
		vec[is.element(as.character(geno.train$DOWN), g)] = 1
		gg = geno.train[vec==1,,drop=FALSE]
		nr = nrow(gg)
		sel.mrks = c(sel.mrks, as.character(gg$MARKER)[1])
		while (nr > 1) {
			X = as.matrix(gg[,-(1:14)]) + 1
			Mcor = genoCorCoef(X, X)
			rmv = (1 -Mcor[,1]) > 0.3
			gg = gg[rmv,,drop=FALSE]
			if (nrow(gg) > 0) sel.mrks = c(sel.mrks, as.character(gg$MARKER)[1])
			nr = nrow(gg)
		}
	} #j
	outfile = paste("outputs/fig4_", phenoi, "_sel.mrks.rda", sep="")
	save(sel.mrks, file=outfile)
}

