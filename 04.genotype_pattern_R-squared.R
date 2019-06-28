### genotype pattern R-squared : Table 1 ###

## libraries ##
source("inhouse.functions.R")

## data ##
load("data/geno.train.rda")
genetic.map = read.csv("data/genetic.map.csv")

## execution ##
geno.train = geno.train[geno.train$CHROM != "SL3.0ch00",]	
geno.train = loessMapPos(genetic.map, geno.train)

chroms = sort(as.character(unique(geno.train$CHROM)))
Res = matrix(NA, nr=length(chroms), nc=3)
rownames(Res) = chroms
colnames(Res) = c("5cM","10cM","20cM")

for (i in 1:length(chroms)) {
	g.in = geno.train[geno.train$CHROM==chroms[i],]
	bin = seq(0,max(g.in$cM),by=20)
	n.bin = length(bin)-1
	y = c()
	x = c()
	for (j in 1:n.bin) {
		print(paste("CHROM = ",chroms[i], " : ", j, "/", n.bin, sep=""))
		g = g.in[(g.in$cM >= bin[j]) & (g.in$cM < bin[j+1]),]
		if (nrow(g)==0) next
		r = as.matrix(g[,-c(1:14)])
		r = genoCorCoef(r,r)
		d = g$cM
		d = abs(outer(d, d, "-"))
		res = as.data.frame(cbind(d[upper.tri(d)],r[upper.tri(r)]))
		res = as.data.frame(res)
		colnames(res) = c("cM", "r2")	
		res$r2 = res$r2^2
		res$cM = round(res$cM)
		y.b = rep(NA, 21)
		for (k in 0:20) y.b[(k+1)] = mean(res$r2[res$cM==k])
		y = c(y, y.b)
		x = c(x, seq(0,20))
		rm(r,g,d,res); gc()
	}
	fit = loess(y ~ x, span=1)
	Res[i,] = predict(fit, c(5,10,20))
}
write.csv(Res, "outputs/Table1.csv")
