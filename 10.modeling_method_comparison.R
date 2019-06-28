### GS modeling method comparison : Table 3 ###

## libraries ##
library(rrBLUP)
library(VIGoR)
library(randomForest)

## data ##
load("data/pheno.test.rda")
load("data/pheno.train.rda")
load("data/geno.test.rda")
load("data/geno.train.rda")

## execution ##
g.train = t(geno.train[,-(1:14)]); colnames(g.train) = geno.train$MARKER
g.test = t(geno.test[,-(1:14)]); colnames(g.test) = geno.test$MARKER

n = nrow(pheno.train)
Z = matrix(0,n,nrow(g.train))
rownames(Z) = pheno.train[,1]
colnames(Z) = rownames(g.train)
Z[cbind(1:n,match(pheno.train[,1],rownames(g.train)))] = 1
	
fixed.item = unique(pheno.train$year)
X = matrix(0, nr=nrow(pheno.train), nc=length(fixed.item))
colnames(X) = fixed.item
for (i in 1:length(fixed.item)) X[pheno.train$year==fixed.item[i],i] = 1

phenos = colnames(pheno.train)[-c(1,2)]

v.na = rep(NA, 1000)
Res = data.frame(
	method=v.na,
	file=v.na,
	exp=v.na,
	trait=v.na,
	r=v.na
)
test.year = as.character(unique(pheno.test$year))
proc = 1

for (p in 1:length(test.year)) {
	print(test.year[p])
	p.test = pheno.test[pheno.test$year==test.year[p],]
	
	for (i in 1:4) {
		set.seed(1)
		phenoi = phenos[i]
		sel.mrk.name = paste("outputs/fig4_",phenoi,"_sel.mrks.rda",sep="")
		load(sel.mrk.name)
		G.train = g.train[,is.element(colnames(g.train), sel.mrks)]
		ZG.train = Z%*%as.matrix(G.train)
		y = pheno.train[[phenoi]]
		use = !is.na(y)
		Xuse = X[use,]; Xuse = Xuse[,apply(Xuse, 2, sum)!=0]			
		fm = randomForest(x=cbind(Xuse, ZG.train[use,]), y=y[use])
		Xtrain = matrix(0, nr=nrow(G.train), nc=ncol(Xuse)); colnames(Xtrain) = colnames(Xuse)
		train = predict(fm, cbind(Xtrain, G.train))
		G.test = g.test[,is.element(colnames(g.test), sel.mrks)]
		Xtest = matrix(0, nr=nrow(G.test), nc=ncol(Xuse)); colnames(Xtest) = colnames(Xuse)
		pred = predict(fm, cbind(Xtest, G.test))
		cmp = data.frame(Line=p.test$gid, obs=rep(NA, nrow(p.test)), train=rep(NA, nrow(p.test)), pred=rep(NA, nrow(p.test)))
		cmp$obs = p.test[[phenoi]]
		for (j in 1:nrow(p.test)) {
			bv = train[as.character(names(train))==as.character(cmp$Line[j])]
			if (length(bv)!=0) cmp$train[j] = bv
		}
		for (j in 1:nrow(p.test)) {
			bv = pred[as.character(names(pred))==as.character(cmp$Line[j])]
			if (length(bv)!=0) cmp$pred[j] = bv
		}
		ct = cor.test(cmp$obs,cmp$train,use="pairwise.complete.obs")
		Res$method[proc] = "RF"
		Res$file[proc] = test.year[p]
		Res$exp[proc] = "train"
		Res$trait[proc] = phenoi
		Res$r[proc] = ct$estimate
		proc = proc + 1

		ct = cor.test(cmp$obs,cmp$pred,use="pairwise.complete.obs")
		Res$method[proc] = "RF"
		Res$file[proc] = test.year[p]
		Res$exp[proc] = "test"
		Res$trait[proc] = phenoi
		Res$r[proc] = ct$estimate
		proc = proc + 1
	}

	hp.BLasso = matrix(c(1,0.001, 1,0.01, 1,0.1, 1,1, 1,2, 1,5),nc=2,byrow=TRUE)
	for (i in 1:4) {
		set.seed(1)
		phenoi = phenos[i]
		sel.mrk.name = paste("outputs/fig4_",phenoi,"_sel.mrks.rda",sep="")
		load(sel.mrk.name)
		G.test = g.test[,is.element(colnames(g.test), sel.mrks)]
		G.train = g.train[,is.element(colnames(g.train), sel.mrks)]
		ZG.train = Z%*%as.matrix(G.train)
		y = pheno.train[[phenoi]]
		use = !is.na(y)
		yuse = y[use]
		ZGuse = ZG.train[use,]
		Xuse = X[use,]; Xuse = Xuse[,apply(Xuse, 2, sum)!=0]			
		res = vigor(Pheno=yuse, Geno=ZGuse, Covariates=Xuse, Method="BL", Hyperparameters=hp.BLasso, Function="tuning")
		Beta = res$Beta
		train = matrix(Beta, nrow=1) %*% t(G.train)
		pred = matrix(Beta, nrow=1) %*% t(G.test)	
		cmp = data.frame(Line=p.test$gid, obs=rep(NA, nrow(p.test)), train=rep(NA, nrow(p.test)), pred=rep(NA, nrow(p.test)))
		cmp$obs = p.test[[phenoi]]
		for (j in 1:nrow(p.test)) {
			bv = train[as.character(colnames(train))==as.character(cmp$Line[j])]
			if (length(bv)!=0) cmp$train[j] = bv
		}
		for (j in 1:nrow(p.test)) {
			bv = pred[as.character(colnames(pred))==as.character(cmp$Line[j])]
			if (length(bv)!=0) cmp$pred[j] = bv
		}

		ct = cor.test(cmp$obs,cmp$train,use="pairwise.complete.obs")
		Res$method[proc] = "BL"
		Res$file[proc] = test.year[p]
		Res$exp[proc] = "train"
		Res$trait[proc] = phenoi
		Res$r[proc] = ct$estimate
		proc = proc + 1

		ct = cor.test(cmp$obs,cmp$pred,use="pairwise.complete.obs")
		Res$method[proc] = "BL"
		Res$file[proc] = test.year[p]
		Res$exp[proc] = "test"
		Res$trait[proc] = phenoi
		Res$r[proc] = ct$estimate
		proc = proc + 1
	} #i

	for (i in 1:4) {
		phenoi = phenos[i]
		sel.mrk.name = paste("outputs/fig4_",phenoi,"_sel.mrks.rda",sep="")
		load(sel.mrk.name)
		G.train = g.train[,is.element(colnames(g.train), sel.mrks)]
		G.test = g.test[,is.element(colnames(g.test), sel.mrks)]
		y = pheno.train[[phenoi]]
		use = !is.na(y)
		y.in = y[use]
		Z.in = Z[use,]
		X.in = X[use,]; X.in = X.in[,apply(X.in, 2, sum)!=0]
		ans = kinship.BLUP(y=y.in, G.train=G.train, G.pred=G.test, Z.train=Z.in, X=X.in, K.method="RR", n.core=10)

		cmp = data.frame(Line=p.test$gid, obs=rep(NA, nrow(p.test)), train=rep(NA, nrow(p.test)), pred=rep(NA, nrow(p.test)))
		cmp$obs = p.test[[phenoi]]
		for (j in 1:nrow(p.test)) {
			bv = ans$g.train[names(ans$g.train)==cmp$Line[j]]
			if (length(bv)!=0) cmp$train[j] = bv
		}
		for (j in 1:nrow(p.test)) {
			bv = ans$g.pred[names(ans$g.pred)==cmp$Line[j]]
			if (length(bv)!=0) cmp$pred[j] = bv
		}

		ct = cor.test(cmp$obs,cmp$train,use="pairwise.complete.obs")
		Res$method[proc] = "RR"
		Res$file[proc] = test.year[p]
		Res$exp[proc] = "train"
		Res$trait[proc] = phenoi
		Res$r[proc] = ct$estimate
		proc = proc + 1

		ct = cor.test(cmp$obs,cmp$pred,use="pairwise.complete.obs")
		Res$method[proc] = "RR"
		Res$file[proc] = test.year[p]
		Res$exp[proc] = "test"
		Res$trait[proc] = phenoi
		Res$r[proc] = ct$estimate
		proc = proc + 1
	}
}

Res = Res[!is.na(Res[,1]),]
v.na = rep(NA, 24)
M = data.frame(
  Trait=v.na,
  Year=v.na,
  Pop=v.na,
  RR=v.na,
  BL=v.na,
  RF=v.na)
ct = 1
years = sort(unique(Res$file))
pops = rev(sort(unique(Res$exp)))
for (pop in pops) {
  for (phenoi in phenos) {
    for (year in years) {
      x = Res
      x = x[x$trait == phenoi,]
      x = x[x$file == year,]
      x = x[x$exp == pop,]
      M$Trait[ct] = phenoi
      M$Year[ct] = year
      M$Pop[ct] = pop
      M$RR[ct] = x$r[x$method=="RR"]
      M$BL[ct] = x$r[x$method=="BL"]
      M$RF[ct] = x$r[x$method=="RF"]
      ct = ct + 1
    }
  }
}
write.csv(M, "outputs/Table3.csv", row.names=F)

