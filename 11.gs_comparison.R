### WGS-based and selected variant-based GS : Figure 6C,D ###

## libraries ##
library(rrBLUP)
library(randomForest)
library(ggplot2)

### data ###
load("data/geno.train.rda")
load("data/geno.test.rda")
load("data/pheno.train.rda")
load("data/pheno.test.rda")

### execution ###
Ncores = 10

phenos = colnames(pheno.train)[-c(1:2)]

g.train = t(geno.train[,-(1:14)])
colnames(g.train) = geno.train$MARKER
g.test = t(geno.test[,-(1:14)])
colnames(g.test) = geno.test$MARKER

n = nrow(pheno.train)
Z = matrix(0,n,nrow(g.train))
rownames(Z) = pheno.train[,1]
colnames(Z) = rownames(g.train)
Z[cbind(1:n,match(pheno.train[,1],rownames(g.train)))] = 1
	
fixed.item = unique(pheno.train$year)
X = matrix(0, nr=nrow(pheno.train), nc=length(fixed.item))
colnames(X) = fixed.item
for (i in 1:length(fixed.item)) X[pheno.train$year==fixed.item[i],i] = 1

v.na = rep(NA, 48)
Res = data.frame(
	variants=v.na,
	file=v.na,
	variants=v.na,
	exp=v.na,
	trait=v.na,
	lower=v.na,
	r=v.na,
	upper=v.na)

proc=1
test.year = as.character(unique(pheno.test$year))

for (p in 1:length(test.year)) {
	print(test.year[p])
	p.test = pheno.test[pheno.test$year==test.year[p],]
	for (i in 1:length(phenos)) {
		print(i)
		phenoi = phenos[i]
		y = pheno.train[[phenoi]]
		use = !is.na(y)
		y.in = y[use]
		Z.in = Z[use,]
		X.in = X[use,]; X.in = X.in[,apply(X.in, 2, sum)!=0]
		ans = kinship.BLUP(y=y.in, G.train=g.train, G.pred=g.test, Z.train=Z.in, X=X.in, K.method="RR", n.core=Ncores)
		v.na = rep(NA, nrow(p.test))
		cmp = data.frame(gid=p.test$gid, obs=v.na, train=v.na, pred=v.na)
		cmp$obs = p.test[[phenoi]]
		for (j in 1:nrow(p.test)) {
			bv = ans$g.train[names(ans$g.train)==cmp$gid[j]]
			if (length(bv)!=0) cmp$train[j] = bv
		}
		for (j in 1:nrow(p.test)) {
			bv = ans$g.pred[names(ans$g.pred)==cmp$gid[j]]
			if (length(bv)!=0) cmp$pred[j] = bv
		}
		ct = cor.test(cmp$obs,cmp$train,use="pairwise.complete.obs")
		Res$variants[proc] = 1
		Res$file[proc] = p
		Res$exp[proc] = 1
		Res$trait[proc] = i
		Res$lower[proc] = ct$conf.int[1]
		Res$r[proc] = ct$estimate
		Res$upper[proc] = ct$conf.int[2]
		proc = proc + 1
		ct = cor.test(cmp$obs,cmp$pred,use="pairwise.complete.obs")
		Res$variants[proc] = 1
		Res$file[proc] = p
		Res$exp[proc] = 2
		Res$trait[proc] = i
		Res$lower[proc] = ct$conf.int[1]
		Res$r[proc] = ct$estimate
		Res$upper[proc] = ct$conf.int[2]
		proc = proc + 1
	}
}

for (p in 1:length(test.year)) {
	print(test.year[p])
	p.test = pheno.test[pheno.test$year==test.year[p],]
	for (i in 1:length(phenos)) {
		print(i)
		phenoi = phenos[i]
		y = pheno.train[[phenoi]]
			sel.mrk.name = paste("outputs/fig4_",phenoi,"_sel.mrks.rda",sep="")
			load(sel.mrk.name)
		G.train = g.train[,is.element(colnames(g.train), sel.mrks)]
		ZG.train = Z%*%as.matrix(G.train)
		y = pheno.train[[phenoi]]
		use = !is.na(y)
		Xuse = X[use,]; Xuse = Xuse[,apply(Xuse, 2, sum)!=0]			
		set.seed(1)
		fm = randomForest(x=cbind(Xuse, ZG.train[use,]), y=y[use])
		Xtrain = matrix(0, nr=nrow(G.train), nc=ncol(Xuse)); colnames(Xtrain) = colnames(Xuse)
		train = predict(fm, cbind(Xtrain, G.train))
		G.test = g.test[,is.element(colnames(g.test), sel.mrks)]
		Xtest = matrix(0, nr=nrow(G.test), nc=ncol(Xuse)); colnames(Xtest) = colnames(Xuse)
		pred = predict(fm, cbind(Xtest, G.test))
		v.na = rep(NA, nrow(p.test))
		cmp = data.frame(gid=p.test$gid, obs=v.na, train=v.na, pred=v.na)
		cmp$obs = p.test[[phenoi]]
		for (j in 1:nrow(p.test)) {
			bv = train[as.character(names(train))==as.character(cmp$gid[j])]
			if (length(bv)!=0) cmp$train[j] = bv
		}
		for (j in 1:nrow(p.test)) {
			bv = pred[as.character(names(pred))==as.character(cmp$gid[j])]
			if (length(bv)!=0) cmp$pred[j] = bv
		}
		ct = cor.test(cmp$obs,cmp$train,use="pairwise.complete.obs")
		Res$variants[proc] = 2
		Res$file[proc] = p
		Res$exp[proc] = 1
		Res$trait[proc] = i
		Res$lower[proc] = ct$conf.int[1]
		Res$r[proc] = ct$estimate
		Res$upper[proc] = ct$conf.int[2]
		proc = proc + 1
		ct = cor.test(cmp$obs,cmp$pred,use="pairwise.complete.obs")
		Res$variants[proc] = 2
		Res$file[proc] = p
		Res$exp[proc] = 2
		Res$trait[proc] = i
		Res$lower[proc] = ct$conf.int[1]
		Res$r[proc] = ct$estimate
		Res$upper[proc] = ct$conf.int[2]
		proc = proc + 1
	}
}

x = Res[Res$exp==1,]
x$r[x$r < 0] = 0.001
x$lower[x$lower < 0] = 0
x$file = as.factor(x$file)
gp = ggplot(x, aes(y=r,x=file, fill=factor(variants))) + 
	geom_bar(stat="identity", position="dodge") +
	labs(x="Phenotype data", y=expression(italic(r))) + 
	scale_fill_manual(name="GS model",labels=c("1"="WGS data","2"="Sel. variants"), values=c("#D6604D","#2166AC")) +
	scale_x_discrete(labels=c("1"="2018","2"="2019","3"="Rep. val.")) + 
	theme(axis.title=element_text(size=12)) + 
	theme(axis.text=element_text(size=9)) + 
	theme(strip.text=element_text(size=12,face="bold")) + 
	facet_grid(~ factor(trait), labeller=as_labeller(c("1"="SSC","2"="Fruit weigth","3"="4th truss height","4"="Yield")))
ggsave(file="outputs/Fig6C.png", plot = gp, dpi = 400, width = 8, height = 2.5)

x = Res[Res$exp==2,]
x$r[x$r < 0] = 0.001
x$lower[x$lower < 0] = 0
x$file = as.factor(x$file)
gp = ggplot(x, aes(y=r,x=file, fill=factor(variants))) + 
	geom_bar(stat="identity", position="dodge") +
	labs(x="Phenotype data", y=expression(italic(r))) + 
	scale_fill_manual(name="GS model",labels=c("1"="WGS data","2"="Sel. variants"), values=c("#D6604D","#2166AC")) +
	scale_x_discrete(labels=c("1"="2018","2"="2019","3"="Rep. val.")) + 
	theme(axis.title=element_text(size=12)) + 
	theme(axis.text=element_text(size=9)) + 
	theme(strip.text=element_text(size=12,face="bold")) + 
	facet_grid(~ factor(trait), labeller=as_labeller(c("1"="SSC","2"="Fruit weigth","3"="4th truss height","4"="Yield")))
ggsave(file="outputs/Fig6D.png", plot = gp, dpi = 400, width = 8, height = 2.5)

