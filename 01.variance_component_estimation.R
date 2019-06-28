### variance components estimation : Table 2 ##

## libraries ###
library(rrBLUP)

### data ###
load("data/geno.train.rda")
load("data/pheno.train.rda")

### execution ###
g = geno.train[,-c(1:14)]
K.train = A.mat(t(g))
save(K.train, file="data/K.train.rda")

Z = matrix(0, nr=nrow(pheno.train), nc=nrow(K.train))
rownames(Z) = pheno.train$gid
colnames(Z) = rownames(K.train)
Z[cbind(1:nrow(Z), match(pheno.train$gid, rownames(K.train)))] = 1
	
X = matrix(0, nr=nrow(pheno.train), nc=length(unique(pheno.train$year)))
rownames(X) = pheno.train$gid
colnames(X) = unique(pheno.train$year)
X[cbind(1:nrow(X), match(pheno.train$year, colnames(X)))] = 1

phenos = colnames(pheno.train)[-c(1,2)]
H2 = numeric(length(phenos))
names(H2) = phenos

for (phenoi in phenos) {
	y = pheno.train[[phenoi]]
	use = !is.na(y)
	y.use = y[use]
	Z.use = Z[use,]
	X.use = X[use,]; X.use = X.use[,apply(X.use, 2, sum)!=0]
	res = mixed.solve(y=y.use, Z=Z.use, K=K.train, X=X.use)
	H2[names(H2)==phenoi] = (var(y.use) - res$Ve) / var(y.use)
}

H2 = t(H2)
rownames(H2) = "H2"	
write.csv(t(H2), "outputs/Table2.csv")
