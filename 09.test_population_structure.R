### test popultion structure : Figure6A ###

## libraries ##
library(ggplot2)
library(rrBLUP)

## data ##
load("data/geno.train.rda")
load("data/geno.test.rda")

### execution ###
g = cbind(geno.train[,-c(1:14)], geno.test[,-c(1:14)])
K = A.mat(t(g), min.MAF=0.05, n.core=10)

pca = prcomp(K)
PoV = pca$sdev^2/sum(pca$sdev^2)*100
PoV = round(PoV,1)
x = as.data.frame(pca$x[,1:2])

Population = rep(1, nrow(x))
Population[grep("PS.", rownames(x))] = 2
x = cbind(x, Population)
x$Population = as.factor(x$Population)

pches = rep(1, nrow(x))
pches[grep("PS.", rownames(x))] = 0
colors <- rep("hotpink1", nrow(x))
colors[grep("PS.", rownames(x))] = "dodgerblue3"

gp = ggplot(x, aes(x=PC2, y=PC1, shape=Population, colour=Population)) + 
	geom_point() + 
	labs(x=paste("PC2 (",PoV[2]," %)",sep=""), y=paste("PC1 (",PoV[1]," %)",sep="")) + 
	theme(axis.title=element_text(size=12)) + 
	theme(axis.text=element_text(size=9)) + 
	scale_colour_manual(labels=c("Training","Test"), values=c("hotpink1", "dodgerblue3")) + 
	scale_shape_manual(labels=c("Training","Test"), values=c(1,0))
ggsave(file = "outputs/Fig6A.png", plot = gp, dpi = 400, width = 3.6, height = 2.5)

