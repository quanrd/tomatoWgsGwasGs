### population structire analysis : Figure 2 ###

## libraries ###
library(rrBLUP)
library(ggplot2)
source("inhouse.functions.R")

## data ##
load("data/K.train.rda")

## execution ##
pca = prcomp(K.train)
PoV = pca$sdev^2/sum(pca$sdev^2)*100
PoV = round(PoV,1)
x = as.data.frame(pca$x[,1:3])

gp = ggplot(x, aes(x=PC2, y=PC1, shape="", colour="")) + 
	geom_point() + 
	labs(x=paste("PC2 (",PoV[2]," %)",sep=""), y=paste("PC1 (",PoV[1]," %)",sep="")) + 
	theme(axis.title=element_text(size=12)) + 
	theme(axis.text=element_text(size=9)) + 
	scale_fill_discrete("") + 
	theme(legend.position="none") +
	scale_colour_manual(values="hotpink") + 
	scale_shape_manual(values=1)
ggsave(file = "outputs/Fig2A.png", plot = gp, dpi = 400, width = 2.5, height = 2.5)
	
gp = ggplot(x, aes(x=PC3, y=PC1, shape="", colour="")) + 
	geom_point() + 
	labs(x=paste("PC3 (",PoV[3]," %)",sep=""), y=paste("PC1 (",PoV[1]," %)",sep="")) + 
	theme(axis.title=element_text(size=12)) + 
	theme(axis.text=element_text(size=9)) + 
	scale_fill_discrete("") + 
	theme(legend.position="none") +
	scale_colour_manual(values="hotpink") + 
	scale_shape_manual(values=1)
ggsave(file = "outputs/Fig2B.png", plot = gp, dpi = 400, width = 2.5, height = 2.5)
