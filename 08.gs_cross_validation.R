### GS cross validation : Figure 5 ##

## libraries ##
library(rrBLUP)
library(ggplot2)
source("inhouse.functions.R")

### data ###
load("data/geno.train.rda")
load("data/pheno.train.rda")

### execution ###
g = t(geno.train[,-c(1:14)])
phenos = colnames(pheno.train)[-c(1,2)]

n.markers = c(10, 100, 1000, 10000, ncol(g))
id.markers = c(" 10 ", " 100 ", " 1000 ", " 10K ", " WGS ")

X = c()
ct = 1
for (phenoi in phenos) {
	res.r = res.p = c()
	for (i in 1:length(n.markers)) {
		n.p = n.markers[i]
		print(paste(phenoi, n.p))
		ans <- CV(pheno.train, phenoi, g=g, n.fold=10, n.rep=10, n.p, seed.numb=1, n.core=10)
		res.r = c(res.r, ans)
		res.p = c(res.p, rep(id.markers[i], 10))
	}
	x = data.frame(r=res.r, n.variants=res.p, t=rep(ct, length(n.markers)))
	X = rbind(X, x)
	ct = ct + 1
}
	
gp = ggplot(X, aes(y=r,x=n.variants, fill=n.variants)) + 
	geom_boxplot() + 
	labs(x="No. of variants", y=expression(italic(r))) + 
	theme(axis.title=element_text(size=12)) + 
	theme(axis.text=element_text(size=9)) + 
	theme(strip.text=element_text(size=12,face="bold")) + 
	scale_fill_discrete("") + 
	theme(legend.position = 'none') +
	facet_wrap(~ t, labeller=as_labeller(c("1"="SSC","2"="Fruit weigth","3"="4th truss height","4"="Yield")))
ggsave(file = "outputs/Fig5.png", plot = gp, dpi = 400, width = 4, height = 4)

