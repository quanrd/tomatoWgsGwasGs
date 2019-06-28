### haplotyping : Figure 1B ###

## library ##
library(ggplot2)
source("inhouse.functions.R")

## data ##
load("data/geno.train.rda")

### execution ###
pro = haplotyping(geno.train,"PROM")
cds = haplotyping(geno.train,"GENE")
dow = haplotyping(geno.train,"DOWN")
haplo = rbind(pro,cds,dow)
haplo = haplo[order(haplo[,2],haplo[,3]),]
save(haplo, file="data/haplo.rda")

H = haplo[,-c(1:4)]
n.hap = apply(H, 1, count.hap)
t.hap = table(n.hap)
x = data.frame(t.hap)	

gp = ggplot(x, aes(x=n.hap, y=Freq, fill="")) + 
	geom_bar(stat="identity") +
	geom_text(aes(label=Freq), vjust=c(2,2,-0.4,-0.4,-0.4,-0.5,-0.5,-0.2,-0.2,-0.2,-0.2), size=3) +
	labs(x="No. of alleles", y="No. of loci") + 
	theme(axis.title=element_text(size=12)) + 
	theme(axis.text=element_text(size=9)) + 
	scale_fill_discrete("") + 
	theme(legend.position="none") +
	scale_fill_manual(values="#00a0e9")
ggsave(file = "outputs/Fig1B.png", plot = gp, dpi = 400, width = 4.5, height = 2)
