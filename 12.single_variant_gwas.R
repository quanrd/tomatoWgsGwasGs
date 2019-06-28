## single variant-based GWAS : Figure S1 ##

### libraries ###
library(rrBLUP)
library(qqman)

### data ###
load("data/geno.train.rda")
load("data/K.train.rda")
load("data/pheno.train.rda")

### execution ###
Ncore = 16

g = geno.train[,-c(4:14)]
scores = GWAS(pheno.train, g, fixed="year", K=K.train, min.MAF=0.05, n.core=20, n.PC=3, P3D=TRUE, plot=FALSE)	
phenos = colnames(pheno.train)[-c(1,2)]
jpeg(filename="outputs/FigS1.jpeg", height = 1600, width=600)
par(mfrow=c(4,1))
for (phenoi in phenos) {
  input = data.frame(SNP = scores$MARKER,
                     CHR = scores$CHROM,
                     BP = scores$POS,
                     P = 10^(-1*scores[[phenoi]])
                     )
  input$CHR = as.numeric(input$CHR)
  manhattan(input, 
            suggestiveline = -log10(0.05/nrow(scores)),
            genomewideline = -log10(0.01/nrow(scores)))
}
dev.off()
