## Load packages
library(affy)   # Affymetrix pre-processing
library(limma)  # two-color pre-processing; differential expression

phenoData <- read.AnnotatedDataFrame(
    system.file("extdata", "pdata.txt", package = "arrays")
)
show(pData(phenoData))

celfiles <- system.file("extdata", package="arrays")
eset <- justRMA(phenoData=phenoData,
                celfile.path=celfiles)


combn <- factor(paste(pData(phenoData)[,1],
                      pData(phenoData)[,2], sep = "_"))
design <- model.matrix(~combn)

fit <- lmFit(eset, design)  # fit each probeset to model
efit <- eBayes(fit)        # empirical Bayes adjustment
topTable(efit, coef=2)      # table of differentially expressed probesets


fit2 <- treat(fit, lfc=0.1)
topTreat(fit2, coef=2)
