# source("http://bioconductor.org/biocLite.R")
# biocLite(c("GEOquery", "affy", "limma"))
library(GEOquery)

getGEOSuppFiles("GSE10940", baseDir = "./GEO_datasets/")
# tar xf GSE10940_RAW.tar
# gunzip [--keep] *.gz
gse10940 <- getGEO(
    "GSE10940", GSEMatrix = FALSE,
    destdir = "./GEO_datasets/"
)
gse10940_pheno_tab <- pData(phenoData(gse10940[[1]]))
gse10940_pheno_tab[, c("title", "type", "source_name_ch1")]
gpl <- getGEO(annotation(gse10940[[1]]), destdir = "./GEO_datasets/")
cat(Meta(gpl)$title)

gsmplatforms <- lapply(GSMList(gse10940), function(x){Meta(x)$platform})
pub_probe_val <- Table(GSMList(gse10940)[[1]])
Columns(GSMList(gse10940)[[1]])
