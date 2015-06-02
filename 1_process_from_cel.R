library(GEOquery)
library(affy)
library(limma)
library(plyr)
library(stringr)
options(stringsAsFactors = FALSE)

gse10940 <- getGEO(
    "GSE10940", GSEMatrix = TRUE,
    destdir = "./GEO_datasets/"
)

pheno_data <- phenoData(gse10940[[1]])
cel_data <- ReadAffy(
    # celfile.path = "./GEO_datasets/GSE10940/",
    sampleNames = sampleNames(pheno_data),
    filenames = list.celfiles("./GEO_datasets/GSE10940/", full.names = TRUE)
)
phenoData(cel_data) <- pheno_data
est <- rma(cel_data)

gpl <- getGEO(annotation(gse10940[[1]]), destdir = "./GEO_datasets/")
cat(Meta(gpl)$title)
probe_tab <- Table(gpl)
fData(est) <- probe_tab[, c("Gene Title", "Gene Symbol", "ENTREZ_GENE_ID", "RefSeq Transcript ID")]



# featureNames(est)[which(featureNames(est) != probe_tab$ID)]
# probe_tab[which(probe_tab$ID != featureNames(est)), c("ID")]

boxplot(cel_data)
boxplot(exprs(est))

# MAplot(cel_data[, 5:8], pairs = TRUE, plot.method = "smoothScatter")
extract_word_after_colon <- function(col) {
    str_match(col, ": ([^:]+)$")[, 2]
}

pheno_df <- pData(pheno_data)
cleaned_pheno_df <- dplyr::transmute(
    pheno_df,
    title = title,
    geo_accession = geo_accession,
    sample_group = factor(extract_word_after_colon(characteristics_ch1)),
    tissue = factor(extract_word_after_colon(characteristics_ch1.1))
)
design <- model.matrix(~ sample_group + tissue, data = cleaned_pheno_df)
fit <- lmFit(est, design)  # fit each probeset to model
efit <- eBayes(fit)        # empirical Bayes adjustment
top_de_probes <- topTable(efit, coef = c("sample_groupLogjam"), n = 20)      # table of differentially expressed probesets


fit2 <- treat(fit, lfc=0.1)
topTreat(fit2, coef=2)
