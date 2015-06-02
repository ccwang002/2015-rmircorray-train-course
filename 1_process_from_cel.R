# Load required libraries
library(GEOquery)
library(affy)
library(limma)
library(plyr)
library(stringr)
options(stringsAsFactors = FALSE)

# read in GSE dataset
gse10940 <- getGEO(
    "GSE10940", GSEMatrix = TRUE,
    destdir = "./GEO_datasets/"
)

# Load the sample information
pheno_data <- phenoData(gse10940[[1]])
# Read CEL files
cel_data <- ReadAffy(
    # celfile.path = "./GEO_datasets/GSE10940/",
    sampleNames = sampleNames(pheno_data),
    filenames = list.celfiles("./GEO_datasets/GSE10940/", full.names = TRUE)
)
phenoData(cel_data) <- pheno_data  # put back the sample information

# Normalization using RMA
eset <- rma(cel_data)

# Read in sample platform annotation (GPL)
gpl <- getGEO(annotation(gse10940[[1]]), destdir = "./GEO_datasets/")
cat(Meta(gpl)$title)
probe_tab <- Table(gpl)
# Load only a part of the annotation with the array data
fData(eset) <- probe_tab[, c("Gene Title", "Gene Symbol", "ENTREZ_GENE_ID", "RefSeq Transcript ID")]

# # Check if two probe names match
# featureNames(est)[which(featureNames(est) != probe_tab$ID)]
# probe_tab[which(probe_tab$ID != featureNames(est)), c("ID")]

# View the probe expression normalization result
boxplot(cel_data)
boxplot(exprs(est))
# # For MAplot see the following code snippet
# MAplot(cel_data[, 5:8], pairs = TRUE, plot.method = "smoothScatter")

# Define a helper function to extract sample information
extract_word_after_colon <- function(col) {str_match(col, ": ([^:]+)$")[, 2]}
pheno_df <- pData(pheno_data)
cleaned_pheno_df <- dplyr::transmute(
    pheno_df,
    title = title,
    geo_accession = geo_accession,
    sample_group = factor(extract_word_after_colon(characteristics_ch1)),
    tissue = factor(extract_word_after_colon(characteristics_ch1.1))
)

# Design our experiment modle
design <- model.matrix(~ sample_group + tissue, data = cleaned_pheno_df)
fit <- lmFit(eset, design)  # fit each probeset to model

# Use ANOVA-like test to find differentially expressed genes across different sample groups
efit <- eBayes(fit)        # empirical Bayes adjustment
top_de_probes_anova <- topTable(efit, coef = c("sample_groupLogjam"), n = 20)      # table of differentially expressed probesets

# Use (log) fold change to find DE genes
fit2 <- treat(fit, lfc = 0.1)
top_de_probes_foldchange <- topTreat(fit2, coef = c("sample_groupLogjam"))

# Export as CSV files
# write.csv(top_de_probes_foldchange, "fold_change_genes.csv")
