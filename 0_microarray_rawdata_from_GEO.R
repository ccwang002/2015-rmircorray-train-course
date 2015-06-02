# source("http://bioconductor.org/biocLite.R")
# biocLite(c("GEOquery", "affy", "limma"))
library(GEOquery)

# Download raw data from GEO database
getGEOSuppFiles("GSE10940", baseDir = "./GEO_datasets/")
# Extract the CEL files by the following commands (already done)
# tar xf GSE10940_RAW.tar
# gunzip [--keep] *.gz

# Read in the published array data (processed probe values)
gse10940 <- getGEO(
    "GSE10940", GSEMatrix = FALSE,
    destdir = "./GEO_datasets/"
)
# Get the description about each samples(GSM files)
gse10940_pheno_tab <- pData(phenoData(gse10940[[1]]))
gse10940_pheno_tab[, c("title", "type", "source_name_ch1")]

# Get the platform of each samples
gsmplatforms <- lapply(GSMList(gse10940), function(x){Meta(x)$platform})

# Get the platform information
gpl <- getGEO(gsmplatforms[[1]], destdir = "./GEO_datasets/")
cat(Meta(gpl)$title)  # Print out the platform name

# Get published probe data
pub_probe_val <- Table(GSMList(gse10940)[[1]])
Columns(GSMList(gse10940)[[1]])  # Get the description for each column

