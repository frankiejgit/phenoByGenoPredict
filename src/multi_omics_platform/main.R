# Import packages
library(data.table)

# Import modules
source("1_data_load.R")

### USER ARGUMENTS ###

args <- commandArgs(trailingOnly = TRUE)
nan.freq <- 0.2                            # NaN threshold frequency
col.env.id <- 11                           # Env ID in phenos file
marker.path <- "../../data/SNPs.rda"       # Covariate matrix file
phenos.path <- "../../data/phenos.csv"     # Phenotype/Environment file

### 1 - Data Load ###
loaded.data <- loadData(marker.path, phenos.path)

markers <- loaded.data$markers
phenos <- loaded.data$phenos


# Create the NaNs csv files
createNaNFiles(phenos, markers, nan.freq)