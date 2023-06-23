# Import packages
library(data.table)
library(ggplot2)
library(gridExtra)

# Import modules
source("1_data_load.R")
source("2_ge_matrices.R")

### USER ARGUMENTS ###

args <- commandArgs(trailingOnly = TRUE)
nan.freq <- 0.2                            # NaN threshold frequency
col.env.id <- 11                           # Env ID in phenos file
marker.path <- "../../data/SNPs.rda"       # Covariate matrix file
phenos.path <- "../../data/phenos.csv"     # Phenotype/Environment file

ctr <- TRUE
std <- TRUE
weighting <- FALSE                      # Remove line
prop.maf.j <- NULL                         # Remove line


### 1 - Data Load ###
loaded.data <- loadData(phenos.path, marker.path)

markers <- loaded.data$markers
phenos <- loaded.data$phenos

rm(loaded.data)

# Create the NaNs csv files
createNaNFiles(phenos, markers, nan.freq)  # These files are the Mod1 outputs

### 2 - G/E Matrices ###

# E matrix
createMatrixPdf("../../output/E/", phenos, markers = NULL, col.env.id = col.env.id )

# G matrix
createMatrixPdf("../../output/G/", phenos=phenos, markers=markers, 
                col.env.id = 2, prop.maf.j =  0.03)

# ZE matrix
createZMatrix(phenos, 11, "../../output/ZE/")

# ZL matrix
createZMatrix(markers, 2, "../../output/ZL/") # TODO Confirm what dataset to use



