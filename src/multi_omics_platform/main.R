# Import packages
library(data.table)
library(ggplot2)
library(gridExtra)

# Import modules
source("1_data_load.R")
source("2_matrices.R")
source("3_cv_prep.R")

### USER ARGUMENTS ###

args <- commandArgs(trailingOnly = TRUE)
nan.freq <- 0.2                            # NaN threshold frequency
col.env.id <- 11                           # Env ID in phenos file
marker.path <- "../../data/SNPs.rda"       # Covariate matrix file
phenos.path <- "../../data/Phenos.csv"     # Phenotype/Environment file

ctr <- TRUE
std <- TRUE
weighting <- FALSE
prop.maf.j <- NULL 

cv.reps <- 5
trait.cols <- 3:11

### 1 - Data Load ###
loaded.data <- loadData(phenos.path, marker.path)

markers <- loaded.data$markers
phenos <- loaded.data$phenos

rm(loaded.data)

# Create the NaNs csv files
createNaNFiles(phenos, markers, nan.freq)  # These files are the Mod1 outputs

### 2 - G/E Matrices ###

# E matrix
generateMatrix("../../output/E/", phenos, markers = NULL, col.env.id = col.env.id )
# G matrix
generateMatrix("../../output/G/", phenos=phenos, markers=markers, 
                col.env.id = 2, prop.maf.j =  NULL)
# ZE matrix
createZMatrix(phenos, 11, "../../output/ZE/")
# ZL matrix
createZMatrix(markers, 2, "../../output/ZL/") # TODO Confirm what dataset to use

# Interaction matrix
g1.file <- '../../output/G/G.rda'          # path to matrix file 1
g2.file <- '../../output/E/G.rda'          # path to matrix file 2

generateIntMatrix(g1.file, g2.file, output.path='../../output/GE/')

### 3 - Phenotype data prep ###
set.seed(1)

# col.id = 3, 4, 5, 6, 7, 8, 9, 10, 11
# col.folds = 14

# Do only CV1 and CV2
phenos.cv <- cvPrep(phenos, paste0("../../output/"), col.id=2, folds = 5, cv0 = FALSE, cv00 = FALSE)

# Do only CV0 and CV00
for (col in trait.cols) {
  col.name <- colnames(phenos.cv)[col]
  phenos.cv <- cvPrep(phenos.cv, paste0("../../output/", col.name,"/cv0/"), col.id = col, col.folds = 14, folds = 5, cv1 = FALSE, cv2 = FALSE, cv00 = FALSE)
  phenos.cv <- cvPrep(phenos.cv, paste0("../../output/", col.name,"/cv00/"), col.id = col, col.folds = 14, folds = 5, cv1 = FALSE, cv2 = FALSE, cv00 = FALSE)
}



