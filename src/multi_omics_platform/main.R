# Import packages
library(data.table)
library(ggplot2)
library(gridExtra)
library(BGLR)

# Import modules
source("1_data_load.R")
source("2_matrices.R")
source("3_cv_prep.R")
source("4_fit_models.R")

### USER ARGUMENTS ###
env.col <- 11                           # column with ENV ID
gid.col <- 2                            # column with genotype ID
trait.col <- 3                          # column with desired trait
phen.col <- 15                          # starting CV0 results
cv.col <- 12

marker.path <- "../../data/SNPs.rda"
phenos.path <- "../../data/Phenos.csv"

nan.freq <- 0.2                         # NaN threshold limit
reps <- 1                               # Number of repetitions for each analysis
folds <- 5                              # Value for k-fold cross-validation

ctr <- TRUE
std <- TRUE
weighting <- FALSE
prop.maf.j <- NULL 

cv1 <- TRUE
cv2 <- TRUE
cv0 <- TRUE
cv00 <- TRUE

nIter <- 12000
burnIn <- 2000
esc <- FALSE

col.cv <- 12

### 1 - Data Load ###
loaded.data <- loadData(phenos.path, marker.path)

markers <- loaded.data$markers
phenos <- loaded.data$phenos

rm(loaded.data)

# Create the NaNs csv files
createNaNFiles(phenos, markers, nan.freq)  # These files are the Mod1 outputs

### 2 - G/E Matrices ###

# E matrix
generateMatrix("../../output/E/", phenos, markers = NULL, col.env.id = env.col)
# G matrix
generateMatrix("../../output/G/", phenos=phenos, markers=markers, 
                col.env.id = gid.col, prop.maf.j =  NULL)
# ZE matrix
createZMatrix(phenos, env.col, "../../output/ZE/")
# ZL matrix
createZMatrix(phenos, gid.col, "../../output/ZL/")

# Interaction matrix
g1.file <- '../../output/G/G.rda'          # path to matrix file 1
g2.file <- '../../output/E/G.rda'          # path to matrix file 2

generateIntMatrix(g1.file, g2.file, output.path='../../output/GE/')

### 3 - Phenotype data prep ###
set.seed(1)
phenos.cv <- phenos

if (cv1) {
  phenos.cv <- cvPrep(phenos.cv, "../../output/CV1/", col.id= gid.col , folds = folds, cv1 = TRUE)
}

if (cv2) {
  phenos.cv <- cvPrep(phenos.cv, "../../output/CV2/", col.id= gid.col, folds = folds, cv2 = TRUE)
}

if (cv0) {
  phenos.cv <- cvPrep(phenos.cv, "../../output/CV0/", col.id = trait.col, folds = folds, cv0 = TRUE)
}

if (cv00) {
  phenos.cv <- cvPrep(phenos.cv, "../../output/CV00/", col.id = trait.col, folds = folds, cv00 = TRUE)
}

### 4 - Fit models ####

# Output structure: trait (e.g height) --> CV --> fold_n --> predictions.csv 

ab.list <- list()
ab.list[[1]] <- '../../output/ZE/Z.rda'
ab.list[[2]] <- '../../output/ZL/Z.rda'

# E + L
runBGLR(phenos.cv, phen.col, gid.col, cv.col, env.col = NULL, file_list = ab.list, 
        folds = folds)

# E + L + G
ab.list[[3]] <- '../../output/G/EVD.rda'  
set.seed(1)

runBGLR(phenos.cv, phen.col, gid.col, cv.col, env.col = NULL, file_list = ab.list, 
        folds = folds)

# E + L + G + GE
ab.list[[4]] <- '../../output/GE/EVD.rda'  
set.seed(1)

runBGLR(phenos.cv, phen.col, gid.col, cv.col, env.col = NULL, file_list = ab.list, 
        folds = folds)