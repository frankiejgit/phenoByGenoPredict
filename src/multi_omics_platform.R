# Import packages
library(data.table)
library(ggplot2)
library(gridExtra)
library(BGLR)
library(matrixStats)

# Import modules
source("modules/1_data_load.R")
source("modules/2_matrices.R")
source("modules/3_cv_prep.R")
source("modules/4_fit_models.R")

### USER ARGUMENTS ###

marker.path <- "../data/SNPs.rda"    # Provided by user as marker file
phenos.path <- "../data/Phenos.csv"  # Provided as phenotype file
nan.freq <- 0.2                         # NaN threshold limit

env.col <- 11                           # column with ENV ID
gid.col <- 2                            # column with genotype ID
trait.col <- 3                          # column with desired trait

ctr <- TRUE
std <- TRUE
weighting <- FALSE
prop.maf.j <- NULL

cv1 <- TRUE
cv2 <- TRUE
cv0 <- TRUE
cv00 <- TRUE

reps <- 1
folds <- 5
nIter <- 12000
burnIn <- 2000
esc <- FALSE

### END USER ARGS ###

### 1 - Data Load ###
loaded.data <- loadData(phenos.path, marker.path)

# Create output directory if it doesn't exist
if (!dir.exists("../output")) { dir.create("../output") }

markers <- loaded.data$markers
phenos <- loaded.data$phenos

rm(loaded.data)

# Create the NaNs csv files
createNaNFiles(phenos, markers, nan.freq)  # These files are the Mod1 outputs

### 2 - G/E Matrices ###

# E matrix
generateMatrix("../output/E/", phenos, markers = NULL, col.env.id = env.col)
# G matrix
generateMatrix("../output/G/", phenos=phenos, markers=markers, 
                col.env.id = gid.col, prop.maf.j =  NULL)
# ZE matrix
createZMatrix(phenos, env.col, "../output/ZE/")
# ZL matrix
createZMatrix(phenos, gid.col, "../output/ZL/")

# Interaction matrix
g1.file <- '../output/G/G.rda'          # path to matrix file 1
g2.file <- '../output/E/G.rda'          # path to matrix file 2

generateIntMatrix(g1.file, g2.file, output.path='../output/GE/')

### 3 - Phenotype data prep ###
set.seed(1)

# Do CV1 and CV2 first
phenos.cv <- phenos

phenos.cv <- cvPrep(phenos.cv, "../output/cv/", col.id = gid.col, folds = folds,
                    cv1 = cv1, cv2 = cv2)

phenos.cv <- cvPrep(phenos.cv, "../output/cv/", col.id = trait.col, folds = folds,
                    cv0 = cv0, cv00 = cv00)

### 4 - Fit models ####

# current structure: cv --> trait --> E+L --> folds --> cols 15:19
# Output structure: trait (e.g height) --> CV --> fold_n --> predictions.csv 

ab.list <- list()
ab.list[1] <- '../output/ZE/Z.rda'
ab.list[2] <- '../output/ZL/Z.rda'
ab.list[3] <- '../../output/G/EVD.rda'  
ab.list[4] <- '../../output/GE/EVD.rda'  


# Find the CV columns
cv.list <- list(
  cv1 = match("CV1", colnames(phenos.cv)),
  cv2 = match("CV2", colnames(phenos.cv)),
  cv0 = grep("CV0_", colnames(phenos.cv)),
  cv00 = grep("CV00_", colnames(phenos.cv))
)

for (i in 1:length(cv.list)) {
  cv <- names(cv.list)[i]
  val <- cv.list[i]$cv
  
  if (get(cv)) {
    if (length(val) > 1) {
      # Call the function for all the columns
      for (v in val) {
        print(v)
      }
    } else {
      # Run E + L
      runBGLR(cv, phenos.cv, trait.col, gid.col, val, env.col = NULL, 
              file_list = ab.list[1:2], folds = folds)
      # Run E + L + G
      runBGLR(cv, phenos.cv, trait.col, gid.col, val, env.col = NULL, 
              file_list = ab.list[1:3], folds = folds)
      # Run E + L + G + GE
      runBGLR(cv, phenos.cv, trait.col, gid.col, val, env.col = NULL, 
              file_list = ab.list, folds = folds)
    }
  }
  
}


# E + L
runBGLR(phenos.cv, phen.col, gid.col, cv1.col, env.col = NULL, file_list = ab.list, 
        folds = folds)

# E + L + G
set.seed(1)

runBGLR(phenos.cv, phen.col, gid.col, cv.col, env.col = NULL, file_list = ab.list, 
        folds = folds)

# E + L + G + GE
set.seed(1)

runBGLR(phenos.cv, phen.col, gid.col, cv1.col, env.col = NULL, file_list = ab.list, 
        folds = folds)

### 5 - Get Results ###

