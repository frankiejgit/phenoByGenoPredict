# Import packages
library(data.table)
library(ggplot2)
library(gridExtra)
library(BGLR)

# Import modules
source("1_data_load.R")
source("2_matrices.R")
source("3_cv_prep.R")

### USER ARGUMENTS ###
env.col <- 11                           # column with ENV ID
gid.col <- 2                            # column with genotype ID
trait.col <- NULL                       # column with desired trait

marker.path <- "../../data/SNPs.rda"
phenos.path <- "../../data/Phenos.csv"

nan.freq <- 0.2                         # NaN threshold limit
reps <- 1                               # Number of repetitions for each analysis

ctr <- TRUE
std <- TRUE
weighting <- FALSE
prop.maf.j <- NULL 

nIter <- 12000
burnIn <- 2000

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
createZMatrix(phenos, 2, "../../output/ZL/") # TODO Confirm what dataset to use

# Interaction matrix
g1.file <- '../../output/G/G.rda'          # path to matrix file 1
g2.file <- '../../output/E/G.rda'          # path to matrix file 2

generateIntMatrix(g1.file, g2.file, output.path='../../output/GE/')

### 3 - Phenotype data prep ###
set.seed(1)

# Do only CV1 and CV2
phenos.cv <- cvPrep(phenos, paste0("../../output/"), col.id=2, folds = 5, cv0 = FALSE, cv00 = FALSE)

# Do only CV0 and CV00
for (col in trait.cols) {
  col.name <- colnames(phenos.cv)[col]
  phenos.cv <- cvPrep(phenos.cv, "../../output/cv/", col.id = col, col.folds = 14,
                      folds = 5, cv1 = FALSE, cv2 = FALSE)
  }

### 4 - Fit models ####

# Output structure: trait (e.g height) --> CV --> fold_n --> predictions.csv 

# CV0 - E+L
# Only thing that changes between folds is col.phen


nIter <- 12000
burnIn <- 2000
folds <- 6
reps <- 10

col.phen <- 15 # starting CV0 results
col.var <- 2  # ID of variety
col.cv <- 12
phen.folds <- 1:5
col.env <- NULL
cv0 <- TRUE
esc <- FALSE

phen.name <- strsplit(colnames(phenos.cv)[col.phen], "_")[[1]][1]
y   <- phenos.cv[, col.phen]
gid <- phenos.cv[, col.var]

if(esc) { y <- scale(y, center=TRUE, scale=TRUE) }

z.list <- list()
z.list[[1]] <- '../../output/ZE/Z.rda'
z.list[[2]] <- '../../output/ZL/Z.rda'

models <- c('FIXED')     # Option to add more

eta <- list()
for (i in seq_along(z.list)) {
  Z <- get(load(z.list[[i]]))
  eta[[i]] <- list(X=Z, model='FIXED')
  rm(Z)
}

for (fold in 1:folds) {
  if (fold != -999) {
    
    output.path <- paste("../../output/", phen.name, "/E+L/fold_", fold, "/", sep = "")
    if (!dir.exists(output.path)) { dir.create(output.path, recursive = TRUE) }
    
    testing <- which(phenos.cv[, col.cv] == fold)
    
    if (cv0) { testing <- intersect(testing, which(gid %in% gid[testing])) }
    
    y.na <- y
    y.na[testing] <- NA
    
    fm <- BGLR(y = y.na, ETA = eta, nIter = nIter, burnIn = burnIn, verbose=TRUE)
    fm$y <- y
    
    predictions <- data.frame(testing = testing, Individual = gid[testing], y = y[testing], yHat = fm$yHat[testing])
    write.table(predictions, file = paste(output.path, "predictions_", fold, ".csv", sep=''), row.names = FALSE, sep = ",")

  } else {
    output.path <- paste("../../output/", phen.name, "/full_data/", sep='')
    if (!dir.exists(output.path)) { dir.create(output.path) }
    
    fm <- BGLR(y = y, ETA = ETA, nIter = nIter, burnIn = burnIn, verbose = TRUE)
    save(fm, file = 'fm_full.RData')
  }
  
  rm(fm)
  file.remove(list.files(pattern = "*.dat"))
}

# E + L + G
evd.path <- '../../output/G/EVD.rda'  
cv0 <- FALSE
ESC <- FALSE 
set.seed(1)

if(esc) { y <- scale(y, center=TRUE, scale=TRUE) }

z.list <- list()
z.list[[1]] <- '../../output/ZE/Z.rda'
z.list[[2]] <- '../../output/ZL/Z.rda'
z.list[[3]] <- evd.path

models <- c('FIXED','FIXED','RKHS')     # Option to add more

eta <- list()
for (i in seq_along(z.list)) {
  
  Z <- get(load(z.list[[i]]))
  eta[[i]] <- list(X=Z, model='FIXED')
  
  if (models[i] == "RKHS") {
    eta[[i]] <- list(V = EVD$vectors, d=EVD$values, model=models[i])
    rm(EVD)
  }
  
  rm(Z)
  
}

for (fold in 1:folds) {
  if (fold != -999) {
    
    output.path <- paste("../../output/", phen.name, "/E+L+G/fold_", fold, "/", sep = "")
    if (!dir.exists(output.path)) { dir.create(output.path, recursive = TRUE) }
    
    testing <- which(phenos.cv[, col.cv] == fold)
    
    if (cv0) { testing <- intersect(testing, which(gid %in% gid[testing])) }
    
    y.na <- y
    y.na[testing] <- NA
    
    fm <- BGLR(y = y.na, ETA = eta, nIter = nIter, burnIn = burnIn, verbose=TRUE)
    fm$y <- y
    
    predictions <- data.frame(testing = testing, Individual = gid[testing], y = y[testing], yHat = fm$yHat[testing])
    write.table(predictions, file = paste(output.path, "predictions_", fold, ".csv", sep=''), row.names = FALSE, sep = ",")
    
  } else {
    output.path <- paste("../../output/", phen.name, "/full_data/", sep='')
    if (!dir.exists(output.path)) { dir.create(output.path) }
    
    fm <- BGLR(y = y, ETA = ETA, nIter = nIter, burnIn = burnIn, verbose = TRUE)
    save(fm, file = 'fm_full.RData')
  }
  
  rm(fm)
  file.remove(list.files(pattern = "*.dat"))
}