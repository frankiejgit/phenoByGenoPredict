library(data.table)

### ARGUMENTS TO PASS ###

args <- commandArgs(trailingOnly = TRUE)
nan.freq <- 0.2 # 1 - NaN threshold frequency

# This may be something created by the program
col.env.id <- 11 # The environment ID in the phenos file

### 1 - Splitting SNPs ###

# Get the molecular marker information 
load('../../data/markers.rda')
assign("markers", SNPs, envir = .GlobalEnv)
rm(SNPs)
phenos <- read.csv('../../data/phenos.csv')

# Use args[1] for missing value threshold
nan.freq.i <- length(unique(phenos$strain)) * nan.freq
NaNs <- colSums(is.na(markers))

index.1 <- order(NaNs)
NaN.freq <- unique(NaNs[index.1])
percentage.NaNs <- NaN.freq * 100 / nrow(markers)

# Write the NaN frequencies and percentages to separate CSV files
write.table(NaN.freq, file = '../../data/NaNs.freq.csv', sep = ',', row.names = FALSE, col.names = FALSE)
write.table(percentage.NaNs, file = '../../data/percentage.NaNs.csv', sep = ',', row.names = FALSE, col.names = FALSE)

# Subset the 'SNPs' matrix based on a condition and save it as 'X'
index.2 <- NaNs <= nan.freq.i
X <- markers[, index.2]
write.csv(X, file = '../../data/X.csv', row.names = FALSE)

### 2 - Create G/E matrices ###

## E - Phenos ##

# Check parameters 
if(is.null(markers)) {
  cat('Are you working with a grouping factor?
       If not, please provide a path of the covariate matrix file')
  
  # Prompt the user for input
  user_input <- readline("Do you want to continue without the matrix file? (y/n): ")
  
  # Check the user's response
  if (tolower(user_input) != "n") {
    stop("Please provide the covariate matrix file and rerun the pipeline.")
  }
}

### Reads data and perform consistency checks
cat('=2===> Reading Data '); cat('\n')

if (is.null(markers)) {
  
  Y <- phenos
  Y[, col.env.id] <- factor(Y[, col.env.id])
  n <- length(levels(Y[, col.env.id]))
  Z <- as.matrix(model.matrix(~Y[, col.env.id] - 1))
  d <- colSums(Z)
  V <- Z
  
  for (i in 1:ncol(Z)) {
    V[, i] <- V[, i] / sqrt(d[i])
  }
  
  EVD <- list(vectors = V, values = d)
  G <- tcrossprod(Z)
  
  save(G, file = '../../output/G.rda')
  save(EVD, file = '../../output/EVD.rda')
  } else {
  ## Reads covariates
  n <- ncol(markers)
  p <- length(colnames(markers))
  
  Y <- matrix(nrow = n, ncol = p, NA)
  id.Y <- rep(NA, n)
  
  for (i in 1:n) {
    if (!is.null(col.env.id)) { id.Y[i] <- markers[1] }
    
    ifelse(!is.null(co))
    
    ifelse(!is.null(colIDy) == FALSE, X[i, ] <- as.numeric(tmp), X[i, ] <- as.numeric(tmp[-1]))
    print(i)
  }
  close(fileIn)
  rownames(X) <- IDx
  colnames(X) <- colNames
  
  if (weighting) {
    weight <- scan(weight.file, skip = 1)
  }
  
  S <- 0
  
  ## Imputing, centering, standardizing and weighting.
  for (i in 1:ncol(X)) {
    meanXi <- mean(X[, i], na.rm = TRUE)
    X[, i] <- ifelse(is.na(X[, i]), meanXi, X[, i])  # Naive imputation
    if (ctr) { X[, i] <- X[, i] - meanXi }  # Centering
    if (std) { X[, i] <- X[, i] / sd(X[, i]) }  # Standard
    