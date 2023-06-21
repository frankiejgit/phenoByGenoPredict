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
loaded.data <- loadData(phenos.path, marker.path)

markers <- loaded.data$markers
phenos <- loaded.data$phenos

# Create the NaNs csv files
createNaNFiles(phenos, markers, nan.freq)  # These files are the Mod1 outputs

### 2 - G/E Matrices ###

## E matrix ##

cat('== 1 ==> Reading and checking parameters\n')
if(is.null(markers)) {
  cat('Are you working with a grouping factor?
       If not, please provide a path of the covariate matrix file\n')
  
  # Prompt the user for input
  user_input <- readline("Do you want to continue without the matrix file? (y/n): ")
  
  # Check the user's response
  if (tolower(user_input) != "n") {
    stop("Please provide the covariate matrix file and rerun the pipeline.\n")
  }
}

# Read data and perform consistency checks
cat("== 2 ==> Reading data\n")

if (!is.null(col.env.id)) {
  env.IDs <- phenos[, col.env.id]
}

if (is.null(markers)) {
  env.IDs <- factor(phenos[, col.env.id])
  n <- length(levels(env.IDs))
  Z <- as.matrix(model.matrix(~env.IDs - 1))
  d <- colSums(Z)
  V <- Z
  
  for (i in 1:ncol(Z)) {
    V[, i] <- V[, i] / sqrt(d[i])
  }
  
  EVD <- list(vectors = V, values = d)
  E.G <- tcrossprod(Z)
  
  save(E.G, file = '../../output/e_G.rda')
  save(EVD, file = '../../output/e_EVD.rda')
  
} else {
  
  # Reads covariates
  colNames <- colnames(markers)
  n <- nrow(markers)
  p <- length(colNames)
  
  ##### TODO Remove mm.matrix, it is identical to `markers` #####
  
  # Create empty matrix with same dimensions as `markers`
  mm.matrix <- matrix(nrow=n, ncol=p, NA)
  # List of NaN values equal to length as rows in `markers`
  id.mm <- rep(NA, n)
  
  for (i in 1:n) {
    if (!is.null(col.env.id)) {
      tmp <- markers[i,]    # Look at each line at a time
      id.mm[i] <- row.names(markers)[i]    # TODO Original code uses markers[1], verify
      mm.matrix[i,] <- as.numeric(tmp)
    }
  }
  
  rownames(mm.matrix) <- id.mm
  colnames(mm.matrix) <- colNames
  
  ###### END of TODO ######
  
  if (weighting) {
    weight <- scan(weight.file, skip =1 )
  }
  
  # Imputing, centering, standarizing and weighting
  S <- 0
  for (i in 1:p) {
    mean.i <- mean(mm.matrix[, i], na.rm=TRUE)
    # Naive imputation
    mm.matrix[,i] <- ifelse(is.na(mm.matrix(X[,i])), mean.i, mm.matrix[,i])
    print(i)
    S <- S + var(mm.matrix[,i])
  }
  
  p2 <- colMeans(mm.matrix, na.rm=TRUE) / 2
  p2 <- ifelse(p <= 0.5, p, 1-p)
  
  
}

