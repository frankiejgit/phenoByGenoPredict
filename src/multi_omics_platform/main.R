# Import packages
library(data.table)
library(ggplot2)
library(gridExtra)

# Import modules
source("1_data_load.R")

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

markers_2 <- NULL
if (is.null(markers_2)) {
  env.IDs <- factor(phenos[, col.env.id])
  n <- length(levels(env.IDs))
  Z <- as.matrix(model.matrix(~env.IDs - 1))
  d <- colSums(Z)
  V <- Z
  
  for (i in 1:ncol(Z)) {
    V[, i] <- V[, i] / sqrt(d[i])
  }
  
  EVD <- list(vectors = V, values = d)
  e.G <- tcrossprod(Z)
  
  save(e.G, file = '../../output/e_G.rda')
  save(EVD, file = '../../output/e_EVD.rda')
  
} else {
  
  # Reads covariates
  mm.matrix <- markers
  s <- 0
  col.count <- ncol(mm.matrix)
  
  if (weighting) {
    weight <- scan(weight.file, skip =1 ) # TODO check with Dr. Jarquin about this
  }
  
  # Naive imputation
  for (i in 1:col.count) {
    mean.i <- mean(mm.matrix[, i], na.rm=TRUE)
    mm.matrix[,i] <- ifelse(is.na(mm.matrix[,i]), mean.i, mm.matrix[,i])
    s <- s + var(mm.matrix[,i])
  }
  
  p <- colMeans(mm.matrix, na.rm=TRUE) / 2
  p <- ifelse(p <= 0.5, p, 1-p)
  
  # Adjusting for proportion of MAF
  if (!is.null(prop.maf.j)) {
    maf.idx <- which( p >= prop.maf.j )
    mm.matrix <- mm.matrix[, maf.idx]
  }
  
  # Imputing, centering, standardizing and weighting
  for (i in 1:col.count) {
    mean.i <- mean(mm.matrix[,i], na.rm=TRUE)
    if (ctr) { mm.matrix[,i] <- mm.matrix[,i] - mean.i }            # Centering
    if (std) { mm.matrix[,i] <- mm.matrix[,i]/sd(mm.matrix[,i]) }   # Standardizing
    if (weighting) { mm.matrix[,i] <- mm.matrix[,i] * weight[i] }   # Weighting
  }
  
  e.G <- tcrossprod(mm.matrix) / s

  # Check all IDs are in e.G
  if(!is.null(col.env.id)){
    stopifnot(all(env.IDs%in%rownames(e.G))) 
    }
  
  if(!is.null(col.env.id)){
    env.IDs <- factor(env.IDs, levels=rownames(e.G))
    e.Z <- as.matrix(model.matrix(~env.IDs-1))
    e.G <- tcrossprod(tcrossprod(e.Z,e.G),e.Z)
  }
  
  EVD <- eigen(e.G)
  rownames(EVD$vectors) <- rownames(e.G)
  
  save(E.G, file = '../../output/e_G.rda')
  save(EVD, file = '../../output/e_EVD.rda')
  
}

plt1 <- ggplot(data = as.data.frame(EVD$vectors), aes(x = EVD$vectors[, 1], y = EVD$vectors[, 2])) +
  geom_point() +
  labs(title = "", xlab = "First Component", ylab = "Second Component")
plt2 <- ggplot(data = as.data.frame(EVD$vectors), aes(x = EVD$vectors[, 1], y = EVD$vectors[, 3])) +
  geom_point() +
  labs(title = "", xlab = "First Component", ylab = "Third Component")
plt3 <- ggplot(data = as.data.frame(EVD$vectors), aes(x = EVD$vectors[, 2], y = EVD$vectors[, 2])) +
  geom_point() +
  labs(title = "", xlab = "Second Component", ylab = "Third Component")

plt4 <- ggplot(data = as.data.frame(EVD$values), aes(x = seq_along(EVD$values), y = EVD$values)) +
  geom_line() +
  labs(title = "Eigen-values", x = "Components", y = "Eigen-value") +
  theme_minimal()

plt5 <- ggplot(data = as.data.frame(EVD$values), aes(x = seq_along(EVD$values), y = cumsum(EVD$values) / sum(EVD$values))) +
  geom_line() +
  labs(title = "Cumulative Variance", x = "Components", y = "Cumulative Variance") +
  theme_minimal()

# Add horizontal and vertical lines to Plot 5
plt5 <- plt5 + geom_hline(yintercept = 0.8, linetype = "dashed") +
  geom_vline(xintercept = which(floor(cumsum(EVD$values) / sum(EVD$values) * 10) == 8)[1], linetype = "dashed")

# Combine the two plots
plots <- grid.arrange(plt1, plt2, plt3, plt4, plt5)
ggsave("../../output/e_Eigen.pdf", plots)
