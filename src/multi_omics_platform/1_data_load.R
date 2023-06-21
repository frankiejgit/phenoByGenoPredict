### 1 - Data Load ###

function loadData(marker.path, phenos.path) {
  # Load the covariate matrix file
  if(grepl("\\.rda$", marker.path)) {
    markers <- load(marker.path)
  } else if (grepl("\\.csv$", marker.path)) {
    markers <- read.csv(markers.path)
  } else {
    stop("Invalid file format for covariate matrix file, Only CSV and RDA files are supported.")
  }
  
  # Load the phenos/environment file
  if(grepl("\\.rda$", phenos.path)) {
    phenos <- load(marker.path)
  } else if (grepl("\\.csv$", phenos.path)) {
    phenos <- read.csv(phenos.path)
  } else {
    stop("Invalid file format for covariate matrix file, Only CSV and RDA files are supported.")
  }
  
}



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