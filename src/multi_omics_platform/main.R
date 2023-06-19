library(data.table)

### ARGUMENTS TO PASS ###

args <- commandArgs(trailingOnly = TRUE)
nan.freq <- 0.2 # 1 - NaN threshold frequency

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

