### 1 - Data Load ###

loadData <- function(phenos.path, marker.path = NULL) {
  markers <- NULL
  phenos <- NULL
  
  if (!is.null(marker.path)) {
    # Load the covariate matrix file
    if(grepl("\\.rda$", marker.path)) {
      markers <- get(load(marker.path))
    } else if (grepl("\\.csv$", marker.path)) {
      markers <- read.csv(markers.path)
    } else {
      stop("Invalid file format for covariate matrix file, Only CSV and RDA files are supported.")
    }
  }
  
  if (!is.null(phenos.path)) {
    # Load the phenos/environment file
    if(grepl("\\.rda$", phenos.path)) {
      phenos <- get(load(phenos.path))
    } else if (grepl("\\.csv$", phenos.path)) {
      phenos <- read.csv(phenos.path)
    } else {
      stop("Invalid file format for covariate matrix file, Only CSV and RDA files are supported.")
    }
  }
  
  return(list(markers = markers, phenos = phenos))  
}

# These files are the output of module 1
createNaNFiles <- function(phenos, markers, nan.freq) {
  # Use args[1] for NaN threshold, # of strains with missing values to keep
  nan.freq.i <- length(unique(phenos[["strain"]])) * nan.freq
  NaNs <- colSums(is.na(markers))
  
  index.1 <- order(NaNs)
  NaN.freq <- unique(NaNs[index.1])
  percentage.NaNs <- NaN.freq * 100 / nrow(markers)
  
  # Write the NaN frequencies and percentages to separate CSV files
  write.table(NaN.freq, file = '../../tmp/NaNs_freq.csv', sep = ',', row.names = FALSE, col.names = FALSE)
  cat("Saved the NaN.freq file in the tmp directory\n")
  
  write.table(percentage.NaNs, file = '../../tmp/percentage_NaNs.csv', sep = ',', row.names = FALSE, col.names = FALSE)
  cat("Saved the percentage.NaNs file in the tmp directory\n")
  
  # Subset the markers matrix based on a condition and save it as 'X'
  index.2 <- NaNs <= nan.freq.i
  no.NaNs.markers <- markers[, index.2]
  write.csv(no.NaNs.markers, file = '../../tmp/no_nans_markers.csv', row.names = FALSE)
  cat("Saved the no.NaNs.markers file in the tmp directory\n")
  
}