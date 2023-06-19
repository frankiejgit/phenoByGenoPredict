#' Author: Francisco Gonzalez
#' Creation Date: 06/07/2023
#' Last Updated: 06/08/2023

# Print out logs 
logMessage <- function(message, status="INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  print(paste(status,"- [", timestamp, "]", message))
}

# Get mean values for each phenotypic traits, classified by genotype
getMeansByGenotype <- function(data, traits, geno.col) {
  col.name <- geno.col
  
  # Check if geno_col is a number or not
  if (is.numeric(geno.col)) {
    col.name <- colnames(data)[geno.col]
    logMessage("Converted geno_col in getMeansByGenoType() from index to name")
  }
  
  trait.names <- colnames(data)[traits]
  unique.geno <- unique(data[[col.name]])
  
  ret.data <- matrix(NA, nrow=length(unique.geno), ncol=length(trait.names),
                     dimnames = list(unique.geno, trait.names))

  for (t in trait.names) {
    col_mean <- data %>%
      group_by(!!sym(col.name)) %>%
      summarize(!!paste0("mean_", t) := mean(!!as.name(t), na.rm = TRUE))
    
    ret.data[, t] <- as.numeric(col_mean[[paste0("mean_", t)]])
  }
  
  return(ret.data)
  
}