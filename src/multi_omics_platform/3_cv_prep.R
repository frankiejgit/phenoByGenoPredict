CV1 = TRUE
CV2 = TRUE

cvPrep <- function(phenos, output.path, col.id=NULL, col.folds=NULL, folds = 10, cv1 = TRUE,
                   cv2 = TRUE, cv0 = TRUE, cv00 = TRUE) {
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output.path)) { dir.create(output.path) }
  
  IDs <- unique(phenos[, col.id])
  
  if (cv1) {
    fold <- rep(1:folds,each=ceiling(length(IDs)/folds))[order(runif(length(IDs)))]
    phenos$CV1 <- NA
    for (f in 1:folds) {
      fold.lines <- IDs[fold == f]
      phenos$CV1[phenos[,col.id] %in% fold.lines] <- f
    }
  }
  
  if (cv2) {
    phenos$CV2 <- NA
    
    for (id in IDs) {
      line.idx <- which(phenos[, col.id] == id)
      num.lines <- length(line.idx)
      fold <- sample(1:folds, size=num.lines, replace = num.lines > folds)
      phenos$CV2[line.idx] <- fold
    }
  }
  
  
  if (cv0) {
    # Specify col for phenotypic trait and for folds
    # Read environmental covariates

    tf <- matrix(NA, nrow=nrow(phenos), ncol=folds)
    tf[,] <- phenos[, col.id]
    colnames(tf) <- paste0(colnames(phenos)[col.id], '_CV0_fold_', 1:folds)
    
    # delete the existing columns to override with new cv values
    phenos <- phenos[, -which(colnames(phenos) %in% colnames(tf))]
    
    for (f in 1:folds) {
      tf[phenos[,col.folds] == f, f] <- NA
    }
    
    phenos <- cbind(phenos, tf)
    
  }
  
  if (cv00) {
    # Specify col for phenotypic trait and for folds
    # Read environmental covariates
    
    tf <- matrix(NA, nrow=nrow(phenos), ncol=folds)
    colnames(tf) <- paste0(colnames(phenos)[col.id], '_CV00_fold_', 1:folds)
    tf[,] <- phenos[,col.id]
    
    for (f in 1:folds) {
      tf[phenos[,col.folds] == f, f] <- NA
    }
    
    phenos <- cbind(phenos, tf)
  }
  
  # Save the latest phenos file
  write.table(phenos, file=paste0(output.path, 'phenos_cv.csv'), sep=',', row.names=F, col.names=T)
  
  return(phenos)
  
}



