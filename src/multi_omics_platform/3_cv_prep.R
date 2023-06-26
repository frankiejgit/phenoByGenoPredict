CV1 = TRUE
CV2 = TRUE

cvPrep <- function(phenos, col.id, folds = 10, cv1 = TRUE, cv2 = TRUE) {
  IDs <- unique(phenos[, col.id])
  fold <- rep(1:folds,each=ceiling(length(IDs)/folds))[order(runif(length(IDs)))]
  
  if (cv1) {
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
  
  # Save the latest phenos file
  write.table(phenos, file='../../output/phenos_cv.csv', sep=',', row.names=F, col.names=T)
  
  return(phenos)
  
}




# Save the latest phenos file
write.table(phenos, file='../../output/phenos_cv.csv', sep=',', row.names=F, col.names=T)




