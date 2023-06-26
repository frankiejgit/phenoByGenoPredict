runBGLR <- function(phenos, phen.col, var.col, cv.col, env.col = NULL, file_list = NULL, 
                    folds = 5, cv0 = FALSE, esc=FALSE, nIter = 5000, burnIn = 500) {
  
  # Check which model to use
  items <- length(file_list)
  if (items == 2) { 
    mod <- "E+L" 
  } else if (items == 3) {
    mod <- "E+L+G"
  } else if (items == 4) {
    mod <- "E+L+G+GE"
  }

  # Get values for model fitting
  phen.name <- strsplit(colnames(phenos)[phen.col], "_")[[1]][1]
  y <- phenos[, phen.col]
  gid <- phenos[, var.col]
  
  if (esc) { y <- scale(y, center = TRUE, scale = TRUE) }
  
  # TODO -- make this option an argument for end user
  models <- c('FIXED','FIXED', 'RKHS', 'RKHS')
  
  eta <- list()
  for (i in seq_along(file_list)) {
    Z <- get(load(file_list[[i]]))
    eta[[i]] <- list(X=Z, model='FIXED')
    
    if (models[i] == "RKHS") {
      eta[[i]] <- list(V = EVD$vectors, d=EVD$values, model=models[i])
      rm(EVD)
    }
    
    rm(Z)
  }
  
  for (fold in 1:folds) {
    if (fold != -999) {
      
      output.path <- paste("../../output/", phen.name, "/", mod, "/fold_", fold, "/", sep = "")
      if (!dir.exists(output.path)) { dir.create(output.path, recursive = TRUE) }
      
      testing <- which(phenos[, cv.col] == fold)
      
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
  
}
