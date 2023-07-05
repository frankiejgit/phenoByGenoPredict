### 4 - Predict trait values using BGLR ###

getPredictions <- function(data, cv, phen.col, gid.col, cv.col, file.list, 
                           folds = 5, esc = FALSE, nIter = 5000, burnIn = 500) {
  
  # Check which model to use
  items <- length(file.list)
  if (items == 2) { 
    mod <- "E+L" 
  } else if (items == 3) {
    mod <- "E+L+G"
  } else if (items == 4) {
    mod <- "E+L+G+GE"
  }
  
  # Initialize return object
  predictions <- data.frame(testing = integer(), fold = integer(), individual = character(),
                            observed = numeric(), predicted = numeric())
  
  # Extract target variable and the id of the variety
  y <- data[, phen.col]               
  gid <- data[, gid.col]
  
  # Apply scaling to y if esc is enabled
  if (esc) { y <- scale(y, center = TRUE, scale = TRUE) }
  
  # Assign model type to each data
  models <- c('FIXED','BRR','RKHS','RKHS')    # TODO: give user choice of model selection
  eta <- list()
  
  for (i in seq_along(file.list)) {
    # Add data and assigned model to ETA
    Z <- get(load(file.list[[i]]))
    eta[[i]] <- list(X=Z, model=models[i])
    
    if (models[i] == "RKHS") {
      eta[[i]] <- list(V = EVD$vectors, d=EVD$values, model=models[i])
      rm(EVD)
    }
    
    rm(Z)
  }
  
  # Perform CV1 and CV2
  if (cv %in% c("cv1","cv2")) {
    predictions <- fitCv1Cv2(data, y, gid, folds, predictions, cv.col, eta, nIter, burnIn)
  }
  
  if (cv %in% c("cv0","cv00")) {
    predictions <- fitCv0Cv00(data, y, gid, cv, predictions, cv.col, eta, nIter, burnIn)
  }
  
  # Clean up workspace
  unlink("*.dat")
  
  # Save predictions as a file
  output.path <- paste("../output/", colnames(data)[phen.col], "/", mod, "/", cv, "/", sep = "")
  if (!dir.exists(output.path)) { dir.create(output.path, recursive = TRUE) }
  
  write.table( predictions, file = paste0(output.path, "predictions.csv"),
               row.names = FALSE, sep = ",")
  
  return(predictions)
  
}

fitCv0Cv00 <- function(data, y, gid, cv, predictions, cv.col, eta, nIter, burnIn) {
  for (i in cv.col) {
    # Get the train-test data (test data is shown with NA value)
    y.na <- data[, i]
    
    # Record indices of all rows used for testing
    testing <- which(is.na(y.na))
      
    if (tolower(cv) == "cv0") { testing <- intersect(testing, which(gid %in% gid[testing])) }
    
    # Get the current fold
    curr.fold <- tail(strsplit(colnames(data)[cv.col], "_")[[1]], 1)
    
    # Fit the BGLR model with the training data
    fm <- BGLR(y = y.na, ETA = eta, nIter = nIter, burnIn = burnIn, verbose = TRUE)
      
    # Generate predictions for testing samples
    predi <- data.frame(testing = testing, fold = curr.fold, 
                        individual = gid[testing], observed = y[testing],
                        predicted = fm$yHat[testing])
      
    # Add predictions to a dataset keeping track of results from all folds
    predictions <- rbind(predictions, predi)
    } 
  
  # Clean up workspace
  rm(fm)
  closeAllConnections()
  return(predictions)
  
}

fitCv1Cv2 <- function(data, y, gid, folds, predictions, cv.col, eta, nIter, burnIn) {
  for (fold in 1:folds) {
    y.na <- y
    
    if (fold != -999) {
      # Get the indices for rows to test on for a given fold
      testing <- which(data[, cv.col] == fold)     # CV0/00 = env.col, CV1 = 13, CV2 = 14
      
      # This is the testing data, all indices in testing have been marked as NA
      y.na[testing] <- NA
      
      # Fit BGLR model with the training data
      fm <- BGLR(y = y.na, ETA = eta, nIter = nIter, burnIn = burnIn, verbose = TRUE)
      
      # Generate predictions for testing samples
      predi <- data.frame(testing = testing, fold = fold, 
                          individual = gid[testing], observed = y[testing],
                          predicted = fm$yHat[testing])
      
      # Add predictions to a dataset keeping track of results from all folds
      predictions <- rbind(predictions, predi)
      
    } else {
      fm <- BGLR(y = y, ETA = ETA, nIter = nIter, burnIn = burnIn, verbose = TRUE)
      save(fm, file = 'fm_full.RData')
      predictions <- NULL
    }
    
  }
  
  # Clean up workspace
  rm(fm)
  closeAllConnections()
  return(predictions)
}