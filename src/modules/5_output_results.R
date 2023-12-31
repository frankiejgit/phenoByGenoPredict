# 5 - Generate final report ###

getCvResults <- function(data, env.col, trait.col, source.dir) {
  # Define models, environment, and trait
  models <- rev(c('E+L', 'E+L+G', 'E+L+G+GE'))
  envs <- 1:max(data[, env.col])
  trait <- colnames(data)[trait.col]
  
  rep <- 1 # TODO: Pass as user argument later
  mean.col <- rep + 1
  sd.col <- rep + 2

  source.path <- file.path(source.dir, "output")
  
  # Loop over models
  for (mod in seq_along(models)) {
    
    for (cv in c("cv1","cv2", "cv0", "cv00")) {
      
      # Initialize matrices
      tf2 <- matrix(NA, nrow = max(envs), ncol = sd.col)
      tf <- matrix(NA, nrow = max(envs), ncol = 3)

      # Initialize temporary matrix 
      tm <- matrix(NA, nrow = nrow(data), ncol = sd.col + 4)
      path <- file.path(source.path, trait, models[mod], cv, "predictions.csv")
      
      # Fill in temp matrix with values
      if(file.exists(path)) {
        pred <- read.csv(path, stringsAsFactors = FALSE)
        print(dim(pred))
        tm[pred$testing, 1:5] <- as.matrix(pred)
        # TODO: think of how to get env name
        tm[pred$testing, 6:7] <- as.matrix(data[, c(1, env.col)])
      } else {
        print(paste("No data for", "rep", rep, "for", toupper(cv), sep = " "))
      }
      
      # Track model used and cross validation
      print(c(models[mod], cv))
      
      # Compute correlations and fill in matrices
      tm2 <- matrix(NA, nrow = max(envs), ncol = 1)
      
      for (env in seq_along(envs)) {
        idx <- as.numeric(tm[, 7]) == envs[env]
        tf[env, 2] <- length(which(idx))
        tf[env, 1] <- tm[idx, 7][1]
        
        if(all(is.na(as.numeric(tm[idx, 7])))) {
          print(paste("No data for Env", env, "for", toupper(cv), sep=" "))
        } else {
          # Put accuracy metrics -obs vs pred- in matrix
          tf[env, 3] <- cor(as.numeric(tm[idx, 4]), as.numeric(tm[idx, 5]), use="complete.obs")
          tm2[env] <- tf[env,3]
        }
        
        # Assign correlation to certain rep
        tf2[,1] <- tm2
        
      }
      
      # Compute summary statistics
      if(rep > 1) {
        tf2[, mean.col] <- apply(tf2[, 1:rep], 1, mean, na.rm = T)
        tf2[, sd.col] <- apply(tf2[, 1:rep], 1, sd, na.rm = T)
      } else {
        tf2[, mean.col] <- tf2[, 1]
        tf2[, sd.col] <- 0
      }
      
      tf3 <- cbind(tf[,-3], tf2)
      colnames(tf3) <- c("environment","sample_size","accuracy","mean","std")
      
      # Mean for all reps 
      ri <- as.numeric(tf2[, mean.col])
      # Variance inflation factor
      vif <- (1 - ri^2) / (as.numeric(tf2[, 2]) - 2)
      # Average VIF
      avc <- sum(ri / vif) / sum(1 / vif)
      
      # Write results to output as CSVs
      out.path <- paste(rep, cv, trait, models[mod], sep="_")
      out.dir <- file.path(source.path, 'report')

      if (!dir.exists(out.dir)) { dir.create(out.dir, recursive = TRUE) }
      
      write.csv(tf3, file = file.path(out.dir, paste('pa_', out.path, '.csv', sep = "")), row.names = FALSE)
      write.csv(avc, file = file.path(out.dir, paste('avc_', out.path, '.csv', sep= "")), row.names = FALSE)
      
    }
    
  }
}