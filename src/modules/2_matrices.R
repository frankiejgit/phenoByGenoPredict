### 2 - Generate Matrices ###

generateMatrix <- function(output.path, phenos, col.env.id, markers=NULL, 
                           weight.file=NULL, prop.maf.j=NULL, ctr=TRUE, std=TRUE) {
  
  env.IDs <- NULL
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output.path)) { dir.create(output.path) }
  
  cat('== 1 ==> Checking parameters\n')
  
  if (is.null(markers)) {
    cat("Continuing without a covariate matrix file\n")
  }
  
  cat("== 2 ==> Reading data\n")
  
  if (!is.null(col.env.id)) {
    env.IDs <- phenos[, col.env.id]
  }
  
  if(is.null(markers)) {
    env.IDs <- factor(phenos[, col.env.id])
    Z <- as.matrix(model.matrix(~env.IDs - 1))
    d <- colSums(Z)
    V <- Z / sqrt(d)
    
    EVD <- list(vectors = V, values = d)
    G <- tcrossprod(Z)
    
    cat('== 3 ==> Saving RDA files\n')
    
    save(G, file = paste0(output.path, "G.rda"))
    save(EVD, file = paste0(output.path, "EVD.rda"))
    
  } else {
    env.IDs <- factor(phenos[, col.env.id])
    
    # Read covariates
    cat("reading covariates\n")
    mm.matrix <- matrix(NA, nrow = nrow(markers), ncol = ncol(markers))
    s <- 0
    col.count <- ncol(markers)
    
    cat("checking weights\n")
    if (weighting) {
      weight <- scan(weight.file, skip =1 ) # TODO check with Dr. Jarquin about this
    }
    
    # Naive imputation
    cat("naive imputation\n")
    for (i in 1:col.count) {
      mean.i <- mean(markers[, i], na.rm=TRUE)
      mm.matrix[,i] <- ifelse(is.na(markers[,i]), mean.i, markers[,i])
      s <- s + var(mm.matrix[,i])
    }
    
    rownames(mm.matrix) <- rownames(markers)
    colnames(mm.matrix) <- colnames(markers)
    
    p <- colMeans(mm.matrix, na.rm=T) / 2
    p <- ifelse(p <= 0.5, p, 1-p)
    
    # Adjusting for MAF proportion
    cat("adjusting for MAF\n")
    if (!is.null(prop.maf.j)) {
      maf.idx <- which(p >= prop.maf.j)
      mm.matrix <- mm.matrix[, maf.idx]
    }
    
    # Imputing, centering, standardizing, and weighting
    cat("Imputing, centering, standardizing, and weighting\n")
    for (i in 1:col.count) {
      mean.i <- mean(mm.matrix[,i] , na.rm=TRUE)
      if (ctr) { mm.matrix[,i] <- mm.matrix[,i] - mean.i }            # Centering
      # TODO - discuss issue with Inf values
      if (std) {                                                      # Standardizing
        sd_i <- sd(mm.matrix[, i])
        if (sd_i == 0) {
          # Handle zero standard deviation
          # Set specific value or omit rows
          # Example: Set the column to zero
          mm.matrix[, i] <- 0
        } else {
          mm.matrix[, i] <- mm.matrix[, i] / sd_i
        }
      }
      if (weighting) { mm.matrix[,i] <- mm.matrix[,i] * weight[i] }   # Weighting
    }
    
    G <- tcrossprod(mm.matrix) / s
    
    # Check that all IDs are in G
    if(!is.null(col.env.id)) {
      stopifnot(all(env.IDs %in% rownames(G)))
    }
    
    if(!is.null(col.env.id)){
      env.IDs <- factor(env.IDs, levels=rownames(G))
      Z <- as.matrix(model.matrix(~env.IDs-1))
      G <- tcrossprod(tcrossprod(Z,G),Z)
    }
    
    EVD <- eigen(G)
    rownames(EVD$vectors) <- rownames(G)
    
    cat('== 3 ==> Saving RDA files\n')
    
    save(G, file = paste0(output.path, "G.rda"))
    save(EVD, file = paste0(output.path, "EVD.rda"))
  }
  
  cat(' == 4 ==> Creating EVD plot file\n')
  createEvdPdf(EVD, output.path)
  
  return(paste0("The results were saved in the following folder: ", output.path))
  
}

createEvdPdf <- function(EVD, output.path) {
  cat("Initiating plot creation")
  
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
  ggsave(paste0(output.path, "Eigen.pdf"), plots)
}

createZMatrix <- function(dataset, col.env.id, output.path) {
  # Create the output directory if it doesn't exist
  if (!dir.exists(output.path)) { dir.create(output.path) }
  
  # Create Z matrix
  env.IDs <- dataset[,col.env.id]
  z0 <- model.matrix(~env.IDs-1)
  z <- as.matrix(z0[1:nrow(z0),1:ncol(z0)])
  
  # Save files in output dir
  save(z, file = paste0(output.path, 'Z.rda'))
  
  pdf(paste0(output.path, 'exp_des.pdf'))
  image(z,col=c('white','black'))
  dev.off()
  
}

generateIntMatrix <- function(g1.file, g2.file, output.path) {
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output.path)) { dir.create(output.path) }
  
  load(g1.file)
  G1 <- G
  load(g2.file)
  G2 <- G
  
  GI <- G1 * G2
  EVD <- eigen(GI)
  
  save(GI, file = paste0(output.path, "G.rda"))
  save(EVD, file = paste0(output.path, "EVD.rda"))
  
  createEvdPdf(EVD, output.path)
  
  return("Interaction matrix has been created")
  
}




