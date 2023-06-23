# 2 - Matrices #

### For E matrix ###

# TODO make output.path full output path
createMatrixPdf <- function(output.path, phenos, col.env.id, 
                            markers = NULL, weight.file = NULL, 
                            prop.maf.j = NULL, ctr = TRUE, std = TRUE) {
  
  # Create the output directory if it doesn't exist
  if (!dir.exists(output.path)) { dir.create(output.path) }
  
  cat('== 1 ==> Reading and checking parameters\n')
  
  if(is.null(markers)) {
    cat('Are you working with a grouping factor?
       If not, please provide a path of the covariate matrix file\n')
    
    # Prompt the user for input
    user_input <- readline("Do you want to continue without the matrix file? (y/n): ")
    
    # Check the user's response
    if (tolower(user_input) == "n") {
      stop("Please provide the covariate matrix file and rerun the pipeline.\n")
    }
  }
  
  # Read data and perform consistency checks
  cat("== 2 ==> Reading data\n")
  
  if (!is.null(col.env.id)) {
    env.IDs <- phenos[, col.env.id]
  }
  
  if (is.null(markers)) {
    env.IDs <- factor(phenos[, col.env.id])
    n <- length(levels(env.IDs))
    Z <- as.matrix(model.matrix(~env.IDs - 1))
    d <- colSums(Z)
    V <- Z
    
    for (i in 1:ncol(Z)) {
      V[, i] <- V[, i] / sqrt(d[i])
    }
    
    EVD <- list(vectors = V, values = d)
    G <- tcrossprod(Z)
    
    save(G, file = paste0(output.path, "G.rda"))
    save(EVD, file = paste0(output.path, "EVD.rda"))
    
  } else {
    # Reads covariates
    mm.matrix <- markers
    s <- 0
    col.count <- ncol(mm.matrix)
    
    if (weighting) {
      weight <- scan(weight.file, skip =1 ) # TODO check with Dr. Jarquin about this
      
      # Naive imputation
      for (i in 1:col.count) {
        mean.i <- mean(mm.matrix[, i], na.rm=TRUE)
        mm.matrix[,i] <- ifelse(is.na(mm.matrix[,i]), mean.i, mm.matrix[,i])
        s <- s + var(mm.matrix[,i])
      }
      
      p <- colMeans(mm.matrix, na.rm=TRUE) / 2
      p <- ifelse(p <= 0.5, p, 1-p)
      
      # Adjusting for proportion of MAF
      if (!is.null(prop.maf.j)) {
        maf.idx <- which( p >= prop.maf.j )
        mm.matrix <- mm.matrix[, maf.idx]
      }
      
      # Imputing, centering, standardizing and weighting
      for (i in 1:col.count) {
        mean.i <- mean(mm.matrix[,i], na.rm=TRUE)
        if (ctr) { mm.matrix[,i] <- mm.matrix[,i] - mean.i }            # Centering
        if (std) { mm.matrix[,i] <- mm.matrix[,i]/sd(mm.matrix[,i]) }   # Standardizing
        if (weighting) { mm.matrix[,i] <- mm.matrix[,i] * weight[i] }   # Weighting
      }
      
      G <- tcrossprod(mm.matrix) / s
      
      # Check all IDs are in G
      if(!is.null(col.env.id)){
        stopifnot(all(env.IDs%in%rownames(G))) 
      }
      
      if(!is.null(col.env.id)){
        env.IDs <- factor(env.IDs, levels=rownames(G))
        Z <- as.matrix(model.matrix(~env.IDs-1))
        G <- tcrossprod(tcrossprod(Z,G),Z)
      }
      
      EVD <- eigen(G)
      rownames(EVD$vectors) <- rownames(G)
      
      save(G, file = paste0(output.path, "G.rda"))
      save(EVD, file = paste0(output.path, "EVD.rda"))
      
    }
    
  }
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
    cat("Eigen.pdf has been created.")
    
    return(paste0("The results were saved in the following folder: ", output.path))
    
}


createZMatrix <- function(dataset, col.env.id, output.path) {
  # Create the output directory if it doesn't exist
  if (!dir.exists(output.path)) { dir.create(output.path) }
  
  # Create Z matrix
  env.IDs <- dataset[,col.env.id]
  z0 <- model.matrix(~as.factor(env.IDs)-1)
  z <- as.matrix(z0[1:nrow(z0),1:ncol(z0)])
  
  # Save files in output dir
  save(z, file = paste0(output.path, 'Z.rda'))
  
  pdf(paste0(output.path, 'exp_des.pdf'))
  image(z,col=c('white','black'))
  dev.off()
  
}



