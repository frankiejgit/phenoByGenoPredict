runBGLR <- function(phenos, col.phen, col.var, col.cv, col.env = NULL,folds = 5,
                    cv0 = FALSE, esc=FALSE, nIter = 5000, burnIn = 500) {
  
  # Get values for model fitting
  phen.name <- strsplit(colnames(phenos)[col.phen], "_")[[1]][1]
  y <- phenos[, col.phen]
  gid <- phenos[, col.var]
  
  if (esc) { y <- scale()}
  
}


phen.name <- strsplit(colnames(phenos.cv)[col.phen], "_")[[1]][1]
y   <- phenos.cv[, col.phen]
gid <- phenos.cv[, col.var]

if(esc) { y <- scale(y, center=TRUE, scale=TRUE) }

z.list <- list()
z.list[[1]] <- '../../output/ZE/Z.rda'
z.list[[2]] <- '../../output/ZL/Z.rda'

models <- c('FIXED')     # Option to add more

eta <- list()
for (i in seq_along(z.list)) {
  Z <- get(load(z.list[[i]]))
  eta[[i]] <- list(X=Z, model='FIXED')
  rm(Z)
}

for (fold in 1:folds) {
  if (fold != -999) {
    
    output.path <- paste("../../output/", phen.name, "/E+L/fold_", fold, "/", sep = "")
    if (!dir.exists(output.path)) { dir.create(output.path, recursive = TRUE) }
    
    testing <- which(phenos.cv[, col.cv] == fold)
    
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

# E + L + G
evd.path <- '../../output/G/EVD.rda'  
cv0 <- FALSE
ESC <- FALSE 
set.seed(1)

if(esc) { y <- scale(y, center=TRUE, scale=TRUE) }

z.list <- list()
z.list[[1]] <- '../../output/ZE/Z.rda'
z.list[[2]] <- '../../output/ZL/Z.rda'
z.list[[3]] <- evd.path

models <- c('FIXED','FIXED','RKHS')     # Option to add more

eta <- list()
for (i in seq_along(z.list)) {
  
  Z <- get(load(z.list[[i]]))
  eta[[i]] <- list(X=Z, model='FIXED')
  
  if (models[i] == "RKHS") {
    eta[[i]] <- list(V = EVD$vectors, d=EVD$values, model=models[i])
    rm(EVD)
  }
  
  rm(Z)
  
}

for (fold in 1:folds) {
  if (fold != -999) {
    
    output.path <- paste("../../output/", phen.name, "/E+L+G/fold_", fold, "/", sep = "")
    if (!dir.exists(output.path)) { dir.create(output.path, recursive = TRUE) }
    
    testing <- which(phenos.cv[, col.cv] == fold)
    
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





nk <-length(AB)
ETA<-list(nk)
for(i in 1:nk){
  
  if(type[i]=='BRR')
  {
    load(AB[[i]])
    ETA[[i]] <- list(X=Z,model='BRR')
    rm(Z)
  }
  
  if(type[i]=='FIXED')
  {
    load(AB[[i]])
    ETA[[i]] <- list(X=Z,model='FIXED')
    rm(Z)
  }
  
  if(type[i]=='RKHS')
  {
    load(AB[[i]])
    ETA[[i]] <- list(V=EVD$vectors,d=EVD$values,model='RKHS')
    rm(EVD)
  }
  
  if(type[i]=='BayesA')
  {
    load(AB[[i]])
    ETA[[i]] <- list(X=X,model='BayesA')
    rm(X)
  }
  
  if(type[i]=='BayesB')
  {
    load(AB[[i]])
    ETA[[i]] <- list(X=X,model='BayesB')
    rm(X)
  }
  
  if(type[i]=='BayesC')
  {
    load(AB[[i]])
    ETA[[i]] <- list(X=X,model='BayesC')
    rm(X)
  }
  
  if(type[i]=='BL')
  {
    load(AB[[i]])
    ETA[[i]] <- list(X=X,model='BL')
    rm(X)
  }  
  print(i)
}


for(fold in folds)
{
  
  yNA<-y
  print(fold)
  
  if(fold != -999)
  {
    
    dir.create(paste('fold_',fold,sep=''))
    setwd(paste('fold_',fold,sep=''))
    
    testing=which(Y[,colCV]==fold)
    
    if(CV0)
    {         
      testing <- which(gid %in% gid[testing])   
    }  
    
    
    yNA=y
    yNA[testing]=NA
    
    
    fm=BGLR(y=yNA,ETA=ETA,nIter=nIter,burnIn=burnIn,verbose=TRUE)
    fm$y=y
    
    predictions=data.frame(testing,Individual=gid[testing], y=y[testing], yHat=fm$yHat[testing])
    
    write.table(predictions,file=paste("predictions_",fold,".csv",sep=""),row.names=FALSE,sep=",") # Change to a unique name?
    
  }else{
    
    dir.create('fullData')
    setwd('fullData')
    
    fm=BGLR(y=y,ETA=ETA,nIter=nIter,burnIn=burnIn,verbose=TRUE)
    save(fm,file='fm_full.RData')
    
  }          
  
  print(str(fm))
  rm(fm)
  
  
  unlink("*.dat")
  
  setwd('..')
  
}







