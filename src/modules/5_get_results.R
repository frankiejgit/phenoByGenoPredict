phenos 
trait.col <- 3

model <- 1
cv <- 2
Env <- 1
trait <- 1
rep <- 5
fold <- 2

models <- rev(c('E+L','E+L+G','E+L+G+GE'))
envs <- 1:length(unique(phenos$EID))

 