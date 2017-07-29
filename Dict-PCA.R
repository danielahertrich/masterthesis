library(tictoc)

## load data
load(file="faces400.rda")
load(file="faces100.rda")
Xtrain <- rbind(test[51:100,],X[1:350,])

################################################################################
## Compute PCA approximation

## set parameter
q <- 100 #300

## learn dictionary and measure time
tic.clearlog()
tic("Total")
sv <- svd(Xtrain)
A <- sv$u[,1:q,drop=FALSE] %*% diag(sv$d[1:q])
H <- t(sv$v[,1:q,drop=FALSE])
toc(log = TRUE)
timing <- tic.log(format = TRUE)
save(timing, file=paste("time_PCA",q,".rda",sep=""))

## save dictionary
spdict <- list(A=A,H=H)
save(spdict,file = paste("dict_PCA",q,".rda",sep = ""))
