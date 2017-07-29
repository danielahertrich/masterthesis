###############################################################################
### Non-negative Matrix Factorization
###############################

library(NNLM)
library(tictoc)
set.seed(42)

## load data
load(file="faces400.rda")
load(file="faces100.rda")
Xtrain <- rbind(test[51:100,],X[1:350,])

## set the parameters
lam_A <- 0.01
q <- 100 #300,600
it <- 50 #150

## learn dictionary and measure time
tic.clearlog()
tic("Total")
spdict <- nnmf(Xtrain, k=q, alpha=c(0,0,lam_A), beta=c(0,0,0), max.iter = it)
toc(log = TRUE)

## store time information
timing <- tic.log(format = TRUE)
save(timing, file=paste("time_NMF",q,"_it",it,".rda",sep=""))

## delete zero rows (and corresponding scores)
sel <- which(apply(spdict$H,1,sd)!=0)
spdict$H <- spdict$H[sel,]
spdict$W <- spdict$W[,sel]

## save dictionary
save(spdict,file=paste("dict_NMF",q,"_it",it,".rda",sep=""))