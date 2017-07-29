###############################################################################
### Online Learning
###############################

library(progress)
library(lars)
library(tictoc)

## dictionary update
updateD <- function(H0,B0,C0,eps=1e-06){
  H <- H0; B <- B0; C <- C0
  u <- numeric(ncol(H))
  if(sum(sqrt(apply(H^2,1,sum))>1+eps)>0){stop("Dictionary not in C")}
  for(i in 1:nrow(C)){
    if((abs(C[i,i])<=1e-08)==TRUE){C[i,i] <- C[i,i]+eps}
  }
  for(k in 1:5){
    for(j in 1:nrow(H)){
      u <- 1/C[j,j]*(B[,j]-t(H)%*%C[,j]) + H[j,]
      H[j,] <- 1/max(sqrt(sum(u^2)),1)*u
    }
    k <- k+1
  }
  return(H)
}
## OL algorithm
onlineL <- function(X, H0, lambda, iter){
  pb <- progress_bar$new(format = "  progress [:bar] :percent",
                         total = iter, clear = FALSE, width= 60)
  pb$tick(0)
  set.seed(133)
  index <- sample(1:nrow(X),replace = FALSE)
  X <- X[index,]
  H <- H0
  B <- matrix(0,nrow=ncol(H),ncol=nrow(H))
  C <- matrix(0,nrow=nrow(H),ncol=nrow(H))
  alpha <- matrix(0, nrow=nrow(H),ncol=iter)
  for(t in 1:iter){
    pb$tick()
    tic(paste("Update scores in iteration",t))
    a.lars <- lars(t(H),X[t,])
    alpha[,t] <- coef(a.lars, s=lambda, mode = 'lambda')
    toc(log=TRUE,quiet = TRUE)
    tic(paste("Update H in iteration",t))
    B <- B + as.matrix(X[t,])%*%alpha[,t]
    C <- C + as.matrix(alpha[,t])%*%alpha[,t]
    H <- updateD(H0=H, B0=B, C0=C)
    toc(log = TRUE,quiet = TRUE)
    t <- t+1
  }
  spdict <- list(H=H, Alpha=alpha)
  return(spdict)
}

###############################################################################
### Use OL Algorithm to Learn Dictionaries
###############################

## load data
load(file="faces400.rda") #contains matrix X
load(file="faces100.rda") #contains matrix test
Xtrain <- rbind(test[51:100,],X[1:350,])

set.seed(42)

## set parameters
lam_A <- 0.01
q <- 100 #300,600
it <-50  # 150

## choose start dictionary
art <- "OLp" # "OLr"
load(file=paste("dict_PCA",q,".rda",sep=""))
Hpca <- spdict$H
# Hrand <- matrix(rnorm(q*ncol(Xtrain)),nrow=q,ncol=ncol(Xtrain))
# for(i in 1:nrow(Hrand)){
#   Hrand[i,] <- 1/max(sqrt(sum(Hrand[i,]^2)),1)*Hrand[i,]
# }

## learn dictionary and measure time
tic.clearlog()
tic("Total")
spdict <- onlineL(X=Xtrain, H0=Hpca,lambda = lam_A, iter = it)
toc(log = TRUE)

## store time information
log.txt <- tic.log(format = TRUE)
n <- length(log.txt)
upA<-upH<-character(length = (n-1)/2)
for(i in seq(1,(n-1)/2,1)){upA[i]<-log.txt[2*i-1];upH[i]<-log.txt[2*i]}
timing <- list(A=upA,H=upH,total=log.txt[n])
save(timing, file=paste("time_OL",q,"p_it",it,".rda",sep=""))

## delete zero rows
sel <- which(apply(spdict$H,1,sd)!=0)
spdict$H <- spdict$H[sel,]

## compute score matrix A
A <- matrix(nrow = nrow(Xtrain),ncol = nrow(spdict$H))
a.lars <- apply(Xtrain,1,lars,x=t(spdict$H))
a <- lapply(a.lars,coef, s=lam_A, mode='lambda')
for(i in 1:nrow(A)){
  A[i,] <- a[[i]]
}
spdict <- list(H=spdict$H,A=A,Alpha=spdict$Alpha)

## save dictionary
save(spdict, file=paste("dict_",art,q,"_it",it,".rda",sep=""))