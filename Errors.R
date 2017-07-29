set.seed(123)

## load data
load(file="faces400.rda")
load(file ="faces100.rda")
Xtest <- rbind(test[1:50,],X[351:400,])
Xtrain <- rbind(test[51:100,],X[1:350,])
ind <- sample(1:400,10,replace = FALSE)
Xtr <- Xtrain[ind,]

## set parameters
art <- "NMF" #"AMr","AMp","OLr","OLp", "PCA"
q <- 100 #300, 600
it <- 50 # 150

## load dictionary
load(file=paste("dict_",art,q,"_it",it,".rda",sep = ""))
# load(file=paste("dict_",art,q,".rda",sep = "")) # for art="PCA"
dict <- spdict$H
Z <- t(dict)
score <- spdict$W

## load bhat_test
load(file=paste("bhat_test",art,q,"it",it,".rda",sep = ""))
## 1=r10, 2=r20, 3=r50, 4=r75, 5=p10, 6=p20, 7=p50, 8=p75, 9=Mp10, 10=Mp20, 11=Mp50, 12=Mp75

####################################################
## check l_2 error and sparsity for test images
####################
bhat <- bhat_test[1,,]
Yhat <- Z%*%bhat
# miss <- rand10  #1
# miss <- rand20  #2
# miss <- rand50  #3
# miss <- rand75  #4
# miss <- patch10  #5
# miss <- patch20  #6
# miss <- patch50  #7
# miss <- patch75  #8
# miss <- Mpatch10  #9
# miss <- Mpatch20  #10
# miss <- Mpatch50  #11
# miss <- Mpatch75  #12

err <- errobs <- errm <- sp1 <- sp2 <- numeric(nrow(Xtest))
for(i in 1:nrow(Xtest)){
  err[i] <- mean((Xtest[i,]-Yhat[,i])^2) ## total error
  errm[i] <- mean((Xtest[i,miss]-Yhat[miss,i])^2)  ## error on missing pixels
  errobs[i] <- mean((Xtest[i,-miss]-Yhat[-miss,i])^2) ## error on observed pixels
  sp1[i] <-sum(bhat[,i]==0)  
  sp2[i] <-sum(abs(bhat[,i])<=0.005)
}

## look at range and mean of errors and sparsity
range(err); round(mean(err),4)
range(errobs); round(mean(errobs),4)
range(errm); round(mean(errm),4)
range(sp1); mean(sp1)
range(sp2); mean(sp2)


####################################################
## Errors for Table 5.14
####################

## load bhat_train
load(file=paste("bhat_train",art,q,"it",it,".rda",sep = ""))

## column1 (approximations)
Xapprox <- score%*%dict
err_app <- numeric(nrow(Xtrain))
for(i in 1:nrow(Xtrain)){err_app[i] <- mean((Xtrain[i,] - Xapprox[i,])^2)}
round(mean(err_app),4)

err_test <- matrix(nrow=nrow(Xtest),ncol = 12); err_train <- matrix(nrow=nrow(Xtr),ncol = 12)
for(j in 1:12){
  btest <- bhat_test[j,,]; btrain <- bhat_train[j,,]
  Ytest <- Z%*%btest; Ytrain <- Z%*%btrain
  for(i in 1:nrow(Xtest)){err_test[i,j] <- mean((Xtest[i,]-Ytest[,i])^2)}
  for(i in 1:nrow(Xtr)){err_train[i,j] <- mean((Xtr[i,]-Ytrain[,i])^2)}
}
round(c(mean(apply(err_test,2,mean)),mean(apply(err_train,2,mean))),4)
