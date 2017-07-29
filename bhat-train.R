library(glmnet)
library(progress)

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

## compute and plot approximation of Xtrain
Xapprox <- score%*%dict
png(filename = paste(art,q,"it50","approxXtrain.png",sep = ""), width = 600, height = 875)       #!!!#
par(mfrow=c(8,5), mar=rep(0,4))
for(i in seq(1,nrow(Xapprox),10)){
  image(matrix(Xapprox[i,],nrow=92), col=gray.colors(100),axes=FALSE)
}
dev.off()

## delete pixels  - create Y^(S)_(train) for S in M_o
set.seed(42)
rand10 <- sample(1:ncol(X),100,replace=FALSE); Ytrainr10 <- Xtr[,-rand10]
rand20 <- sample(1:ncol(X),400,replace=FALSE); Ytrainr20 <- Xtr[,-rand20]
rand50 <- sample(1:ncol(X),2500,replace=FALSE); Ytrainr50 <- Xtr[,-rand50]
rand75 <- sample(1:ncol(X),5625,replace=FALSE); Ytrainr75 <- Xtr[,-rand75]

m10 <- matrix(0, nrow=10, ncol=10)
for(i in 55:64){m10[i-54,] <- (i*92+25):(i*92+34)}
patch10 <- as.vector(m10); Ytrainp10 <- Xtr[,-patch10]
m20 <- matrix(0, nrow=20, ncol=20)
for(i in 55:74){m20[i-54,] <- (i*92+18):(i*92+37)}
patch20 <- as.vector(m20); Ytrainp20 <- Xtr[,-patch20]
m50 <- matrix(0, nrow=50, ncol=50)
for(i in 45:94){m50[i-44,] <- (i*92+5):(i*92+54)}
patch50 <- as.vector(m50); Ytrainp50 <- Xtr[,-patch50]
m75 <- matrix(0, nrow=75, ncol=75)
for(i in 25:99){m75[i-24,] <- (i*92+3):(i*92+77)}
patch75 <- as.vector(m75); Ytrainp75 <- Xtr[,-patch75]

Mm10 <- matrix(0, nrow=10, ncol=10)
for(i in 20:29){Mm10[i-19,] <- (i*92+45):(i*92+54)}
Mpatch10 <- as.vector(Mm10); YtrainMp10 <- Xtr[,-Mpatch10]
Mm20 <- matrix(0, nrow=20, ncol=20)
for(i in 15:34){Mm20[i-14,] <- (i*92+40):(i*92+59)}
Mpatch20 <- as.vector(Mm20); YtrainMp20 <- Xtr[,-Mpatch20]
Mm50 <- matrix(0, nrow=50, ncol=50)
for(i in 2:51){Mm50[i-1,] <- (i*92+25):(i*92+74)}
Mpatch50 <- as.vector(Mm50); YtrainMp50 <- Xtr[,-Mpatch50]
Mm75 <- matrix(0, nrow=75, ncol=75)
for(i in 2:76){Mm75[i-1,] <- (i*92+10):(i*92+84)}
Mpatch75 <- as.vector(Mm75); YtrainMp75 <- Xtr[,-Mpatch75]

#######################################################################################################
## Learn bhat for 10 randomly selected Training Samples
###############

lambdB <- 0.001

## initialize array to save matrices bhat for each S in M_o 
bhat_train <- array(0, dim=c(12,nrow(dict),nrow(Xtr)))

######################
miss <- rand10; Y <- Ytrainr10
## compute bhat
bhat<- matrix(0,nrow = nrow(dict),ncol = nrow(Xtr))
pb <- progress_bar$new(format = "progress: [:bar] :percent",total = nrow(Xtr), 
                       clear = FALSE, width= 60)
for(i in 1:nrow(Xtr)){
  pb$tick()
  lmax <- max(abs(drop(crossprod(Y[i,]-mean(Y[i,]),
                                 scale(Z[-miss,],colMeans(Z[-miss,]),FALSE)))))
  lamda <- exp(seq(log(1e-05), log(lmax), length.out = 100)) ## define lambda sequence
  b.lasso <- glmnet(Z[-miss,],Y[i,],lambda = lamda, alpha = 1)
  ind <- which.min(abs(b.lasso$lambda-lambdB))
  bhat[,i] <- b.lasso$beta[,ind]
}
bhat_train[1,,] <- bhat                                         

## reconstruct training images
Yhat <- Z%*%bhat
png(filename = paste("Pictures/Yr10-recon",art,q,"it",it,"_train.png",sep = ""), width = 770, height = 380)       #!!!#
par(mfrow=c(2,5), mar=rep(0,4))
for(i in 1:10){image(matrix(Yhat[,i],nrow=92), col=gray.colors(100),axes=FALSE)}
dev.off()

######################
miss <- rand20; Y <- Ytrainr20
## compute bhat
bhat<- matrix(0,nrow = nrow(dict),ncol = nrow(Xtr))
pb <- progress_bar$new(format = "progress: [:bar] :percent",total = nrow(Xtr), 
                       clear = FALSE, width= 60)
for(i in 1:nrow(Xtr)){
  pb$tick()
  lmax <- max(abs(drop(crossprod(Y[i,]-mean(Y[i,]),
                                 scale(Z[-miss,],colMeans(Z[-miss,]),FALSE)))))
  lamda <- exp(seq(log(1e-05), log(lmax), length.out = 100)) ## define lambda sequence
  b.lasso <- glmnet(Z[-miss,],Y[i,],lambda = lamda, alpha = 1)
  ind <- which.min(abs(b.lasso$lambda-lambdB))
  bhat[,i] <- b.lasso$beta[,ind]
}
bhat_train[2,,] <- bhat                                         

## reconstruct training images
Yhat <- Z%*%bhat
png(filename = paste("Pictures/Yr20-recon",art,q,"it",it,"_train.png",sep = ""), width = 770, height = 380)       #!!!#
par(mfrow=c(2,5), mar=rep(0,4))
for(i in 1:10){image(matrix(Yhat[,i],nrow=92), col=gray.colors(100),axes=FALSE)}
dev.off()

######################
miss <- rand50; Y <- Ytrainr50
## compute bhat
bhat<- matrix(0,nrow = nrow(dict),ncol = nrow(Xtr))
pb <- progress_bar$new(format = "progress: [:bar] :percent",total = nrow(Xtr), 
                       clear = FALSE, width= 60)
for(i in 1:nrow(Xtr)){
  pb$tick()
  lmax <- max(abs(drop(crossprod(Y[i,]-mean(Y[i,]),
                                 scale(Z[-miss,],colMeans(Z[-miss,]),FALSE)))))
  lamda <- exp(seq(log(1e-05), log(lmax), length.out = 100)) ## define lambda sequence
  b.lasso <- glmnet(Z[-miss,],Y[i,],lambda = lamda, alpha = 1)
  ind <- which.min(abs(b.lasso$lambda-lambdB))
  bhat[,i] <- b.lasso$beta[,ind]
}
bhat_train[3,,] <- bhat                                         

## reconstruct training images
Yhat <- Z%*%bhat
png(filename = paste("Pictures/Yr50-recon",art,q,"it",it,"_train.png",sep = ""), width = 770, height = 380)       #!!!#
par(mfrow=c(2,5), mar=rep(0,4))
for(i in 1:10){image(matrix(Yhat[,i],nrow=92), col=gray.colors(100),axes=FALSE)}
dev.off()

######################
miss <- rand75; Y <- Ytrainr75
## compute bhat
bhat<- matrix(0,nrow = nrow(dict),ncol = nrow(Xtr))
pb <- progress_bar$new(format = "progress: [:bar] :percent",total = nrow(Xtr), 
                       clear = FALSE, width= 60)
for(i in 1:nrow(Xtr)){
  pb$tick()
  lmax <- max(abs(drop(crossprod(Y[i,]-mean(Y[i,]),
                                 scale(Z[-miss,],colMeans(Z[-miss,]),FALSE)))))
  lamda <- exp(seq(log(1e-05), log(lmax), length.out = 100)) ## define lambda sequence
  b.lasso <- glmnet(Z[-miss,],Y[i,],lambda = lamda, alpha = 1)
  ind <- which.min(abs(b.lasso$lambda-lambdB))
  bhat[,i] <- b.lasso$beta[,ind]
}
bhat_train[4,,] <- bhat                                         

## reconstruct training images
Yhat <- Z%*%bhat
png(filename = paste("Pictures/Yr75-recon",art,q,"it",it,"_train.png",sep = ""), width = 770, height = 380)       #!!!#
par(mfrow=c(2,5), mar=rep(0,4))
for(i in 1:10){image(matrix(Yhat[,i],nrow=92), col=gray.colors(100),axes=FALSE)}
dev.off()


######################
miss <- patch10; Y <- Ytrainp10
## compute bhat
bhat<- matrix(0,nrow = nrow(dict),ncol = nrow(Xtr))
pb <- progress_bar$new(format = "progress: [:bar] :percent",total = nrow(Xtr), 
                       clear = FALSE, width= 60)
for(i in 1:nrow(Xtr)){
  pb$tick()
  lmax <- max(abs(drop(crossprod(Y[i,]-mean(Y[i,]),
                                 scale(Z[-miss,],colMeans(Z[-miss,]),FALSE)))))
  lamda <- exp(seq(log(1e-05), log(lmax), length.out = 100)) ## define lambda sequence
  b.lasso <- glmnet(Z[-miss,],Y[i,],lambda = lamda, alpha = 1)
  ind <- which.min(abs(b.lasso$lambda-lambdB))
  bhat[,i] <- b.lasso$beta[,ind]
}
bhat_train[5,,] <- bhat                                         

## reconstruct training images
Yhat <- Z%*%bhat
png(filename = paste("Pictures/Yp10-recon",art,q,"it",it,"_train.png",sep = ""), width = 770, height = 380)       #!!!#
par(mfrow=c(2,5), mar=rep(0,4))
for(i in 1:10){image(matrix(Yhat[,i],nrow=92), col=gray.colors(100),axes=FALSE)}
dev.off()

######################
miss <- patch20; Y <- Ytrainp20
## compute bhat
bhat<- matrix(0,nrow = nrow(dict),ncol = nrow(Xtr))
pb <- progress_bar$new(format = "progress: [:bar] :percent",total = nrow(Xtr), 
                       clear = FALSE, width= 60)
for(i in 1:nrow(Xtr)){
  pb$tick()
  lmax <- max(abs(drop(crossprod(Y[i,]-mean(Y[i,]),
                                 scale(Z[-miss,],colMeans(Z[-miss,]),FALSE)))))
  lamda <- exp(seq(log(1e-05), log(lmax), length.out = 100)) ## define lambda sequence
  b.lasso <- glmnet(Z[-miss,],Y[i,],lambda = lamda, alpha = 1)
  ind <- which.min(abs(b.lasso$lambda-lambdB))
  bhat[,i] <- b.lasso$beta[,ind]
}
bhat_train[6,,] <- bhat                                         

## reconstruct training images
Yhat <- Z%*%bhat
png(filename = paste("Pictures/Yp20-recon",art,q,"it",it,"_train.png",sep = ""), width = 770, height = 380)       #!!!#
par(mfrow=c(2,5), mar=rep(0,4))
for(i in 1:10){image(matrix(Yhat[,i],nrow=92), col=gray.colors(100),axes=FALSE)}
dev.off()

######################
miss <- patch50; Y <- Ytrainp50
## compute bhat
bhat<- matrix(0,nrow = nrow(dict),ncol = nrow(Xtr))
pb <- progress_bar$new(format = "progress: [:bar] :percent",total = nrow(Xtr), 
                       clear = FALSE, width= 60)
for(i in 1:nrow(Xtr)){
  pb$tick()
  lmax <- max(abs(drop(crossprod(Y[i,]-mean(Y[i,]),
                                 scale(Z[-miss,],colMeans(Z[-miss,]),FALSE)))))
  lamda <- exp(seq(log(1e-05), log(lmax), length.out = 100)) ## define lambda sequence
  b.lasso <- glmnet(Z[-miss,],Y[i,],lambda = lamda, alpha = 1)
  ind <- which.min(abs(b.lasso$lambda-lambdB))
  bhat[,i] <- b.lasso$beta[,ind]
}
bhat_train[7,,] <- bhat                                         

## reconstruct training images
Yhat <- Z%*%bhat
png(filename = paste("Pictures/Yp50-recon",art,q,"it",it,"_train.png",sep = ""), width = 770, height = 380)       #!!!#
par(mfrow=c(2,5), mar=rep(0,4))
for(i in 1:10){image(matrix(Yhat[,i],nrow=92), col=gray.colors(100),axes=FALSE)}
dev.off()

######################
miss <- patch75; Y <- Ytrainp75
## compute bhat
bhat<- matrix(0,nrow = nrow(dict),ncol = nrow(Xtr))
pb <- progress_bar$new(format = "progress: [:bar] :percent",total = nrow(Xtr), 
                       clear = FALSE, width= 60)
for(i in 1:nrow(Xtr)){
  pb$tick()
  lmax <- max(abs(drop(crossprod(Y[i,]-mean(Y[i,]),
                                 scale(Z[-miss,],colMeans(Z[-miss,]),FALSE)))))
  lamda <- exp(seq(log(1e-05), log(lmax), length.out = 100)) ## define lambda sequence
  b.lasso <- glmnet(Z[-miss,],Y[i,],lambda = lamda, alpha = 1)
  ind <- which.min(abs(b.lasso$lambda-lambdB))
  bhat[,i] <- b.lasso$beta[,ind]
}
bhat_train[8,,] <- bhat                                         

## reconstruct training images
Yhat <- Z%*%bhat
png(filename = paste("Pictures/Yp75-recon",art,q,"it",it,"_train.png",sep = ""), width = 770, height = 380)       #!!!#
par(mfrow=c(2,5), mar=rep(0,4))
for(i in 1:10){image(matrix(Yhat[,i],nrow=92), col=gray.colors(100),axes=FALSE)}
dev.off()


######################
miss <- Mpatch10; Y <- YtrainMp10
## compute bhat
bhat<- matrix(0,nrow = nrow(dict),ncol = nrow(Xtr))
pb <- progress_bar$new(format = "progress: [:bar] :percent",total = nrow(Xtr), 
                       clear = FALSE, width= 60)
for(i in 1:nrow(Xtr)){
  pb$tick()
  lmax <- max(abs(drop(crossprod(Y[i,]-mean(Y[i,]),
                                 scale(Z[-miss,],colMeans(Z[-miss,]),FALSE)))))
  lamda <- exp(seq(log(1e-05), log(lmax), length.out = 100)) ## define lambda sequence
  b.lasso <- glmnet(Z[-miss,],Y[i,],lambda = lamda, alpha = 1)
  ind <- which.min(abs(b.lasso$lambda-lambdB))
  bhat[,i] <- b.lasso$beta[,ind]
}
bhat_train[9,,] <- bhat                                         

## reconstruct training images
Yhat <- Z%*%bhat
png(filename = paste("Pictures/YMp10-recon",art,q,"it",it,"_train.png",sep = ""), width = 770, height = 380)       #!!!#
par(mfrow=c(2,5), mar=rep(0,4))
for(i in 1:10){image(matrix(Yhat[,i],nrow=92), col=gray.colors(100),axes=FALSE)}
dev.off()

######################
miss <- Mpatch20; Y <- YtrainMp20
## compute bhat
bhat<- matrix(0,nrow = nrow(dict),ncol = nrow(Xtr))
pb <- progress_bar$new(format = "progress: [:bar] :percent",total = nrow(Xtr), 
                       clear = FALSE, width= 60)
for(i in 1:nrow(Xtr)){
  pb$tick()
  lmax <- max(abs(drop(crossprod(Y[i,]-mean(Y[i,]),
                                 scale(Z[-miss,],colMeans(Z[-miss,]),FALSE)))))
  lamda <- exp(seq(log(1e-05), log(lmax), length.out = 100)) ## define lambda sequence
  b.lasso <- glmnet(Z[-miss,],Y[i,],lambda = lamda, alpha = 1)
  ind <- which.min(abs(b.lasso$lambda-lambdB))
  bhat[,i] <- b.lasso$beta[,ind]
}
bhat_train[10,,] <- bhat                                         

## reconstruct training images
Yhat <- Z%*%bhat
png(filename = paste("Pictures/YMp20-recon",art,q,"it",it,"_train.png",sep = ""), width = 770, height = 380)       #!!!#
par(mfrow=c(2,5), mar=rep(0,4))
for(i in 1:10){image(matrix(Yhat[,i],nrow=92), col=gray.colors(100),axes=FALSE)}
dev.off()

######################
miss <- Mpatch50; Y <- YtrainMp50
## compute bhat
bhat<- matrix(0,nrow = nrow(dict),ncol = nrow(Xtr))
pb <- progress_bar$new(format = "progress: [:bar] :percent",total = nrow(Xtr), 
                       clear = FALSE, width= 60)
for(i in 1:nrow(Xtr)){
  pb$tick()
  lmax <- max(abs(drop(crossprod(Y[i,]-mean(Y[i,]),
                                 scale(Z[-miss,],colMeans(Z[-miss,]),FALSE)))))
  lamda <- exp(seq(log(1e-05), log(lmax), length.out = 100)) ## define lambda sequence
  b.lasso <- glmnet(Z[-miss,],Y[i,],lambda = lamda, alpha = 1)
  ind <- which.min(abs(b.lasso$lambda-lambdB))
  bhat[,i] <- b.lasso$beta[,ind]
}
bhat_train[11,,] <- bhat                                         

## reconstruct training images
Yhat <- Z%*%bhat
png(filename = paste("Pictures/YMp50-recon",art,q,"it",it,"_train.png",sep = ""), width = 770, height = 380)       #!!!#
par(mfrow=c(2,5), mar=rep(0,4))
for(i in 1:10){image(matrix(Yhat[,i],nrow=92), col=gray.colors(100),axes=FALSE)}
dev.off()

######################
miss <- Mpatch75; Y <- YtrainMp75
## compute bhat
bhat<- matrix(0,nrow = nrow(dict),ncol = nrow(Xtr))
pb <- progress_bar$new(format = "progress: [:bar] :percent",total = nrow(Xtr), 
                       clear = FALSE, width= 60)
for(i in 1:nrow(Xtr)){
  pb$tick()
  lmax <- max(abs(drop(crossprod(Y[i,]-mean(Y[i,]),
                                 scale(Z[-miss,],colMeans(Z[-miss,]),FALSE)))))
  lamda <- exp(seq(log(1e-05), log(lmax), length.out = 100)) ## define lambda sequence
  b.lasso <- glmnet(Z[-miss,],Y[i,],lambda = lamda, alpha = 1)
  ind <- which.min(abs(b.lasso$lambda-lambdB))
  bhat[,i] <- b.lasso$beta[,ind]
}
bhat_train[12,,] <- bhat                                         

## reconstruct training images
Yhat <- Z%*%bhat
png(filename = paste("Pictures/YMp75-recon",art,q,"it",it,"_train.png",sep = ""), width = 770, height = 380)       #!!!#
par(mfrow=c(2,5), mar=rep(0,4))
for(i in 1:10){image(matrix(Yhat[,i],nrow=92), col=gray.colors(100),axes=FALSE)}
dev.off()

########################
## save  bhat
save(bhat_train, file=paste("bhat_train",art,q,"it",it,".rda",sep = ""))