library(glmnet)
library(progress)

## load data
load(file="faces400.rda")
load(file ="faces100.rda")
Xtest <- rbind(test[1:50,],X[351:400,])

## set parameters
art <- "NMF" #"AMr","AMp","OLr","OLp", "PCA"
q <- 100 #300, 600
it <- 50 # 150

## load dictionary
load(file=paste("dict_",art,q,"_it",it,".rda",sep = ""))
# load(file=paste("dict_",art,q,".rda",sep = "")) # for art="PCA"
dict <- spdict$H
Z <- t(dict)

## delete pixels  - create Y^(S) for S in M_o
set.seed(42)
rand10 <- sample(1:ncol(X),100,replace=FALSE); Yr10 <- Xtest[,-rand10]
rand20 <- sample(1:ncol(X),400,replace=FALSE); Yr20 <- Xtest[,-rand20]
rand50 <- sample(1:ncol(X),2500,replace=FALSE); Yr50 <- Xtest[,-rand50]
rand75 <- sample(1:ncol(X),5625,replace=FALSE); Yr75 <- Xtest[,-rand75]

m10 <- matrix(0, nrow=10, ncol=10)
for(i in 55:64){m10[i-54,] <- (i*92+25):(i*92+34)}
patch10 <- as.vector(m10); Yp10 <- Xtest[,-patch10]
m20 <- matrix(0, nrow=20, ncol=20)
for(i in 55:74){m20[i-54,] <- (i*92+18):(i*92+37)}
patch20 <- as.vector(m20); Yp20 <- Xtest[,-patch20]
m50 <- matrix(0, nrow=50, ncol=50)
for(i in 45:94){m50[i-44,] <- (i*92+5):(i*92+54)}
patch50 <- as.vector(m50); Yp50 <- Xtest[,-patch50]
m75 <- matrix(0, nrow=75, ncol=75)
for(i in 25:99){m75[i-24,] <- (i*92+3):(i*92+77)}
patch75 <- as.vector(m75); Yp75 <- Xtest[,-patch75]

Mm10 <- matrix(0, nrow=10, ncol=10)
for(i in 20:29){Mm10[i-19,] <- (i*92+45):(i*92+54)}
Mpatch10 <- as.vector(Mm10); YMp10 <- Xtest[,-Mpatch10]
Mm20 <- matrix(0, nrow=20, ncol=20)
for(i in 15:34){Mm20[i-14,] <- (i*92+40):(i*92+59)}
Mpatch20 <- as.vector(Mm20); YMp20 <- Xtest[,-Mpatch20]
Mm50 <- matrix(0, nrow=50, ncol=50)
for(i in 2:51){Mm50[i-1,] <- (i*92+25):(i*92+74)}
Mpatch50 <- as.vector(Mm50); YMp50 <- Xtest[,-Mpatch50]
Mm75 <- matrix(0, nrow=75, ncol=75)
for(i in 2:76){Mm75[i-1,] <- (i*92+10):(i*92+84)}
Mpatch75 <- as.vector(Mm75); YMp75 <- Xtest[,-Mpatch75]


#######################################################################################################
## Learn bhat for 10 randomly selected Training Samples
###############

lambdB <- 0.001

## initialize array to save matrices bhat for each S in M_o 
bhat_test <- array(0, dim=c(12,nrow(dict),nrow(Xtest)))

######################
miss <- rand10; Y <- Yr10
## compute bhat
bhat<- matrix(0,nrow = nrow(dict),ncol = nrow(Xtest))
pb <- progress_bar$new(format = "progress: [:bar] :percent",total = nrow(Xtest), 
                       clear = FALSE, width= 60)
for(i in 1:nrow(Xtest)){
  pb$tick()
  lmax <- max(abs(drop(crossprod(Y[i,]-mean(Y[i,]),
                                 scale(Z[-miss,],colMeans(Z[-miss,]),FALSE)))))
  lamda <- exp(seq(log(1e-05), log(lmax), length.out = 100)) ## define lambda sequence
  b.lasso <- glmnet(Z[-miss,],Y[i,],lambda = lamda, alpha = 1)
  ind <- which.min(abs(b.lasso$lambda-lambdB))
  bhat[,i] <- b.lasso$beta[,ind]
}

## save bhat
bhat_test[1,,] <- bhat

## reconstruct test images
Yhat <- Z%*%bhat
png(filename = paste("Pictures/Yr10-recon",art,q,"it",it,"_test.png",sep = ""), width = 770, height = 380)       #!!!#
par(mfrow=c(2,5), mar=rep(0,4))
for(i in seq(1,nrow(Xtest),10)){image(matrix(Yhat[,i],nrow=92), col=gray.colors(100),axes=FALSE)}
dev.off()

######################
miss <- rand20; Y <- Yr20
## compute bhat
bhat<- matrix(0,nrow = nrow(dict),ncol = nrow(Xtest))
pb <- progress_bar$new(format = "progress: [:bar] :percent",total = nrow(Xtest), 
                       clear = FALSE, width= 60)
for(i in 1:nrow(Xtest)){
  pb$tick()
  lmax <- max(abs(drop(crossprod(Y[i,]-mean(Y[i,]),
                                 scale(Z[-miss,],colMeans(Z[-miss,]),FALSE)))))
  lamda <- exp(seq(log(1e-05), log(lmax), length.out = 100)) ## define lambda sequence
  b.lasso <- glmnet(Z[-miss,],Y[i,],lambda = lamda, alpha = 1)
  ind <- which.min(abs(b.lasso$lambda-lambdB))
  bhat[,i] <- b.lasso$beta[,ind]
}

## save bhat
bhat_test[2,,] <- bhat

## reconstruct test images
Yhat <- Z%*%bhat
png(filename = paste("Pictures/Yr20-recon",art,q,"it",it,"_test.png",sep = ""), width = 770, height = 380)       #!!!#
par(mfrow=c(2,5), mar=rep(0,4))
for(i in seq(1,nrow(Xtest),10)){image(matrix(Yhat[,i],nrow=92), col=gray.colors(100),axes=FALSE)}
dev.off()

######################
miss <- rand50; Y <- Yr50
## compute bhat
bhat<- matrix(0,nrow = nrow(dict),ncol = nrow(Xtest))
pb <- progress_bar$new(format = "progress: [:bar] :percent",total = nrow(Xtest), 
                       clear = FALSE, width= 60)
for(i in 1:nrow(Xtest)){
  pb$tick()
  lmax <- max(abs(drop(crossprod(Y[i,]-mean(Y[i,]),
                                 scale(Z[-miss,],colMeans(Z[-miss,]),FALSE)))))
  lamda <- exp(seq(log(1e-05), log(lmax), length.out = 100)) ## define lambda sequence
  b.lasso <- glmnet(Z[-miss,],Y[i,],lambda = lamda, alpha = 1)
  ind <- which.min(abs(b.lasso$lambda-lambdB))
  bhat[,i] <- b.lasso$beta[,ind]
}

## save bhat
bhat_test[3,,] <- bhat

## reconstruct test images
Yhat <- Z%*%bhat
png(filename = paste("Pictures/Yr50-recon",art,q,"it",it,"_test.png",sep = ""), width = 770, height = 380)       #!!!#
par(mfrow=c(2,5), mar=rep(0,4))
for(i in seq(1,nrow(Xtest),10)){image(matrix(Yhat[,i],nrow=92), col=gray.colors(100),axes=FALSE)}
dev.off()

######################
miss <- rand75; Y <- Yr75
## compute bhat
bhat<- matrix(0,nrow = nrow(dict),ncol = nrow(Xtest))
pb <- progress_bar$new(format = "progress: [:bar] :percent",total = nrow(Xtest), 
                       clear = FALSE, width= 60)
for(i in 1:nrow(Xtest)){
  pb$tick()
  lmax <- max(abs(drop(crossprod(Y[i,]-mean(Y[i,]),
                                 scale(Z[-miss,],colMeans(Z[-miss,]),FALSE)))))
  lamda <- exp(seq(log(1e-05), log(lmax), length.out = 100)) ## define lambda sequence
  b.lasso <- glmnet(Z[-miss,],Y[i,],lambda = lamda, alpha = 1)
  ind <- which.min(abs(b.lasso$lambda-lambdB))
  bhat[,i] <- b.lasso$beta[,ind]
}

## save bhat
bhat_test[4,,] <- bhat

## reconstruct test images
Yhat <- Z%*%bhat
png(filename = paste("Pictures/Yr75-recon",art,q,"it",it,"_test.png",sep = ""), width = 770, height = 380)       #!!!#
par(mfrow=c(2,5), mar=rep(0,4))
for(i in seq(1,nrow(Xtest),10)){image(matrix(Yhat[,i],nrow=92), col=gray.colors(100),axes=FALSE)}
dev.off()

######################
miss <- patch10; Y <- Yp10
## compute bhat
bhat<- matrix(0,nrow = nrow(dict),ncol = nrow(Xtest))
pb <- progress_bar$new(format = "progress: [:bar] :percent",total = nrow(Xtest), 
                       clear = FALSE, width= 60)
for(i in 1:nrow(Xtest)){
  pb$tick()
  lmax <- max(abs(drop(crossprod(Y[i,]-mean(Y[i,]),
                                 scale(Z[-miss,],colMeans(Z[-miss,]),FALSE)))))
  lamda <- exp(seq(log(1e-05), log(lmax), length.out = 100)) ## define lambda sequence
  b.lasso <- glmnet(Z[-miss,],Y[i,],lambda = lamda, alpha = 1)
  ind <- which.min(abs(b.lasso$lambda-lambdB))
  bhat[,i] <- b.lasso$beta[,ind]
}

## save bhat
bhat_test[5,,] <- bhat

## reconstruct test images
Yhat <- Z%*%bhat
png(filename = paste("Pictures/Yp10-recon",art,q,"it",it,"_test.png",sep = ""), width = 770, height = 380)       #!!!#
par(mfrow=c(2,5), mar=rep(0,4))
for(i in seq(1,nrow(Xtest),10)){image(matrix(Yhat[,i],nrow=92), col=gray.colors(100),axes=FALSE)}
dev.off()

######################
miss <- patch20; Y <- Yp20
## compute bhat
bhat<- matrix(0,nrow = nrow(dict),ncol = nrow(Xtest))
pb <- progress_bar$new(format = "progress: [:bar] :percent",total = nrow(Xtest), 
                       clear = FALSE, width= 60)
for(i in 1:nrow(Xtest)){
  pb$tick()
  lmax <- max(abs(drop(crossprod(Y[i,]-mean(Y[i,]),
                                 scale(Z[-miss,],colMeans(Z[-miss,]),FALSE)))))
  lamda <- exp(seq(log(1e-05), log(lmax), length.out = 100)) ## define lambda sequence
  b.lasso <- glmnet(Z[-miss,],Y[i,],lambda = lamda, alpha = 1)
  ind <- which.min(abs(b.lasso$lambda-lambdB))
  bhat[,i] <- b.lasso$beta[,ind]
}

## save bhat
bhat_test[6,,] <- bhat

## reconstruct test images
Yhat <- Z%*%bhat
png(filename = paste("Pictures/Yp20-recon",art,q,"it",it,"_test.png",sep = ""), width = 770, height = 380)       #!!!#
par(mfrow=c(2,5), mar=rep(0,4))
for(i in seq(1,nrow(Xtest),10)){image(matrix(Yhat[,i],nrow=92), col=gray.colors(100),axes=FALSE)}
dev.off()

######################
miss <- patch50; Y <- Yp50
## compute bhat
bhat<- matrix(0,nrow = nrow(dict),ncol = nrow(Xtest))
pb <- progress_bar$new(format = "progress: [:bar] :percent",total = nrow(Xtest), 
                       clear = FALSE, width= 60)
for(i in 1:nrow(Xtest)){
  pb$tick()
  lmax <- max(abs(drop(crossprod(Y[i,]-mean(Y[i,]),
                                 scale(Z[-miss,],colMeans(Z[-miss,]),FALSE)))))
  lamda <- exp(seq(log(1e-05), log(lmax), length.out = 100)) ## define lambda sequence
  b.lasso <- glmnet(Z[-miss,],Y[i,],lambda = lamda, alpha = 1)
  ind <- which.min(abs(b.lasso$lambda-lambdB))
  bhat[,i] <- b.lasso$beta[,ind]
}

## save bhat
bhat_test[7,,] <- bhat

## reconstruct test images
Yhat <- Z%*%bhat
png(filename = paste("Pictures/Yp50-recon",art,q,"it",it,"_test.png",sep = ""), width = 770, height = 380)       #!!!#
par(mfrow=c(2,5), mar=rep(0,4))
for(i in seq(1,nrow(Xtest),10)){image(matrix(Yhat[,i],nrow=92), col=gray.colors(100),axes=FALSE)}
dev.off()

######################
miss <- patch75; Y <- Yp75
## compute bhat
bhat<- matrix(0,nrow = nrow(dict),ncol = nrow(Xtest))
pb <- progress_bar$new(format = "progress: [:bar] :percent",total = nrow(Xtest), 
                       clear = FALSE, width= 60)
for(i in 1:nrow(Xtest)){
  pb$tick()
  lmax <- max(abs(drop(crossprod(Y[i,]-mean(Y[i,]),
                                 scale(Z[-miss,],colMeans(Z[-miss,]),FALSE)))))
  lamda <- exp(seq(log(1e-05), log(lmax), length.out = 100)) ## define lambda sequence
  b.lasso <- glmnet(Z[-miss,],Y[i,],lambda = lamda, alpha = 1)
  ind <- which.min(abs(b.lasso$lambda-lambdB))
  bhat[,i] <- b.lasso$beta[,ind]
}

## save bhat
bhat_test[8,,] <- bhat

## reconstruct test images
Yhat <- Z%*%bhat
png(filename = paste("Pictures/Yp75-recon",art,q,"it",it,"_test.png",sep = ""), width = 770, height = 380)       #!!!#
par(mfrow=c(2,5), mar=rep(0,4))
for(i in seq(1,nrow(Xtest),10)){image(matrix(Yhat[,i],nrow=92), col=gray.colors(100),axes=FALSE)}
dev.off()

######################
miss <- Mpatch10; Y <- YMp10
## compute bhat
bhat<- matrix(0,nrow = nrow(dict),ncol = nrow(Xtest))
pb <- progress_bar$new(format = "progress: [:bar] :percent",total = nrow(Xtest), 
                       clear = FALSE, width= 60)
for(i in 1:nrow(Xtest)){
  pb$tick()
  lmax <- max(abs(drop(crossprod(Y[i,]-mean(Y[i,]),
                                 scale(Z[-miss,],colMeans(Z[-miss,]),FALSE)))))
  lamda <- exp(seq(log(1e-05), log(lmax), length.out = 100)) ## define lambda sequence
  b.lasso <- glmnet(Z[-miss,],Y[i,],lambda = lamda, alpha = 1)
  ind <- which.min(abs(b.lasso$lambda-lambdB))
  bhat[,i] <- b.lasso$beta[,ind]
}

## save bhat
bhat_test[9,,] <- bhat

## reconstruct test images
Yhat <- Z%*%bhat
png(filename = paste("Pictures/YMp10-recon",art,q,"it",it,"_test.png",sep = ""), width = 770, height = 380)       #!!!#
par(mfrow=c(2,5), mar=rep(0,4))
for(i in seq(1,nrow(Xtest),10)){image(matrix(Yhat[,i],nrow=92), col=gray.colors(100),axes=FALSE)}
dev.off()

######################
miss <- Mpatch20; Y <- YMp20
## compute bhat
bhat<- matrix(0,nrow = nrow(dict),ncol = nrow(Xtest))
pb <- progress_bar$new(format = "progress: [:bar] :percent",total = nrow(Xtest), 
                       clear = FALSE, width= 60)
for(i in 1:nrow(Xtest)){
  pb$tick()
  lmax <- max(abs(drop(crossprod(Y[i,]-mean(Y[i,]),
                                 scale(Z[-miss,],colMeans(Z[-miss,]),FALSE)))))
  lamda <- exp(seq(log(1e-05), log(lmax), length.out = 100)) ## define lambda sequence
  b.lasso <- glmnet(Z[-miss,],Y[i,],lambda = lamda, alpha = 1)
  ind <- which.min(abs(b.lasso$lambda-lambdB))
  bhat[,i] <- b.lasso$beta[,ind]
}

## save bhat
bhat_test[10,,] <- bhat

## reconstruct test images
Yhat <- Z%*%bhat
png(filename = paste("Pictures/YMp20-recon",art,q,"it",it,"_test.png",sep = ""), width = 770, height = 380)       #!!!#
par(mfrow=c(2,5), mar=rep(0,4))
for(i in seq(1,nrow(Xtest),10)){image(matrix(Yhat[,i],nrow=92), col=gray.colors(100),axes=FALSE)}
dev.off()

######################
miss <- Mpatch50; Y <- YMp50
## compute bhat
bhat<- matrix(0,nrow = nrow(dict),ncol = nrow(Xtest))
pb <- progress_bar$new(format = "progress: [:bar] :percent",total = nrow(Xtest), 
                       clear = FALSE, width= 60)
for(i in 1:nrow(Xtest)){
  pb$tick()
  lmax <- max(abs(drop(crossprod(Y[i,]-mean(Y[i,]),
                                 scale(Z[-miss,],colMeans(Z[-miss,]),FALSE)))))
  lamda <- exp(seq(log(1e-05), log(lmax), length.out = 100)) ## define lambda sequence
  b.lasso <- glmnet(Z[-miss,],Y[i,],lambda = lamda, alpha = 1)
  ind <- which.min(abs(b.lasso$lambda-lambdB))
  bhat[,i] <- b.lasso$beta[,ind]
}

## save bhat
bhat_test[11,,] <- bhat

## reconstruct test images
Yhat <- Z%*%bhat
png(filename = paste("Pictures/YMp50-recon",art,q,"it",it,"_test.png",sep = ""), width = 770, height = 380)       #!!!#
par(mfrow=c(2,5), mar=rep(0,4))
for(i in seq(1,nrow(Xtest),10)){image(matrix(Yhat[,i],nrow=92), col=gray.colors(100),axes=FALSE)}
dev.off()

######################
miss <- Mpatch75; Y <- YMp75
## compute bhat
bhat<- matrix(0,nrow = nrow(dict),ncol = nrow(Xtest))
pb <- progress_bar$new(format = "progress: [:bar] :percent",total = nrow(Xtest), 
                       clear = FALSE, width= 60)
for(i in 1:nrow(Xtest)){
  pb$tick()
  lmax <- max(abs(drop(crossprod(Y[i,]-mean(Y[i,]),
                                 scale(Z[-miss,],colMeans(Z[-miss,]),FALSE)))))
  lamda <- exp(seq(log(1e-05), log(lmax), length.out = 100)) ## define lambda sequence
  b.lasso <- glmnet(Z[-miss,],Y[i,],lambda = lamda, alpha = 1)
  ind <- which.min(abs(b.lasso$lambda-lambdB))
  bhat[,i] <- b.lasso$beta[,ind]
}

## save bhat
bhat_test[12,,] <- bhat

## reconstruct test images
Yhat <- Z%*%bhat
png(filename = paste("Pictures/YMp75-recon",art,q,"it",it,"_test.png",sep = ""), width = 770, height = 380)       #!!!#
par(mfrow=c(2,5), mar=rep(0,4))
for(i in seq(1,nrow(Xtest),10)){image(matrix(Yhat[,i],nrow=92), col=gray.colors(100),axes=FALSE)}
dev.off()


########################
## save  bhat
save(bhat_test, file=paste("bhat_test",art,q,"it",it,".rda",sep = ""))