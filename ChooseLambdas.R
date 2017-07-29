############################################################################################
### Error plots
############################################################################################
library(progress)
library(glmnet)

## load data
load(file="faces400.rda")
load(file="faces100.rda")
Xtest <- rbind(test[1:50,],X[351:400,])

## load dictionaries for same q, different lambdas
art <- "NMF"  # "AMp","AMr","OLp"
q <- 100 # 300
it <- 50 # 150
## comment out all variables (whole scirpt) ending with numbers 5,6,7 when using art="AMp","AMr","OL"
load(file=paste("dict_",art,q,"_lA0.001_it",it,".rda",sep = "")); dict1 <- spdict$H
load(file=paste("dict_",art,q,"_lA0.01_it",it,".rda",sep = "")); dict2 <- spdict$H
load(file=paste("dict_",art,q,"_lA0.05_it",it,".rda",sep = "")); dict3 <- spdict$H
load(file=paste("dict_",art,q,"_lA0.1_it",it,".rda",sep = "")); dict4 <- spdict$H
load(file=paste("dict_",art,q,"_lA1_it",it,".rda",sep = "")); dict5 <- spdict$H
load(file=paste("dict_",art,q,"_lA5_it",it,".rda",sep = "")); dict6 <- spdict$H
load(file=paste("dict_",art,q,"_lA10_it",it,".rda",sep = "")); dict7 <- spdict$H

## delete pixels - create Y^(S) for S in M_o
set.seed(42)
rand10 <- sample(1:ncol(X),100,replace=FALSE); Yr10 <- Xtest2[,-rand10]
rand20 <- sample(1:ncol(X),400,replace=FALSE); Yr20 <- Xtest2[,-rand20]
rand50 <- sample(1:ncol(X),2500,replace=FALSE); Yr50 <- Xtest2[,-rand50]
rand75 <- sample(1:ncol(X),5625,replace=FALSE); Yr75 <- Xtest2[,-rand75]

m10 <- matrix(0, nrow=10, ncol=10)
for(i in 55:64){m10[i-54,] <- (i*92+25):(i*92+34)}
patch10 <- as.vector(m10); Yp10 <- Xtest2[,-patch10]
m20 <- matrix(0, nrow=20, ncol=20)
for(i in 55:74){m20[i-54,] <- (i*92+18):(i*92+37)}
patch20 <- as.vector(m20); Yp20 <- Xtest2[,-patch20]
m50 <- matrix(0, nrow=50, ncol=50)
for(i in 45:94){m50[i-44,] <- (i*92+5):(i*92+54)}
patch50 <- as.vector(m50); Yp50 <- Xtest2[,-patch50]
m75 <- matrix(0, nrow=75, ncol=75)
for(i in 25:99){m75[i-24,] <- (i*92+3):(i*92+77)}
patch75 <- as.vector(m75); Yp75 <- Xtest2[,-patch75]

Mm10 <- matrix(0, nrow=10, ncol=10)
for(i in 20:29){Mm10[i-19,] <- (i*92+45):(i*92+54)}
Mpatch10 <- as.vector(Mm10); YMp10 <- Xtest2[,-Mpatch10]
Mm20 <- matrix(0, nrow=20, ncol=20)
for(i in 15:34){Mm20[i-14,] <- (i*92+40):(i*92+59)}
Mpatch20 <- as.vector(Mm20); YMp20 <- Xtest2[,-Mpatch20]
Mm50 <- matrix(0, nrow=50, ncol=50)
for(i in 2:51){Mm50[i-1,] <- (i*92+25):(i*92+74)}
Mpatch50 <- as.vector(Mm50); YMp50 <- Xtest2[,-Mpatch50]
Mm75 <- matrix(0, nrow=75, ncol=75)
for(i in 2:76){Mm75[i-1,] <- (i*92+10):(i*92+84)}
Mpatch75 <- as.vector(Mm75); YMp75 <- Xtest2[,-Mpatch75]

## choose set S^c and Y^(S)  (change name of plot at the end accordingly)
miss <- rand10;Y <- Yr10
# miss <- rand20;Y <- Yr20
# miss <- rand50;Y <- Yr50
# miss <- rand75;Y <- Yr75
# miss <- patch10;Y <- Yp10
# miss <- patch20;Y <- Yp20
# miss <- patch50;Y <- Yp50
# miss <- patch75;Y <- Yp75
# miss <- Mpatch10;Y <- YMp10
# miss <- Mpatch20;Y <- YMp20
# miss <- Mpatch50;Y <- YMp50
# miss <- Mpatch75;Y <- YMp75

## compute lambda_max^star
lam.max1 <- lam.max2 <- lam.max3 <- lam.max4 <- lam.max5 <- lam.max6 <- lam.max7 <- numeric(nrow(Xtest))
for(i in 1:nrow(Xtest)){
  lam.max1[i] <- max(abs(crossprod(Y[i,]-mean(Y[i,]),
                                  scale(t(dict1)[-miss,],colMeans(t(dict1)[-miss,])))))/length(Y[i,])
  lam.max2[i] <- max(abs(crossprod(Y[i,]-mean(Y[i,]),
                                  scale(t(dict2)[-miss,],colMeans(t(dict2)[-miss,])))))/length(Y[i,])
  lam.max3[i] <- max(abs(crossprod(Y[i,]-mean(Y[i,]),
                                  scale(t(dict3)[-miss,],colMeans(t(dict3)[-miss,])))))/length(Y[i,])
  lam.max4[i] <- max(abs(crossprod(Y[i,]-mean(Y[i,]),
                                   scale(t(dict4)[-miss,],colMeans(t(dict4)[-miss,])))))/length(Y[i,])
  lam.max5[i] <- max(abs(crossprod(Y[i,]-mean(Y[i,]),
                                   scale(t(dict5)[-miss,],colMeans(t(dict5)[-miss,])))))/length(Y[i,])
  lam.max6[i] <- max(abs(crossprod(Y[i,]-mean(Y[i,]),
                                   scale(t(dict6)[-miss,],colMeans(t(dict6)[-miss,])))))/length(Y[i,])
  lam.max7[i] <- max(abs(crossprod(Y[i,]-mean(Y[i,]),
                                   scale(t(dict7)[-miss,],colMeans(t(dict7)[-miss,])))))/length(Y[i,])
}

l.max.min <- min(lam.max1,lam.max2,lam.max3,lam.max4,lam.max5,lam.max6,lam.max7)


## define lambda sequence for which to compute reconstruction
lam <- exp(seq(log(1e-04),log(l.max.min),length.out = 10))


## check for which lambda in lam we want to learn b for incomplete test images
errmiss1 <- errmiss2 <- errmiss3 <- errmiss4 <- errmiss5 <- errmiss6 <- errmiss7 <- matrix(0,nrow=nrow(Xtest),ncol=length(lam))
lamA1 <- lamA2 <- lamA3 <- lamA4 <- lamA5 <- lamA6 <- lamA7 <- numeric(length(lam))
# lamA <- numeric(length(lam))
pb <- progress_bar$new(format = "progress: [:bar] :percent",total = nrow(Xtest), 
                       clear = FALSE, width= 60)
for(i in 1:nrow(Xtest)){
  pb$tick()
  ## for each test image & dictionary get smallest lambda value lmax for which all coefficients are zero
  lmax1 <- lam.max1[i];lmax2 <- lam.max2[i];lmax3 <- lam.max3[i];lmax4 <- lam.max4[i]
  lmax5 <- lam.max5[i];lmax6 <- lam.max6[i];lmax7 <- lam.max7[i]
  lamda1 <- exp(seq(log(1e-05), log(lmax1), length.out = 100))
  lamda2 <- exp(seq(log(1e-05), log(lmax2), length.out = 100))
  lamda3 <- exp(seq(log(1e-05), log(lmax3), length.out = 100))
  lamda4 <- exp(seq(log(1e-05), log(lmax4), length.out = 100))
  lamda5 <- exp(seq(log(1e-05), log(lmax5), length.out = 100))
  lamda6 <- exp(seq(log(1e-05), log(lmax6), length.out = 100))
  lamda7 <- exp(seq(log(1e-05), log(lmax7), length.out = 100))
  ## compute bhat
  b.lasso1 <- glmnet(t(dict1)[-miss,],Y[i,],lambda = lamda1, alpha = 1)
  b.lasso2 <- glmnet(t(dict2)[-miss,],Y[i,],lambda = lamda2, alpha = 1)
  b.lasso3 <- glmnet(t(dict3)[-miss,],Y[i,],lambda = lamda3, alpha = 1)
  b.lasso4 <- glmnet(t(dict4)[-miss,],Y[i,],lambda = lamda4, alpha = 1)
  b.lasso5 <- glmnet(t(dict5)[-miss,],Y[i,],lambda = lamda5, alpha = 1)
  b.lasso6 <- glmnet(t(dict6)[-miss,],Y[i,],lambda = lamda6, alpha = 1)
  b.lasso7 <- glmnet(t(dict7)[-miss,],Y[i,],lambda = lamda7, alpha = 1)
  ## Compute errors for the 10 values of lambda_b in sequence l
  for(k in 1:length(lam)){
    ind1 <- which.min(abs(b.lasso1$lambda-lam[k]));ind2 <- which.min(abs(b.lasso2$lambda-lam[k]))
    ind3 <- which.min(abs(b.lasso3$lambda-lam[k]));ind4 <- which.min(abs(b.lasso4$lambda-lam[k]))
    ind5 <- which.min(abs(b.lasso5$lambda-lam[k]));ind6 <- which.min(abs(b.lasso6$lambda-lam[k]))
    ind7 <- which.min(abs(b.lasso7$lambda-lam[k]))
    lamA1[k] <- b.lasso1$lambda[ind1];lamA2[k] <- b.lasso2$lambda[ind2]
    lamA3[k] <- b.lasso3$lambda[ind3];lamA4[k] <- b.lasso4$lambda[ind4]
    lamA5[k] <- b.lasso5$lambda[ind5];lamA6[k] <- b.lasso6$lambda[ind6]
    lamA7[k] <- b.lasso7$lambda[ind7]
    bh1 <- b.lasso1$beta[,ind1];bh2 <- b.lasso2$beta[,ind2]
    bh3 <- b.lasso3$beta[,ind3];bh4 <- b.lasso4$beta[,ind4]
    bh5 <- b.lasso5$beta[,ind5];bh6 <- b.lasso6$beta[,ind6]
    bh7 <- b.lasso7$beta[,ind7]
    Yh1 <- t(dict1)%*%bh1;Yh2 <- t(dict2)%*%bh2
    Yh3 <- t(dict3)%*%bh3;Yh4 <- t(dict4)%*%bh4
    Yh5 <- t(dict5)%*%bh5;Yh6 <- t(dict6)%*%bh6
    Yh7 <- t(dict7)%*%bh7
    errmiss1[i,k] <-  mean((Xtest[i,miss]-Yh1[miss])^2);errmiss2[i,k] <-  mean((Xtest[i,miss]-Yh2[miss])^2)
    errmiss3[i,k] <-  mean((Xtest[i,miss]-Yh3[miss])^2);errmiss4[i,k] <-  mean((Xtest[i,miss]-Yh4[miss])^2)
    errmiss5[i,k] <-  mean((Xtest[i,miss]-Yh5[miss])^2);errmiss6[i,k] <-  mean((Xtest[i,miss]-Yh6[miss])^2)
    errmiss7[i,k] <-  mean((Xtest[i,miss]-Yh7[miss])^2)
  }
}

## for every value in lam compute the mean error (on missing pixels) over all test images
a1 <- apply(errmiss1,2,mean);a2 <- apply(errmiss2,2,mean);a3 <- apply(errmiss3,2,mean);a4 <- apply(errmiss4,2,mean)
a5 <- apply(errmiss5,2,mean);a6 <- apply(errmiss6,2,mean);a7 <- apply(errmiss7,2,mean)

## plot errors
png(filename = paste("Pictures/Errors-",art,q,"it",it,"-Yr10.png"), width = 770, height = 380)
par(mfrow=c(1,1), mar=rep(5,4))
plot(lamA1,a1, xlab = expression(lambda[b]),ylab = expression(paste(l[2],'-error')), ylim = c(0,0.3),lwd=2)  ## plot mean errors against the lambda values 
lines(lamA1,a1)
points(lamA2,a2,col='red'); lines(lamA2,a2,col='red',lwd=2)
points(lamA3,a3,col='green'); lines(lamA3,a3,col='green',lwd=2)
points(lamA4,a4,col='blue'); lines(lamA4,a4,col='blue',lwd=2)
points(lamA5,a5,col='magenta'); lines(lamA5,a5,col='magenta',lwd=2)
points(lamA6,a6,col='cyan'); lines(lamA6,a6,col='cyan',lwd=2)
points(lamA7,a7,col='yellow'); lines(lamA7,a7,col='yellow',lwd=2)
legend('topright',c(expression(paste(lambda[A],'=0.001')),expression(paste(lambda[A],'=0.01')),
                    expression(paste(lambda[A],'=0.05')),expression(paste(lambda[A],'=0.1')),
                    expression(paste(lambda[A],'=1')),expression(paste(lambda[A],'=5')),
                    expression(paste(lambda[A],'=10'))),lty = rep(1,5), 
                    col=c('black','red','green','blue','magenta','cyan','yellow'), ncol=2)
dev.off()