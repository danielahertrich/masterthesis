library(fields)

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

## Look at range of values in H and A
range(dict)
range(score)

###########################################
## Print first 20 atoms of dictionary
###############
## Prepare color palette for blue(neg)-white(0)-red(pos) depiction of dictonary elements
## (adjusted code from 'https://stackoverflow.com/questions/29262824/r-center-color-palette-on-0' for our purpose)
nHalf <- round(nrow(dict)/2)
Min <- -0.01                  #!# choose according to dictionary range      
Max <-  0.5                   #!#choose according to dictionary range 
Thresh <-0
## Make vector of colors for values below threshold
rc1 <- colorRampPalette(colors = c("blue", "white"), space="Lab")(nHalf)    
## Make vector of colors for values above threshold
rc2 <- colorRampPalette(colors = c("white", "red"), space="Lab")(nHalf)
rampcols <- c(rc1, rc2)
rb1 <- seq(Min, Thresh, length.out=nHalf+1)
rb2 <- seq(Thresh, Max, length.out=nHalf+1)[-1]
rampbreaks <- c(rb1, rb2)

## dictionary elements (red and blue) (print first 20 atoms)
png(filename = paste("Dictionary_elements_",art,q,"it",it,".png",sep = ""), width = 770, height = 530)
par(mfrow=c(4,5),mar=rep(0,4),oma=c(0.5,rep(0,2),4.7))
for (qc in 1:20){
  image(matrix(dict[qc,], nrow=92),col = rampcols,
                                    breaks = rampbreaks, axes=FALSE)                     
}
## add color scale
image.plot(matrix(dict[20,], nrow=92),col = rampcols,
           breaks = rampbreaks, axes=FALSE, axis.args = list(cex.axis = 2), legend.width = 3.5, legend.only = TRUE,
           smallplot= c(0.95,1,0.04,1))
dev.off()

## some info about dictionaries
spH1 <- spH2 <- spH3 <- numeric(nrow(score))
for(i in 1:nrow(dict)){
  spH1[i] <- sum(dict[i,]==0)
  spH2[i] <- sum(abs(dict[i,])<=0.005)
  spH3[i] <- sum(abs(dict[i,]) <=1)
}
sum(spH1);sum(spH1)/length(dict);mean(spH1);range(spH1)
sum(spH2);sum(spH2)/length(dict);mean(spH2);range(spH2)
sum(spH3);sum(spH3)/length(dict);mean(spH3);range(spH3)


###########################################
## Print image of score matrix
###############
png(filename = paste("Scores",art,q,"it",it,".png",sep = ""), width = 520, height = 310)
par(mfrow=c(1,1),mar=c(2,2,1,0),oma=c(0,0,0,1.5))
image.plot(t(score)[,nrow(score):1],axis.args = list(cex.axis=1.5))
dev.off()



## some info about dictionaries
spA1 <- spA2 <- spA3 <- numeric(nrow(score))
for(i in 1:nrow(score)){
  spA1[i] <- sum(score[i,]==0)
  spA2[i] <- sum(abs(score[i,])<=0.005)
  spA3[i] <- sum(abs(score[i,]) <=1)
}
sum(spA1);sum(spA1)/length(score);mean(spA1);range(spA1)
sum(spA2);sum(spA2)/length(score);mean(spA2);range(spA2)
sum(spA3);sum(spA3)/length(score);mean(spA3);range(spA3)