rm(list=ls())
graphics.off()
cat('\014')

library('zoo')
library('tidyverse')
library('patchwork')

thetasDF= rbind.data.frame(
  c(-Inf,seq(1.5,6.5,1),Inf),
  c(-Inf,seq(1.5,6.5,1),Inf),
  c(-Inf,1.5,3.1,3.7,4.3,4.9,6.5,Inf),
  c(-Inf,1.5,2.25,3,5,5.75,6.5,Inf))
colnames(thetasDF)= 0:7

centersDF= matrix(NA,nrow=dim(thetasDF)[1],ncol=dim(thetasDF)[2]-1)
#get the "centers" between the cutoffs
for(i in 1:dim(thetasDF)[1]){
  centers= c(1, 
    rollapply(as.numeric(thetasDF[i,2:(length(thetasDF[i,])-1)]), FUN= mean, width=2),
    ceiling(thetasDF[i,dim(thetasDF)[2]-1])
  )
  centersDF[i,]= centers
}

#smooth= seq(-1,8,.01)
#rollapply(thetas,FUN= function(x){x[2]-x[1]},width=2,align='right')

GaussDF= cbind.data.frame(mu=c(4,1,4,4),sig=c(1.5,2.5,1,3))
p= list()
for(i in 1:dim(GaussDF)[1]){
  print(paste(i))
  
  thetas= as.numeric(thetasDF[i,])
  mu= GaussDF$mu[i]; sig= GaussDF$sig[i]
  
  bars= rollapply(thetas,FUN= function(x){pnorm(mean=mu,sd=sig,x[2])-pnorm(mean=mu,sd=sig,x[1])},
            width=2,align='right')
  # curve(dnorm(x,mu,sig),xlim=c(-1,8))
  # barplot(x=thetas,y=bars)
  
  p[[i]] <- ggplot()+
    geom_function(fun= dnorm,args=list(mean=mu,sd=sig),xlim=c(-1,10))+
    geom_bar(aes(x=centersDF[i,],y=bars),stat='identity',fill='skyblue')
}

p[[1]] / p[[2]] / p[[3]] / p[[4]]

### !!!!!! There is a problem with saving the "geom_bar" part of the plot: apparently 
 #all list entries are overwritten by the last entry
 #it does not matter if in the ggplot block I call "geom_function" or "geom_bar" first, it is
 #always the bars that overwrite previous entries

## B
mu= 3.6
sig= 1.6
thr= c(-Inf, 1.5, 3.4, 3.9, 4.3, 4.9, 6.5, Inf)

prob= rollapply(thr, FUN= function(x){pnorm(x[2],mu,sig)-pnorm(x[1],mu,sig)},width=2)
barplot(height=prob)

## C
b0= 3.16
b1= 4.14E-6
sig= .854

thr= c(-Inf, 1.5, 2.25, 3.1, 4.5,Inf)

x1= 1.6e5
x2= 4.9e5

mu1= b0+b1*x1
prob1= rollapply(thr, FUN= function(t){pnorm(t[2],mu1,sig)-pnorm(t[1],mu,sig)},width=2)
barplot(height=prob1,horiz=TRUE)

mu2= b0+b1*x2
prob2= rollapply(thr,FUN= function(t){pnorm(t[2],mu2,sig)-pnorm(t[1],mu2,sig)},width=2)
barplot(height=prob2,horiz=TRUE)

#####
# 2,
# Ordinal response variable with metric predictors 

rm(list=ls())
graphics.off()
cat('\014')
library('rjags')
source('G:/My Drive/learning/20230220_BayesianStatistics/20230504_BayesianDogsBook/DBDA2Eprograms/DBDA2E-utilities.R')

modStr= 
'
data{
   # rescale
   x1_mean= mean(x1)
   x1_sd= sd(x1)
   x2_mean= mean(x2)
   x2_sd= sd(x2)
   for(i in 1:Ntotal){
      z1[i]= (x1[i]-x1_mean)/x1_sd
      z2[i]= (x2[i]-x2_mean)/x2_sd
   }
}
model{
   for (i in 1:Ntotal){
    y[i] ~ dcat(pr[i,1:L])
    
    pr[i,1]= pnorm(theta[1],mu[i],1/sig^2)
    pr[i,L]= 1 - pnorm(theta[L-1],mu[i],1/sig^2)
    for(l in 2:(L-1)){
      pr[i,l]= max(0,
                pnorm(theta[l],mu[i],1/sig^2) - pnorm(theta[l-1],mu[i],1/sig^2)
               )
    }
    
    mu[i]= a0 + a1*z1[i] + a2*z2[i]
   }
    
    #Priors
    a0 ~ dnorm((1+L)/2, 1/L^2)
    a1 ~ dnorm(0, 1/L^2)
    a2 ~ dnorm(0, 1/L^2)
    sig ~ dunif(L/1000, L*10)
    #theta[1]= 1.5
    #theta[L-1]= L-1+0.5
    for(l in 2:(L-2)){
      theta[l] ~ dnorm(l+0.5,1/2^2)
    }
  
  #rescale
  b0= a0 - a1*x1_mean/x1_sd - a2*x2_mean/x2_sd
  b1= a1/x1_sd
  b2= a2/x2_sd
}
'

data= read.csv('G:/My Drive/learning/20230220_BayesianStatistics/20230504_BayesianDogsBook/DBDA2Eprograms/Movies.csv')
head(data)
theta= rep(NA, 6)
theta[1]=1.5
theta[6]= 6.5
myDat= list(L=7, x1= data$Year, x2=data$Length, y= data$Rating, Ntotal= dim(data)[1],theta=theta)

jagMod= jags.model(file=textConnection(modStr),data=myDat,n.chains=3)
cS= coda.samples(model=jagMod, variable.names=c('b0','b1','b2','theta','sig'),n.iter=10000)

for(i in 1:length(varnames(cS))){
   diagMCMC(codaObject= cS,parName=varnames(cS)[i]  )
}

graphics.off()
DF1= data.frame(cS[[1]])
head(DF1)

hist(DF1$b0,breaks=40,main='Intercept',col='skyblue')
hist(DF1$b1,breaks=40,main='Years',col='skyblue')
hist(DF1$b2,breaks=40,main='Length',col='skyblue')

plot(type='n',x=1:7, y=1:7)
for(i in 1:dim(DF1)[1]/10){ #showing all points would be too much
  points(x=DF1[i,5:10],y=rep(mean(as.numeric(DF1[i,5:10])),6),col='skyblue')
}

#for each value in the data, sample one step from the chain
set.seed(21092023)
samp_idx= sample(1:dim(DF1)[1],dim(data)[1])

mu_pred= DF1$b0[samp_idx] + DF1$b1[samp_idx]*data$Year + DF1$b2[samp_idx]*data$Length
#find level into which mu falls
lvl_pred= rep(NA,length(mu_pred))
for (i in 1:length(mu_pred)){
  lvl_pred[i]= min(which(mu_pred[i] < DF1[samp_idx[i],5:10]))
}

cols= palette('ggplot2')[1:7]

plot(type='n',x=data$Year,y=data$Length)
text(x=data$Year,y=data$Length, label= lvl_pred, col=cols[lvl_pred])
