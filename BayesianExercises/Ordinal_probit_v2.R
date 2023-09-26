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

# a
modStr= 
'
data{
   # rescale
   for(j in 1:Nx){
     xm[j]= mean(x[1:Ntotal,j])
     xsd[j]= sd(x[1:Ntotal,j])
     for(i in 1:Ntotal){
        zx[i,j]= (x[i,j]-xm[j])/xsd[j]
     }
   }
}
model{
   for (i in 1:Ntotal){
    y[i] ~ dcat(pr[i,1:nYlevels])
    pr[i,1]= pnorm(thresh[1],mu[i],1/sigma^2)
    for(l in 2:(nYlevels-1)){
      pr[i,l]= max(0,
                pnorm(thresh[l],mu[i],1/sigma^2) - pnorm(thresh[l-1],mu[i],1/sigma^2)
               )
    }
    pr[i,nYlevels]= 1 - pnorm(thresh[nYlevels-1],mu[i],1/sigma^2)
    
    mu[i]= zbeta0 + sum(zbeta[1:Nx]*zx[i,1:Nx])
   }

    #Priors
    zbeta0 ~ dnorm((1+nYlevels)/2, 1/nYlevels^2)
    for(j in 1:Nx){
      zbeta[j] ~ dnorm(0, 1/nYlevels^2)
    }
    
    zsigma ~ dunif(nYlevels/1000, nYlevels*10)
    #thresh[1]= 1.5
    #thresh[nYlevels-1]= nYlevels-1+0.5
    
  #rescale
  beta0= zbeta0 - sum(zbeta[1:Nx] * xm[1:Nx] / xsd[1:Nx])
  beta[1:Nx]= zbeta[1:Nx]/xsd[1:Nx]

  sigma= zsigma
  for(l in 2:(nYlevels-2)){
      thresh[l] ~ dnorm(l+0.5,1/2^2)
    }
}
'

data= read.csv('G:/My Drive/learning/20230220_BayesianStatistics/20230504_BayesianDogsBook/DBDA2Eprograms/Movies.csv')
head(data)
data$Rating_float= data$Rating
rating_float_sort= unique(sort(data$Rating_float))
unique(rating_float_sort)
#data$Rating= as.numeric(as.factor(data$Rating))
data$Rating= match(data$Rating_float, rating_float_sort)  #convert to int
cbind.data.frame(data$Rating_float,data$Rating)

thresh= rep(NA, 6)
thresh[1]=1.5
thresh[6]= 6.5
myDat= list(nYlevels=7, x= cbind.data.frame(data$Year, data$Length), Nx= 2,
            y= data$Rating, Ntotal= dim(data)[1],thresh=thresh)

jagMod= jags.model(file=textConnection(modStr),data=myDat,n.chains=3)
cS= coda.samples(model=jagMod, variable.names=c('beta0','beta','thresh','sigma'),n.iter=10000)

####################
# RUN THE CHAINS
# parameters = c( "beta0" ,  "beta" ,  "sigma", "thresh" ,
#                 "zbeta0" , "zbeta" , "zsigma" )
# adaptSteps = 500  # Number of steps to "tune" the samplers
# burnInSteps = 1000
# numSavedSteps=10000 
# thinSteps=1
# runJagsOut <- run.jags( method=runjagsMethodDefault ,
#                         model= textConnection(modStr) , 
#                         monitor=parameters , 
#                         data=myDat ,  
#                         #inits=initsList , 
#                         n.chains=3 ,
#                         adapt=adaptSteps ,
#                         burnin=burnInSteps , 
#                         sample=ceiling(numSavedSteps/3) ,
#                         thin=thinSteps ,
#                         summarise=FALSE ,
#                         plots=FALSE )
# codaSamples = as.mcmc.list( runJagsOut )
###################
# myData = read.csv( file="Movies.csv" ) # Real data
# myData$Rating = as.numeric(as.factor(myData$Rating))
# data= myData
# yName = "Rating" ; xName = c("Year","Length")
# 
# y = data[,yName]
# x = as.matrix(data[,xName],ncol=length(xName))
# # Do some checking that data make sense:
# if ( any( !is.finite(x) ) ) { stop("All x values must be finite.") }
# if ( any( y!=round(y) ) ) { stop("All y values must be integers (whole numbers).") }
# if ( any( y < 1 ) ) { stop("All y values must be 1 or larger.") }
# # COMPRESS OUT ANY EMPTY VALUES OF Y:
# yOrig=y
# y=as.numeric(factor(y,levels=names(table(y))))
# if ( any(y != yOrig) ) { 
#   warning("*** WARNING: Y RE-CODED TO REMOVE EMPTY LEVELS ***")
# }
# nYlevels = max(y)  
# thresh = rep(NA,nYlevels-1)
# thresh[1] = 1 + 0.5
# thresh[nYlevels-1] = nYlevels-1 + 0.5
# 
# dataList = list(
#   x = x ,
#   y = y ,
#   nYlevels = nYlevels ,
#   thresh = thresh ,
#   Nx = dim(x)[2] ,
#   Ntotal = dim(x)[1]
# )
# 
# 
# jagMod= jags.model(file="TEMPmodel.txt",data=dataList,n.chains=3)
# cS= coda.samples(model=jagMod, variable.names=c('beta0','beta','thresh','sigma'),n.iter=10000)

##################

for(i in 1:length(varnames(cS))){
   diagMCMC(codaObject= cS,parName=varnames(cS)[i]  )
}

graphics.off()
DF1= data.frame(cS[[1]])
head(DF1)

hist(DF1$beta0,breaks=40,main='Intercept',col='skyblue')
hist(DF1$beta.1.,breaks=40,main='Years',col='skyblue')
hist(DF1$beta.2.,breaks=40,main='Length',col='skyblue')

plot(type='n',x=1:7, y=1:7)
for(i in 1:dim(DF1)[1]/10){ #showing all points would be too much
  points(x=DF1[i,5:10],y=rep(mean(as.numeric(DF1[i,5:10])),6),col='skyblue')
}

#for each value in the data, sample one step from the chain
set.seed(21092023)
samp_idx= sample(1:dim(DF1)[1],dim(data)[1])

mu_pred= DF1$beta0[samp_idx] + DF1$beta.1.[samp_idx]*data$Year + DF1$beta.2.[samp_idx]*data$Length
#find level into which mu falls
lvl_pred= rep(NA,length(mu_pred))
for (i in 1:length(mu_pred)){
  lvl_pred[i]= min( c(7, which(mu_pred[i] < DF1[samp_idx[i],5:10]) ))
}

#convert to the original ratings
rating_pred= rating_float_sort[lvl_pred]

cols= palette('ggplot2')[1:7]

plot(type='n',x=data$Year,y=data$Length)
text(x=data$Year,y=data$Length, label= rating_pred, col=cols[lvl_pred])

######## Robustify
modStr_robu= 
  '
data{
   # rescale
   for(j in 1:Nx){
     xm[j]= mean(x[1:Ntotal,j])
     xsd[j]= sd(x[1:Ntotal,j])
     for(i in 1:Ntotal){
        zx[i,j]= (x[i,j]-xm[j])/xsd[j]
     }
   }
   for (l in 1:nYlevels){
      guess_vec[l]= 1/nYlevels
   }
}
model{
   for (i in 1:Ntotal){
    y[i] ~ dcat(pr_adj[i,1:nYlevels])
    pr[i,1]= pnorm(thresh[1],mu[i],1/sigma^2)
    for(l in 2:(nYlevels-1)){
      pr[i,l]= max(0,
                pnorm(thresh[l],mu[i],1/sigma^2) - pnorm(thresh[l-1],mu[i],1/sigma^2)
               )
    }
    pr[i,nYlevels]= 1 - pnorm(thresh[nYlevels-1],mu[i],1/sigma^2)
    
    pr_adj[i,1:nYlevels]= (1-alpha)*pr[i,1:nYlevels]+ alpha*guess_vec[1:nYlevels]  #robustify
    
    mu[i]= zbeta0 + sum(zbeta[1:Nx]*zx[i,1:Nx])
   }

    #Priors
    zbeta0 ~ dnorm((1+nYlevels)/2, 1/nYlevels^2)
    for(j in 1:Nx){
      zbeta[j] ~ dnorm(0, 1/nYlevels^2)
    }
    
    zsigma ~ dunif(nYlevels/1000, nYlevels*10)
    #thresh[1]= 1.5
    #thresh[nYlevels-1]= nYlevels-1+0.5
    
    alpha~ dbeta(1,9)   ## robustify
    
  #rescale
  beta0= zbeta0 - sum(zbeta[1:Nx] * xm[1:Nx] / xsd[1:Nx])
  beta[1:Nx]= zbeta[1:Nx]/xsd[1:Nx]

  sigma= zsigma
  for(l in 2:(nYlevels-2)){
      thresh[l] ~ dnorm(l+0.5,1/2^2)
    }
}
'

jagMod_robu= jags.model(file=textConnection(modStr_robu),data=myDat,n.chains=3)
cS_robu= coda.samples(model=jagMod_robu, variable.names=c('beta0','beta','thresh','sigma','alpha'),n.iter=10000)

for(i in 1:length(varnames(cS_robu))){
  diagMCMC(codaObject= cS_robu,parName=varnames(cS_robu)[i]  )
}

graphics.off()
DF1_robu= data.frame(cS_robu[[1]])
head(DF1_robu)

hist(DF1_robu$beta0,breaks=40,main='Intercept',col='skyblue')
hist(DF1_robu$beta.1.,breaks=40,main='Years',col='skyblue')
hist(DF1_robu$beta.2.,breaks=40,main='Length',col='skyblue')
hist(DF1_robu$alpha,breaks=40,main='Length',col='skyblue')
hist(DF1_robu$sigma,breaks=40,main='Length',col='skyblue')

plot(type='n',x=1:7, y=1:7)
for(i in 1:dim(DF1_robu)[1]/10){ #showing all points would be too much
  points(x=DF1_robu[i,6:11],y=rep(mean(as.numeric(DF1_robu[i,6:11])),6),col='skyblue')
}

#for each value in the data, sample one step from the chain
set.seed(21092023)
samp_idx= sample(1:dim(DF1_robu)[1],dim(data)[1])

mu_pred_robu= DF1_robu$beta0[samp_idx] + DF1_robu$beta.1.[samp_idx]*data$Year + 
  DF1_robu$beta.2.[samp_idx]*data$Length
#find level into which mu falls
lvl_pred_robu= rep(NA,length(mu_pred_robu))
for (i in 1:length(mu_pred_robu)){
  lvl_pred_robu[i]= min( c(7, which(mu_pred_robu[i] < DF1_robu[samp_idx[i],6:11]) ))
}

#convert to the original ratings
rating_pred_robu= rating_float_sort[lvl_pred_robu]

cols= palette('ggplot2')[1:7]

plot(type='n',x=data$Year,y=data$Length)
text(x=data$Year,y=data$Length, label= rating_pred_robu, col=cols[lvl_pred_robu])


# b: Robust via t-dist
I don`t know what the arguments for 'pt' in rjags are; for dt: dt(mu, df), but for pt you need a
 way to specify the quantile (here: threshold value), don`t you?
I tried dt(t,df) as in R, but then I get error "Incorrect number of arguments in function pt"