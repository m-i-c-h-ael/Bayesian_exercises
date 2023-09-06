rm(list=ls())
cat('\014')
graphics.off()

library('rjags')
source('G:/My Drive/learning/20230220_BayesianStatistics/20230504_BayesianDogsBook/DBDA2Eprograms/DBDA2E-utilities.R')

logis= function(x){
  y= 1/(1+exp(-x))
  return(y)
}

logis(2)

#simple logistic regression with one metric explanatory variable
modStr= '
data{
  #re-scale x
  mean_x= mean(x)
  sd_x= sd(x)
  for(i in 1:Ntotal){
    z[i]= (x[i]-mean(x))/sd_x
  }
}
model{
  for(i in 1:Ntotal){
    y[i] ~ dbern(theta[i])
    theta[i]= ilogit(a0 + a1*z[i])
  }
  
  #prior
  a0 ~ dnorm(0,1/2^2)     #sd=2 -> .12 < theta < .88
  a1 ~ dnorm(0,1/2^2)
  
  #re-transform
  b0= a0+ a1/sd_x*mean_x
  b1= a1/sd_x
}'

#robust logistic regression with one metric explanatory variable
modStr_rob= '
data{
  #re-scale x
  mean_x= mean(x)
  sd_x= sd(x)
  for(i in 1:Ntotal){
    z[i]= (x[i]-mean(x))/sd_x
  }
}
model{
  for(i in 1:Ntotal){
    y[i] ~ dbern(theta[i])
    theta[i]= 1/2*guess + (1-guess)*ilogit(a0 + a1*z[i])
  }
  
  #prior
  a0 ~ dnorm(0,1/2^2)     #sd=2 -> .12 < theta < .88
  a1 ~ dnorm(0,1/2^2)
  guess ~ dbeta(1,9)
  
  #re-transform
  b0= a0+ a1/sd_x*mean_x
  b1= a1/sd_x
}'

##########
## 21.1
set.seed(06092023)

# A, Simulate noisy dichotomous data
N=500
#x= runif(n=N,min=0,max=1)
x= rnorm(n=N,mean=0,sd=1)
x= (x-mean(x))/sd(x)  #center, scale

b0=0
b1=4
mu= 1/(1+exp(-(b0+b1*x)))

alpha= 0.1
mu_robust= 1/2*alpha + mu*(1-alpha)

y= rbinom(n=length(x),size=1,prob= mu_robust)

plot(x,y)

# B, Non-robust robust logistic regression
dataList= list(x=x,y=y,Ntotal= length(y))

jagMod= jags.model(file=textConnection(modStr),data=dataList,n.chains=3)
cS= coda.samples(model=jagMod,variable.names=c('b0','b1'),n.iter=10000)

for(i in 1:length(varnames(cS))){
  diagMCMC(codaObject=cS , parName=varnames(cS)[i])
}
graphics.off()

means1= apply(data.frame(cS[[1]]),2,mean)
means2= apply(data.frame(cS[[2]]),2,mean)
means3= apply(data.frame(cS[[3]]),2,mean)
x1= seq(-1.5,1.5,.01)
plot(x,y)
#lines(curve(1/(1+exp(-(means1[1]+means1[2]*x)))))
lines( x1, 1/(1+exp(-(means1[1]+means1[2]*x1))),col='blue' )
lines( x1, 1/(1+exp(-(means2[1]+means2[2]*x1))),col='violet' )
lines( x1, 1/(1+exp(-(means3[1]+means3[2]*x1))),col='pink' )

DF1= data.frame(cS[[1]])
hist(DF1$b0,breaks=50,col='skyblue')
hist(DF1$b1,breaks=50,col='skyblue')
DF2= data.frame(cS[[1]])
hist(DF2$b0,breaks=50,col='skyblue')
hist(DF2$b1,breaks=50,col='skyblue')
 #the estimates are off: posterior mode of b0 is about -0.3 (instead of 0) and
                                  #     of b1 is about 2.7 (instead of 4)
  #when x are generated from normal distribution: mode of b0 is 0, but b1 is
     #around 2.5

#C, Robust logistic regression
jagMod_rob= jags.model(file=textConnection(modStr_rob),data=dataList,n.chains=3)
cS_rob= coda.samples(model=jagMod_rob,variable.names=c('b0','b1','guess'),n.iter=10000)

for(i in 1:length(varnames(cS_rob))){
  diagMCMC(codaObject=cS_rob , parName=varnames(cS_rob)[i])
}
graphics.off()

means1_rob= apply(data.frame(cS_rob[[1]]),2,mean)
means2_rob= apply(data.frame(cS_rob[[2]]),2,mean)
means3_rob= apply(data.frame(cS_rob[[3]]),2,mean)
x1= seq(-1.5,1.5,.01)
plot(x,y)
#lines(curve(1/(1+exp(-(means1[1]+means1[2]*x)))))
logis1_rob= 1/(1+exp(-(means1_rob[1]+means1_rob[2]*x1)))
lines( x1, 1/2*means1_rob[3]+(1-means1_rob[3])*logis1_rob ,col='blue' )

logis2_rob= 1/(1+exp(-(means2_rob[1]+means2_rob[2]*x1)))
lines( x1, 1/2*means2_rob[3]+(1-means2_rob[3])*logis2_rob ,col='violet' )

logis3_rob= 1/(1+exp(-(means3_rob[1]+means3_rob[2]*x1)))
lines( x1, 1/2*means3_rob[3]+(1-means3_rob[3])*logis3_rob ,col='pink' )

DF1_rob= data.frame(cS_rob[[1]])
hist(DF1_rob$b0,breaks=50,col='skyblue')
hist(DF1_rob$b1,breaks=50,col='skyblue')
DF2_rob= data.frame(cS_rob[[1]])
hist(DF2_rob$b0,breaks=50,col='skyblue')
hist(DF2_rob$b1,breaks=50,col='skyblue')
hist(DF2_rob$guess,breaks=50,col='skyblue')
 #the estimates are still off:
  #median of the marginal posterior for b0 is around -0.25 (instead of 0) and
                                  # for b1 around 3 (instead of 4)
  #guess is underestimated with modes at 0 and 0.07 (instead of 0.1)

  #with data from normal distribution, the estimates are pretty good; mode of b1
   #is around 3.7