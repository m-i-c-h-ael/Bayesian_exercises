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

###########
#21.2. Nominal predictor: Modeled as logistic function with normal prior vs.
 #beta distribution w/o transformation

#A, Specify model
rm(list=ls())
graphics.off()
cat('\014')

library('rjags')
source('G:/My Drive/learning/20230220_BayesianStatistics/20230504_BayesianDogsBook/DBDA2Eprograms/DBDA2E-utilities.R')

modStr=
'model{
  for(i in 1:Ntotal){
    x[i] ~ dbern(theta)
  }
  theta= ilogit(beta)
  
  #Prior
  beta ~ dnorm(0,1/2^2)
}'

modStr_beta=
'model{
  for(i in 1:Ntotal){
    x[i] ~ dbern(theta)
  }
  #Prior
  theta ~ dbeta(1,1)
}'

data= read.csv("../DBDA2Eprograms/z15N50.csv")
head(data)
dataList= list(x=data$y, Ntotal= dim(data)[1])

jagMod= jags.model(file=textConnection(modStr),data=dataList,n.chains=3)

#B, Plot prior
 #Prior from dedicated model
modStrPr=
  'model{
  theta= ilogit(beta)
  
  #Prior
  beta ~ dnorm(0,1/2^2)
}'

#dataListPr= list(x=data$y,Ntotal=0)
jagModPr= jags.model(file=textConnection(modStrPr),data=dataList,n.chains=3)
cS_Pr= coda.samples(model=jagModPr,variable.names = c('theta'),n.iter=10000)

diagMCMC(cS_Pr)
dev.new()

  #Prior from full model w/o data
dataListPr= list(x=0,Ntotal=0)

jagModPr= jags.model(file=textConnection(modStr),data=dataListPr,n.chains=3)
cS_Pr= coda.samples(model=jagModPr,variable.names = c('theta'),n.iter=10000)
diagMCMC(cS_Pr)
b_1_Pr= data.frame(cS_Pr[[1]])[,1]

jagModPr_beta= jags.model(file=textConnection(modStr_beta),data=dataListPr,n.chains=3)
cS_Pr_beta= coda.samples(model=jagModPr_beta,variable.names = c('theta'),n.iter=10000)
diagMCMC(cS_Pr_beta)
b_1_beta_Pr= data.frame(cS_Pr_beta[[1]])[,1]

#C One flip, one head
graphics.off()
oneHead= list(x=1,Ntotal=1)

jagMod1H= jags.model(file=textConnection(modStr),data=oneHead,n.chains=3)
cS_1H= coda.samples(model=jagMod1H,variable.names = c('theta'),n.iter=10000)
diagMCMC(cS_1H)
b_1_1H= data.frame(cS_1H[[1]])[,1]
#dev.new()
par(mfrow=c(2,1))
hist(b_1_Pr,col='skyblue',main='Prior (logistic)')
hist(b_1_1H,col='skyblue',main='1 head (logistic)')


jagMod_beta_1H= jags.model(file=textConnection(modStr_beta),data=oneHead,n.chains=3)
cS_beta_1H= coda.samples(model=jagMod_beta_1H,variable.names = c('theta'),n.iter=10000)
diagMCMC(cS_beta_1H)
b_1_beta_1H= data.frame(cS_beta_1H[[1]])[,1]
par(mfrow=c(2,1))
hist(b_1_beta_Pr,col='skyblue',main='Prior (beta)')
hist(b_1_beta_1H,col='skyblue',main='1 head (beta)')

#Posteriors look similar for the two models; logistic model puts even a bit
 #more credibility on H, as prior had more weight on extreme theta-values

#C 40 flips with 30 heads
graphics.off()
H30= list(x=c(rep(1,30),rep(0,40-30)), Ntotal=40)

jagMod30H= jags.model(file=textConnection(modStr),data=H30,n.chains=3)
cS_30H= coda.samples(model=jagMod30H,variable.names = c('theta'),n.iter=10000)
diagMCMC(cS_30H)
dev.new()
b_1_30H= data.frame(cS_30H[[1]])[,1]
hist(b_1_30H,col='skyblue',main='30 heads (logistic)')

jagMod_beta_30H= jags.model(file=textConnection(modStr_beta),data=H30,n.chains=3)
cS_beta_30H= coda.samples(model=jagMod_beta_30H,variable.names = c('theta'),n.iter=10000)
diagMCMC(cS_beta_30H)
dev.new()
b_1_beta_30H= data.frame(cS_beta_30H[[1]])[,1]
hist(b_1_beta_30H,col='skyblue',main='30 heads (beta)')
 #Posterior distributions look virtually identical