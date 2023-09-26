rm(list=ls())
graphics.off()
cat('\014')

library('rjags')
source('H:/My Drive/learning/20230220_BayesianStatistics/20230504_BayesianDogsBook/DBDA2Eprograms/DBDA2E-utilities.R')

#22.1
#A Intuition about coefficients
set.seed(07092023)

# Reference: 2
x= seq(0,2,.01)

coefDF= cbind.data.frame(
  b0= c(1,0,-2),
  b1= c(-2,0,2)
)
coefDF

linfun= function(x,b){
  y= exp(b[1]+b[2]*x)
  return(y)
}

p= matrix(nrow=dim(coefDF)[1],ncol=length(x))
for(r in 1:nrow(coefDF)){
   p[r,]= linfun(x=x,b=as.numeric(coefDF[r,]))/ 
       apply( apply(coefDF,1,function(b){linfun(x=x,as.numeric(b))}),1,sum )
}

colors= c('black','red','gold')
plot(type='n',x=0,y=0,xlim=range(x),ylim=c(0,1),xlab='x',ylab='prob')
for(r in 1:nrow(p)){
  lines(x,p[r,],col=colors[r])
}

samp= apply(p,2,function(x){sample(c(1,2,3),size=1,prob=x)})
plot(type='n',x=0,y=0,xlim=range(x),ylim=c(0,1),xlab='x',ylab='')
text(x,0.5+rnorm(length(samp),mean=0,sd=.1),labels=samp,col=colors[samp])


# Reference: 1
x= seq(0,2,.01)

coefDF= cbind.data.frame(
  b0= c(0,-0.8,-2.5),
  b1= c(0,1.5,2.7)
)
coefDF

linfun= function(x,b){
  y= exp(b[1]+b[2]*x)
  return(y)
}

p= matrix(nrow=dim(coefDF)[1],ncol=length(x))
for(r in 1:nrow(coefDF)){
  p[r,]= linfun(x=x,b=as.numeric(coefDF[r,]))/ 
    apply( apply(coefDF,1,function(b){linfun(x=x,as.numeric(b))}),1,sum )
}

colors= c('black','red','gold')
plot(type='n',x=0,y=0,xlim=range(x),ylim=c(0,1),xlab='x',ylab='prob')
for(r in 1:nrow(p)){
  lines(x,p[r,],col=colors[r])
}

samp= apply(p,2,function(x){sample(c(1,2,3),size=1,prob=x)})
plot(type='n',x=0,y=0,xlim=range(x),ylim=c(0,1),xlab='x',ylab='')
text(x,0.5+rnorm(length(samp),mean=0,sd=.1),labels=samp,col=colors[samp])

#B, Softmax model
softM= '
data{
  mean_x1= mean(x1)
  sd_x1= sd(x1)
  mean_x2= mean(x2)
  sd_x2= sd(x2)
  for(i in 1:Ntotal){
    z1[i]= (x1[i]-mean_x1)/sd_x1
    z2[i]= (x2[i]-mean_x2)/sd_x2
  }
}
model{
  for(i in 1:Ntotal){
    y[i] ~ dcat(mu[1:Nout,i]) #vectorized over k
    mu[1:Nout,i]= explambda[1:Nout,i]/expall[i]
  
    for(k in 1:Nout){
      explambda[k,i]= exp( b0[k] + b1[k]*z1[i] + b2[k]*z2[i] )
    }
    expall[i]= sum( explambda[1:Nout,i] )        #b0,b1,b2 vectorized over k
  }
  
  for(k in 2:Nout){
    b0[k] ~ dnorm(0,1/20^2)
    b1[k] ~ dnorm(0,1/20^2)
    b2[k] ~ dnorm(0,1/20^2)
  }
  b0[1]= 0
  b1[1]= 0
  b2[1]= 0
  
  # Transform to original scale:
  c1 <- b1 / sd_x1 
  c2 <- b2 / sd_x2 
  c0 <- b0 - b1*mean_x1/sd_x1 - b2*mean_x2/sd_x2
}
'

# dataList= list(
#   x1=x, x2=x,
#   y= samp,   #pretend it also depends on x2, not only x1
#   Nout= length(unique(samp)),
#   Ntotal= length(samp)
# )
dSoft= read.csv( file="H:/My Drive/learning/20230220_BayesianStatistics/20230504_BayesianDogsBook/DBDA2Eprograms/SoftmaxRegData2.csv" )
  #assuming regData1 is the left and regData2 the right side of fig. 22.1
head(dSoft)
softList= list(
  x1= dSoft$X1, x2= dSoft$X2,
  y= dSoft$Y,
  Nout= length(unique(dSoft$Y)), Ntotal= dim(dSoft)[1]
)

jagSoft= jags.model(file=textConnection(softM),data=softList,n.chains=3)
cSoft= coda.samples(model=jagSoft,variable.names = c('c0','c1','c2'),n.iter=10000)

for (i in 1:length(varnames(cSoft))){
  diagMCMC(cSoft,parName=varnames(cSoft)[i])
}
 #warnings OK: references have fixed values
graphics.off()

DF1= data.frame(cSoft[[1]])
hist(DF1[,'c0.1.'])
hist(DF1[,'c0.2.'])
hist(DF1[,'c0.3.'])
hist(DF1[,'c0.4.'])

hist(DF1[,'c1.1.'])
hist(DF1[,'c1.2.'])
hist(DF1[,'c1.3.'])
hist(DF1[,'c1.4.'])

hist(DF1[,'c2.1.'])
hist(DF1[,'c2.2.'])
hist(DF1[,'c2.3.'])
hist(DF1[,'c2.4.'])

# Posterior predictive
predFun= function(x1,x2,DF){
  #subset estimated coefficients to fit data length
  DF= DF[1:length(x1),]
  
  y1= rep(0,dim(DF)[1])
  y2= DF[,'c0.2.'] + DF[,'c1.2.']*x1 + DF[,'c2.2.']*x2
  y3= DF[,'c0.3.'] + DF[,'c1.3.']*x1 + DF[,'c2.3.']*x2
  y4= DF[,'c0.4.'] + DF[,'c1.4.']*x1 + DF[,'c2.4.']*x2
  
  all= exp(y1) + exp(y2) + exp(y3) + exp(y4)
  
  p1= exp(y1)/all
  p2= exp(y2)/all
  p3= exp(y3)/all
  p4= exp(y4)/all
  
  return(list(p1=p1,p2=p2,p3=p3,p4=p4))
}

predDF= predFun(x1= dSoft$X1, x2= dSoft$X2, DF=DF1)
predDF= data.frame(predDF)

sampOutcome= apply(predDF,1,function(x) {sample(x=1:4,size=1,prob=x)} )

colors= c('black','red','gold','blue')
plot(x=x1,y=x2,type='n')
for(i in 1:length(sampOutcome)){
  text(x=x1[i],y=x2[i],label=sampOutcome[i],col=colors[sampOutcome[i]])
}


###-------------------------------###
# 22.5  - no working solution found
rm(list=ls())
graphics.off()
cat('\014')

library('rstan')

modStr='
data{
  int<lower=1> Ntotal;
  real x1[Ntotal];
  real x2[Ntotal];
  int<lower=1> y[Ntotal];
  int<lower=1> Nout;
}
transformed data{
  real mean_x1;
  real<lower=0> sd_x1;
  real mean_x2;
  real<lower=0> sd_x2;
  real z1;
  real z2;
  real lambda;
  
  mean_x1= mean(x1);
  sd_x1= sd(x1);
  mean_x2= mean(x2);
  sd_x2= sd(x2);
  z1= (x1-mean_x1)/sd_x1;
  z2= (x2-mean_x2)/sd_x2;
}

parameters{
  real a0[Nout];
  real a1[Nout];
  real a2[Nout];
}
transformed parameters{
  real b0[Nout];
  real b1[Nout];
  real b2[Nout];
  
  b0= a0- a1*mean_x1/sd_x1 - a2*mean_x2/sd_x2;    //vectorized over output categories
  b1= a1/sd_x1;    //vectorized over output categories
  b2= a2/sd_x2;    //vectorized over output categories
}
model{
  # Priors
  a0[1]= 0;
  a1[1]= 0;
  a2[1]= 0;
  
  for(k in 2:Nout){
    a0[k] ~ normal(0,20);
    a1[k] ~ normal(0,20);
    a2[k] ~ normal(0,20);
  }

  // Softmax regression model
  for(i in 1:Ntotal){  
    lambda[i]= a0 + a1*z1[i] + a2*z2[i];    //vectorized over output categories
    y[i] ~ categorical( softmax(lambda[i]) );
  }
}'

myData = read.csv( file="H:/My Drive/learning/20230220_BayesianStatistics/20230504_BayesianDogsBook/DBDA2Eprograms/SoftmaxRegData1.csv" )
head(myData)
datList= list(x1=myData$X1, x2=myData$X2, Ntotal=dim(myData)[1], Nout=length(unique(myData$Y)),
              y=myData$Y)

sM= stan_model(model_code= modStr)



#####-----------------------
modStr2=
'
data{
  int<lower=1> Ntotal;
  int<lower=1> Nout;
  int<lower=1,upper=Nout> y[Ntotal];
}

parameters{
  real<lower=0,upper=1> a[Nout];
}

model{
  for(k in 1:Nout){
    a[k] ~ uniform(0,1);
  }

  for(i in 1:Ntotal){
    y[i] ~ categorical( a[1:Nout] ) ;
  }
}
'

datList2= list(Ntotal=dim(myData)[1], Nout=length(unique(myData$Y)),
              y=myData$Y)

sM2= stan_model(model_code= modStr2)
