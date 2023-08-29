#17.1 Compare models for family income as a function of family size via linear model
 #w/o and with quadratic term
  #based on median income from different states

 #next: try if it works with quadratic term

rm(list=ls())
cat('\014')
graphics.off()

library('rjags')
library('tidyverse')

modStr=' 
data{
  #transform data
  mean_x= mean(x)
  sd_x= sd(x)
  mean_y= mean(y)
  sd_y= sd(y)
  for(i in 1:Ntotal){
    zx[i]= (x[i]-mean_x)/sd_x
    zy[i]= (y[i]-mean_y)/sd_y
  }
}

model{
  for(i in 1:Ntotal){
    zy[i] ~ dt(mu[i],1/zsig^2,nu)   #common sig and nu for all states
    mu[i]= zbeta0[s[i]] + zbeta1[s[i]] * zx[i] + zbeta2[s[i]] * zx[i]^2  ##
  }
  
  for(j in 1:Nstates){
    zbeta0[j] ~ dnorm(mu0,1/sig0^2)
    zbeta1[j] ~ dnorm(mu1,1/sig1^2)
    zbeta2[j] ~ dnorm(mu2,1/sig2^2) ##
  }
  
  mu0 ~ dnorm(M,1/S^2)
  sig0 ~ dunif(L,H)
  mu1 ~ dnorm(M,1/S^2)
  sig1 ~ dunif(L,H)
  mu2 ~ dnorm(M,1/S^2)  ##
  sig2 ~ dunif(L,H)    ##
  
  nuMin1 ~ dexp(1/29.0)
  nu= 1 + nuMin1
  zsig ~ dunif(L,H)
  
  M= 0
  S= 10
  L= 1/1000
  H= 1000
  
  #Retransform
  beta0= zbeta0*sd_y + mean_y - zbeta1*mean_x*sd_y/sd_x + zbeta2*mean_x^2 *sd_y/sd_x^2 ##
  beta1= zbeta1*sd_y/sd_x - 2*zbeta2*mean_x*sd_y/sd_x^2 ##
  beta2= zbeta2*sd_y/sd_x^2  ##
  sig= zsig*sd_y
}'

csv= read.csv('H:/My Drive/learning/20230220_BayesianStatistics/20230504_BayesianDogsBook/DBDA2Eprograms/IncomeFamszState.csv')
csv$s= as.numeric(as.factor(csv$State))
head(csv)
dataList= list(
  y= csv$Income,
  x= csv$Famsz,
  s= csv$s,
  Ntotal= dim(csv)[1],
  Nstates= length(unique(csv$State))
)

jagMod= jags.model(file= textConnection(modStr),data=dataList,n.chains=3)
cS= coda.samples(model= jagMod,variable.names= c('beta0','beta1','sig','nu','beta2'),n.iter=10000)

chain1= data.frame(cS[[1]])
chain2= data.frame(cS[[2]])
chain3= data.frame(cS[[3]])
head(chain1)

#plot state 20
state20= levels(as.factor(csv$State))[20]
set.seed(11082023)
rand_idx= sample(1:dim(chain1)[1],size=50)
quad_fun= function(x,b0,b1,b2){
  y= b0+b1*x+b2*x^2
  return(y)
}
ggplot(data= csv[csv$State==state20,],aes(x=Famsz,y=Income))+
  geom_point()+
  coord_cartesian(ylim=c(0,max(csv$Income[csv$State==state20])))+
  #geom_abline(intercept=chain1$beta0.20.[rand_idx],slope=chain1$beta1.20.[rand_idx],col='blue')
  geom_function(fun=quad_fun,args=list(b0=chain1$beta0.20.[1],
                                       b1=chain1$beta1.20.[1],
                                       b2=chain1$beta2.20.[1] ))

#generate posterior predictive of credible values (variation not considered)
chain1_50= chain1[rand_idx,]
head(csv)
postPredMtx= matrix(NA,nrow=dim(csv)[1],ncol=length(rand_idx))  #50 simulations for each datapoint
for(i in 1:dim(csv)[1]){
  postPredMtx[i,]= chain1_50[,csv$s[i]] + chain1_50[,csv$s[i]+52] * csv$Famsz[i]
}
head(postPredMtx)

#calculate mean difference of each datapoint from prediction (absolute difference!)
csv$meanDiff= rep(NA,dim(csv)[1])
for(i in 1:dim(csv)[1]){
  csv$meanDiff[i]= mean(csv$Income[i]-postPredMtx[i,])
}

ggplot(data=csv,aes(x=State,y=meanDiff,group=as.factor(Famsz),fill=as.factor(Famsz)))+
  geom_bar(stat='identity',position='dodge')

#the data show that the true value is below the predicted value for small and big families
 #the true value is above predicted value for intermediate-size families

####### ///////////////////////////////////// ###########
# using Stan

   #not working: try scaling

library('rstan')

stanMod= 
'data{
  int<lower=1> Nstates;
  int<lower=1> Ntotal;
  int s[Ntotal];
  real x[Ntotal];
  real y[Ntotal];
}
transformed data{
  real M;
  real<lower=0> S;
  real<lower=0> L;
  real<lower=0> H;
  M= mean(y)
  S= sd(y)
  L= M/1000;
  H= M*1000;
}
parameters{
  real<lower=0> nuMin1;
  real <lower=0> mu0;
  real mu1;
  real<lower=0> sig;
  real<lower=0> sig0;
  real<lower=0> sig1;
  real beta0[Nstates];
  real beta1[Nstates];
  //real mu[Ntotal];
}
transformed parameters{
    real<lower=1> nu;
    nu= nuMin1+1;
}
model{
  mu0 ~ normal(M,S);
  mu1 ~ normal(M,S);
  sig0 ~ uniform(L,H);
  sig1 ~ uniform(L,H);
  sig ~ uniform(L,H);
  nuMin1 ~ exponential(1/29.0);
  for(j in 1:Nstates){
    beta0[j] ~ normal(mu0,sig0);
    beta1[j] ~ normal(mu1,sig1);
  }
  for(i in 1:Ntotal){
    //mu[i] = beta0[s[i]] + beta1[s[i]]*x[i];
    y[i] ~ student_t(nu,beta0[s[i]] + beta1[s[i]]*x[i],sig);
  }
}
'

dataList2= dataList
dataList2$M= mean(dataList$y)
dataList2$S= sd(dataList$y)

DSO= stan_model(model_code= stanMod)
# stanOUT= sampling(DSO,data=dataList2,iter=1000,warmup=500,chains=1,
#     init=list(list(nuMin1= 10,mu0=mean(csv$Income),mu1=0,sig=mean(csv$Income)/10,sig0=mean(csv$Income)/10,sig1=1,
#               beta0=rep(mean(csv$Income),dataList2$Nstates),beta1=rep(0,dataList2$Nstates))))
#with the values you got from rjags
stanOUT= sampling(DSO,data=dataList2,iter=20000,warmup=5000,chains=1,
  init=list(list(nuMin1= 16,mu0=50000,mu1=2000,
                 sig=7000,sig0=10000,sig1=100,
                 beta0=rep(50000,dataList2$Nstates),beta1=rep(2000,dataList2$Nstates))))

################### ////////////////////////////////////// ###############
#Stan model with scaling
library('rstan')






stanScalMod= '
data{
  int<lower=1> Ntotal;
  real<lower=0> y[Ntotal];
  real<lower=0> x[Ntotal];
  int<lower=1> Nstates;
  int<lower=0> s[Ntotal];
  //real mu;
}
transformed data{
  real ymean;
  real <lower=0> ysd;
  real xmean;
  real<lower=0> xsd;
  real zx[Ntotal];
  real zy[Ntotal];
  ymean= mean(y);
  ysd= sd(y);
  xmean= mean(x);
  xsd= sd(x);
  for(i in 1:Ntotal){  //apparently not vectorizable
    zx[i]= (x[i]-xmean)/xsd; //x; //
    zy[i]= (y[i]-ymean)/ysd; //y; //
  }
}
parameters{
  real muzbeta0;
  real<lower=0> sdzbeta0;
  real muzbeta1;
  real<lower=0> sdzbeta1;
  real<lower=0> zsigma;
  real<lower=0> nuMin1;
  real zbeta0[Nstates];
  real zbeta1[Nstates];
}
transformed parameters{
  real<lower=1> nu;
  real beta0[Nstates];
  real beta1[Nstates];
  real sigma;
  nu= nuMin1+1;
  for(j in 1:Nstates){ //apparently not vectorizeable
    beta0[j]= ymean + zbeta0[j]*ysd + zbeta1[j]*xmean*ysd/xsd; //zbeta0; //
    beta1[j]= zbeta1[j]*ysd/xsd; //zbeta1; //
  }
  sigma= zsigma * ysd;
}
model{
  muzbeta0 ~ normal(0,1);
  sdzbeta0 ~ uniform(1.0/1000,1000);
  muzbeta1 ~ normal(0,1);
  sdzbeta1 ~ uniform(1.0/1000,1000);
  zbeta0 ~ normal(muzbeta0, sdzbeta0);
  zbeta1 ~ normal(muzbeta1, sdzbeta1);
  zsigma ~ uniform(1.0/1000, 1000);
  nuMin1 ~ exponential(1/29.0);
  for(i in 1:Ntotal){
    //mu[i]= zbeta0[s[i]] + zbeta1[s[i]]*x[i];
    zy[i] ~ student_t(nu, zbeta0[s[i]] + zbeta1[s[i]]*x[i], zsigma);
  }
}'

DSO_scale= stan_model(model_code= stanScalMod)
stanOUT_scale= sampling(object= DSO_scale,data=dataList, iter=1000, warmup=500, chains=3)

chain1= (as.matrix(stanOUT_scale))
hist(chain1$beta1.1.,col='skyblue')
