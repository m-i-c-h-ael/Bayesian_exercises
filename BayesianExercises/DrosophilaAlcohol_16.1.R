# Compare alcohol consumption of sexually frustrated vs. mated male fruit flies

library('rstan')

#using normal distribution

modStr= '
data{
  int<lower=1> Ntotal;
  real<lower=0> y[Ntotal];
  int s[Ntotal];
  real y_mean;
  real y_sd;
}
transformed data{
  real <lower =0> unif_low;
  real <lower=0> unif_up;
  real <lower=0> normal_sigma;
  unif_low= y_sd/1000;
  unif_up= y_sd*1000;
  normal_sigma= y_sd*100;
}
parameters{
  real mu[2];
  real <lower=0> sigma[2];
}
model{
  mu[2] ~ normal(y_mean,normal_sigma);
  sigma[2] ~ uniform(unif_low, unif_up);
  for(i in 1:Ntotal){
    y[i] ~ normal(mu[s[i]],sigma[s[i]]);
  }
}
'

data= read.csv('H:/My Drive/learning/20230220_BayesianStatistics/20230504_BayesianDogsBook/DBDA2Eprograms/ShohatOphirKAMH2012dataReduced.csv')
head(data)
data2= list(
  y= data$GrandTotal,
  s= as.numeric(as.factor(data$Group)),
  Ntotal= dim(data)[1],
  y_mean= mean(data$GrandTotal),
  y_sd= sd(data$GrandTotal)
)

DSO= stan_model(model_code= modStr)
stanOUT= sampling(DSO, data= data2, chains=3, iter=1000, warmup=500)

head( as.array(stanOUT)[,1,] )
chain1= as.array(stanOUT)[,1,]
chain2= as.array(stanOUT)[,2,]
chain3= as.array(stanOUT)[,3,]
plot(chain1[,1],type='l')
lines(chain2[,1],col='red')
lines(chain3[,1],col='green')

par(mfrow=c(2,1))
plot(density(chain1[,1]),main= 'mu1')
plot(density(chain1[,2]),main= 'mu2')

# Difference of means
diff= chain1[,1]-chain1[,2]
par(mfrow=c(1,1))
plot(density(diff))

# Posterior predictive
par(mfrow=c(2,1))
postPred1= rnorm(dim(chain1)[1],mean=chain1[,1],sd=chain1[,3])
plot(density(data$GrandTotal[data$Group=='MatedGrouped']),xlim=c(100,220))
lines(density(postPred1),col='red')
postPred2= rnorm(dim(chain1)[1],mean=chain1[,2],sd=chain1[,4])
plot(density(data$GrandTotal[data$Group=='RejectedIsolated']),xlim=c(100,220))
lines(density(postPred2),col='red')


## using student t-distribution
modt= '
data{
  int<lower=1> Ntotal;
  real<lower=0> y_mean;
  real<lower=0> y_sd;
  int x[Ntotal];
  real<lower=0> y[Ntotal];
}
transformed data{
  real<lower=0> sd_norm;
  real<lower=0> u_low;
  real<lower=0> u_up;
  real lambda;
  sd_norm= y_sd*100;
  u_low= y_sd/1000;
  u_up= y_sd*1000;
  lambda= 1/29.;
}
parameters{
  real <lower=0> nu_min1;
  real <lower=0> mu[2];
  real <lower=0> sigma[2];
}
transformed parameters{
  real<lower=1> nu;
  nu= nu_min1+1;
}
model{
  mu ~ normal(y_mean,sd_norm);
  sigma ~ uniform(u_low,u_up);
  nu_min1 ~ exponential(lambda);
  for(i in 1:Ntotal){
    y[i] ~ student_t(nu,mu[x[i]],sigma[x[i]]);
  }
}
'
data3= list(
  y= data$GrandTotal,
  x= as.numeric(as.factor(data$Group)),
  Ntotal= dim(data)[1],
  y_mean= mean(data$GrandTotal),
  y_sd= sd(data$GrandTotal)
)

DSOt= stan_model(model_code= modt)
stanFit_t= sampling(object= DSOt, data= data3, chains=3, iter=1000, warmup=500)

t_chain1= as.array(stanFit_t)[,1,]
t_chain2= as.array(stanFit_t)[,2,]
t_chain3= as.array(stanFit_t)[,3,]
plot(t_chain1[,2],main='mu1',type='l',ylim=c(130,190))
lines(t_chain2[,2],col='red')
lines(t_chain3[,2],col='green')
plot(t_chain1[,3],main='mu2',type='l',ylim=c(130,190))
lines(t_chain2[,3],col='red')
lines(t_chain3[,3],col='green')

plot(density(t_chain1[,2]),main='mu1',xlim=c(130,190))
plot(density(t_chain1[,3]),main='mu2',xlim=c(130,190))

plot(density(t_chain1[,4]),main='sigma1',xlim=c(0,70))
plot(density(t_chain1[,5]),main='sigma2',xlim=c(0,70))

t_diff= t_chain1[,3]-t_chain1[,2]
plot(density(t_diff))

#find HDI
# HDI_fun= function(x,vec) {   #x is lower density boundary
#   mean1_vec= vec
#   q1= quantile(mean1_vec,x)
#   q2= quantile(mean1_vec,0.95+x)
#   width_mean1= q2-q1
#   return(width_mean1)
# }
# 
# nlm(f=HDI_fun,p=c('x'=0.01),vec=t_chain1[,2])

# Posterior predictive
plot(density(t_chain1[,2]+rt(dim(t_chain1)[1],df= t_chain1[,6],ncp= t_chain1[,4])),main='PostPr mu1')
plot(density(t_chain1[,3]+rt(dim(t_chain1)[1],df= t_chain1[,6],ncp= t_chain1[,5])),main='PostPr mu1')
