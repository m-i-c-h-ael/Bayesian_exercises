# 16.2. Caloric restriction extends the lifespan of rats; but does it also increase the
 #variability in lifespan?

#Two-group analysis describing the two groups with student t-distributions with
 #individual shape parameters, but common normality parameter

rm(list=ls())
cat('\014')
graphics.off()

library('rstan')

#specify model
stanCode= '
data {
  int<lower=1> Ntotal;       //number of rats
  real<lower=0> y[Ntotal];  //rat lifespan
  int x[Ntotal];            //category: ad libidum or CR
  real<lower=0> y_mean;
  real<lower=0> y_sd;
}
transformed data{
  real<lower=0> sd_norm;
  real<lower=0> unif_lo;
  real<lower=1> unif_hi;
  sd_norm= y_sd*100;       //sd for prior on mu  
  unif_lo= y_sd/1000;
  unif_hi= y_sd*1000;
}
parameters{
  real<lower=0> nuMin1;
  real<lower=0> mu[2];    //mean for both groups
  real<lower=0> sigma[2];
}
transformed parameters{
  real<lower=0> lambda;
  real<lower=1> nu;
  lambda= 1/29.0;
  nu= nuMin1+1;
}
model{
  nuMin1 ~ exponential(lambda);
  mu ~ normal(y_mean, sd_norm);
  sigma ~ uniform(unif_lo, unif_hi);
  for (i in 1:Ntotal) {
    y[i] ~ student_t(nu, mu[x[i]], sigma[x[i]]);
  }
}
'
# prepare data
data= read.csv('H:/My Drive/learning/20230220_BayesianStatistics/20230504_BayesianDogsBook/DBDA2Eprograms/RatLives.csv')
head(data)
dataList= list(
  x= as.numeric(as.factor(data$Group)),
  y= data$DaysLive,
  Ntotal= dim(data)[1],
  y_mean= mean(data$DaysLive),
  y_sd= sd(data$DaysLive)
)

#run MCMC
DSO= stan_model(model_code= stanCode)
stanOUT= sampling(DSO, data=dataList, chains=3, iter=1000, warmup=500)

# isolate chains
chain1= as.array(stanOUT)[,1,]
chain2= as.array(stanOUT)[,2,]
chain3= as.array(stanOUT)[,3,]
head(chain1)

#check runs
plot(chain1[,2],main='mu1',type='l')
lines(chain2[,2],col='red')
