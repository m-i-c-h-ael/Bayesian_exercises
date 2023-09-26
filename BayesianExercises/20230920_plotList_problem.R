co=1:5

l= list()

for(i in 1:length(co)){
  l[[i]]= co[i]
}


#######################
library('patchwork')
rm(list=ls())

GaussDF= cbind.data.frame(mu=c(4,1,4,4),sig=c(1.5,2.5,1,3))
thetasDF= rbind.data.frame(1:5,2:6,3:7,4:8)
centersDF= rbind.data.frame(seq(1.5,4.5,1),seq(2.5,5.5,1),seq(3.5,6.5,1),seq(4.5,7.5,1))

p= list()
for(i in 1:dim(GaussDF)[1]){
  print(paste(i))
  
  thetas= as.numeric(thetasDF[i,])
  mu= GaussDF$mu[i]; sig= GaussDF$sig[i]
  
  bars= rollapply(thetas,FUN= function(x){pnorm(mean=mu,sd=sig,x[2])-pnorm(mean=mu,sd=sig,x[1])},
                  width=2,align='right')

  p[[i]] <- ggplot()+
    geom_bar(aes(x=as.numeric(centersDF[i,]),y=bars),stat='identity',fill='skyblue')+
    geom_function(fun= dnorm,args=list(mean=mu,sd=sig),xlim=c(-1,10))+
    coord_cartesian(ylim=c(0,0.4))
}

p[[1]] / p[[2]] / p[[3]] / p[[4]]


