getHDI= function(chain, HDI_width){
  #function to obtain highest-density interval (i.e. the region that covers HDI_width (e.g. 95%)
    #of indices where interval of the values at these indices is smallest)
  
  chain_sort= sort(chain)
  idxAsFrac= 1:length(chain_sort)/length(chain_sort)  #put indices of sorted chain on 0-1 scale
  upperLimIdx_lowTail=  max(which(idxAsFrac<(1-HDI_width)))  #max index where lower tail can end
  #so that HDI can still be accommodated
  idxIntervalMtx= matrix(NA,nrow=upperLimIdx_lowTail,ncol=2) #stores lower and upper interval
  #limits; row number is lower idx
  for(loIdx in 1:upperLimIdx_lowTail){ #move through allowed lower indices
    upperLimIdx= max(which( idxAsFrac < idxAsFrac[loIdx]+HDI_width))
    intervalWidth= chain_sort[upperLimIdx]-chain_sort[loIdx]
    idxIntervalMtx[loIdx,]= c(chain_sort[loIdx],chain_sort[upperLimIdx])
  }
  
  #find idx for narrowest interval
  intervalWidths= idxIntervalMtx[,2]-idxIntervalMtx[,1]
  HDI_loIdx= which( intervalWidths==min(intervalWidths))[1]  #take first in case of tie
  HDI= idxIntervalMtx[HDI_loIdx,] 
  
  return(HDI)
}