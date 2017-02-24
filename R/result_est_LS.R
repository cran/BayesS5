result_est_LS <- function(res, X, y,verbose=TRUE){
  hppm = res$hppm
  hppm.prob = res$hppm.prob 
  marg.prob = res$marg.gam
  gam = res$gam
  obj = res$obj 
  post = res$post
  tuning = res$tuning
  p = nrow(gam)
  
  ind.LS = which(gam[,which.max(obj)] == 1)
  p.LS = length(ind.LS)
  if(p.LS>0){
      beta.LS = rep(0,p+1)
      beta.LS[c(1,1+ind.LS)] = stats::lm(y~X[,ind.LS])$coefficients
  }else{
    beta.LS = rep(0,p+1)
    beta.LS[1] = mean(y)
  }
  
  beta.BMA.LS = rep(0,p+1)
  
  for(i in 1:length(post)){
    ind.BMA = which(gam[,i] == 1)
    p.BMA = length(ind.BMA)
    if(p.BMA>0){
      beta.BMA = rep(0,p+1)
      beta.BMA[c(1,1+ind.BMA)] = stats::lm(y~X[,ind.BMA])$coefficients
      beta.BMA.LS = beta.BMA.LS + beta.BMA*post[i]
    }else{
      beta.BMA = rep(0,p+1)
      beta.BMA[1] = mean(y)
      beta.BMA.LS = beta.BMA.LS + beta.BMA*post[i]
    }
    if(verbose==TRUE&&i%%1000==0){
      print(paste("The number of evaluated models: ", i))
    }
  }
  
  return(list(intercept.MAP = beta.LS[1], beta.MAP = beta.LS[-1],intercept.BMA = beta.BMA.LS[1] ,beta.BMA = beta.BMA.LS[-1]))
}
