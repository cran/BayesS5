result_est_MAP <- function(res, X, y, obj_fun,verbose=TRUE){
  hppm = res$hppm
  hppm.prob = res$hppm.prob 
  marg.prob = res$marg.gam
  gam = res$gam
  obj = res$obj 
  post = res$post
  tuning = res$tuning
  p = nrow(gam)
  n = nrow(X)
  if(missing(obj_fun)){
    print("The prior on regression coefficients is unspecified. The default is piMoM")
    ind_fun = BayesS5::obj_fun_pimom
    tuning <- BayesS5::hyper_par(type="pimom",X,y,thre = p^-0.5)  # tuning parameter selection for nonlocal priors
    print("The choosen hyperparameter tau")
    print(tuning)
    #assign("tuning", tuning, .GlobalEnv)
  }
  
  ind.MAP = which(gam[,which.max(obj)] == 1)
  o = obj_fun(ind.MAP,X=X,y=y,n=n,p=p,tuning=tuning)
  beta.MAP = o$beta
  sig.MAP = o$sig
  
  beta.BMA.MAP = rep(0,p+1)
  
  for(i in 1:length(post)){
    ind.BMA = which(gam[,i] == 1)
    p.BMA = length(ind.BMA)
    if(p.BMA>0){
      o = obj_fun(ind.BMA,X=X,y=y,n=n,p=p,tuning=tuning)
      beta.BMA = o$beta
      beta.BMA.MAP = beta.BMA.MAP + beta.BMA*post[i]
    }else{
      beta.BMA = rep(0,p+1)
      beta.BMA[1] = mean(y)
      beta.BMA.MAP = beta.BMA.MAP + beta.BMA*post[i]
    }
    if(verbose==TRUE&&i%%500==0){
      print(paste("The number of evaluated models: ", i))
    }
  }
  
  return(list(intercept.MAP = beta.MAP[1], beta.MAP = beta.MAP[-1],sig.MAP = sig.MAP, intercept.BMA = beta.BMA.MAP[1] ,beta.BMA = beta.BMA.MAP[-1]))
}
