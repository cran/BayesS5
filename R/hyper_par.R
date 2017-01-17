hyper_par <-
function(type,X,y,thre){
  n =nrow(X)
  p =ncol(X)
  if(missing(thre)){thre = p^-0.5}
  
  if(type=="pimom"){
    
    betas = matrix(0,3,50000)
    for(k in 1:50000){
      sam = sample(1:p,3)
      ind = sample(1:n,n)
      betas[,k] = as.vector(solve(crossprod(X[ind,sam]))%*%crossprod(X[ind,sam],y))
    }
    res=y
    corr = as.vector(cor(res,X))
    ind.ix = sort.int(abs(corr),decreasing=TRUE,index.return=TRUE)$ix
    s = ind.ix[1:3]
    #p = p+1
    beta.hat =solve(crossprod(X[,s]))%*%crossprod(X[,s],y)
    sig.hat = crossprod(y - X[,s]%*%beta.hat)/n
    betas=as.vector(betas)
    tau.cand = seq(0.1,(sd(y)+0.1),length.out=300)^2
    pro = rep(0,300)
    for(k in 1:300){
      tau = tau.cand[k]
      den = function(x){tau^0.5*x^-2*exp(-1*tau/(x^2) )/gamma(1/2)}
      den.null1 = density(betas)
      data = list(x=den.null1$x,y=den.null1$y)
      den.null = approxfun(data[[1]], data[[2]], rule=1,method = "linear")
      f = function(x){den(x) - den.null(x)}
      tryCatch({
        th=1
        
        a = uniroot(f,interval = c(0.001,max(betas))) 
        th = a$root
        loc = integrate(den.null,lower = th, upper =max(betas)-0.001)$value
        nonloc =  integrate(den,lower = 0, upper = th)$value
        pro[k] = loc + nonloc}, error=function(e){})
    }
    
    tau=1
    B = tau.cand[which.min((pro-thre)^2)]
  }
  
  if(type=="pemom"){
    
    betas = matrix(0,3,50000)
    for(k in 1:50000){
      sam = sample(1:p,3)
      ind = sample(1:n,n)
      betas[,k] = as.vector(solve(crossprod(X[ind,sam]))%*%crossprod(X[ind,sam],y))
    }
    
    res=y
    corr = as.vector(cor(res,X))
    ind.ix = sort.int(abs(corr),decreasing=TRUE,index.return=TRUE)$ix
    s = ind.ix[1:3]
    #p = p+1
    beta.hat =solve(crossprod(X[,s]))%*%crossprod(X[,s],y)
    sig.hat = crossprod(y - X[,s]%*%beta.hat)/n
    betas=as.vector(betas)
    tau.cand = seq(0.1,(sd(y)+0.1),length.out=300)^2
    pro = rep(1,300)
    for(k in 1:300){
      tau = tau.cand[k]
      den = function(x){ sqrt(2*pi*sig.hat)^-1*exp( -1*tau/(x^2) - x^2/(sig.hat*tau) + sqrt(2/sig.hat) ) }
      den.null1 = density(betas)
      data = list(x=den.null1$x,y=den.null1$y)
      den.null = approxfun(data[[1]], data[[2]], rule=1,method = "linear")
      #curve(den.null,-5,5) 
      #curve(den,add=T,col="red") 
      f = function(x){den(x) - den.null(x)}
      th=1
      tryCatch({
        a = uniroot(f,interval = c(0.001,max(betas))) 
        th = a$root
        loc = integrate(den.null,lower = th, upper =max(betas)-0.001)$value
        nonloc =  integrate(den,lower = 0, upper = th)$value
        pro[k] = loc + nonloc},  error=function(e){})
    }
    
    B = tau.cand[which.min(abs(pro-thre))]
  }
  
  return(B)
}
