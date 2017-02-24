obj_fun_g <- function(ind,X,y,n,p,tuning){ 
  p.g = length(ind)
  X0 = cbind(rep(1,n),X[,ind])
  beta= rep(0,p+1)
  if(p.g >0){
    #X0 = X[,ind2];QR=qr(X0)
    #ress = crossprod(qr.resid(QR, y))
    
    beta[c(1,1+ind)] = tuning*solve(crossprod(X0))%*%crossprod(X0,y)/(1+tuning)
    ress = crossprod(y-X0%*%beta[c(1,1+ind)])
    sig = ress/(n + 2)
    
  }else{
    beta[1] = mean(y)
    sig = crossprod(y-mean(y))/(n + 2)
    }
  
  return(list(beta=beta,sig=sig))
}
