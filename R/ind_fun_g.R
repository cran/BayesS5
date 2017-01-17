ind_fun_g <-
function(X.ind,y,n,p,tuning){ 
  g= tuning
  a0=0.01;b0=0.01
  tau = 1
  p.g = ncol(X.ind)
  if(length(p.g)==0){p.g=0}
  if(p.g >1){
    #X0 = X[,ind2];QR=qr(X0)
    #ress = crossprod(qr.resid(QR, y))
    fit = solve(crossprod(X.ind))%*%crossprod(X.ind,y)
    ress = crossprod(y-X.ind%*%fit)
    v = crossprod(y-mean(y))
    int = -0.5*p.g*log(1+g)-0.5*(n-1)*log(1+g*(ress/v))
    
  }else{
    if(p.g==1){
      fit = as.numeric(as.numeric(sum(X.ind*y))/as.numeric(crossprod(X.ind)))
      ress = crossprod(y-X.ind*fit)
      v = crossprod(y-mean(y))
      int = -0.5*p.g*log(1+g)-0.5*(n-1)*log(1+g*(ress/v))
      
    }else{if(p.g==0){
      int = -0.5*(n-1)*log(1+g)
    }
    }
  }
  return(int)
}
