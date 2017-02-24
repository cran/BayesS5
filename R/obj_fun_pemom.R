obj_fun_pemom <- function(ind,X,y,n,p,tuning){ 
  p.g = length(ind)
  X0 = cbind(rep(1,n),X[,ind])
  
  obj01_pimom <-
    function(x2,B,y,X0,a0,b0){
      a0=0.01;b0=0.01;r0=1
      tau = 1
      
      p.g = length(x2)-1
      
      if(p.g>1){
        beta = x2[2:p.g] 
        beta0 = x2[1]
        prec = x2[(1+p.g)]
        if(prec<0){prec=10^-5}
        a = 0.5*prec*crossprod(y-X0%*%c(beta0,beta) ) + 0.5*prec*tau*crossprod(beta) 
        b = sum(B/beta^2)- (a0-1)*log(prec) + b0*prec - (p.g-1)*sqrt(2*prec*tau*B)-0.5*(n+p.g-1)*log(prec)
        c = a+b
      }else{
        prec = x2[(1+p.g)]
        if(prec<0){prec=10^-5}
        a = 0.5*prec*crossprod(y-mean(y)) 
        b = -1*(a0-1)*log(prec) + b0*prec -0.5*n*log(prec)
        c = a+b                       
      }
      return(c)
    }
  
  r0=1
  B = tuning
  a0=0.01;b0=0.01
  tau = 1
  if(p.g >0){
    fit = solve(crossprod(X0)+diag(p.g+1)/B)%*%crossprod(X0,y)
    ress = crossprod(y-X0%*%fit)
    prec0 = (n+p.g+2*a0)/(ress+2*b0)
    initial_x = c(fit,prec0)
    wrapper <- function(theta) obj01_pimom(theta,tuning,y,X0,a0,b0)
    o = nlm(wrapper, initial_x)
    f_x = o$minimum  
    x0 = o$estimate
    beta = rep(0,p+1)
    beta[c(1,1+ind)] = x0[1:(p.g+1)]
    sig = x0[p.g+2]
  }else{
    beta = rep(0,p+1)
    beta[1] = mean(y)
    sig = crossprod(y-mean(y))/(n + a0)
    }
    
  return(list(beta=beta,sig=sig))
}
