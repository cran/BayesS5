
ind_fun_pemom <- function(X.ind,y,n,p,tuning){ 
  B = tuning
  obj01 = function(x2,B,y,X.ind,a0,b0){
    a0=0.01;b0=0.01
    p.g = ncol(X.ind)
    if(length(p.g)==0){p.g=0}
    tau=1
    if(p.g>1){
      beta = x2[1:p.g] 
      prec = x2[(1+p.g)]
      if(prec<0){prec=10^-5}
      a = 0.5*prec*crossprod(y-X.ind%*%beta)+0.5*prec*tau*crossprod(beta) 
      b = sum(B/beta^2)- (a0-1)*log(prec) + b0*prec - p.g*sqrt(2*prec*tau*B)-0.5*(n+p.g)*log(prec)
      c = a+b
    }else{if(p.g==1){
      beta = x2[1:p.g] 
      prec = x2[(1+p.g)]
      if(prec<0){prec=10^-5}
      a = 0.5*prec*crossprod(y-X.ind*beta)+0.5*prec*tau*beta^2 
      b = B/beta^2 - (a0-1)*log(prec) + b0*prec - p.g*sqrt(2*prec*tau*B)-0.5*(n+p.g)*log(prec)
      c = a+b
    }else{
      prec = x2[(1+p.g)]
      if(prec<0){prec=10^-5}
      a = 0.5*prec*crossprod(y) 
      b = -1*(a0-1)*log(prec) + b0*prec -0.5*n*log(prec)
      c = a+b                       
    }
      
    }
    return(c)
  }
  
  J = function(x1,B,y,X.ind,a0,b0){
    a0=0.01;b0=0.01
    p.g1 = ncol(X.ind)
    if(length(p.g1)==0){p.g1=0}
    tau=1
    if(p.g1>1){
      beta1 = x1[1:p.g1]
      prec1 = x1[(p.g1+1)]
      D = matrix(0,(p.g1+1),(p.g1+1))
      D[1:p.g1,1:p.g1] =prec1*crossprod(X.ind)+prec1*tau*diag(p.g1)+diag(6*B/beta1^4)
      
      FF = as.vector(crossprod(X.ind,y-X.ind%*%beta1)-tau*beta1)*prec1^2
      D[(p.g1+1),1:p.g1]=FF
      D[1:p.g1,(p.g1+1)]=FF
      D[(p.g1+1),(p.g1+1)] = -1*(n/2+p.g1/2+a0+1)*prec1^2+prec1^3*(crossprod(y-X.ind%*%beta1)-tau*crossprod(beta1))-0.75*p.g1*sqrt(2*B*tau)*prec1^2.5+prec1^3*2*b0
      #0.5*(n+p.g1-2+2*a0)/prec1^2+p.g1*0.25*sqrt(2*tau*B)*prec1^(-1.5)
      
    }else{
      if(p.g1==1){
        beta1 = x1[1:p.g1]
        prec1 = x1[(p.g1+1)]
        D = matrix(0,(p.g1+1),(p.g1+1))
        D[1:p.g1,1:p.g1] =prec1*crossprod(X.ind)+prec1*tau*diag(p.g1)+6*B/beta1^4
        FF = as.vector(crossprod(X.ind,y-X.ind*beta1)-tau*beta1)*prec1^2
        D[(p.g1+1),1:p.g1]=FF
        D[1:p.g1,(p.g1+1)]=FF
        D[(p.g1+1),(p.g1+1)] = -1*(n/2+p.g1+a0+1)*prec1^2+prec1^3*(crossprod(y-X.ind*beta1)-tau*beta1^2)-0.75*p.g1*sqrt(2*B*tau)*prec1^2.5+prec1^3*2*b0
        
        #0.5*(n+p.g1-2+2*a0)/prec1^2+p.g1*0.25*sqrt(2*tau*B)*prec1^(-1.5)
      }else{
        prec1 = x1[(p.g1+1)]
        D =abs(-0.5*(n+2*a0+2)*prec1^2+prec1^3*2*b0+prec1^3*crossprod(y))
      }
    }
    
    return(D)
  }
  
  
  
  r0=1
  B = tuning
  a0=0.01;b0=0.01
  tau = 1
  p.g = ncol(X.ind)
  if(length(p.g)==0){p.g=0}
  if(p.g >1){
    fit = solve(crossprod(X.ind)+diag(p.g))%*%crossprod(X.ind,y)
    ress = crossprod(y-X.ind%*%fit)
    prec0 = n/ress
    initial_x = fit
    initial_x = c(initial_x,prec0)
    wrapper <- function(theta) obj01(theta,B,y,X.ind,a0,b0)
    #o <- optim(initial_x, wrapper) 
    #f_x = o$value
    #x0 = o$par
    o = nlm(wrapper, initial_x)
    f_x = o$minimum  
    x0 = o$estimate
    ccc = 0.5*log(2*pi)+0.5*p.g*log(tau)
    #int = ccc-0.5*as.numeric(determinant(J(x0,ind2,B,XtX,y,X,tau,a0,b0),logarithm=TRUE)$modulus)-f_x
    int = ccc-0.5*log(det(J(x0,B,y,X.ind,a0,b0)))-f_x
    
  }else{
    if(p.g==1){
      fit = solve(crossprod(X.ind)+diag(p.g))%*%crossprod(X.ind,y)
      ress = crossprod(y-X.ind%*%fit)
      prec0 = n/ress
      initial_x = fit
      initial_x = c(initial_x,prec0)
      wrapper <- function(theta) obj01(theta,B,y,X.ind,a0,b0)
      #o <- optim(initial_x, wrapper) 
      #f_x = o$value
      #x0 = o$par
      o = nlm(wrapper, initial_x)
      f_x = o$minimum  
      x0 = o$estimate
      ccc = 0.5*log(2*pi)+0.5*p.g*log(tau)
      #int = ccc-0.5*as.numeric(determinant(J(x0,ind2,B,XtX,y,X,tau,a0,b0),logarithm=TRUE)$modulus)-f_x
      int = ccc-0.5*log(det(J(x0,B,y,X.ind,a0,b0)))-f_x
      
    }else{if(p.g==0){
      int = lgamma(0.5*n+a0) -(0.5*n+a0)*log(crossprod(y)+2*b0)
    }
    }
  }
  return(int)
}
