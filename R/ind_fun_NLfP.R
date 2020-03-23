ind_fun_NLfP =function(ind2, y, phi, n, p, K, IP.phi, C.prior1, tuning){
  tau = g = tuning
  a0 = b0 = 1
#  assign("g", g, .GlobalEnv)
#  assign("tau", tau, .GlobalEnv)
  #assign("a0", a0, .GlobalEnv)
  #assign("b0", b0, .GlobalEnv)
  
  index = function(j){
    a = (K*(j-1)+2):(K*j+1)
    return(a)
  }
  
  index.tot = function(ind2){
    #a = rep(0,K*p+1)
    ind = sapply(ind2,index)#;a[ind] = 1
    return(as.vector(ind))
  }
  
  nonlocal = function(j,beta){
    ind = index(ind2[j])
    return(1/crossprod(IP.phi[,ind]%*%beta[(K*(j-1)+1):(K*j)]))
  }

  obj01 = function(x2,ind2,tau,g,a0,b0){
    p.g = length(ind2)
    if(p.g>0){
      beta0 = x2[1]
      beta = x2[2:(p.g*K+1)] 
      prec = x2[p.g*K+2]
      ind3 = index.tot(ind2)
      if(prec<0){prec=10^-5}
      
      a = 0.5*prec*crossprod(y-phi[,c(1,ind3)]%*%c(beta0,beta)) + 0.5*prec*crossprod(beta)/g + tau*sum(sapply(1:p.g,nonlocal,beta = beta))/prec*n
      b = -1*(a0-1)*log(prec) + b0*prec -0.5*(n+K*p.g)*log(prec) 
      c = a+b
    }else{
      prec = x2[(1+p.g)]
      if(prec<0){prec=10^-5}
      a = 0.5*prec*crossprod(y) 
      b = -1*(a0-1)*log(prec) + b0*prec -0.5*n*log(prec)
      c = a+b                       
    }
    return(c)
  }
  
  
  J = function(x0,ind2,tau,g,a0,b0){
    p.g = length(ind2) 
    beta0 = x0[1]
    beta = x0[2:(p.g*K+1)] 
    prec = x0[p.g*K+2]
    ind3 = index.tot(ind2)
    
    if(p.g>0){
      a22 = matrix(0,p.g*K,p.g*K)
      a23 = rep(0,p.g*K)
      a33 = 0
      for(k in 1:p.g){
        j = ind2[k]
        ind = index.tot(j)
        ind1 = (K*(k-1)+1):(K*k)
        phi.beta = IP.phi[,ind]%*%beta[ind1]
        a22.1 = 8*tau*t(phi[,ind])%*%tcrossprod(phi.beta)%*%(phi[,ind])*as.numeric(crossprod(phi.beta))^-3/prec
        a22.2 = -2*tau*crossprod(IP.phi[,ind])*as.numeric(crossprod(phi.beta))^-2/prec
        a22[ind1,ind1] = a22.1 + a22.2
        a23[ind1] = 2*tau*(crossprod(IP.phi[,ind])*as.numeric(crossprod(phi.beta))^-1)%*%beta[ind1]/(prec^2)
        a33 = a33 +  2*tau*as.numeric(crossprod(phi.beta))^-1/(prec^3)
      }
      u = matrix(0,K*p.g+2,K*p.g+2)
      b11 =  prec*n
      b12 = prec*as.vector(prec*crossprod(phi[,ind3],rep(1,n)))
      b13 = -1*sum(y-phi[,c(1,ind3)]%*%c(beta0,beta))
      b22 = prec*(crossprod(phi[,ind3]) + diag(K*p.g)/g) + a22
      b23 = -1*crossprod(phi[,ind3],y-phi[,c(1,ind3)]%*%c(beta0,beta)) + beta/g + a23 
      b33 = (n/2+p.g*K/2+a0-1)*prec^-2 + a33
      
      u[1,1] = b11
      u[1,2:(p.g*K+1)] = u[2:(p.g*K+1),1] = b12
      u[1,(p.g*K+2)] = u[(p.g*K+2),1] = b13
      u[2:(p.g*K+1),(p.g*K+2)] = u[(p.g*K+2),2:(p.g*K+1)] = b23
      u[2:(p.g*K+1),2:(p.g*K+1)] = b22
      u[(p.g*K+2),(p.g*K+2)] = b33
    }else{
      u = diag(2)
    }
    return(u)
  } 
  
  
  ind_fun1=function(ind2,tau,g,a0,b0){ 
    #a0=0.01;b0=0.01
    p.g=length(ind2)
    ind3 = index.tot(ind2)
    if(p.g >0){
      fit = solve(crossprod(phi[,c(1,ind3)])+0.001*diag(p.g*K+1))%*%crossprod(phi[,c(1,ind3)],y)
      ress = crossprod(y-phi[,c(1,ind3)]%*%fit)
      prec0 = n/ress
      initial_x = c(fit,prec0)
      wrapper <- function(theta) obj01(theta,ind2,tau,g,a0,b0)
      #o <- optim(initial_x, wrapper,hessian=TRUE) 
      o <- optim(initial_x, wrapper) 
      f_x = o$value
      x0 = o$par
      
      Hess = J(x0,ind2,tau,g,a0,b0) 
      det.J = 10^100
      tryCatch({
        det.J = determinant(Hess,logarithm = TRUE)$modulus
      },error=function(e){})
      #; if(det.J<0){det.J = 100^100}
      int = -1*p.g*C.prior1 - 0.5*det.J - f_x - 0.5*K*p.g*log(g)
    }else{
      int = -(0.5*n+a0)*log(crossprod(y-mean(y))/2+b0)
    }
    
    return(as.numeric(int))
  }
  
  return(ind_fun1(ind2,tau,g,a0,b0))
}
