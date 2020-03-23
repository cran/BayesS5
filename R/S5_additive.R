S5_additive <- function(X,y, K = 5, model, tuning = 10,  tem, ITER=20,S=20,C0=5, verbose=TRUE){
  requireNamespace("splines2")
  requireNamespace("Matrix")
  n = nrow(X)
  p = ncol(X)
  y = y -mean(y)

  tau = g = tuning
  Matrix = Matrix::Matrix
  if(missing(tem)){tem = seq(0.4,1,length.out=20)^2}
  #assign("K", K, .GlobalEnv)
  #assign("n", n, .GlobalEnv)
  #assign("p", p, .GlobalEnv)
  #assign("tau", tau, .GlobalEnv)
  #assign("g", g, .GlobalEnv)
  
  ##################################################
  ind_fun = BayesS5::ind_fun_NLfP
  #ind_fun = ind_fun_NLfP
  
  index = function(j){
    a = (K*(j-1)+2):(K*j+1)
    return(a)
  }  
  #assign("index", index, .GlobalEnv)
  index.tot = function(ind2){
    ind = sapply(ind2,index)#;a[ind] = 1
    return(as.vector(ind))
  }
  
  screening = function(j,phi,res){
    ind3 = index.tot(j)
    fit = solve(crossprod(phi[,c(1,ind3)])+0.0001*diag(K+1))%*%crossprod(phi[,c(1,ind3)],res)
    fit.f = phi[,c(1,ind3)]%*%fit
    a = crossprod(fit.f - mean(fit.f))
    return(a)
  }
  ###########################################################
  
  if(missing(model)){
    print("The model prior is unspecified. The default is Bernoulli_Uniform")
    model = BayesS5::Bernoulli_Uniform
  }
  A3 = S;r0=1
  verb = verbose
  P0 = tcrossprod(rep(1,n))/n
  phi0 = matrix(0,n,K*p)
  Knots = matrix(0,p,2)
  colnames(Knots) = c("Lower","Upper")
  for(j in 1:p){
    Knots[j, ] = c(min(X[,j])-1.0,max(X[,j])+1.0)
    phi0[,(K*(j-1)+1):(K*j)] = splines2::bSpline(X[,j], df = K, Boundary.knots = Knots[j, ])
  }
  phi = cbind(rep(1,n),phi0)
  IP = diag(n) - tcrossprod(rep(1,n))/n
  IP.phi = IP%*%phi
  
  #assign("IP", IP, .GlobalEnv)
  #assign("IP.phi", IP.phi, .GlobalEnv)
  
  ind2 = sample(1:p,2)
  #ind2 = true
  gam = rep(0,p);
  gam[ind2]=1
  ind2 = which(gam==1)
  GAM.screen = Matrix(0,p,50000,sparse=TRUE)
  ID.screen = rep(-100000000,50000)
  
  save = rep(0,p)
  p.g = length(ind2)
  ind3 = index.tot(ind2)
  
  if(p.g>0){
    fit = solve(crossprod(phi[,c(1,ind3)])+0.001*diag(p.g*K+1))%*%crossprod(phi[,c(1,ind3)],y)
    res = y-phi[,c(1,ind3)]%*%fit}else{res=y}
  save = sapply(1:p,screening,phi,res)
  ind.ix = sort.int(save,decreasing=TRUE,index.return=TRUE)$ix
  corr = as.vector(cor(res,X))
  ind.l = sort.int(abs(corr),decreasing=TRUE,index.return=TRUE)$ix
  
  IND = c(ind2,union(ind.ix[1:S],ind.l[1:5]))
  IND = unique(IND)
  p.ind = length(IND)
  
  
  ID.screen[1] = sum(5^(log(ind2)))
  GAM.screen[IND,1] = 1
  #####
  
  j = 1; NNN = 10000
  kk = stats::rchisq(NNN,K-1)
  aa = log(mean(exp(-1/kk)))
  C.prior3 = rep(0,p)
  C.g2 = rep(0,p)
  for(j in 1:p){
    C.g2[j] = -0.5*log(det(crossprod(phi[,(K*(j-1)+2):(K*j+1)])))
    C.prior3[j] = C.g2[j] + aa
  }
   # assign("C.prior3", C.prior3, .GlobalEnv)
  #  assign("C.g2", C.prior3, .GlobalEnv)
  
  
    #assign("tau", tau, .GlobalEnv)
    #assign("g", g, .GlobalEnv)
    aa = 0; j = 1; NNN = 10000
    for(h in 1:NNN){
      kk = stats::rnorm(K)*sqrt(g)
      aa = aa + exp(-tau*n/crossprod(IP.phi[,(K*(j-1)+2):(K*j+1)]%*%kk))
    }
    
    
    C.prior1 = log(aa/NNN)
    C.prior1 = as.numeric(C.prior1)
   # assign("C.prior1", C.prior1, .GlobalEnv)
    
    aa = 0; j = 1; NNN = 10000
    for(h in 1:NNN){
      kk = stats::rcauchy(K)
      aa = aa + exp(-tau*n/crossprod(IP.phi[,(K*(j-1)+2):(K*j+1)]%*%kk))
    }
    
    
    C.prior2 = log(aa/NNN)
    C.prior2 = as.numeric(C.prior2)
    #assign("C.prior2", C.prior2, .GlobalEnv)
    
    
    pmt = proc.time()
    print("#################################")
    print("S5 starts running")
    
    IT = length(tem)
    IT.seq = rep(ITER,IT)
    
    curr = -1000000000
    tryCatch({
      curr  = ind_fun(ind2, y, phi, n, p, K, IP.phi, C.prior1, tuning) + model(ind2, p )    
    },error=function(e){})
    p.g=sum(gam)
    GAM.fin0 = NULL
    OBJ.fin0 = NULL
    
    for(uu in 1:C0){
      
      C.p = rep(-1000000000,p)
      C.m = rep(-1000000000,p)

      GAM = gam
      OBJ = curr
      obj = OBJ
      p.g=sum(gam)
      
      C.p = rep(-100000000,p)
      for(i in (p.g+1):p.ind){
        j=IND[i]
        gam.p = gam;gam.p[j]=1;ind.p=which(gam.p==1)
        int = -10000000
        int  = ind_fun(ind.p, y, phi, n, p, K, IP.phi, C.prior1, tuning) + model(ind.p, p )    
        obj.p =  c(int)
        if(is.na(obj.p)==TRUE){obj.p = -100000000}
        C.p[j] = obj.p
      }
      C.m = rep(-100000000,p)
      IND.m = ind2
      p.ind.m = length(IND.m)
      for(i in 1:p.g){
        j=ind2[i]
        gam.m = gam;gam.m[j]=0;ind.m=which(gam.m==1)     
        int = -10000000
        int  = ind_fun(ind.m, y, phi, n, p, K, IP.phi, C.prior1, tuning) + model(ind.m, p)    
        obj.m =  c(int)
        if(is.na(obj.m)==TRUE){obj.m = -100000000}
        C.m[j] = obj.m  
      }
      
      p.g = sum(gam)
      
      OBJ.m0 = matrix(C.m,p,1)
      OBJ.p0 = matrix(C.p,p,1)
      ID = sum(5^(log(ind2)))
      ID.obj = ID
      it=1
      #GAM.total = matrix(0,p,50000)
      GAM.total = Matrix(0,p,50000,sparse=TRUE)
      OBJ.total = rep(-100000000,50000)
      GAM.total[,1] = gam
      OBJ.total[1] = obj
      time.total = rep(0,50000)
      it=1
      INT = NULL
      
      pmt0 = proc.time()
      for(it in 1:IT){   
        IT0 = IT.seq[it] 
        pq=0
        for(iter in 1:IT0){
          
          id = sum(5^(log(ind2)))
          id.ind = which(id==ID)
          leng = length(id.ind)
          
          if(leng==0){
            ID = c(ID,id)
            C.p = rep(-100000000,p)
            for(i in (p.g+1):p.ind){
              j=IND[i]
              gam.p = gam;gam.p[j]=1;ind.p=which(gam.p==1)
              int = -10000000
              #tryCatch({
              int  = ind_fun(ind.p,y, phi, n, p, K, IP.phi, C.prior1, tuning) + model(ind.p, p)    
              #},error=function(e){})
              obj.p =  c(int)
              if(is.na(obj.p)==TRUE){obj.p = -100000000}
              C.p[j] = obj.p
              ind.total = which(OBJ.total< -90000000)[1]
              OBJ.total[ind.total] = obj.p
              GAM.total[,ind.total] =  gam.p
              time.total[ind.total] = (proc.time()-pmt0)[3]
            }
            p.g = sum(gam)
            C.m = rep(-100000000,p)
            IND.m = ind2
            p.ind.m = length(IND.m)
            for(i in 1:p.g){
              j=ind2[i]
              gam.m = gam;gam.m[j]=0;ind.m=which(gam.m==1)     
              int = -10000000
              #tryCatch({
              int  = ind_fun(ind.m, y, phi, n, p, K, IP.phi, C.prior1, tuning) + model(ind.m,p)    
              #},error=function(e){})
              obj.m =  c(int)
              if(is.na(obj.m)==TRUE){obj.m = -100000000}
              C.m[j] = obj.m  
              ind.total = which(OBJ.total< -90000000)[1]
              OBJ.total[ind.total] = obj.m
              GAM.total[,ind.total] =  gam.m
              time.total[ind.total] = (proc.time()-pmt0)[3]
            }
            
            OBJ.p0 = cbind(OBJ.p0,C.p)
            OBJ.m0 = cbind(OBJ.m0,C.m)  
            
          }else{
            pq= pq+1
            C.p = OBJ.p0[,(id.ind[1])];C.m = OBJ.m0[,(id.ind[1])]
          }
          
          prop = exp(tem[it]*(C.p-max(C.p)))
          sample.p = sample(1:length(prop),1,prob=prop)
          obj.p = C.p[sample.p]
          #obj.p
          prop = exp(tem[it]*(C.m-max(C.m)))
          sample.m = sample(1:length(prop),1,prob=prop)
          obj.m = C.m[sample.m]
          #obj.m
          
          l = 1/(1+exp(tem[it]*obj.m-tem[it]*obj.p))
          if(l>runif(1)){ gam[sample.p]=1;obj = obj.p;curr=obj.p
          }else{
            gam[sample.m]=0;obj = obj.m;curr=obj.m
          } 
          ind2 = which(gam==1)
          p.g = sum(gam)
          #int = -100000000
          #tryCatch({
          #  int  = ind_fun(ind2) + model(ind2)    
          #},error=function(e){})
          #curr = int
          #jjj = sample(1:3,1)
          #if(jjj==1){
          #pmt0 = proc.time()
          id = sum(5^(log(ind2)))
          id.ind = which(id==ID.screen)
          leng = length(id.ind)
          if(leng==0){
            jjj = sample(1:2,1)
            if(jjj==1){
              save = rep(0,p)
              ind3 = index.tot(ind2)
              if(p.g>0){
                fit = solve(crossprod(phi[,c(1,ind3)])+0.01*diag(p.g*K+1))%*%crossprod(phi[,c(1,ind3)],y)
                res = y-phi[,c(1,ind3)]%*%fit}else{res=y}
              save = sapply(1:p,screening,phi,res)
              ind.ix = sort.int(save,decreasing=TRUE,index.return=TRUE)$ix
              
              corr = as.vector(cor(res,X))
              ind.l = sort.int(abs(corr),decreasing=TRUE,index.return=TRUE)$ix
              
              IND = c(ind2,union(ind.ix[1:S],ind.l[1:5]))
              p.ind = length(IND)
              ind.id = which(ID.screen< 0)[1]
              ID.screen[ind.id] = id
              GAM.screen[IND,ind.id] =  1}
          }else{
            IND = which(GAM.screen[,id.ind[1]]==1)
            p.ind = length(IND)
          }
          #print(proc.time()-pmt0)
          #}
          id = sum(5^(log(ind2)))
          id.ind = which(id==ID.obj)
          
          leng = length(id.ind)
          if(leng==0){
            ID.obj = c(ID.obj,id)
            OBJ = c(OBJ,curr)
            GAM= cbind(GAM,gam)
          }
        }
        
        if(verbose==TRUE){
          print("#################################")
          gam.pr = GAM.total[,which.max(OBJ.total)]
          obj.pr = max(OBJ.total)
          ind2.pr = which(gam.pr==1)
          print("Inverse Temperature");print(tem[it]);print("The Selected Variables in the Searched MAP Model");
          print(ind2.pr);print("The Evaluated Object Value at the Searched MAP Model");print(obj.pr);
          print("Current Model");print(ind2);  
          print("The Evaluated Object Value at the Current Model");print(curr);
          print("Total Searched Variables");
          print(IND)
          print("The Number of Total Searched Models");
          print(length(unique(OBJ.total)))
          #print(length(which(OBJ.total> -10000))) 
          print(paste("tuning parameter = ", tau));
        }    
      }
      time0  = proc.time()-pmt0
      print(time0)
      rm(OBJ.p0);rm(C.p)
      rm(OBJ.m0);rm(C.m)
      
      gam = GAM.total[,which.max(OBJ.total)]
      ind2 = which(gam==1)
      
      ind.total = which(OBJ.total> -100000000)
      OBJ.fin = unique(OBJ.total[ind.total])
      
      w = length(OBJ.fin)
      time.fin = rep(0,w)
      GAM.fin = matrix(0,p,w);GAM.fin[,1] = GAM.total[,which(OBJ.total==OBJ.fin[1])[1]]
      if(w>1){
        for(i in 2:length(OBJ.fin)){
          GAM.fin[,i] = GAM.total[,which(OBJ.total==OBJ.fin[i])[1]]
        }
      }
      
      rm(GAM.total);rm(OBJ.total)
      const = sum(exp(OBJ.fin-max(OBJ.fin)))
      posterior = exp(OBJ.fin-max(OBJ.fin))/const
      total.size = length(OBJ.fin)
      m = max(OBJ.fin)
      ind.m0 = which.max(OBJ.fin)
      gam = GAM.fin[,ind.m0]
      ind2 = which(gam==1);p.g = sum(gam)
      GAM.fin0 = cbind(GAM.fin0,GAM.fin)
      OBJ.fin0 = c(OBJ.fin0,OBJ.fin)
      
      #ind2 = true
      #gam = rep(0,p);
      #gam[ind2]=1
    }
    print("#################################")
    print("Post-process starts")
    print("#################################")
    
    OBJ.fin1 = unique(OBJ.fin0)
    
    w = length(OBJ.fin1)
    time.fin = rep(0,w)
    GAM.fin1 = Matrix(0,p,w,sparse=TRUE);GAM.fin1[,1] = GAM.fin0[,which(OBJ.fin0==OBJ.fin1[1])[1]]
    if(w>1){  
      for(i in 2:w){
        GAM.fin1[,i] = GAM.fin0[,which(OBJ.fin0==OBJ.fin1[i])[1]]
      }
    }
    rm(GAM.fin0)
    GAM = GAM.fin1
    OBJ = OBJ.fin1
    print("Done!")
    
    gam.map = GAM[, which.max(OBJ)]
    ind.map = which(gam.map==1);p.map = length(ind.map)
    POST_model = exp(OBJ - max(OBJ))/sum(exp(OBJ - max(OBJ)))
    POST_incl_prob = GAM%*%POST_model
    hppm = 1/sum(exp(OBJ - max(OBJ)))
    ind.MAP = which(gam.map == 1)
    print(ind.MAP)
    print("# of Searched Models by S5")
    print(length(OBJ))
    ind.marg=which(POST_incl_prob>0.5)
  return(list(GAM = GAM, OBJ = OBJ, phi = phi, Knots= Knots, K = K, post = POST_model, marg.inc = POST_incl_prob, 
              ind.MAP = ind.MAP, ind.marg = ind.marg, hppm.prob = hppm, tuning=tau ))
}
