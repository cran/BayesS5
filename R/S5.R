S5 <- function(X,y,ind_fun,model,tuning,tem,ITER=20,S=20,C0=5,verbose=TRUE){
  n = nrow(X)
  p = ncol(X)
   y = y -mean(y)
  #requireNamespace()
  requireNamespace("Matrix")
  Matrix = Matrix::Matrix
  if(missing(tem)){tem = seq(0.4,1,length.out=20)^2}
  
  if(missing(ind_fun)){
    print("The prior on regression coefficients is unspecified. The default is piMoM")
    ind_fun = BayesS5::ind_fun_pimom
    tuning <- BayesS5::hyper_par(type="pimom",X,y,thre = p^-0.5)  # tuning parameter selection for nonlocal priors
    print("The choosen hyperparameter tau")
    print(tuning)
    #assign("tuning", tuning, .GlobalEnv)
  }
  
  #else{
  #  a = 0
  #  if(ind_fun == "pimom"){ind_fun = BayesS5::ind_fun_pimom; a = 1}
  #  if(ind_fun == "pemom"){ind_fun = BayesS5::ind_fun_pemom; a = 1}
  #  if(ind_fun == "g-prior"){ind_fun = BayesS5::ind_fun_g; a = 1}
  #  if(a == 0){stop("The ind_fun is not in the list!")}
  #  if(missing(tuning)){ stop("The tuning parameter is missing!")  }
  #}
  
  if(missing(model)){
    print("The model prior is unspecified. The default is Bernoulli_Uniform")
    model = BayesS5::Bernoulli_Uniform
  }
  print("#################################")
  print("S5 starts running")
  
  A3 = S;r0=1
  verb = verbose
  
  a0=0.01;b0=0.01
  tau = 1
  IT = length(tem)
  IT.seq = rep(ITER,IT)
  #require(Matrix)
  g = B = tuning
  sam = sample(1:p,3)
  gam = rep(0,p);
  gam[sam]=1
  curr = ind_fun(X[,sam],y,n,p,tuning) + model(sam,p)
  p.g=sum(gam)
  ind2= which(gam==1)
  GAM.fin0 = NULL
  OBJ.fin0 = NULL
  
  for(uu in 1:C0){
    
    C.p = rep(-1000000,p)
    C.m = rep(-1000000,p)
    
    
    
    #prior based on model1
    curr = ind_fun(X[,ind2],y,n,p,tuning) + model(ind2,p)
    GAM = gam
    OBJ = curr
    obj = OBJ
    
    if(p.g>0){
      fit = solve(crossprod(X[,ind2])+diag(p.g))%*%crossprod(X[,ind2],y)
      res = y-X[,ind2]%*%fit
      corr = as.vector(cor(res,X))
      ind.ix = sort.int(abs(corr),decreasing=TRUE,index.return=TRUE)$ix
      s = c(ind2,ind.ix)
    }else{res=y
    corr = as.vector(cor(res,X))
    ind.ix = sort.int(abs(corr),decreasing=TRUE,index.return=TRUE)$ix
    s = ind.ix
    }
    
    if(p<50){p00=10}else{p00=round(p/10)}
    
    
    size = A3
    IND = s[1:(size+p.g)]
    p.ind = length(IND)
    
    C.p = rep(-100000,p)
    for(i in (p.g+1):p.ind){
      j=IND[i]
      gam.p = gam;gam.p[j]=1;ind.p=which(gam.p==1)
      int  = ind_fun(X[,ind.p],y,n,p,tuning) + model(ind.p,p)    
      obj.p =  c(int)
      if(is.na(obj.p)==TRUE){obj.p = -100000}
      C.p[j] = obj.p
    }
    
    C.m = rep(-1000000,p)
    IND.m = ind2
    p.ind.m = length(IND.m)
    for(i in 1:p.g){
      j=ind2[i]
      gam.m = gam;gam.m[j]=0;ind.m=which(gam.m==1)     
      int  = ind_fun(X[,ind.m],y,n,p,tuning) +model(ind.m,p)
      obj.m =  c(int)
      if(is.na(obj.m)==TRUE){obj.m = -1000000}
      C.m[j] = obj.m  
    }
    
    p.g = sum(gam)
    
    OBJ.m0 = matrix(C.m,p,1)
    OBJ.p0 = matrix(C.p,p,1)
    ID = sum(2^(3*log(ind2)))
    ID.obj = ID
    it=1
    #GAM.total = matrix(0,p,50000)
    GAM.total = Matrix(0,p,50000,sparse=TRUE)
    OBJ.total = rep(-100000,50000)
    GAM.total[,1] = gam
    OBJ.total[1] = obj
    time.total = rep(0,50000)
    it=1
    ID0 = NULL
    INT = NULL
    
    pmt0 = proc.time()
    for(it in 1:IT){   
      IT0 = IT.seq[it] 
      pq=0
      for(iter in 1:IT0){
        
        
        id = sum(2^(3*log(ind2)))
        id.ind = which(id==ID)
        leng = length(id.ind)
        
        
        if(leng==0){
          ID = c(ID,id)
          C.p = rep(-100000,p)
          for(i in (p.g+1):p.ind){
            j=IND[i]
            gam.p = gam;gam.p[j]=1;ind.p=which(gam.p==1)
            int  = ind_fun(X[,ind.p],y,n,p,tuning) + model(ind.p,p)    
            obj.p =  c(int)
            if(is.na(obj.p)==TRUE){obj.p = -100000}
            C.p[j] = obj.p
            ind.total = which(OBJ.total< -10000)[1]
            OBJ.total[ind.total] = obj.p
            GAM.total[,ind.total] =  gam.p
            time.total[ind.total] = (proc.time()-pmt0)[3]
          }
          p.g = sum(gam)
          C.m = rep(-100000,p)
          IND.m = ind2
          p.ind.m = length(IND.m)
          for(i in 1:p.g){
            j=ind2[i]
            gam.m = gam;gam.m[j]=0;ind.m=which(gam.m==1)     
            int  = ind_fun(X[,ind.m],y,n,p,tuning) + model(ind.m,p)    
            obj.m =  c(int)
            if(is.na(obj.m)==TRUE){obj.m = -100000}
            C.m[j] = obj.m  
            ind.total = which(OBJ.total< -10000)[1]
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
        curr = ind_fun(X[,ind2],y,n,p,tuning)  +model(ind2,p)
        
        
        if(p.g>0){
          fit = solve(crossprod(X[,ind2])+diag(p.g))%*%crossprod(X[,ind2],y)
          res = y-X[,ind2]%*%fit
          corr = as.vector(crossprod(res,X))
          ind.ix = sort.int(abs(corr),decreasing=TRUE,index.return=TRUE)$ix
          s = c(ind2,ind.ix)
        }else{res=y
        corr = as.vector(crossprod(res,X))
        ind.ix = sort.int(abs(corr),decreasing=TRUE,index.return=TRUE)$ix
        s = ind.ix
        }   
        
        size = A3
        IND = s[1:(size+p.g)]
        p.ind = length(IND)
        
        id = sum(2^(3*log(ind2)))
        id.ind = which(id==ID.obj)
        
        leng = length(id.ind)
        if(leng==0){
          ID.obj = c(ID.obj,id)
          OBJ = c(OBJ,curr)
          GAM= cbind(GAM,gam)
        }
      }
      
      gam.pr = GAM.total[,which.max(OBJ.total)]
      obj.pr = max(OBJ.total)
      ind2.pr = which(gam.pr==1)
      if(verb==TRUE){
        print("#################################")
        curr = ind_fun(X[,ind2],y,n,p,tuning) + model(ind2,p)  
        print("Inverse Temperature");print(tem[it]);print("The Selected Variables in the Searched MAP Model");
        print(ind2.pr);print("The Evaluated Object Value at the Searched MAP Model");print(obj.pr);
        print("Current Model");print(ind2);  
        print("The Evaluated Object Value at the Current Model");print(curr);
        print("The Number of Total Searched Models");
        print(length(unique(OBJ.total))) 
      }    
    }
    time0  =proc.time()-pmt0
    print(time0)
    rm(OBJ.p0);rm(C.p)
    rm(OBJ.m0);rm(C.m)
    
    gam = GAM.total[,which.max(OBJ.total)]
    ind2 = which(gam==1)
    
    ind.total = which(OBJ.total> -100000)
    OBJ.fin = unique(OBJ.total[ind.total])
    
    w = length(OBJ.fin)
    time.fin = rep(0,w)
    GAM.fin = matrix(0,p,w);GAM.fin[,1] = GAM.total[,which(OBJ.total==OBJ.fin[1])[1]]
    for(i in 2:length(OBJ.fin)){
      GAM.fin[,i] = GAM.total[,which(OBJ.total==OBJ.fin[i])[1]]
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
  }
  print("#################################")
  print("Post-process starts")
  print("#################################")
  
  OBJ.fin1 = unique(OBJ.fin0)
  
  w = length(OBJ.fin1)
  time.fin = rep(0,w)
  GAM.fin1 = Matrix(0,p,w,sparse=TRUE);GAM.fin1[,1] = GAM.fin0[,which(OBJ.fin0==OBJ.fin1[1])[1]]
  for(i in 2:w){
    GAM.fin1[,i] = GAM.fin0[,which(OBJ.fin0==OBJ.fin1[i])[1]]
    #  time.fin[i] = time.total[which(OBJ.total==OBJ.fin[i])[1]]
  }
  rm(GAM.fin0)
  GAM = GAM.fin1
  OBJ = OBJ.fin1
  print("Done!")
  return(list(GAM = GAM,OBJ = OBJ, tuning=tuning))
}
