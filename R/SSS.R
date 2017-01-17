SSS = function(X,y,ind_fun,model,tuning,N=1000,C0=1,verbose=TRUE){
  n = nrow(X)
  p = ncol(X)
#  assign("p",p,.GlobalEnv)
#  assign("n",n,.GlobalEnv)
  y = y -mean(y)
  requireNamespace("Matrix")
  requireNamespace("abind")
  abind = abind::abind
  
  if(missing(ind_fun)){
    print("The prior on regression coefficients is unspecified. The default is piMoM")
    ind_fun = BayesS5::ind_fun_pimom
    tuning <- BayesS5::hyper_par(type="pimom",X,y,thre = p^-0.5)  # tuning parameter selection for nonlocal priors
    print("The choosen hyperparameter tau")
    print(tuning)
    #assign("tuning", tuning, .GlobalEnv)
  }
  
  if(missing(model)){
    print("The model prior is unspecified. The default is Bernoulli_Uniform")
    model = BayesS5::Bernoulli_Uniform
  }
  
  r0=1 # auxiilary value
  verb = verbose # set verbose
  a0=0.01;b0=0.01 # hyperparameter for inverse gamma prior on sigma^2
  tau = 1 # tau
  #require(Matrix)
  #library(abind)

  ### initialize the model
  sam = sample(1:p,3)
  gam = rep(0,p)#;gam[sam]=1;curr = ind_fun(sam) + model(sam)
  ind2= which(gam==1)
  p.g=sum(gam) 
  ###
  
  curr = ind_fun(X[,ind2],y,n,p,tuning)  + model (ind2,p) # evaluate the objective value on the initial model
  GAM.fin0 = NULL # to save the binary vector of the searched models 
  OBJ.fin0 = NULL # to save the objective values of the searched models 
  g = tuning # set g for Zellner's g-prior
  B = tuning # set tau for piMoM or peMoM prior
  
  GAM = gam 
  OBJ = curr
  obj = OBJ
  
  ### search the neighborhood of the initial model
  C.p = rep(-1000000,p)
  IND = (1:p)[-ind2]
  for(i in 1:(p-p.g)){
    j=IND[i]
    gam.p = gam;gam.p[j]=1;ind.p=which(gam.p==1)
    int  = ind_fun(X[,ind.p],y,n,p,tuning) +model(ind.p,p)  
    obj.p =  c(int)
    if(is.na(obj.p)==TRUE){obj.p = -1000000}
    C.p[j] = obj.p
  }
  p.g = sum(gam)
  
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
  
  C.s = matrix(-1000000,round(n/2),p)
  
  for(i in 1:p.g){
    for(w in 1:(p-p.g)){
      j=ind2[i]
      u = (1:p)[-ind2][w]
      gam.s = gam;gam.s[j]=0;gam.s[u]=1;ind.s=which(gam.s==1)     
      int  = ind_fun(X[,ind.s],y,n,p,tuning) +model(ind.s,p)
      obj.s =  c(int)
      if(is.na(obj.m)==TRUE){obj.m = -1000000}
      C.s[i,u] = obj.s  
    }
  }
  
  print("#################################")
  print("SSS starts running")
  
  ### creat the arrays to save the neighborhood of visited models 
  
  d0 = 200
  OBJ.s0 = array(-1000000,dim=c(round(n/2),p,d0));OBJ.s0[,,1] = C.s
  OBJ.m0 = matrix(-1000000,p,d0);OBJ.m0[,1] = C.m
  OBJ.p0 = matrix(-1000000,p,d0);OBJ.p0[,1] = C.p
  
  ID = sum(2^(3*log(ind2))) # set the id of the initial model
  
  ############################################# Start !!!
  for(uu in 1:C0){ # repeat the SSS algorithm C0 times  
    
    GAM.total = Matrix(0,p,1000000,sparse=TRUE)
    OBJ.total = rep(-1000000,1000000)
    GAM.total[,1] = gam
    OBJ.total[1] = obj
    time.total = rep(0,1000000)
    
    pmt0 = proc.time()
    for(iter in 1:N){ 
        
        id = sum(2^(3*log(ind2))) # calculate the id of the current model
        id.ind = which(id==ID) # check whether the current model has already been visited using the id  
        leng = length(id.ind) # how many times visited?
        
        if(leng==0){ 
          ### if the current model has not been visited, search the neighborhood of the current model
          ID = c(ID,id)
          C.p = rep(-1000000,p)
          IND = (1:p)[-ind2]
          for(i in 1:(p-p.g)){
            j=IND[i]
            gam.p = gam;gam.p[j]=1;ind.p=which(gam.p==1)
            int  = ind_fun(X[,ind.p],y,n,p,tuning) +model(ind.p,p)    
            obj.p =  c(int)
            if(is.na(obj.p)==TRUE){obj.p = -100000}
            C.p[j] = obj.p
            ind.total = which(OBJ.total< -100000)[1]
            OBJ.total[ind.total] = obj.p
            GAM.total[,ind.total] =  gam.p
            time.total[ind.total] = (proc.time()-pmt0)[3]
          }
          p.g = sum(gam)
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
            ind.total = which(OBJ.total< -100000)[1]
            OBJ.total[ind.total] = obj.m
            GAM.total[,ind.total] =  gam.m
            time.total[ind.total] = (proc.time()-pmt0)[3]
          }
          
          C.s = matrix(-1000000,round(n/2),p)
          IND.s = ind2
          
          for(i in 1:p.g){
            for(w in 1:(p-p.g)){
              j=ind2[i]
              u = (1:p)[-ind2][w]
              gam.s = gam;gam.s[j]=0;gam.s[u]=1;ind.s=which(gam.s==1)     
              int  = ind_fun(X[,ind.s],y,n,p,tuning) +model(ind.s,p)
              obj.s =  c(int)
              if(is.na(obj.m)==TRUE){obj.m = -1000000}
              C.s[i,u] = obj.s 
              ind.total = which(OBJ.total< -100000)[1]
              OBJ.total[ind.total] = obj.s
              GAM.total[,ind.total] =  gam.s
              time.total[ind.total] = (proc.time()-pmt0)[3]
            }
          }
          d =  which(apply(OBJ.p0,2,mean)< (min(OBJ.p0[,1])+1))[1]
          if(is.na(d)==TRUE){
          OBJ.s0 = abind(OBJ.s0,C.s)
          OBJ.p0 = cbind(OBJ.p0,C.p)
          OBJ.m0 = cbind(OBJ.m0,C.m)  
          }else{
            OBJ.s0[,,d] = C.s
            OBJ.p0[,d] = C.p
            OBJ.m0[,d] = C.m
          }
        ###  
        }else{
        ### if the current model has already been visited, call the saved neighborhood 
          C.p = OBJ.p0[,(id.ind[1])];C.m = OBJ.m0[,(id.ind[1])];C.s = OBJ.s0[,,(id.ind[1])]
        }
        
        ### based on the neighborhood, choose a next model to move by randomly sampling proportion to the posterior probability 
        prop = exp(C.s-max(C.s))
        prop[which(is.na(prop))] = 0
        if(sum(prop)<0.1){prop[1]=1}
        sample.s = sample(1:length(prop),1,prob=prop)
        
        obj.s = C.s[sample.s]
        
        prop = exp(C.p-max(C.p))
        prop[which(is.na(prop))] = 0
        if(sum(prop)<0.1){prop[1]=1}
        sample.p = sample(1:length(prop),1,prob=prop)
        
        obj.p = C.p[sample.p]
        prop = exp(C.m-max(C.m))
        prop[which(is.na(prop))] = 0
        if(sum(prop)<0.1){prop[1]=1}
        sample.m = sample(1:length(prop),1,prob=prop)
        obj.m = C.m[sample.m]
        
        l1 = 1/(1+exp(obj.m-obj.p)+exp(obj.s-obj.p))
        l2 = 1/(1+exp(obj.p-obj.m)+exp(obj.s-obj.m))
        l3 = 1-l1 - l2
        if(l3<0){l3=0}
        z = sample(1:3,1,prob=c(l1,l2,l3))
        
        
        if(z==1){gam[sample.p]=1;obj = obj.p;curr=obj.p}
        if(z==2){gam[sample.m]=0;obj = obj.m;curr=obj.m}
        if(z==3){ wh = which(obj.s==C.s,arr.ind=TRUE)
        gam[ind2[wh[1]]]=0;gam[wh[2]] = 1
        obj = obj.s;curr=obj.s
        }
        ###
        
        ind2 = which(gam==1) # generate the binary variable for the new model
        p.g = sum(gam) # the size of the new model
        
        
      gam.pr = GAM.total[,which.max(OBJ.total)] # the MAP model among the searched models 
      obj.pr = max(OBJ.total) # its objective value
      ind2.pr = which(gam.pr==1) # its size
      if(verb==TRUE&iter%%50==0){
        print("#################################")
        curr = ind_fun(X[,ind2],y,n,p,tuning) + model(ind2,p) # print the objective value of the current model
        print("# of iterations");print(iter);print("The Selected Variables in the Searched MAP Model");print(ind2.pr);print("The Evaluated Object Value at the Searched MAP Model");print(obj.pr);
        print("Current Model");print(ind2);  print("The Evaluated Object Value at the Current Model");print(curr);
        print("The Number of Total Searched Models");print(length(unique(OBJ.total)))     }
  }
    time0  =proc.time()-pmt0
    print(time0)
  
   ### after a single SSS algorithm runs, summarize the searched models and their posterior model probability
    ind.total = which(OBJ.total> -100000) 
    OBJ.fin = unique(OBJ.total[ind.total]) # generate the unique objective values from models that have been visited multiple times
  
    w = length(OBJ.fin)
    time.fin = rep(0,w)
    GAM.fin = matrix(0,p,w)
    for(i in 1:length(OBJ.fin)){
      GAM.fin[,i] = GAM.total[,which(OBJ.total==OBJ.fin[i])[1]] # generate the unique binary variables of models from models that have been visited multiple times
    }
  
    const = sum(exp(OBJ.fin-max(OBJ.fin))) # calculate the normalizing constant
    posterior = exp(OBJ.fin-max(OBJ.fin))/const # calculate the posterior model probability
    total.size = length(OBJ.fin) # the total number of models that are searched by the SSS
    m = max(OBJ.fin) # the objective value of the MAP model
    ind.m0 = which.max(OBJ.fin) 
    gam = GAM.fin[,ind.m0] # the binary variable of the MAP model
    ind2 = which(gam==1);p.g = sum(gam) # its index and size
    GAM.fin0 = cbind(GAM.fin0,GAM.fin) # save the binary variables
    OBJ.fin0 = c(OBJ.fin0,OBJ.fin) # save the the objective values corresponding to the binary variables
  }
  print("#################################")
  print("Post-process starts")
  print("#################################")
  
  ### after C0 number of single SSS algorithm runs, summarize the searched models and their posterior model probability
  OBJ.fin1 = unique(OBJ.fin0)
  w = length(OBJ.fin1)
  time.fin = rep(0,w)
  GAM.fin1 = Matrix(0,p,w,sparse=TRUE);GAM.fin1[,1] = GAM.fin0[,which(OBJ.fin0==OBJ.fin1[1])[1]]
  for(i in 2:w){
    GAM.fin1[,i] = GAM.fin0[,which(OBJ.fin0==OBJ.fin1[i])[1]]
    #  time.fin[i] = time.total[which(OBJ.total==OBJ.fin[i])[1]]
  }
  rm(GAM.fin0)
  GAM = GAM.fin1 # the binary vaiables of searched models
  OBJ = OBJ.fin1 # the the objective values corresponding to the binary variables
  print("Done!")
  return(list(GAM = GAM,OBJ = OBJ) )
  }