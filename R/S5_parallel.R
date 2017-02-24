S5_parallel = function(NC,X,y,ind_fun,model,tuning,tem,ITER=20,S=20,C0=2){
  requireNamespace("snowfall")
  requireNamespace("Matrix")
  #require(snowfall)
  #require(Matrix)
  n = nrow(X)
  p = ncol(X)
  y = y - mean(y)
  sfInit = snowfall::sfInit
  sfLibrary = snowfall::sfLibrary
  sfExportAll = snowfall::sfExportAll
  sfStop = snowfall::sfStop
  
  if(missing(tem)){tem = seq(0.4,1,length.out=20)^2}
  
  if(missing(ind_fun)){
    print("The prior on regression coefficietns is unspecified. The default is piMoM")
    ind_fun = BayesS5::ind_fun_pimom
    tuning <- BayesS5::hyper_par(type="pimom",X,y,thre = p^-0.5)  # tuning parameter selection for nonlocal priors
    print("The choosen hyperparameter tau")
    print(tuning)
  }
  
  if(missing(model)){
    print("The model prior is unspecified. The default is Bernoulli_Uniform")
    model = BayesS5::Bernoulli_Uniform
  }
  
  sfInit(parallel=TRUE, cpus=NC)
  sfLibrary(Matrix)
  sfExportAll()
  #sfExport( "S5" )
  pmt=proc.time()
  wrapper  = function(i){
    fit = BayesS5::S5(X=X,y=y,ind_fun=ind_fun, model = model,tuning=tuning,tem=tem,ITER=ITER,S=S,C0=C0,verbose=FALSE)
    return(fit)
  }
  out = sfLapply(1:NC,wrapper)
  print(proc.time()-pmt)
  sfStop()
  
  OBJ = NULL
  IND = NULL
  for(i in 1:NC){
    OBJ = c(OBJ,out[[i]]$OBJ)
    IND = c(IND,length(out[[i]]$OBJ))
  }
  IND = c(0,IND,0)
GAM = Matrix(0,p,length(OBJ),sparse= TRUE)
    for(i in 1:NC){
      gam = out[[i]]$GAM
      ind = (sum(IND[1:i])+1):(sum(IND[1:(i+1)]))
      GAM[,ind] = gam
    } 
GAM.fin0 = GAM; OBJ.fin0 = OBJ
OBJ.fin1 = unique(OBJ.fin0)

w = length(OBJ.fin1)
time.fin = rep(0,w)
GAM.fin1 = Matrix(0,p,w,sparse=TRUE);GAM.fin1[,1] = GAM.fin0[,which(OBJ.fin0==OBJ.fin1[1])[1]]
for(i in 2:w){
  GAM.fin1[,i] = GAM.fin0[,which(OBJ.fin0==OBJ.fin1[i])[1]]
  #  time.fin[i] = time.total[which(OBJ.total==OBJ.fin[i])[1]]
}
return(list(GAM=GAM.fin1,OBJ = OBJ.fin1, tuning = tuning))
}