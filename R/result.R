result <-
function(fit){
#GAM = fit[-1,]; OBJ = fit[1,]

  GAM = fit$GAM; OBJ = fit$OBJ; tuning= fit$tuning
p = nrow(GAM)
  marg.gam = rep(0,p)
for(u in 1:ncol(GAM)){
  marg.gam = marg.gam + GAM[,u]*exp(OBJ[u]-max(OBJ))
}
marg.gam = marg.gam / sum(exp(OBJ-max(OBJ)))
gam0 = GAM[,which.max(OBJ)]
ind2 = which(gam0==1)
post = exp(OBJ-max(OBJ))/sum(exp(OBJ-max(OBJ)))
hppm = 1/sum(exp(OBJ-max(OBJ)))
print("# of Searched Models by S5");print(length(OBJ))
print("The MAP model is ")
print(which(gam0==1))
print(paste("with posterior probability",round(hppm,3) )) 
return(list(hppm = which(gam0==1), hppm.prob = hppm, marg.prob = marg.gam,gam = GAM, obj = OBJ, post = post, tuning = tuning) )
}
