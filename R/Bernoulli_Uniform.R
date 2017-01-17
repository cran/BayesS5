Bernoulli_Uniform <-
function(ind,p){
  p.g=length(ind)
  sb = lbeta(1+p.g,1+p-p.g)
  return(sb) }
