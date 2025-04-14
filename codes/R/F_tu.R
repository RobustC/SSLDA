F_tu = function(x,t,del){
  freq = sum(x<=t)/length(x)
  if(freq<del) {T=del
  }else if(freq>=del & freq<= 1-del) {T=freq
  }else if(freq>=1-del){T=1-del}
  F = qnorm(T)
  
  return(F)
}