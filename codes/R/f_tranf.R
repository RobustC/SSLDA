f_tranf <- function(x,xy,n1,S,mx,my){
  p=dim(xy)[2]
  n=dim(xy)[1]
  n2=n-n1
  f_h = numeric(p)
  for(j in 1:p){
    f_1 = mx[j]+S[j,j]^(-1/2)*F_tu(x=xy[1:n1,j],t=x[j],del=1/(2*n1))
    f_2 = my[j]+S[j,j]^(-1/2)*F_tu(x=xy[(n1+1):n,j],t=x[j],del=1/(2*n2))
    f_h[j] = (n1*f_1+ n2*f_2)/n
  }
  
  return(f_h)
  
}



