MixN <- function(n,p,mu,Sigma,Kappa=0.8) {
  
  euf <- runif(n)
  X <- matrix(0, n, p)
  index1 <- (euf < Kappa)
  index2 <- !index1
  
  X[index1, ] <- mvtnorm::rmvnorm(sum(index1), mu, Sigma)
  X[index2, ] <- mvtnorm::rmvnorm(sum(index2), mu, 100*Sigma)
  
  X <- X / (sqrt(Kappa+100*(1-Kappa)))
  
  return(X)
  
}
