# This follows the methods used by Zhu et al (2014)
DatGen <- function(p0,p,d,rank,n){
  
  # Generate the covariances
  Cogen <- function(P,P0){
    A = matrix(NA,P,P)
    for(j in 1:P){
      A[j,] = runif(P)*rbinom(P,1,P0)
      A[j,j] = 1
    }
    out = tcrossprod(A,A)
#     for(j in 1:P){
#       out[j,] = out[j,]/out[j,j]
#     }
    # out = abs(cor(out))
    out
  }
  
  SigX = Cogen(p,p0)
  Sig = Cogen(d,p0)
  
  
  # Generate the covariates
  require(mvtnorm)
  X = rmvnorm(n=n,sigma = SigX)
  
  # Generate the coefficients
  if(missing(rank)){
    B <- matrix(rbinom(p*d,1,0.5),p,d)
  }else{
    U = matrix(rnorm(p*rank),p,rank)
    V = matrix(rnorm(d*rank),d,rank)
    L = diag(seq(100,20,-20))
    B <- U%*%L%*%t(V)
  }
  
  # Generate the responses
  Y = t(apply(X,1,function(x){x%*%B + rmvnorm(1,sigma = Sig)}))
  
  # Spit it back out
  list(Y=Y,X=X)
}
