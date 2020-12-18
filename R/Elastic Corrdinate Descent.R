
#' Elastic Coordiante Descent Algorithm
#' @param lambda penalty term
#' @param X input dataset
#' @param Y output dataset
#' @param alpha constant multiplying penalty term
#' @param maxit maxium number of iterations
#' @return a sequence of coefficients beta by fitting X to Y
#' @export

elnet_coord=function(lambda,X,Y,alpha,maxit=1e5){
  # compute the solution on a grid of lambdas
  nl<-length(lambda)
  p<-ncol(X)
  beta_cd <- matrix(0,nl,p)

  for(l in 1:nl){
    # Initial values
    b<-rep(0,p)
    r <- Y - X%*%b
    # Coordiante descent
    for(step in 1:maxit){
      for(j in 1:p){

        # Partial Residuals
        r <- r + X[,j]*b[j]

        # Soft-Threshold
        xr <- sum(X[,j]*r)
        xx <- sum(X[,j]^2)
        b[j] <-(abs(xr)-lambda[l]*alpha)/(xx+lambda[l]*(1-alpha))
        b[j] <- sign(xr)*ifelse(b[j]>0,b[j],0)

        # Update Residuals
        r <- r - X[,j]*b[j]

      }
    }

    beta_cd[l,]<-b
  }
  return(list(beta_cd=beta_cd))
}
