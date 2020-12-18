#' Gauss Seidel and Jacobi Algorithm for solving linear system of equations
#' @param A matrix input
#' @param b output vector
#' @param n number of cores users can specify
#' @param maxiter number of iterations
#' @param type specify method's type, "GS" denotes "Gauss-Seidel", "Jc" represents Jacobi,"PJC" denotes Parallel Jacobi
#' @return a vector of length n by solving AX=b
#' @export

solve_ols<-function(A,b,c=2,maxit=1e4,type="GS"){

  n = length(b)
  D = diag(diag(A))
  L = matrix(0, nrow = n, ncol = n)
  U = matrix(0, nrow = n, ncol = n)

  ## Upper and Low triangle matrix
  L[lower.tri(L)] = A[lower.tri(A)]
  U[upper.tri(U)] = A[upper.tri(A)]
  ## Initial matrix
  x=rep(0,n)

  ## Number of iterations
  iter=1


  ## Gauss-Seidel algorithm
  if(type=="GS")
   {
      while(iter<=maxit){
        x=solve(L+D)%*%(b-U%*%x)
        iter=iter+1
      }
    return(x)
  }

  ## Jacobi algorithm
  if(type=="JC"){
    while(iter<=maxit){
      x=solve(D)%*%(b-(L+U)%*%x)
      iter=iter+1
    }
    return(x)
  }

  ## Parallel Jacobi Algorithm
  if(type=="PJC"){
    require(doParallel)
    core <- makeCluster(c)
    registerDoParallel(core)
    while(iter<=maxit){
      outlist = foreach (j=1:n) %dopar%{
        tempM=matrix(A[j,],nrow=1,ncol=n)
        temp=tempM[,-j]
        x_1=x[-j]
        x[j]=(b[j]-temp%*%x_1)/A[j,j]
      }
      x=unlist(outlist)
    }
    return(x)
  }
}












