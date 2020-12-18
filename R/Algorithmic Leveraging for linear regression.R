#' Algorithmic Leveraging for linear linear regression
#' @param X input vector
#' @param r a list of subsample
#' @return simulated sampled based on the value of leverage
#' @export


algo_leverage<-function(X,r){
  n<-length(X)
  x_list<-X
  y_list<--x_list+rnorm(n,mean=0,sd=1)
  z_list<-cbind(x_list,y_list)
  ## Fit linear regression
  lm_1<-lm(y_list~x_list)

  ## Leverage
  h_mat<-diag(x_list%*%solve(t(x_list)%*%x_list)%*%t(x_list))

  ## length of each subsample
  n_1<-length(r)

  ## Randomly draw a subsample from data with uniform probability
  b_unif <- matrix(0,nrow =n_1 ,ncol = n)
  b_blev <- matrix(0,nrow =n_1,ncol = n)

  for (i in 1:n_1){
    for (j in 1:n){
      unif <- sample(n,r[i],replace = TRUE)
      x_unif <- x_list[unif]
      y_unif <- y_list[unif]
      b_unif[i,j]= solve(t(x_unif)%*%x_unif)%*%t(x_unif)%*%y_unif

      blev <- sample(n,r[i],replace = TRUE,prob =h_mat/sum(h_mat))
      x_blev <- x_list[blev]
      y_blev <- y_list[blev]
      b_blev[i,j]= solve(t(x_blev)%*%x_blev)%*%t(x_blev)%*%y_blev
    }

  }
  return(list(b_unif=b_unif,b_blev=b_blev))
}

