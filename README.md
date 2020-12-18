# Efficient_Elastic_Solver
Install the package by using the following commands in R:
```{r}
install_github("Luffy-Yao/Efficient_Elastic_Solver")
library(Efficient_Elastic_Solver)

```

##  A. Ordinary least sqaure solver by iterative methods

We are going to solve $x$ for a linear system $Ax = b$
$A$ is required to be a $n \times n$ square invertible matrix, $b$ is a vector with the same dimension $n$.

An example is shown here.

```{r}
A=matrix(rnorm(10*10),nrow=10,ncol=10)
b = runif(10)
x_GS = solve_ols(A, b, c=1, maxit=10000, type = "GS")
```


## B. Leverage subsampling regression

$Y = X \beta + E$. 

An example is given as:
```{r}
  r=c(10,20,30)
  X= rt(500,df=3)
  b_unif = algo_leverage(X, r)$b_unif
  b_blev = algo_leverage(X, r)$b_blev
  
```

## C. Elastic net regression
A similar example is given as:
```{r}
  n = 200; p = 200
  X = matrix(rnorm(n * p), nrow = n)
  t_beta = 2 * runif(p)
  Y = X %*% t_beta + rnorm(n)
  
  b_elnet = elnet_coord(lambda = 5,X, Y, a = 0.5, maxit)$beta_cd
```
