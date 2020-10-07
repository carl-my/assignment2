#### OLS
olsFun <- function(data){
  ### set column names and add intercept column to X
  Y <- data[,1]
  N <- nrow(data)
  X <- cbind(rep(1, N), data[,2])
  
  ### calculate the formel and extract the Beta coefficient 
  beta_ols=(solve(t(X)%*% X) %*% (t(X) %*% Y))[2,1]
  
  return(beta_ols)
}


testData <- cbind( c(0.62, 0.18, 3.92, 0.80, -5.15),
c(0.44, 1.49, 0.69, 0.13, 1.90) )

olsFun(data = testData)

##### weighted least squares
wlsFun <- function(data, lambda){
  if (is.numeric(lambda)==FALSE){print("lambda is not numeric") } 
  else{
    ### create variable for the number of observations in the dataset
    N <- nrow(data)
    
    ### set column names and add intercept column for X
    Y <- data[,1]
    X <- cbind(rep(1,N), data[,2])
    
    ### create a zero matrix N x N
    Z <- matrix(0, N, N)
    
    ## make a forloop to put in the error terms on the diagonal 
    ## to create the error covariance matrix
    er <- NULL
    for (i in 1:N) {
      er[i] <- exp(X[i,2]*lambda)
      Z[i,i] <- er[i] 
    }
    ### calculate and extract the Beta coefficient
    beta_wls = ((solve(t(X)%*%(solve(Z))%*%X)%*%t(X)%*%
                   (solve(Z))%*%Y))[2,1]
  return(beta_wls)
  }
}

wlsFun(data = testData, lambda = 2)

#### FWLS
fwlsFun <- function(data, trueVar){
  y = data[,1]
  N <- nrow(data)
  X = cbind(rep(1,N), data[,2:ncol(data)])
  
  mod = lm(y ~ -1 +X)
  res = mod$residuals
  res2 = res^2
  
  ln_res2 = lm(log(res2) ~  -1 +X)
  lambda_hat = ln_res2$coefficients[2]
  
  if(trueVar == TRUE){
    error_cov = matrix(0, N, N)
    for(i in 1:N){
      error_cov[i,i] <- exp(X[i,2]*lambda_hat)
    }
    
  }
  else if(trueVar == FALSE)
  {
    error_cov = matrix(0, N, N)
    for(i in 1:N){
      error_cov[i,i] <- 1 + X[i,2]*lambda_hat
    }
    
  }
beta_fwls = (solve(t(X)%*%solve(error_cov)%*%X)%*%t(X)%*%solve(error_cov)%*%y)[2,1]

return(beta_fwls)
}
fwlsFun(data = testData, trueVar = TRUE)
fwlsFun(data = testData, trueVar = FALSE)


#### Data simulation
DataFun <- function(n, lambda) {

    # independent variable
    x <- runif(n, min = 0, max = 2)

    # standard deviation in epsilons normal distribution
    s <- NULL
    for (i in seq_len(n)) {
        s[i] <- exp(x[i]*lambda)
    }

    # error term
    epsilon <- rnorm(n, mean = rep(0, n), sd = s)

    beta <- 2

    # dependent variable
    y <- NULL
    for (i in seq_len(n)) {
        y[i] <- beta*x[i] + epsilon[i]
    }

    # container matrix
    mat <- matrix(data = 0, ncol = 2, nrow = n)

    # creating matrix of indepedent and dependent variables
    for (i in seq_len(n)) {
        mat[i, 1] <- y[i]
        mat[i, 2] <- x[i]
    }
    return(mat)
}

#####

SimFun <- function(n, sim_reps, seed, lambda) {
    set.seed(seed)
    R <- sim_reps

    # saving betas
    mat <- matrix(0, nrow = R, ncol = 4)

    for (i in seq_len(R)) {
    # data sim
    dat <- DataFun(n, lambda)
    # estimate sim
    mat[i,1] <- olsFun(dat)
    mat[i,2] <- wlsFun(dat, lambda)
    mat[i,3]<- fwlsFun(dat, trueVar = FALSE)
    mat[i,4]<- fwlsFun(dat, trueVar = TRUE)
    }
    betas <- apply(mat, 2, var)

    return(betas)
}



SimFun(25, 4, 2020, 2)

