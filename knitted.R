## ----echo=FALSE--------------------------------------------------------------------------------
library(StatProg)
library(dplyr) 
library(ggplot2)

#### OLS
olsFun <- function(data){
  ### create variable for the number of observations in the dataset
  N <- nrow(data)
  
  ### set column names and add intercept column to X
  Y <- data[,1]
  X <- cbind(rep(1, N), data[,2])
  
  ### calculate the formel and extract the Beta coefficient 
  beta_ols=(solve(t(X)%*% X) %*% (t(X) %*% Y))[2,1]
  
  return(beta_ols)
}


## ----echo=FALSE--------------------------------------------------------------------------------
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


## ----echo=FALSE--------------------------------------------------------------------------------
fwlsFun <- function(data, trueVar){
  y = data[,1]                             # y is given the value of the first column in the data  
  N <- nrow(data)                          # N is given the value of the number of rows in the data
  X = cbind(rep(1,N), data[,2:ncol(data)]) # X is given the value of the remaining columns in the data 
                                           # and an intercept is added
  
  mod = lm(y ~ -1 +X)         # We use lm() to estimate a linear regression model, mod, with an intercept 
  res = mod$residuals         # res is assigned the value of the residuals of mod
  res2 = res^2
  
  ln_res2 = lm(log(res2) ~  -1 +X)        # We use the lm() function to estimate a new model
                                          # In this model, the dependent variable is ln(res2), 
                                          # which is explained by X*lambda + v
  lambda_hat = ln_res2$coefficients[2]    # Thereby, the value of lambda_hat is the value of 
                                          # the second estimated coefficient, 
                                          # since the first is the intercept 
  
  # We check if trueVar is TRUE, if so, we assume the structure of the error variance to be exp(x*lambda)
  if(trueVar == TRUE){
    error_cov = matrix(0, N, N)                 # We create an empty square matrix 
    for(i in 1:N){                              # The loop repeats itself for every row of X
      error_cov[i,i] <- exp(X[i,2]*lambda_hat)  # We fill the diagonal in the empty square matrix 
                                                # with the error variance, exp(x*lambda_hat) for each value of x
    }
    
  }
  # We check if trueVar is TRUE, if not, we assume the structure of the error variance to be 1 + x*lambda
  else if(trueVar == FALSE)
  {
    error_cov = matrix(0, N, N)                 # We create an empty square matrix.
    for(i in 1:N){                              # The loop repeats itself for every row of X.
      error_cov[i,i] <- 1 + X[i,2]*lambda_hat   # We fill the diagonal in the empty square matrix 
                                                # with the error variance, 1 + x*lambda_hat, for each value of x
    }
    
  }
beta_fwls = (solve(t(X)%*%solve(error_cov)%*%X)%*%t(X)%*%solve(error_cov)%*%y)[2,1]
# beta_fwls is calculated and returned 

return(beta_fwls)
}


## ----echo=FALSE--------------------------------------------------------------------------------
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


## ----echo=FALSE--------------------------------------------------------------------------------
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
    mat[i,4]<- fwlsFun(dat, trueVar = TRUE)
    mat[i,3]<- fwlsFun(dat, trueVar = FALSE)
    }

    # calculate variance of betas and returns the values as a matrix
    betas <- apply(mat, 2, var)
    return(betas)
}


## ----echo=FALSE--------------------------------------------------------------------------------
##### 4. Plot variance estimates
## create variable for the different values on n
x <- c(25, 50, 100, 200, 400)

## make zero matrix to use in for loop
var_obs <- matrix(0, ncol = 4, nrow = 5)

for (i in seq_along(x)) {
    var_obs[i,] <- (SimFun(x[i], 100, 2020, 2))
}
### pute variance estimates in a data frame
var_obs = as.data.frame(var_obs)
rownames(var_obs) = c(25, 50, 100, 200, 400)

### add column and row names
colnames(var_obs) = c("OLS","WLS","FWLST","FWLSF")
rownames(var_obs) = c("25","50","100","200","400")


## ----------------------------------------------------------------------------------------------
### create line plot with sample size on x-axis and variance estimate on y-axis
ggplot(as.data.frame(var_obs),aes(x=as.numeric(rownames(var_obs)))) +
  geom_line(aes(y = OLS, colour = "OLS")) + 
  geom_line(aes(y = WLS, colour = "WLS")) +
  geom_line(aes(y = FWLST, colour = "FWLST")) + 
  geom_line(aes(y = FWLSF, colour = "FWLSF")) +
  labs(x="Number of obsvervations",y="Variance") +
  theme(legend.title = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks=c(25,50,100,200,400),limits=c(25, 400))


## ----echo=FALSE--------------------------------------------------------------------------------
##################################################
# start of part2
##################################################

# importing the galaxy dataset as a manageable data frame
galaxies <- as.data.frame(galaxies)         
names(galaxies) <- "km"


## ----------------------------------------------------------------------------------------------
# plotting the density
ggplot(galaxies, aes(x = km)) +
  geom_density()


## ----echo=FALSE--------------------------------------------------------------------------------
# turning the dataframe into a vector
galaxies = galaxies$km


## ----echo=FALSE--------------------------------------------------------------------------------
# Function for gamma update 
gammaUpdate = function(x, mu, sigma, pi){
    pdf = t(sapply(x, dnorm, mean = mu, sd = sigma))    
    # pdf is assigned the value of the pdf of the normal distribution, 
    # with mean=mu and standard deviation=sigma evaluated at each x
    # this creates a n*K matrix where n = #observations in data and 
    # K = #components in the mixture

    # we create a for loop that for each row...  
    for(n in 1:length(x)){
        # ... goes through every component 
        for(k in 1:length(mu)){
            # each value of pdf is multiplied by pi  
            pdf[n, k] = pi[k]*pdf[n, k]      
        }
    }
    # an empty matrix created 
    gamma = as.data.frame(matrix(NA, ncol = length(pi), nrow = length(x)))

    # we create a loop that for each row...
    for(n in 1:length(x)){
        # ... goes through every component 
        for(k in 1:length(mu)){
            # gamma[n,k] is assigned the value of pdf[n,k] divided by 
            # the sum of all values on the same row
            gamma[n, k] = pdf[n, k]/sum(pdf[n,])
        }
    }
    # This creates a n*K matrix filled with the gamma-values for each x
    return(as.data.frame(gamma))
}


## ----echo=FALSE--------------------------------------------------------------------------------
# mu
muUpdate = function(x, gamma){
    # extracting number of components
    K <- ncol(gamma)
    # container
    mu <- NULL
    for (i in seq_len(K)) {
        # implementing function
        mu[i] <- sum(gamma[,i]*x)/sum(gamma[,i])
    }
    return(mu)
}


## ----echo=FALSE--------------------------------------------------------------------------------
### Sigma
sigmaUpdate = function(x, gamma, mu){
    N = matrix(0, ncol= ncol(gamma))          # A vector of zeros to be filled with the sum of gamma values in one cluster 
    sigma = matrix(0, ncol = ncol(gamma))     # A vector of zeros to be filled with the K values of sigma
    # a for loop that for every cluster...
    for(k in 1:ncol(gamma)){
        # ... goes through every row in x
        for(n in 1:length(x)){
            sigma[k] = sigma[k] + gamma[n,k]*(x[n]-mu[k])^2
            N[k] = N[k] + gamma[n,k]
        }
        sigma[k] = sqrt(sigma[k]/N[k])
    }
    return(sigma)
}


## ----echo=FALSE--------------------------------------------------------------------------------
### Pi
piUpdate = function(gamma){
  pi <- NULL
  ## a for loop for every cluster that takes the sum of each column gamma
  ## divided with the sum of gamma
  for (i in 1:ncol(gamma)) {
    pi[i] <- sum(gamma[,i])/sum(gamma)
  }
  return(pi)
}


## ----echo=FALSE--------------------------------------------------------------------------------
### A function for log-likelihood 
loglik = function(x, pi, mu, sigma){
    sum_pdf = matrix(0, nrow = length(x))     # A column vector of zeros
    loglike = 0     # To make sure loglike is 0

    # A for loop that for each row...
    for(n in 1:length(x)){
        # ... goes through every cluster 
        for(k in 1:length(pi)){
            sum_pdf[n] = sum_pdf[n] + (pi[k] * dnorm(x[n], mu[k], sigma[k]))
            # for this row, n, we add the value of the pdf of the normal distribution, 
            # with mean=mu[k] and standard deviation=sigma[k] evaluated at x[n], 
            # multiplied by pi[k], to sum_pdf[n]
        }
        # we take the log of this sum and add it to loglike 
        loglike = loglike + log(sum_pdf[n])
        # the loop is then started again 
    }
    # When the loop is finished, loglike has the value of the log-likelihood 
    return(loglike)
}


## ----echo=FALSE--------------------------------------------------------------------------------
initialValues = function(x, K, reps = 100){
    mu = rnorm(K, mean(x), 5)
    sigma = sqrt(rgamma(K, 5))
    p = runif(K)
    p = p/sum(p)
    currentLogLik = loglik(x, p, mu, sigma)
    for(i in 1:reps){
        mu_temp = rnorm(K, mean(x), 10)
        sigma_temp = sqrt(rgamma(10, 5))
        p_temp = runif(K)
        p_temp = p_temp/sum(p_temp)
        tempLogLik = loglik(x, p_temp, mu_temp, sigma_temp)
        if(tempLogLik > currentLogLik){
            mu = mu_temp
            sigma = sigma_temp
            p = p_temp
            currentLogLik = tempLogLik
        }
    }
    return(list("mu" = mu, "sigma" = sigma, "p" = p))
}


## ----echo=FALSE--------------------------------------------------------------------------------
##################################################
# EM algo
##################################################
EM = function(x, K, tol = 0.001){
    # initializing values
    inits = initialValues(x, K, 1000)
    mu = inits$mu
    sigma = inits$sigma
    prob = inits$p

    # containers
    prevLoglik <- 0
    # needs to have a value > 0.001 for the while loop to start
    loglikDiff<- 1

    # updating and saving parameters until the tolorance condition is met
    while(loglikDiff > tol){

        # updating gamma aswell as the other parameters
        gamma <- gammaUpdate(x, mu, sigma, prob)
        mu <- muUpdate(x, gamma)
        sigma <- sigmaUpdate(x, gamma, mu)
        prob <- piUpdate(gamma)

        # calculating new likelihood
        currentLogLik <- loglik(x, prob, mu, sigma)

        # calculaitng the difference which the while loop tests
        loglikDiff <- abs(prevLoglik - currentLogLik)

        # saving loglik for the next iteration
        prevLoglik <- currentLogLik 

    }

    # returns all parameters as a list
    return(list('loglik' = currentLogLik, 'mu' = mu, 'sigma' = sigma, 'prob' = prob))
}


## ----echo=FALSE--------------------------------------------------------------------------------
# difference starting values produce different local optimums, there we need to set a seed
set.seed(1996)

final_plot = matrix(0, ncol = 4, nrow= length(galaxies))    # an empty matrix is created 
loglik_values = NULL      # loglik_values is made sure to be empty

# a for loop that goes through the different values of k, from 2 to 5
for(k in 2:5){
    # for each value of k, the EM-functions generates a new z
    z = EM(galaxies, k)

    # loglik_values is assigned the value of the log-likelihood of the function
    loglik_values[(k-1)] = z$loglik

    # a for loop that repeats k times and sum the value of the pdf of the normal distribution, 
    # with mean=mu[i] and standard deviation=sigma[i] evaluated at x, multiplied by pi[k].
    for(i in 1:k){
        final_plot[,(k-1)] = final_plot[,(k-1)] + z$prob[i] * dnorm(galaxies, z$mu[i], z$sigma[i])
    }
}

### Plot log likelihood for different values on K
loglik_values <- as.data.frame(loglik_values) %>%  ## convert loglik values to data frame and add one column "K"
  mutate("K" = (2:5))


## ----------------------------------------------------------------------------------------------
### create density plot with number of K:s on x-axis and log likelihood on y-axis
ggplot(data = loglik_values) + 
  geom_line(aes(x=K, y=loglik_values), color = "blue") +
  labs(x="Number of components",y= "Log likelihood")


## ----echo=FALSE--------------------------------------------------------------------------------
final_plot = as.data.frame(final_plot)
colnames(final_plot) = c("K = 2", "K = 3", "K = 4", "K = 5")
final_plot = cbind(final_plot, galaxies)


## ----------------------------------------------------------------------------------------------
ggplot(data = final_plot, aes(x = galaxies)) +
    geom_line(aes(y = `K = 2`, color = "K = 2"), size = 1) + 
    geom_line(aes(y = `K = 3`, color = "K = 3"), size = 1) + 
    geom_line(aes(y = `K = 4`, color = "K = 4"), size = 1) + 
    geom_line(aes(y = `K = 5`, color = "K = 5"), size = 1) + 
    geom_density(aes(fill = "Density plot"), color = "pink", alpha = 0.2, size = 0) + 
    labs(x = "km", y = "Density",title="Plot of different values of K and the density plot of galaxies") +
    theme(legend.title = element_blank(),legend.position = c(.95, .95),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(6, 6, 6, 6)) 

