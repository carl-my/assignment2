library(StatProg)
library(tidyverse)

galaxies <- as.data.frame(galaxies)
names(galaxies) <- "km"

ggplot(galaxies, aes(x = km)) +
  geom_density()

galaxies = galaxies$km

gammaUpdate = function(x, mu, sigma, pi){
  pdf = t(sapply(x, dnorm, mean = mu, sd = sigma))
  for(n in 1:length(x)){
    for(k in 1:length(mu)){
      pdf[n, k] = pi[k]*pdf[n, k]
    }
  }
  gamma = as.data.frame(matrix(NA, ncol = length(pi), nrow = length(x)))
  for(n in 1:length(x)){
    for(k in 1:length(mu)){
      gamma[n, k] = pdf[n, k]/sum(pdf[n,])
    }
  }
  return(as.data.frame(gamma))
}

# mu
muUpdate = function(x, gamma){
K <- ncol(gamma)
mu <- NULL
for (i in seq_len(K)) {
    mu[i] <- sum(gamma[,i]*x)/sum(gamma[,i])
}
return(mu)
}

### Sigma
sigmaUpdate = function(x, gamma, mu){
  N = matrix(0, ncol= ncol(gamma))
  sigma = matrix(0, ncol = ncol(gamma))
  for(k in 1:ncol(gamma)){
    for(n in 1:length(x)){
      sigma[k] = sigma[k] + gamma[n,k]*(x[n]-mu[k])^2
      N[k] = N[k] + gamma[n,k]
    }
    sigma[k] = sqrt(sigma[k]/N[k])
  }
  return(sigma)
}

piUpdate = function(gamma){
  pi <- NULL
  for (i in 1:ncol(gamma)) {
      pi[i] <- sum(gamma[,i])/sum(gamma)
  }
  return(pi)
}


#### Testkod
mu = c(10, 20, 30)
sigma = c(2, 2, 2)
probs = c(1/3, 1/3, 1/3)

resp = gammaUpdate(galaxies, mu, sigma, probs)
mu = muUpdate(galaxies, resp)
sigma = sigmaUpdate(galaxies, resp, mu)
probs = piUpdate(resp)

cat("mu:", mu,
"\nsigma:", sigma,
"\nprobs:", probs,
"\nresp[1,]", resp[1,])

### Log-likelihood 
loglik = function(x, pi, mu, sigma){
  sum_pdf = matrix(0, nrow = length(x))
  loglike = 0
  for(n in 1:length(x)){
    for(k in 1:length(pi)){
      sum_pdf[n] = sum_pdf[n] + (pi[k] * dnorm(x[n], mu[k], sigma[k]))
    }
    loglike = loglike + log(sum_pdf[n])
  }
  return(loglike)
}

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
