library(StatProg)
library(ggplot2)

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

##################################################
# EM algo
##################################################
EM = function(x, K, tol = 0.001){
  inits = initialValues(x, K, 1000)
  mu = inits$mu
  sigma = inits$sigma
  prob = inits$p
  
  prevLoglik <- 0
  loglikDiff<- 1

  # while loop
  while(loglikDiff > tol){
    
    gamma <- gammaUpdate(x, mu, sigma, prob)
    mu <- muUpdate(x, gamma)
    sigma <- sigmaUpdate(x, gamma, mu)
    prob <- piUpdate(gamma)
    
    currentLogLik <- loglik(x, prob, mu, sigma)
    
    loglikDiff <- abs(prevLoglik - currentLogLik)
    
    prevLoglik <- currentLogLik 

  }
  
  return(list('loglik' = currentLogLik, 'mu' = mu, 'sigma' = sigma, 'prob' = prob))
}

# se hur loglikelihood förändras för varje iteration
z = EM(galaxies, 5)

#################################################
## H�R �R MITT F�RSLAG P� HUR VI KAN G�RA
#################################################

test = matrix(0, ncol = 4, nrow= length(galaxies))
for(k in 2:5){
  z = EM(galaxies, k)
  temp = 0
  for(i in 1:k){
    test[,(k-1)] = test[,(k-1)] + z$prob[i] * dnorm(galaxies, z$mu[i], z$sigma[i])
  }
}
test = as.data.frame(test)

ggplot(data = test) +
  geom_density(aes(x = V1), color = "white", fill = "red", alpha = 0.2) + 
  geom_density(aes(x = V2), color = "white",fill = "blue", alpha = 0.2) + 
  geom_density(aes(x = V3), color = "white",fill = "green", alpha = 0.2) + 
  geom_density(aes(x = V4), color = "white",fill = "purple", alpha = 0.2) 

inits = initialValues(x, K, 1000)
mu = inits$mu
sigma = inits$sigma
prob = inits$p

i <- 1
logliktest <- 0

prevLoglik <- 0
loglikDiff<- 0
# while loop
while(loglikDiff > tol){

gamma <- gammaUpdate(x, mu, sigma, pi)
mu <- muUpdate(x, gamma)
sigma <- sigmaUpdate(x, gamma, mu)
pi <- piUpdate(gamma)

currentLogLik <- loglik(x, probs, mu, sigma)

loglikDiff <- abs(prevLoglik - currentLogLik)

prevLoglik <- currentLogLik 

logliktest[i]  <- currentLogLik
i = i + 1

}

return(list('loglik' = currentLogLik, 'loglikt' = logliktest, 'mu' = mu, 'sigma' = sigma, 'prob' = prob))
}

# se hur loglikelihood förändras för varje iteration
EM(galaxies, 3)



