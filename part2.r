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
set.seed(1996)
final_plot = matrix(0, ncol = 4, nrow= length(galaxies))
loglik_values = NULL
for(k in 2:5){
  z = EM(galaxies, k)
  loglik_values[(k-1)] = z$loglik
  for(i in 1:k){
    final_plot[,(k-1)] = final_plot[,(k-1)] + z$prob[i] * dnorm(galaxies, z$mu[i], z$sigma[i])
  }
}
loglik_values
final_plot = as.data.frame(final_plot)
colnames(final_plot) = c("K = 2", "K = 3", "K = 4", "K = 5")
final_plot = cbind(final_plot, galaxies)

ggplot(data = final_plot, aes(x = galaxies)) +
  geom_line(aes(y = `K = 2`, color = "K = 2"), size = 1) + 
  geom_line(aes(y = `K = 3`, color = "K = 3"), size = 1) + 
  geom_line(aes(y = `K = 4`, color = "K = 4"), size = 1) + 
  geom_line(aes(y = `K = 5`, color = "K = 5"), size = 1) + 
  geom_density(aes(fill = "Density plot"), color = "pink", alpha = 0.2, size = 0) + 
  labs(x = "km", y = "Density",title="TITLE") +
  theme(legend.title = element_blank(),legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6))


