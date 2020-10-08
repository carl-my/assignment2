library(StatProg)
library(tidyverse)

data("galaxies")
x = as.numeric(galaxies)

gammaUpdate = function(x, mu, sigma, pi){
  pdf = t(sapply(x, dnorm, mean = mu, sd = sigma))
  for(n in 1:length(x)){
    for(k in 1:length(mu)){
      pdf[1, 2] = pi[k]*pdf[n, k]
    }
  }
  gam = as.data.frame(matrix(NA, ncol = length(pi), nrow = length(x)))
  for(n in 1:length(x)){
    for(k in 1:length(mu)){
      gam[n, k] = pdf[n, k]/sum(pdf[n,])
    }
  }
  return(as.data.frame(gam))
}

mu = c(10000, 20000, 32000)
sigma = c(1000, 2000, 3000)
pi = c(1/3, 1/3, 1/3)
test = gammaUpdate(x, mu, sigma, pi)
