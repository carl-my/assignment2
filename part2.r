library(StatProg)
library(tidyverse)

galaxies <- as.data.frame(galaxies)
names(galaxies) <- "km"

ggplot(galaxies, aes(x = km)) +
  geom_density()

x = galaxies$km

gammaUpdate = function(x, mu, sigma, pi){
  pdf = t(sapply(x, dnorm, mean = mu, sd = sigma))
  for(n in 1:length(x)){
    for(k in 1:length(mu)){
      pdf[n, k] = pi[k]*pdf[n, k]
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

mu = c(10, 20, 30)
sigma = c(2, 2, 2)
pi = c(1/3, 1/3, 1/3)
resp = gammaUpdate(x, mu, sigma, pi)
resp[1,]
