library(StatProg)


DataFun <- function(n, lambda) {

    # independent variable
    x <- runif(n, min = 0, max = 2)

    # standard deviation
    for (i in seq_len(n)) {
        s[i] <- exp(x[i]*lambda)
    }

    # error term
    epsilon <- rnorm(n, mean = rep(0, 5), sd = s)

    beta <- 2

    # dependent variable
    for (i in seq_len(n)) {
        y[i] <- beta*x[i] + epsilon[i]
    }

    # container matrix
    mat <- matrix(data = 0, ncol = 3, nrow = n)

    # creating matrix of indepedent and dependent variables
    for (i in seq_len(n)) {
        mat[i, 1] <- y[i]
        mat[i, 2] <- x[i]
        mat[i, 3] <- epsilon[i]
    }
    return(mat)
}

DataFun(5, 2)

rnorm(2, mean = c(0,10), sd = c(1000, 1))

# OLS

testData <- cbind( c(0.62, 0.18, 3.92, 0.80, -5.15),
c(0.44, 1.49, 0.69, 0.13, 1.90) )

x <- cbind(rep(1, 5), testData[,2])
y <- testData[,1]

ols <- solve(t(x) %*% x) %*% t(x) %*% y
ols[2]

y_hat <- rep(0,5)

for (i in 1:5) {
    y_hat[i] <- 2.94 + x[i,2]*(-3.09)
}



#ZZZZZ