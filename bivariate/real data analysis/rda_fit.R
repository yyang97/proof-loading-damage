library(ggplot2)

source("real_data_preprocessing.R")
source("sub_func.R")

mu <- c(45.13,5.51)
sigma <- c(12.92,1.08)

rho <- 0.7

theta0 <- c(mu,sigma,rho)
data <- list(R20_data,R40_data,R60_data,
             T20_data,T40_data,T60_data,
             R100_data,T100_data)

optimresult <- optim(theta0,nlogpost_full,
                     data = data)

optimresult 



hess <- numDeriv::hessian(nlogpost_full,optimresult$par,
                         method.args=list(r = 4),
                         data = data)

theta_sd <- sqrt(diag(chol2inv(chol(hess))))
theta_est <- optimresult$par

theta_sd 
