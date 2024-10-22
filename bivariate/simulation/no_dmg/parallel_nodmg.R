
require(MASS)
require(optimCheck)
library(doParallel)
doParallel::registerDoParallel(cores=8)
source("sub_func.R")


mu <- c(45.13,5.51)
sigma <- c(12.92,1.08)

rho <- 0.7


N <- 87


alpha <- rep(0,6)
# eta is thresh 
eta <- 0.7

theta0 <- c(mu,sigma,rho,eta,alpha)

##------proof loading-----####

R_pf <- qnorm(c(.2,.4,.6),mu[1],sigma[1])
T_pf <- qnorm(c(.2,.4,.6),mu[2],sigma[2])
theta_init <- c(mu,sigma,logit(rho))


N_rep <- 200
param_res <- foreach(ii=1:N_rep, .combine=rbind) %dopar% {

set.seed(ii)
# R group
R20_data <- pl_gen(mu,sigma,rho,eta,alpha[1],R_pf[1],N)
R40_data <- pl_gen(mu,sigma,rho,eta,alpha[2],R_pf[2],N)
R60_data <- pl_gen(mu,sigma,rho,eta,alpha[3],R_pf[3],N)
# T group
T20_data <- pl_gen(c(mu[2],mu[1]),c(sigma[2],sigma[1]),
                   rho,eta,alpha[4],T_pf[1],N)
T40_data <- pl_gen(c(mu[2],mu[1]),c(sigma[2],sigma[1]),
                   rho,eta,alpha[5],T_pf[2],N)
T60_data <- pl_gen(c(mu[2],mu[1]),c(sigma[2],sigma[1]),
                   rho,eta,alpha[6],T_pf[3],N)

##-----T100-----######

R100_data <- rnorm(2*N,mean = mu[1],sd = sigma[1])
T100_data <- rnorm(2*N,mean = mu[2],sd = sigma[2])


data <- list(R20_data,R40_data,R60_data,
             T20_data,T40_data,T60_data,
             R100_data,T100_data)

optimresult <- optim(theta_init,nlogpost_full,
                     data = data)

hess = numDeriv::hessian(nlogpost_full,optimresult$par,
                         method.args=list(r = 4),
                         data = data)

theta_sd <- sqrt(diag(chol2inv(chol(hess))))
theta_est <- optimresult$par

c(theta_est,theta_sd)

# signif(cbind(theta_init,theta_est,theta_sd, 
#              theta_est - 1.96 * theta_sd,
#              theta_est + 1.96 * theta_sd),3)
}


theta_est <- param_res[,1:5]
theta_sd <- param_res[,6:10]
theta_true <- theta_init


upper <- theta_est + 1.96 * theta_sd
lower <- theta_est - 1.96 * theta_sd
coverage <- matrix(NA, nrow = dim(theta_est)[1], ncol = dim(theta_est)[2]) 
for (icov in 1:dim(coverage)[1]){   
  coverage[icov,] <- (upper[icov,] - theta_init> 0)& 
    (theta_init - lower[icov,] > 0)
}
colMeans(coverage)

