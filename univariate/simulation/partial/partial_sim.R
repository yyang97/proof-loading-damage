# data generation
#source("dmg_func.R")
logit <- function(x){
  return(log(x/(1-x)))
}

expit <- function(x){
  1/(1+exp(-x))
}
dmg_model <- function(y, c, alpha, l) {
  return(ifelse(y > l/c,
                y,
                y - alpha*(l/c - y)
                #- alpha*l/c
  ))
  #return ((alpha+1)*y- alpha*l/c)
  #return ((alpha+1)*y - l )
  
}

dmg_inverse <- function(y,c,alpha,l){
  return ((y + alpha*l/c)/(alpha+1))
  #return ((y +l)/(alpha+1))
}
trunc_pdf_y <- function(x,a,b,mu,sigma){
  truncnorm::dtruncnorm(x,a = a, b = b,mean = mu,
                        sd = sigma)
}



lik_pl_dmg <- function(mu,sigma,alpha, c,l, data){
  l_star <- dmg_inverse(l,c,alpha,l)
  group3 <- data$g3
  ystar_dmg <- group3[group3 < l/c]
  y_dmg <- dmg_inverse(ystar_dmg,c,alpha,l)
  
  lik_dmg <- truncnorm::dtruncnorm(y_dmg,a = l_star,
                                   mean = mu, sd = sigma)
  lik_dmg <- lik_dmg /(alpha+1)
  lik_dmg <- sum(log(lik_dmg))
  # undamaged pieces
  ystar_undmg <- group3[group3 > l/c]
  if(length(ystar_undmg) ==0){
    lik_undmg <- 0
  }
  else{
    lik_undmg <- sum(log(truncnorm::dtruncnorm(ystar_undmg,a = l_star,
                                               mean = mu, sd = sigma)))
  }
  loglik_group3 <- lik_dmg + lik_undmg


  loglik <- loglik_group3 
  return (loglik)
}

lik_pl_undmg <- function(mu,sigma, l, data){
  loglik <- sum(log(truncnorm::dtruncnorm(data$g1,b = l,
                                              mean = mu, sd = sigma)))
  loglik <- loglik + sum(log(truncnorm::dtruncnorm(data$g3,a = l,
                                                           mean = mu, sd = sigma)))
  return (loglik)
}



g3_neglik_fixc <- function(theta,c,R20_1Y,R20_4Y,R5_4Y,R100){
  mu <- theta[1]
  sigma <- theta[2]
  c <- c
  alpha <- theta[3]
  loglik_R20_1Y <- lik_pl_dmg(mu,sigma,alpha, c,l1, R20_1Y)
  loglik_R20_4Y <- lik_pl_undmg(mu,sigma,l1, R20_4Y)
  loglik_R5_4Y <- lik_pl_undmg(mu,sigma,l2, R5_4Y)
  
  # R100 lik
  loglik_R100 <- sum(dnorm(R100, mean = mu, sd = sigma,log = TRUE))
  loglik <- loglik_R20_1Y + loglik_R20_4Y + loglik_R100 + loglik_R5_4Y
  return(-loglik)
}

g3_neglik_evalc <- function(c,theta,R20_1Y,R20_4Y,R5_4Y,R100){
  optimout <- optim(theta,g3_neglik_fixc,c = c,
                    R20_1Y = R20_1Y, R20_4Y = R20_4Y,
                    R5_4Y = R5_4Y,R100 = R100)
  optim_val <- optimout$value
  return(optim_val)
}

g3_neglik_2par <- function(theta,group3,R100){
  mu <- theta[1]
  sigma <- theta[2]
  c <- 1
  alpha <- 0
  l_star <- dmg_inverse(l,c,alpha,l)
  # damaged pieces 
  ystar_dmg <- group3[group3 < l/c]
  y_dmg <- dmg_inverse(ystar_dmg,c,alpha,l)
  
  lik_dmg <- truncnorm::dtruncnorm(y_dmg,a = l_star,
                                   mean = mu, sd = sigma)
  lik_dmg <- lik_dmg /(alpha+1)
  loglik <- sum(log(lik_dmg))
  # undamaged pieces
  ystar_undmg <- group3[group3 > l/c]
  lik_undmg <- sum(log(truncnorm::dtruncnorm(ystar_undmg,a = l_star,
                                             mean = mu, sd = sigma)))
  loglik_R100 <- sum(dnorm(R100,mean = mu,sd = sigma, log = TRUE))


  loglik <- loglik + lik_undmg + loglik_R100 
  
  return(-loglik)
}
library(doParallel)
doParallel::registerDoParallel(cores=8)
N_R100 <- 140
N_R20_1Y <- 300
N_R20_4Y <- 101
N_R5_4Y <- 198

mu <- 51
sigma <-  18


c <-  0.76
alpha <-  5


l1 <- 31.02750
l2 <- 20.685
l <- l2
theta0 <- c(mu,sigma,alpha)

N_rep <- 200

param_res<- foreach(ii=1:N_rep,.combine = 'rbind' ) %dopar% {
  
  set.seed(ii)
  R100 <- rnorm(N_R100, mean = mu, sd = sigma)
  ####### generate R20_1Y
  y <- rnorm(N_R20_1Y, mean = mu, sd = sigma)
  y <- y[order(y)]
  R20_1Y <- list()
  R20_1Y$g1 <- y[y<l1]
  y_star <- dmg_model(y, c, alpha, l1)
  R20_1Y$g3 <- y_star[y_star>l1]
  ######## generate R20_4Y
  y <- rnorm(N_R20_4Y, mean = mu, sd = sigma)
  y <- y[order(y)]
  R20_4Y <- list()
  R20_4Y$g1 <- y[y<l1]
  #y_star <- dmg_model(y, c, alpha, l1)
  R20_4Y$g3 <- y[y>l1]
  ######### generate R5_4Y
  y <- rnorm(N_R5_4Y, mean = mu, sd = sigma)
  y <- y[order(y)]
  R5_4Y <- list()
  R5_4Y$g1 <- y[y<l2]
  #y_star <- dmg_model(y, c, alpha, l1)
  R5_4Y$g3 <- y[y>l2]

  
  #g3_neglik_fixc(theta0,0.77,R20_1Y,R20_4Y,R5_4Y,R100)
 
  
  c_proposed <- seq(from = 0.4, to = 0.9, by = 0.01)
  c_fitvalue <- sapply(c_proposed,g3_neglik_evalc,
                       theta = theta0,
                       R20_1Y = R20_1Y, R20_4Y = R20_4Y,
                       R5_4Y = R5_4Y,R100 = R100)
  c_est <- c_proposed[which.min(c_fitvalue)]
  #c_est
  
  
  optimout <- optim(theta0,g3_neglik_fixc,c = c_est,
                    R20_1Y = R20_1Y, R20_4Y = R20_4Y,
                    R5_4Y = R5_4Y,R100 = R100)
  theta_est <- optimout$par
  hess <- numDeriv::hessian(g3_neglik_fixc,optimout$par,c = c_est,
                            R20_1Y = R20_1Y, R20_4Y = R20_4Y,
                            R5_4Y = R5_4Y,R100 = R100)
  theta_se <- sqrt(diag(solve(hess)))
  c(c_est,theta_est,theta_se)
  #llr1 <- 2*optimout$value
  
  # optimout <- optim(theta0[1:2],g3_neglik_2par,
  #                   group3 = group3, R100 = R100)
  # llr2 <- 2*optimout$value
  # llr_stat <- llr2 - llr1
  # llr_stat
#}
}
c_est <- param_res[,1]
mean(c_est)

theta_est <- param_res[,2:4]
theta_sd <-  param_res[,5:7]
upper <- theta_est + 1.96 * theta_sd
lower <- theta_est - 1.96 * theta_sd 
coverage <- matrix(NA, nrow = dim(theta_est)[1], ncol = dim(theta_est)[2]) 
for (icov in 1:dim(coverage)[1]){   
  coverage[icov,] <- (upper[icov,] - theta0> 0)& (theta0 - lower[icov,] > 0)
}
colMeans(theta_est)
colMeans(coverage)




alpha_id <- theta_est[,3] < 15
colMeans(theta_sd[alpha_id,])
