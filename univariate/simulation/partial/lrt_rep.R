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
  ))

  
}

dmg_inverse <- function(y,c,alpha,l){
  return ((y + alpha*l/c)/(alpha+1))
}
trunc_pdf_y <- function(x,a,b,mu,sigma){
  truncnorm::dtruncnorm(x,a = a, b = b,mean = mu,
                        sd = sigma)
}





g3_neglik_fixc <- function(theta,group3,R100){
  mu <- theta[1]
  sigma <- theta[2]
  c <-  theta[3]
  alpha <- theta[4]
  l_star <- dmg_inverse(l,c,alpha,l)
  # damaged pieces 
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
  
  # R100 lik
  loglik_R100 <- sum(dnorm(R100, mean = mu, sd = sigma,log = TRUE))
  loglik <- loglik_group3 + loglik_R100
  return(-loglik)
}

g3_neglik_evalc <- function(c,theta,group3,R100){
  optimout <- optim(theta,g3_neglik_fixc,c = c,
                    group3 = group3, R100 = R100)
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
N_R20_4Y <- 300
N_R5_4Y <- 300

mu <- 51
sigma <-  18


c <-  0.76
alpha <-  5




l1 <- 31.02750
l2 <- 20.685
l <- l1
theta0 <- c(mu,sigma,logit(c),alpha)

N_rep <- 200

lrt_res <- matrix(NA, ncol = 2,nrow = N_rep)
for(ii in 1:N_rep){
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
  # y <- rnorm(N_R20_4Y, mean = mu, sd = sigma)
  # y <- y[order(y)]
  # R20_4Y <- list()
  # R20_4Y$g1 <- y[y<l1]
  # y_star <- dmg_model(y, c, alpha[2], l1)
  # R20_4Y$g3 <-  y_star[y_star>l1]
  ######### generate R5_4Y
  # y <- rnorm(N_R5_4Y, mean = mu, sd = sigma)
  # y <- y[order(y)]
  # R5_4Y <- list()
  # R5_4Y$g1 <- y[y<l2]
  # y_star <- dmg_model(y, c, alpha[3], l2)
  # R5_4Y$g3 <- y_star[y_star>l2]
  
  
  group3 <- R20_1Y$g3
  #c_est
  
  
  optimout <- optim(theta0,g3_neglik_fixc,
                    group3 = group3, R100 = R100)
  llr1 <- 2*optimout$value
  theta_est <- optimout$par
  optimout <- optim(theta0[1:2],g3_neglik_2par,
                    group3 = group3, R100 = R100)
  llr2 <- 2*optimout$value
  llr_stat <- llr2 - llr1
  theta_sim <- optimout$par
  
  lrt_sim<- foreach(jj=1:1000,.combine = 'rbind' ) %dopar% {
    set.seed(jj)
    R100 <- rnorm(N_R100, mean = mu, sd = sigma)
    
    y <- rnorm(N_R20_1Y, mean = theta_sim[1], sd = theta_sim[2])
    y <- y[order(y)]
    R20_1Y <- list()
    R20_1Y$g1 <- y[y<l1]
    y_star <- dmg_model(y, expit(theta_sim[3]), theta_sim[4], l1)
    R20_1Y$g3 <- y_star[y_star>l1]
    
    group3 <- R20_1Y$g3
    optimout <- optim(theta0,g3_neglik_fixc,
                      group3 = group3, R100 = R100)
    llr1 <- 2*optimout$value
    theta_est <- optimout$par
    optimout <- optim(theta0[1:2],g3_neglik_2par,
                      group3 = group3, R100 = R100)
    llr2 <- 2*optimout$value
    llr_stat <- llr2 - llr1
    llr_stat
  }
  pval <- 1-ecdf(lrt_sim)(llr_stat)
  lrt_res[ii,] <- c(pval,theta_est)
}


mean(lrt_res[ii,1] < 0.05)


