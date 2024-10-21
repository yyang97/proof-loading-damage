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





g3_neglik_fixc <- function(theta,c,group3,R100){
  mu <- theta[1]
  sigma <- theta[2]
  c <- c
  alpha <- theta[3]
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
  # group1 lik

  
  # R100 lik
  loglik_R100 <- sum(dnorm(R100, mean = mu, sd = sigma,log = TRUE))
  loglik <- loglik_group3  + loglik_R100
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
  # group1 lik

  loglik <- loglik + lik_undmg + loglik_R100 
  
  return(-loglik)
}
library(doParallel)
doParallel::registerDoParallel(cores=8)
N_R100 <- 198
N <- 198
mu <- 48.4
sigma <-  18.29


c <-  0.76
alpha <-  0


l1 <- 31.02750
l2 <- 20.685
l <- l2
theta0 <- c(mu,sigma,alpha)

N_rep <- 1000

lrt_res<- foreach(ii=1:N_rep,.combine = 'c' ) %dopar% {
  
  set.seed(ii)
  y <- rnorm(N, mean = mu, sd = sigma)
  y <- y[order(y)]
  group1 <- y[y<l]
  y_star <- dmg_model(y, c, alpha, l)
  group3 <- y_star[y_star>l]
  R100 <- rnorm(N_R100, mean = mu, sd = sigma)
      
      
  c_proposed <- seq(from = 0.2, to = 0.9, by = 0.01)
  c_fitvalue <- sapply(c_proposed,g3_neglik_evalc,
                       theta = theta0, group3 = group3,
                       R100 = R100,group1 = group1)
  c_est <- c_proposed[which.min(c_fitvalue)]
  #c_est
    
    
  optimout <- optim(theta0,g3_neglik_fixc,c = c_est,
                    group3 = group3, R100 = R100,
                    group1 = group1)
  llr1 <- 2*optimout$value
    
  optimout <- optim(theta0[1:2],g3_neglik_2par,
                    group3 = group3, R100 = R100,group1 = group1)
  llr2 <- 2*optimout$value
  llr_stat <- llr2 - llr1
  llr_stat
}


ecdf(lrt_res)


1-ecdf(lrt_res)(4.346443)
