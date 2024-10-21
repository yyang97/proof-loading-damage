# data generation
logit <- function(x){
  return(log(x/(1-x)))
}

expit <- function(x){
  1/(1+exp(-x))
}

myqq <- function(x, y, ...) {
  rg <- range(x, y, na.rm = T)
  qqplot(x, y, xlim = rg, ylim = rg, ...)
  abline(0,1,col = 'red')
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


g3_neglik_fixc <- function(theta,c,group3,R100){
  mu <- theta[1]
  sigma <- theta[2]
  c <- c
  alpha <- theta[3]
  l_star <- dmg_inverse(l,c,alpha,l)
  # damaged pieces 
  ystar_dmg <-group3[group3 < l/c]
  y_dmg <- dmg_inverse(ystar_dmg,c,alpha,l)
  
  lik_dmg <- truncnorm::dtruncnorm(y_dmg,a = l_star,
                                   mean = mu, sd = sigma)
  lik_dmg <- lik_dmg /(alpha+1)
  lik_dmg <- sum(log(lik_dmg))
  # undamaged pieces
  ystar_undmg <- group3[group3 > l/c]
  if (length(ystar_undmg) == 0){
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
  loglik <- loglik_group3 +  loglik_R100
  
  lik_R20_4Y <- sum(log(truncnorm::dtruncnorm(R20_4Y,b = l1,
                                              mean = mu, sd = sigma)))
  lik_R20_4Y <- lik_R20_4Y + sum(log(truncnorm::dtruncnorm(R20_4Y_R100,a = l1,
                                              mean = mu, sd = sigma)))
  lik_R5_4Y <- sum(log(truncnorm::dtruncnorm(R5,b = l2,
                                              mean = mu, sd = sigma)))
  lik_R5_4Y <- lik_R5_4Y + sum(log(truncnorm::dtruncnorm(R5R100,a = l2,
                                                           mean = mu, sd = sigma)))
  loglik <- loglik + lik_R20_4Y + lik_R5_4Y
  return(-loglik)
}

g3_neglik_evalc <- function(c,theta,group3,R100){
  optimout <- optim(theta,g3_neglik_fixc,c = c,
                    group3 = group3, R100 = R100)
  optim_val <- optimout$value
  return(optim_val)
}


bending <- read.csv("bending-pl-4groups.csv", header = F)
group_names <- table(bending[, 2])
names(group_names)
bending[, 1] <- bending[, 1] / 1000 * 6.895
R20 <- bending[which(bending[, 2] == "R20"), 1]
R20R100 <- bending[which(bending[, 2] == "R20R100"), 1]
R20_4Y <- bending[which(bending[, 2] == "R20_4Y"), 1]
R20_4Y_R100 <- bending[which(bending[, 2] == "R20_4Y_R100"), 1]
R5R100 <- bending[which(bending[, 2] == "R5R100"), 1]
R5 <- bending[which(bending[, 2] == "R5"), 1]
R100 <- bending[which(bending[, 2] == "R100"), 1]

group3 <- R20R100
group1 <- R20
mu <- 48
sigma <-  19


c <-  0.62
alpha <-  3.65
l1 <- 31.02750
l2 <- 20.685

l <- l1


theta_init <- c(50,19,2)
g3_neglik_fixc(theta_init,c,group3,R100)



g3_neglik_evalc(c,theta_init,group3,R100)

c_proposed <- seq(from = 0.4, to = 0.9, by = 0.01)
c_fitvalue <- sapply(c_proposed,g3_neglik_evalc,
                     theta = theta_init, group3 = group3,
                     R100 = R100)
c_proposed[which.min(c_fitvalue)]

c_merge <- cbind(c_proposed,c_fitvalue)

c_id <- order(c_merge[,2]) 
c_merge[c_id,]




g3_neglik_3par <- function(theta){
  mu <- theta[1]
  sigma <- theta[2]
  c <- 0.76
  alpha <- theta[3]
  
  
  group3 <- R20R100

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
  lik_undmg <- sum(log(truncnorm::dtruncnorm(ystar_undmg,a = l_star,
                                             mean = mu, sd = sigma)))
  
  loglik_group3 <- lik_dmg + lik_undmg
  # group1 lik
  
  # R100 lik
  loglik_R100 <- sum(dnorm(R100, mean = mu, sd = sigma,log = TRUE))
  loglik <- loglik_group3 + loglik_R100
  lik_R20_4Y <- sum(log(truncnorm::dtruncnorm(R20_4Y,b = l1,
                                              mean = mu, sd = sigma)))
  lik_R20_4Y <- lik_R20_4Y + sum(log(truncnorm::dtruncnorm(R20_4Y_R100,a = l1,
                                                           mean = mu, sd = sigma)))
  lik_R5_4Y <- sum(log(truncnorm::dtruncnorm(R5,b = l2,
                                             mean = mu, sd = sigma)))
  lik_R5_4Y <- lik_R5_4Y + sum(log(truncnorm::dtruncnorm(R5R100,a = l2,
                                                         mean = mu, sd = sigma)))
  loglik <- loglik + lik_R20_4Y + lik_R5_4Y
  return(-loglik)
}

theta_init <- c(50,20,2)
g3_neglik_3par(theta_init)
optimout <- optim(theta_init,g3_neglik_3par)
optimout$par
optimCheck::optim_proj(optimout$par,g3_neglik_3par,xrng = .2)

hess <- numDeriv::hessian(g3_neglik_3par,optimout$par)
sqrt(diag(solve(hess)))
  

optimout$value
