logit <- function(x){
  return(log(x/(1-x)))
}

expit <- function(x){
  1/(1+exp(-x))
}


dmg_model <- function(samples,eta,alpha,l){
  #samples[,2] - alpha*(l/eta - samples[,1])
  ifelse(l > samples[,1]*eta, 
         samples[,2] - alpha*(l/eta - samples[,1]), # damage
         samples[,2]) # undamage
}
biv_pdf <- function(x,y,mu,sigma,rho, log = FALSE){
  
  rho_part <- (1- rho^2)
  if (rho_part <=0){
    rho_part <- 0.000001
  }
  
  factor_part <- 2*pi*sigma[1]*sigma[2]*sqrt(rho_part)
  exp_part <- ((x - mu[1])/sigma[1])^2 + ((y - mu[2])/sigma[2])^2
  exp_part <- exp_part - 2*rho*(x - mu[1])*(y - mu[2])/(sigma[1]*sigma[2])
  exp_part <- -1*exp_part/(2*rho_part )
  return(exp(exp_part)/factor_part)
}



# pdf of multivaraite normal 
dmg_int_inner <- function(x,y,mu,sigma,rho,eta,alpha,l){
  sd_x <- sigma[1]
  sd_y <- sigma[2]
  Sigma <- matrix(c(sd_x^2, sd_x*sd_y*rho, sd_x*sd_y*rho, sd_y^2),2)
  biv_pdf(x,y + alpha*(l/eta - x),mu,sigma,rho)
}
dmg_int <- function(y,mu,sigma,rho,eta,alpha,l){
  integrate(f = dmg_int_inner, lower = l, upper = l/eta, y =y,
            mu = mu, sigma = sigma, rho= rho, eta = eta,
            alpha = alpha, l = l,abs.tol = 1e-5 )$value
}


int_function_R <- function(ystar, mu,sigma, rho,l){
  a_l <- (l - mu[1]-rho*(sigma[1]/sigma[2])*(ystar-mu[2]))/(sigma[1]*sqrt(1-rho^2));
  dnorm(ystar, mean = mu[2], sd = sigma[2],log = F)*pnorm(a_l,lower.tail = F, log.p = F)
}


PFY_lik_R <- function(mu, sigma, rho,eta,alpha,l,data){
  
  (sum(dnorm(data[data[,3] == 1,1],
             mean = mu[1],
             sd = sigma[1],
             log = T)
  ) +
    sum(log(sapply(data[data[,3] == 0,2],dmg_int,mu = mu, sigma = sigma, rho = rho,
                   eta = eta,alpha = alpha,l = l) +  
              int_function_R(data[data[,3] == 0,2],mu,sigma,rho,1/eta*l)
    ))
  )
}


pl_gen <- function(mu,sigma,rho,eta,alpha,l,N){
  sd_x <- sigma[1]
  sd_y <- sigma[2]
  Sigma <- matrix(c(sd_x^2, sd_x*sd_y*rho, sd_x*sd_y*rho, sd_y^2),2)
  bvn1 <- mvrnorm(N, mu = c(mu[1],mu[2]), Sigma = Sigma ) # from MASS package
  samples <- bvn1
  res <- matrix(0, nrow = dim(samples)[1],
                ncol = dim(samples)[2]+1)
  # fail in the proof loading
  id <- which(samples[,1] < l)
  res[id,1] <- samples[id,1]
  res[id,3] <- 1
  # survive in the proof loading 
  id <- which(samples[,1] >= l)
  res[id,2] <- dmg_model(samples[id,],eta,alpha,l)
  res[id,3] <- 0
  return(res)
}

nlogpost_full <- function(theta,data){
  mu <- theta[1:2]
  sigma <- theta[3:4]
  rho <- expit(theta[5])
  eta <- 1
  alpha <- rep(0,6)
  lik <- PFY_lik_R(mu,sigma,rho,eta,alpha[1],R_pf[1],data[[1]]) +
    PFY_lik_R(mu,sigma,rho,eta,alpha[2],R_pf[2],data[[2]]) +
    PFY_lik_R(mu,sigma,rho,eta,alpha[3],R_pf[3],data[[3]]) +
    PFY_lik_R(c(mu[2],mu[1]),c(sigma[2],sigma[1]),rho,eta,alpha[4],T_pf[1],data[[4]]) +
    PFY_lik_R(c(mu[2],mu[1]),c(sigma[2],sigma[1]),rho,eta,alpha[5],T_pf[2],data[[5]]) +
    PFY_lik_R(c(mu[2],mu[1]),c(sigma[2],sigma[1]),rho,eta,alpha[6],T_pf[3],data[[6]]) +
    sum(dnorm(data[[7]],
              mean = mu[1],
              sd = sigma[1], log = T)) + 
    sum(dnorm(data[[8]],
              mean = mu[2],
              sd = sigma[2], log = T))
  if (is.infinite(lik)){
    lik <- -10000
  }
  return(-1*lik)
}