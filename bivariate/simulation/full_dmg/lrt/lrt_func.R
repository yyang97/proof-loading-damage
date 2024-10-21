logit <- function(x){
  return(log(x/(1-x)))
}

expit <- function(x){
  1/(1+exp(-x))
}


dmg_model <- function(samples,eta,alpha,l,mu){
  #samples[,2] - alpha*(l/eta - samples[,1])
  ifelse(l > samples[,1]*eta, 
         samples[,2] - mu[2]/mu[1]*alpha*(l/eta - samples[,1]), # damage
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
  biv_pdf(x,y + mu[2]/mu[1]*alpha*(l/eta - x),mu,sigma,rho)
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
  res[id,2] <- dmg_model(samples[id,],eta,alpha,l,mu)
  res[id,3] <- 0
  return(res)
}


nlogpost_plot <- function(theta){
  nlogpost(theta,data)
}


R100_data <- rnorm(2*N,mean = mu[1],sd = sigma[1])
T100_data <- rnorm(2*N,mean = mu[2],sd = sigma[2])

shoulder_group <- list(R100_data,T100_data)
alpha_checkR <- function(theta,R_group,T_group,
                          R_pl,T_pl,shoulder_group){
  mu <- theta[1:2]
  sigma <- theta[3:4]
  rho <- expit(theta[5])
  eta <- expit(theta[6])
  alpha <- theta[7]
  lik <- PFY_lik_R(mu,sigma,rho,1,0,R_pl,R_group) +
    PFY_lik_R(c(mu[2],mu[1]),c(sigma[2],sigma[1]),rho,eta,
              alpha,T_pl,T_group) +
    sum(dnorm(shoulder_group[[1]],
              mean = mu[1],
              sd = sigma[1], log = T)) + 
    sum(dnorm(shoulder_group[[2]],
              mean = mu[2],
              sd = sigma[2], log = T))
  if (is.infinite(lik)){
    lik <- -10000
  }
  return(-1*lik)
}

alpha_checkT <- function(theta,R_group,T_group,
                         R_pl,T_pl,shoulder_group){
  mu <- theta[1:2]
  sigma <- theta[3:4]
  rho <- expit(theta[5])
  eta <- expit(theta[6])
  alpha <- theta[7]
  lik <- PFY_lik_R(mu,sigma,rho,eta,alpha,R_pl,R_group) +
    PFY_lik_R(c(mu[2],mu[1]),c(sigma[2],sigma[1]),rho,1,
              0,T_pl,T_group) +
    sum(dnorm(shoulder_group[[1]],
              mean = mu[1],
              sd = sigma[1], log = T)) + 
    sum(dnorm(shoulder_group[[2]],
              mean = mu[2],
              sd = sigma[2], log = T))
  if (is.infinite(lik)){
    lik <- -10000
  }
  return(-1*lik)
}



single_fitalpha <- function(theta,group,group_pl,
                            group_name,shoulder_group){
  mu <- theta[1:2]
  sigma <- theta[3:4]
  # rho <- expit(theta[5])
  # eta <- expit(theta[6])
  rho <-theta[5]
  eta <- theta[6]
  alpha <- theta[7]
  if (group_name == "R"){
    lik <- PFY_lik_R(mu,sigma,rho,eta,
                     alpha,group_pl,group)
  }
  if (group_name == "T"){
    lik <- PFY_lik_R(c(mu[2],mu[1]),c(sigma[2],sigma[1]),
                     rho,eta,alpha,group_pl,group)
  }
  lik <- lik +  sum(dnorm(shoulder_group[[1]],
              mean = mu[1],
              sd = sigma[1], log = T)) + 
    sum(dnorm(shoulder_group[[2]],
              mean = mu[2],
              sd = sigma[2], log = T))
  if (is.infinite(lik)){
    lik <- -10000
  }
  return(-1*lik)
}





single_alpha0 <- function(theta,group,group_pl,
                            group_name,shoulder_group){
  mu <- theta[1:2]
  sigma <- theta[3:4]
  rho <- expit(theta[5])
  # eta <- expit(theta[6])
  # alpha <- theta[7]
  if (group_name == "R"){
    lik <- PFY_lik_R(mu,sigma,rho,1,
                     0,group_pl,group)
  }
  if (group_name == "T"){
    lik <- PFY_lik_R(c(mu[2],mu[1]),c(sigma[2],sigma[1]),
                     rho,1,0,group_pl,group)
  }
  lik <- lik +  sum(dnorm(shoulder_group[[1]],
                          mean = mu[1],
                          sd = sigma[1], log = T)) + 
    sum(dnorm(shoulder_group[[2]],
              mean = mu[2],
              sd = sigma[2], log = T))
  if (is.infinite(lik)){
    lik <- -10000
  }
  return(-1*lik)
}

llr_sim <- function(data,shoulder_group){
  optim_checkR <- optim(theta0[1:5],single_alpha0,
                        group = data,group_pl = group_pl,
                        group_name = "R",
                        shoulder_group = shoulder_group)
  
  llr_checkR <- 2*optim_checkR$value
  
  optimout <- optim(theta0,single_fitalpha,
                    method = "L-BFGS-B",
                    lower = c(30,0.1,0.1,0.1,0.03,0.03,-2),
                    upper = c(70,40,30,10,0.97,0.97,10),
                    group = data,group_pl = group_pl,
                    group_name = "R",
                    shoulder_group = shoulder_group)
  llr_full <- 2*optimout$value
  llr_stat <- llr_checkR- llr_full
  llr_stat
}
