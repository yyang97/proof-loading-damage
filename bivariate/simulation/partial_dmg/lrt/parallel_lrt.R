require(MASS)
library(doParallel)
source("lrt_func.R")
doParallel::registerDoParallel(cores=150)

mu <- c(45,5.5)
sigma <- c(13,1)

rho <- 0.7


N <- 87


alpha <- 0

# eta is thresh 
eta <- 0.7

theta0 <- c(mu,sigma,logit(rho),logit(eta),alpha)

##------proof loading-----####

R_pf <- qnorm(c(.2,.4,.6),mu[1],sigma[1])
T_pf <- qnorm(c(.2,.4,.6),mu[2],sigma[2])
N_rep <- 200
p_value_res <- rep(0,N_rep)
for (jj in 1:N_rep){

  set.seed(jj)
  # R group
  R40_data <- pl_gen(mu,sigma,rho,eta,alpha,R_pf[2],N)
  group_pl <- R_pf[2]
  # T group
  # T40_data <- pl_gen(c(mu[2],mu[1]),c(sigma[2],sigma[1]),
  #                    rho,eta,alpha,T_pf[3],N)
  # group_pl <- T_pf[3]
  
  ##-----T100-----######
  
  R100_data <- rnorm(2*N,mean = mu[1],sd = sigma[1])
  T100_data <- rnorm(2*N,mean = mu[2],sd = sigma[2])
  
  shoulder_group <- list(R100_data,T100_data)
  # full data
  
  
  
  
  # check R group
  
  
  
  optim_checkR <- optim(theta0[1:5],single_alpha0,
                        group = R40_data,group_pl = group_pl,
                        group_name = "R",
                        shoulder_group = shoulder_group)
  
  llr_checkR <- 2*optim_checkR$value
  
  
  optimout <- optim(theta0,single_fitalpha,
                    method = "L-BFGS-B",
                    lower = c(30,0.1,0.1,0.1,0.01,0.01,-2),
                    upper = c(70,40,30,10,0.99,0.99,10),
                    group = R40_data,group_pl = group_pl,
                    group_name = "R",
                    shoulder_group = shoulder_group)
  # optimout <- optim(optimout$par,single_fitalpha,
  #                   group = R40_data,group_pl = group_pl,
  #                   group_name = "R",
  #                   shoulder_group = shoulder_group)
  theta_init <- optimout$par
  
  llr_full <- 2*optimout$value
  llr_stat <- llr_checkR- llr_full
  
  
  
  
  
  param_res <- foreach(ii=1:1000, .combine=rbind) %dopar% {
  
    set.seed(ii+1000)
    R40_data <- pl_gen(optim_checkR$par[1:2],
                       optim_checkR$par[3:4],expit(optim_checkR$par[5]),
                       eta,alpha,R_pf[2],N)
    
    llr_res <- llr_sim(R40_data,shoulder_group)
    llr_res
  }
  
  
  p_value_res[jj] <- ecdf(param_res)(llr_stat)
  if (jj %% 10 == 0){
    print(jj)
  }
}

saveRDS(p_value_res,file = "R20_alpha0_lrt.rds")