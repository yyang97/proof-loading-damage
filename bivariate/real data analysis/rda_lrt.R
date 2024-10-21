library(ggplot2)
library(doParallel)
source("real_data_preprocessing.R")
source("lrt_func.R")
doParallel::registerDoParallel(cores=8)

mu <- c(45.13,5.51)
sigma <- c(12.92,1.08)

rho <- 0.7

alpha <- 0

# eta is thresh 
eta <- .8
theta0 <- c(mu,sigma,logit(rho),logit(eta),alpha)

group <- R40_data
group_pl <- R_pf[2]
group_name <- "R"
shoulder_group <- list(R100_data,T100_data)


optim_checkR <- optim(theta0[1:5],single_alpha0,
                      group = group,group_pl = group_pl,
                      group_name = group_name,
                      shoulder_group = shoulder_group)

llr_checkR <- 2*optim_checkR$value


optimout <- optim(theta0,single_fitalpha,
                  method = "L-BFGS-B",
                  lower = c(30,0.1,0.1,0.1,0.01,0.01,-2),
                  upper = c(70,40,30,10,0.99,0.99,10),
                  group = group,group_pl = group_pl,
                  group_name = group_name,
                  shoulder_group = shoulder_group)
theta_init <- optimout$par

llr_full <- 2*optimout$value
llr_stat <- llr_checkR- llr_full
llr_stat


optim_checkR$par


param_res <- foreach(ii=1:1000, .combine=rbind) %dopar% {
  
  set.seed(ii+1000)
  R40_data <- pl_gen(optim_checkR$par[1:2],
                     optim_checkR$par[3:4],expit(optim_checkR$par[5]),
                     eta,alpha,R_pf[2],N)
  R100_data <- rnorm(2*N,mean = optim_checkR$par[1],
                     sd = optim_checkR$par[3])
  T100_data <- rnorm(2*N,mean = optim_checkR$par[2],
                     sd = optim_checkR$par[4])
  
  shoulder_group <- list(R100_data,T100_data)
  llr_res <- llr_sim(R40_data,shoulder_group,group_pl,group_name)
  llr_res
}


1-ecdf(param_res)(llr_stat)

