# in the bivariate real data analysis
# this file calculate the p-values for each proof-loading group
library(ggplot2)
library(doParallel)
library(MASS)
source("real_data_preprocessing.R")
source("lrt_func.R")
doParallel::registerDoParallel(cores=180)

mu <- c(45.13,5.51)
sigma <- c(12.92,1.08)

rho <- 0.7

alpha <- 0

# eta is thresh 
eta <- .7



# specify the proof-loading group
group <- R20_data
group_pl <- R_pf[1]
group_name <- "R"
N <- dim(R20_data)[1]
N_shoulder <- length(R100_data)
shoulder_group <- list(R100_data,T100_data)

# calculate H_0 
optim_checkR <- optim(c(theta0[1:5]),single_alpha0,
                      group = group,group_pl = group_pl,
                      group_name = group_name,
                      shoulder_group = shoulder_group)

llr_checkR <- 2*optim_checkR$value
theta_h0  <- optim_checkR$par

# calculate H_A
optimout <- optim(c(mu,sigma,0.5,0.5,.1),single_fitalpha,
                  method = "L-BFGS-B",
                  lower = c(30,0.1,0.1,0.1,0.01,0.01,-2),
                  upper = c(70,40,30,10,0.98,0.98,10),
                  group = group,group_pl = group_pl,
                  group_name = group_name,
                  shoulder_group = shoulder_group)
theta_init <- optimout$par

llr_full <- 2*optimout$value
llr_stat <- llr_checkR- llr_full

# theta_h0 <- c(45.726222, 5.479856, 12.878392, 1.054855, 1.247164)
# llr_stat <- 0.1267




# find the numerical distribution under H_0
param_res <- foreach(ii=1:1000, .combine=rbind) %dopar% {
  
  set.seed(ii+1000)
  if (group_name == "R"){
    pl_data <- pl_gen(theta_h0[1:2],
                       theta_h0[3:4],
                       expit(theta_h0[5]),
                       1,0,group_pl,N)
  }
  if (group_name == "T"){
    pl_data <- pl_gen(theta_h0[2:1],
                      theta_h0[4:3],
                      expit(theta_h0[5]),
                      1,0,group_pl,N)
  }
  R100_data <- rnorm(N_shoulder,mean = theta_h0[1],
                     sd = theta_h0[3])
  T100_data <- rnorm(N_shoulder,mean = theta_h0[2],
                     sd = theta_h0[4])
  
  shoulder_group <- list(R100_data,T100_data)
  #llr_res <- llr_sim(pl_data,shoulder_group,group_pl,group_name)
  llr_res <- tryCatch({
    llr_sim(pl_data, shoulder_group, group_pl, group_name)
  }, error = function(e) {
    20  # Return 20 if an error occurs
  })
  llr_res
}

print(1-ecdf(param_res)(llr_stat))
