library(ggplot2)
library(gridExtra)

source("real_data_preprocessing.R")


myqq <- function(x, y, ...) {
  rg <- range(x, y, na.rm = T)
  qqplot(x, y, xlim = rg, ylim = rg, ...)
  abline(0,1,col = 'black')
}


dmg_model <- function(samples,eta,alpha,l,mu){
  #samples[,2] - alpha*(l/eta - samples[,1])
  ifelse(l > samples[,1]*eta, 
         samples[,2] - mu[2]/mu[1]*alpha*(l/eta - samples[,1]), # damage
         samples[,2]) # undamage
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
data_to_plot <- function(data){
  return(data[data[,3] == 0,2 ])
}

mu <- c(45.13,5.51)
sigma <- c(12.92,1.08)

rho <- 0.7


N <- 87
eta <- 1
alpha <- 0
theta0 <- c(mu,sigma,rho,alpha)

pdf("qqplots_combined.pdf", width = 8.27, height = 11.69)  # A4 size in inches
par(mfrow = c(3, 2))
# R20
l <- R_pf[1]
set.seed(4)
R20_fit <- pl_gen(mu,sigma,rho,eta,alpha,l,N)
p1 <- myqq(data_to_plot(R20_data),data_to_plot(R20_fit),
     xlab = "empirical quantile",
     ylab = "predicted quantile",
     main = "R20: Q-Q plot")

# R40
l <- R_pf[2]
R40_fit <- pl_gen(mu,sigma,rho,eta,alpha,l,N)
p2 <- myqq(data_to_plot(R40_data),data_to_plot(R40_fit),
     xlab = "empirical quantile",
     ylab = "predicted quantile",
     main = "R40: Q-Q plot")

# R60
l <- R_pf[3]
R60_fit <- pl_gen(mu,sigma,rho,eta,alpha,l,N)
p3 <- myqq(data_to_plot(R60_data),data_to_plot(R60_fit),
           xlab = "empirical quantile",
           ylab = "predicted quantile",
           main = "R60: Q-Q plot")


# T20
l <- T_pf[1]
T20_fit <- pl_gen(c(mu[2],mu[1]),c(sigma[2],sigma[1]),
                  rho,eta,alpha,l,N)
p4 <- myqq(data_to_plot(T20_data),data_to_plot(T20_fit),
     xlab = "empirical quantile",
     ylab = "predicted quantile",
     main = "T20: Q-Q plot")


# T40
l <- T_pf[2]
T40_fit <- pl_gen(c(mu[2],mu[1]),c(sigma[2],sigma[1]),
                  rho,eta,alpha,l,N)
p5 <- myqq(data_to_plot(T40_data),data_to_plot(T40_fit),
     xlab = "empirical quantile",
     ylab = "predicted quantile",
     main = "T40: Q-Q plot")



# T60
l <- T_pf[3]
T60_fit <- pl_gen(c(mu[2],mu[1]),c(sigma[2],sigma[1]),
                  rho,eta,alpha,l,N)
p6 <- myqq(data_to_plot(T60_data),data_to_plot(T60_fit),
     xlab = "empirical quantile",
     ylab = "predicted quantile",
     main = "T60: Q-Q plot")
par(mfrow = c(1, 1))

dev.off()
