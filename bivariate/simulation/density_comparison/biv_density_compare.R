library(ggplot2)
library(gridExtra)
library(MASS)



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
  
  (sapply(data,dmg_int,mu = mu, sigma = sigma, rho = rho,
                   eta = eta,alpha = alpha,l = l) +  
              int_function_R(data,mu,sigma,rho,1/eta*l)
    )
  
}


mu <- c(45,5.5)
sigma <- c(13,1)

rho <- 0.7


N <- 100000
eta <- 0.7
alpha <- 0
theta0 <- c(mu,sigma,rho,alpha)
R_pf <- qnorm(c(.2,.4,.6),mu[1],sigma[1])
T_pf <- qnorm(c(.2,.4,.6),mu[2],sigma[2])


# plot begins
col_nodmg <- rgb(0, 0, 1, 1/4)
col_dmg <- rgb(1, 0, 0, 1/4)



l <- R_pf[1]
set.seed(4)
nodmg<- pl_gen(mu,sigma,rho,eta,0,l,N)
dmg<- pl_gen(mu,sigma,rho,eta,1,l,N)
p1 <- density(data_to_plot(nodmg))     
dmg_data <- data_to_plot(dmg)
p2 <- density(dmg_data[dmg_data>0])   
xlim_range <- range(c(p1$x, p2$x))
ylim_range <- range(c(p1$y, p2$y))
plot( p1, col=col_nodmg , 
      xlim = xlim_range, ylim = ylim_range,
      main = "R20: Density plot ",xlab = "")  # first histogram
lines( p2, col=col_dmg )



x_seq <- seq(from = -5, to = 15, length = 300)
y <- PFY_lik_R(mu, sigma, rho,eta,alpha,l,x_seq)



# normalization to 1:
density_norm <- function(x,y){
  area <- density_calc(x,y)
  factor <- 1/area
  return (y*factor)
}


density_calc <- function(x,y){
  size <- x[2]- x[1]
  
  # na_id <- is.na(undmg_den)
  # x <- x[!na_id]
  # y <- y[!na_id]
  len <- length(y)
  area <- sum(size * (y[1:(len-1)] + y[2:(len)])/2)
  return(area)
}

plot_den <- function(l,dmg_alpha,group_name,RorT,tol = 1e-3){
  if (RorT == "R"){
    x_seq <- seq(from = -1, to = 12, length = 300)
    undmg_y <- PFY_lik_R(mu, sigma, rho,eta,0,l,x_seq)
    dmg_y <- PFY_lik_R(mu, sigma, rho,eta,dmg_alpha,l,x_seq)
  }
  if (RorT == "T"){
    x_seq <- seq(from = -1, to = 90, length = 600)
    undmg_y <- PFY_lik_R(c(mu[2],mu[1]), 
                         c(sigma[2],sigma[1]), rho,eta,0,l,x_seq)
    dmg_y <- PFY_lik_R(c(mu[2],mu[1]), 
                       c(sigma[2],sigma[1]), rho,eta,dmg_alpha,l,x_seq)
  }
  undmg_den <- density_norm(x_seq,undmg_y)
  dmg_den <- density_norm(x_seq,dmg_y)
  undmg_mark <- (undmg_den > tol)
  dmg_mark <- (dmg_den > tol)
  
  x_combined <- c(x_seq[undmg_mark], x_seq[dmg_mark])
  y_combined <- c(undmg_den[undmg_mark], dmg_den[dmg_mark])
  
  xlim <- range(x_combined)
  ylim <- range(y_combined)
  plot(x_seq[undmg_mark],undmg_den[undmg_mark],
       type = 'l',col = col_nodmg,xlab = 'survivors',
       ylab = 'density', main = group_name, 
       xlim = xlim, ylim = ylim)
  lines(x_seq[dmg_mark],dmg_den[dmg_mark],
        type = 'l', col = col_dmg)
}

pdf("density_comparsion.pdf", width = 9.27, height = 11.69)  # A4 size in inches
par(mfrow = c(3, 2), oma = c(6, 1, 2, 1), 
    mar = c(3, 5, 2, 1))  # Adjusted outer and inner margins

# R20 group
l <- R_pf[1]
dmg_alpha <- 0.6
plot_den(l,dmg_alpha,"R20","R")


# T20 group
l <- T_pf[1]
dmg_alpha <- 0.6
plot_den(l,dmg_alpha,"T20","T")

# R40 group
l <- R_pf[2]
dmg_alpha <- 1
plot_den(l,dmg_alpha,"R40","R")

# T40 group
l <- T_pf[2]
dmg_alpha <- 1
plot_den(l,dmg_alpha,"T40","T")

# R60 group
l <- R_pf[3]
dmg_alpha <- 1.4
plot_den(l,dmg_alpha,"R60","R")


# T60 group
l <- T_pf[3]
dmg_alpha <- 1.4
plot_den(l,dmg_alpha,"T60","T")



# Add a common legend at the bottom of all plots
par(fig = c(0, 1, 0, 0.1), oma = c(0, 0, 0, 0), new = TRUE)  # Use the lower part of the figure space for the legend
plot.new()  # Create an empty plot for the legend
legend("center", 
       legend = c(expression("Density of survivor strength generated with true " * bold(alpha)),
                  expression("Density of survivor strength generated with " * bold(alpha) == bold(0))),
       col = c(col_dmg,col_nodmg), pch = 15, pt.cex = 2, 
       bty = "n", horiz = TRUE, inset = 0,
       text.width = 0.4, xjust = 0.5,cex = 1.2) 
dev.off()
