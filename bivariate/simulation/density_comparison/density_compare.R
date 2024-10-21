library(ggplot2)
library(gridExtra)




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
eta <- 0.7
alpha <- 0
theta0 <- c(mu,sigma,rho,alpha)
R_pf <- qnorm(c(.2,.4,.6),mu[1],sigma[1])
T_pf <- qnorm(c(.2,.4,.6),mu[2],sigma[2])


# plot begins
col_nodmg <- rgb(0, 0, 1, 1/4)
col_dmg <- rgb(1, 0, 0, 1/4)

pdf("density_comparsion.pdf", width = 8.27, height = 11.69)  # A4 size in inches
par(mfrow = c(3, 2), oma = c(6, 1, 2, 1), mar = c(3, 3, 2, 1))  # Adjusted outer and inner margins

# R20
l <- R_pf[1]
set.seed(4)
R20_nodmg<- pl_gen(mu,sigma,rho,eta,0,l,N)
R20_dmg<- pl_gen(mu,sigma,rho,eta,1,l,N)
p1 <- density(data_to_plot(R20_nodmg))     
p2 <- density(data_to_plot(R20_dmg))  
xlim_range <- range(c(p1$x, p2$x))
ylim_range <- range(c(p1$y, p2$y))
plot( p1, col=col_nodmg , 
      xlim = xlim_range, ylim = ylim_range,
      main = "R20: Density plot ",xlab = "")  # first histogram
lines( p2, col=col_dmg )  # second
# legend("topright", legend = c("nodmg", "dmg"),
#        col = c(col_nodmg, col_dmg ),pch = 15, pt.cex = 2, bty = "n")


# T20
l <- T_pf[1]
nodmg <- pl_gen(c(mu[2],mu[1]),c(sigma[2],sigma[1]),
                rho,eta,0,l,N)
dmg <- pl_gen(c(mu[2],mu[1]),c(sigma[2],sigma[1]),
              rho,eta,1,l,N)
p1 <- density(data_to_plot(nodmg))     
p2 <- density(data_to_plot(dmg))  
xlim_range <- range(c(p1$x, p2$x))
ylim_range <- range(c(p1$y, p2$y))
plot( p1, col=col_nodmg , 
      xlim = xlim_range, ylim = ylim_range,
      main = "T20: Density plot ",xlab = "")  # first histogram
lines( p2, col=col_dmg )  # second
# legend("topright", legend = c("nodmg", "dmg"),
#        col = c(col_nodmg, col_dmg ),pch = 15, pt.cex = 2, bty = "n")


# R40
l <- R_pf[2]
R40_nodmg<- pl_gen(mu,sigma,rho,eta,0,l,N)
R40_dmg<- pl_gen(mu,sigma,rho,eta,1.5,l,N)
p1 <- density(data_to_plot(R40_nodmg))     
p2 <- density(data_to_plot(R40_dmg))  
xlim_range <- range(c(p1$x, p2$x))
ylim_range <- range(c(p1$y, p2$y))
plot( p1, col=col_nodmg , 
      xlim = xlim_range, ylim = ylim_range,
      main = "R40: Density plot ",xlab = "")  # first histogram
lines( p2, col=col_dmg )  # second
# legend("topright", legend = c("nodmg", "dmg"),
#        col = c(col_nodmg, col_dmg ),pch = 15, pt.cex = 2, bty = "n")






# T40
l <- T_pf[2]
nodmg <- pl_gen(c(mu[2],mu[1]),c(sigma[2],sigma[1]),
                rho,eta,0,l,N)
dmg <- pl_gen(c(mu[2],mu[1]),c(sigma[2],sigma[1]),
              rho,eta,1,l,N)
p1 <- density(data_to_plot(nodmg))     
p2 <- density(data_to_plot(dmg))  
xlim_range <- range(c(p1$x, p2$x))
ylim_range <- range(c(p1$y, p2$y))
plot( p1, col=col_nodmg , 
      xlim = xlim_range, ylim = ylim_range,
      main = "T40: Density plot ",xlab = "")  # first histogram
lines( p2, col=col_dmg )  # second
# legend("topright", legend = c("nodmg", "dmg"),
#        col = c(col_nodmg, col_dmg ),pch = 15, pt.cex = 2, bty = "n")


# R60
l <- R_pf[3]
R60_nodmg<- pl_gen(mu,sigma,rho,eta,0,l,N)
R60_dmg<- pl_gen(mu,sigma,rho,eta,2,l,N)
p1 <- density(data_to_plot(R60_nodmg))     
p2 <- density(data_to_plot(R60_dmg))  
xlim_range <- range(c(p1$x, p2$x))
ylim_range <- range(c(p1$y, p2$y))
plot( p1, col=col_nodmg , 
      xlim = xlim_range, ylim = ylim_range,
      main = "R60: Density plot ",xlab = "")  # first histogram
lines( p2, col=col_dmg )  # second
# legend("topright", legend = c("nodmg", "dmg"),
#        col = c(col_nodmg, col_dmg ),pch = 15, pt.cex = 2, bty = "n")


# T60
l <- T_pf[3]
nodmg <- pl_gen(c(mu[2],mu[1]),c(sigma[2],sigma[1]),
                rho,eta,0,l,N)
dmg <- pl_gen(c(mu[2],mu[1]),c(sigma[2],sigma[1]),
              rho,eta,1,l,N)
p1 <- density(data_to_plot(nodmg))     
p2 <- density(data_to_plot(dmg))  
xlim_range <- range(c(p1$x, p2$x))
ylim_range <- range(c(p1$y, p2$y))
plot( p1, col=col_nodmg , 
      xlim = xlim_range, ylim = ylim_range,
      main = "T60: Density plot ",xlab = "")  # first histogram
lines( p2, col=col_dmg )  # second
# legend("topright", legend = c("nodmg", "dmg"),
#        col = c(col_nodmg, col_dmg ),pch = 15, pt.cex = 2, bty = "n")

# Add a common legend at the bottom of all plots
par(fig = c(0, 1, 0, 0.1), oma = c(0, 0, 0, 0), new = TRUE)  # Use the lower part of the figure space for the legend
plot.new()  # Create an empty plot for the legend
legend("center", legend = c("non-damaged survivors", "damaged survivors"),
       col = c(col_nodmg, col_dmg), pch = 15, pt.cex = 2, bty = "n", horiz = TRUE, inset = 0,
       text.width = 0.25, xjust = 0.5) 
dev.off()
