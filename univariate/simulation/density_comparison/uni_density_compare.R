library(ggplot2)
library(gridExtra)



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
  ))
  
}

dmg_inverse <- function(y,c,alpha,l){
  return ((y + alpha*l/c)/(alpha+1))
}
trunc_pdf_y <- function(x,a,b,mu,sigma){
  truncnorm::dtruncnorm(x,a = a, b = b,mean = mu,
                        sd = sigma)
}


lik_pl_dmg <- function(mu, sigma, alpha, c, l, data){
  l_star <- dmg_inverse(l, c, alpha, l)
  group3 <- data
  
  # Damaged data
  ystar_dmg <- group3[group3 < l/c]
  y_dmg <- dmg_inverse(ystar_dmg, c, alpha, l)
  
  lik_dmg <- truncnorm::dtruncnorm(y_dmg, a = l_star, mean = mu, sd = sigma)
  lik_dmg <- lik_dmg / (alpha + 1)
  
  # Undamaged pieces
  ystar_undmg <- group3[group3 > l/c]
  
  lik_undmg <- numeric(length(ystar_undmg))
  if(length(ystar_undmg) != 0){
    lik_undmg <- truncnorm::dtruncnorm(ystar_undmg, a = l_star, mean = mu, sd = sigma)
  }
  
  # Combine likelihoods for each data point
  likelihoods <- numeric(length(group3))
  likelihoods[group3 < l/c] <- lik_dmg
  likelihoods[group3 > l/c] <- lik_undmg
  likelihoods[group3 < l] <- 0
  return (likelihoods)
}


# normalization to 1:
density_norm <- function(x,y){
  area <- density_calc(x,y)
  factor <- 1/area
  return (y*factor)
}


density_calc <- function(x,y){
  size <- x[2]- x[1]
  
  na_id <- is.na(y)
  x <- x[!na_id]
  y <- y[!na_id]
  len <- length(y)
  area <- sum(size * (y[1:(len-1)] + y[2:(len)])/2)
  return(area)
}



plot_den_biv <- function(l,dmg_alpha,group_name,RorT,tol = 1e-3){

  x_seq <- seq(from = -1, to = 12, length = 300)
  undmg_y <- PFY_lik_R(mu, sigma, rho,eta,0,l,x_seq)
  dmg_y <- PFY_lik_R(mu, sigma, rho,eta,dmg_alpha,l,x_seq)


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



mu <- 51
sigma <-  18


c <-  0.76



l1 <- 31.02750
l2 <- 20.685
l <- l2
theta0 <- c(mu,sigma,alpha)








# plot begins
col_nodmg <- rgb(0, 0, 1, 1/4)
col_dmg <- rgb(1, 0, 0, 1/4)



plot_den <- function(l,dmg_alpha,group_name,tol = -1){
  
  x_seq <- seq(from = l-5, to = 140, length = 300)
  undmg_y <- lik_pl_dmg(mu,sigma,0, c,l, x_seq)
  dmg_y <- lik_pl_dmg(mu,sigma,dmg_alpha, c,l, x_seq)
  

  undmg_den <- density_norm(x_seq,undmg_y)
  dmg_den <- density_norm(x_seq,dmg_y)
  undmg_mark <- (undmg_den > tol)
  dmg_mark <- (dmg_den > tol)

  x_combined <- c(x_seq[undmg_mark], x_seq[dmg_mark])
  y_combined <- c(undmg_den[undmg_mark], dmg_den[dmg_mark])
  # 
  xlim <- range(x_combined)
  ylim <- range(y_combined)
  # plot(x_seq[dmg_mark],dmg_den[dmg_mark],
  #      type = 'l',col = col_dmg,xlab = 'survivors',
  #      ylab = 'density', main = group_name,
  #      xlim = xlim, ylim = ylim)
  plot(x_seq[undmg_mark],undmg_den[undmg_mark],
       type = 'l',col = col_nodmg,xlab = 'survivors',
       ylab = 'density', 
       main = bquote(.(group_name) * ", " ~ alpha == .(dmg_alpha)),
       xlim = xlim, ylim = ylim)
  lines(x_seq[dmg_mark],dmg_den[dmg_mark],
        type = 'l', col = col_dmg)
  return (dmg_y)
}
pdf("density_comparsion_uni.pdf", width = 9.27, height = 9.69)  # A4 size in inches

par(mfrow = c(2, 2), oma = c(6, 1, 2, 1),
    mar = c(3, 5, 2, 1))


plot_den(l1,5,"R20_1Y")
plot_den(l1,2,"R20_4Y")

plot_den(l2,2,"R5_4Y")

# Add a common legend at the bottom of all plots
par(fig = c(0, 1, 0, 0.1), oma = c(0, 0, 0, 0), new = TRUE)  # Use the lower part of the figure space for the legend
plot.new()  # Create an empty plot for the legend
legend("center", 
       legend = c(expression("Density of survivor strength generated with true " * bold(alpha)),
         expression("Density of survivor strength generated with " * bold(alpha) == bold(0))),
                  
       col = c(col_dmg,col_nodmg), pch = 15, 
       pt.cex = 2, bty = "n", horiz = TRUE, inset = 0,
       text.width = 0.5, xjust = 0.5,cex = 1.1) 
dev.off()
