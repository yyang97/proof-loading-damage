myqq <- function(x, y, ...) {
  rg <- range(x, y, na.rm = T)
  qqplot(x, y, xlim = rg, ylim = rg, ...)
  abline(0,1)
}

dmg_model <- function(y, c, alpha, l) {
  return(ifelse(y > l/c,
                y,
                y - alpha*(l/c-y)^(2)
  )
  )
}
bending <- read.csv("bending-pl-4groups.csv", header = F)
group_names <- table(bending[, 2])
names(group_names)
bending[, 1] <- bending[, 1] / 1000 * 6.895
R20R100 <- bending[which(bending[, 2] == "R20R100"), 1]
R20_4Y_R100 <- bending[which(bending[, 2] == "R20_4Y_R100"), 1]
R5R100 <- bending[which(bending[, 2] == "R5R100"), 1]
R100 <- bending[which(bending[, 2] == "R100"), 1]
#hist(R20R100)
# data generation
#source("dmg_func.R")
N_R20_1Y <- 300
N_R20_4Y <- 101
N_R5_4Y <- 198

N_R100 <- 140
N <- N_R20_1Y

# mu <- 48
# sigma <-  19


mu <- 50.816
sigma <-  17.912
c <-  0.76
pdf("qqplots_univariate.pdf", width = 8.27, 
    height = 9.19)  # A4 size in inches
par(mfrow = c(2, 2))





####---------------------- R5_4Y------------------
alpha <-  0


l1 <- 31.02750
l2 <- 20.685
l <- l2
set.seed(9)
N <- N_R5_4Y
y <- rnorm(N_R5_4Y, mean = mu, sd = sigma)
y <- y[order(y)]
# truncated data
y_star <- y[y > l1 / c]
# damage data
# group 2, y^*<l<y
y_obs <- list()
y_obs$group1 <- y[which(y < l)]
group23 <- y[which(y > l)]
group23_star <- dmg_model(group23, c, alpha, l)
y_obs$group2 <-
  length(group23[which(group23 > l & group23_star < l)])
y_obs$group3 <- group23_star[which(group23_star > l)]

#hist(y_obs$group3)



myqq(R5R100 ,y_obs$group3,xlab = "empirical quantile",
     ylab = "simulated quantile",
     main = bquote("R5_4Y: Q-Q plot (" * hat(alpha) == 0 * ")"))




####---------------------- R20_4Y------------------
alpha <-  0


l1 <- 31.02750
l2 <- 20.685
l <- l1
set.seed(8)
N <- N_R20_4Y
y <- rnorm(N_R20_4Y, mean = mu, sd = sigma)
y <- y[order(y)]
# truncated data
y_star <- y[y > l1 / c]
# damage data
# group 2, y^*<l<y
y_obs <- list()
y_obs$group1 <- y[which(y < l)]
group23 <- y[which(y > l)]
group23_star <- dmg_model(group23, c, alpha, l)
y_obs$group2 <-
  length(group23[which(group23 > l & group23_star < l)])
y_obs$group3 <- group23_star[which(group23_star > l)]

#hist(y_obs$group3)



myqq(R20_4Y_R100,y_obs$group3,xlab = "empirical quantile",
     ylab = "simulated quantile",
     main = bquote("R20_4Y: Q-Q plot (" * hat(alpha) == 0 * ")"))
####---------------------- R20_1Y------------------
alpha <-  4.908


l1 <- 31.02750
l2 <- 20.685
l <- l1
set.seed(27)
N <- N_R20_1Y
y <- rnorm(N_R20_1Y, mean = mu, sd = sigma)
y <- y[order(y)]
# truncated data
y_star <- y[y > l1 / c]
# damage data
# group 2, y^*<l<y
y_obs <- list()
y_obs$group1 <- y[which(y < l)]
group23 <- y[which(y > l)]
group23_star <- dmg_model(group23, c, alpha, l)
y_obs$group2 <-
  length(group23[which(group23 > l & group23_star < l)])
y_obs$group3 <- group23_star[which(group23_star > l)]

#hist(y_obs$group3)



myqq(R20R100,y_obs$group3,xlab = "empirical quantile",
     ylab = "simulated quantile",
     main = bquote("R20_1Y: Q-Q plot (" * hat(alpha) == 4.908 * ")"))

####---------------------- R5_4Y------------------
alpha <-  0


l1 <- 31.02750
l2 <- 20.685
l <- l2
set.seed(27)
N <- N_R20_1Y
y <- rnorm(N_R20_1Y, mean = mu, sd = sigma)
y <- y[order(y)]
# truncated data
y_star <- y[y > l1 / c]
# damage data
# group 2, y^*<l<y
y_obs <- list()
y_obs$group1 <- y[which(y < l)]
group23 <- y[which(y > l)]
group23_star <- dmg_model(group23, c, alpha, l)
y_obs$group2 <-
  length(group23[which(group23 > l & group23_star < l)])
y_obs$group3 <- group23_star[which(group23_star > l)]

#hist(y_obs$group3)



myqq(R20R100 ,y_obs$group3,xlab = "empirical quantile",
     ylab = "simulated quantile",
     main = bquote("R20_1Y: Q-Q plot (" * alpha == 0 * ")"))

par(mfrow = c(1, 1))

dev.off()