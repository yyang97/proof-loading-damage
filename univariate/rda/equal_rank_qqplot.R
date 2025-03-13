##--------------get data ----------------
bending <- read.csv("bending-pl-4groups.csv", header = F)
group_names <- table(bending[, 2])
names(group_names)


bending[, 1] <- bending[, 1] / 1000 * 6.895
R20 <- bending[which(bending[, 2] == "R20"), 1]
R20R100 <- bending[which(bending[, 2] == "R20R100"), 1]
R20_4Y <- bending[which(bending[, 2] == "R20_4Y"), 1]
R20_4Y_R100 <- bending[which(bending[, 2] == "R20_4Y_R100"), 1]
R5R100 <- bending[which(bending[, 2] == "R5R100"), 1]
R5 <- bending[which(bending[, 2] == "R5"), 1]
R100 <- bending[which(bending[, 2] == "R100"), 1]
N_R20_1Y_g3 <- length(R20R100)
N_R20_4Y_g3 <- length(R20_4Y_R100)
N_R5_4Y_g3 <- length(R5R100)

N_R20_1Y_g2 <- 97
N_R20_4Y_g2 <- 41
N_R5_4Y_g2 <- 42

l1 <- 31.02750
l2 <- 20.685


# define the 1-to-1 relationship between survivors and control, using quantile function
control_qfun <- function(quant){return(quantile(R100,prob = quant,type = 3))}

##-----------calcuate strength degradation (SD) for each group-----------
R20_1Y_q <- seq(from = (300 - N_R20_1Y_g3 + 1)/300,
                to =1 , 
                length = N_R20_1Y_g3)
R20_1Y_sd <- 1- R20R100/control_qfun(R20_1Y_q)


R20_4Y_q <- seq(from = (101 - N_R20_4Y_g3 + 1)/101,
                to =1 , 
                length = N_R20_4Y_g3)
R20_4Y_sd <- 1- R20_4Y_R100/control_qfun(R20_4Y_q)


R5_4Y_q <- seq(from = (198 - N_R5_4Y_g3 + 1)/198,
               to =1 , 
               length = N_R5_4Y_g3)
R5_4Y_sd <- 1- R5R100/control_qfun(R5_4Y_q)

## -------- implement the linear regression --------
Y1 <- 12
Y4 <- 48

R20_1Y_df <- cbind(R20_1Y_sd,l1/mean(R100),Y1,R20R100,control_qfun(R20_1Y_q))
R20_4Y_df <- cbind(R20_4Y_sd,l1/mean(R100),Y4,R20_4Y_R100,control_qfun(R20_4Y_q))
R5_4Y_df <- cbind(R5_4Y_sd,l2/mean(R100),Y4,R5R100,control_qfun(R5_4Y_q))


dmg_df <- data.frame(rbind(R20_1Y_df,R20_4Y_df,R5_4Y_df))
colnames(dmg_df) <- c("SD","SL","TD","pl_MOR","control_MOR")
dmg_df <- dmg_df[dmg_df$SD > 0.05,]


lm_model <- lm(SD ~ SL + TD, data = dmg_df)
summary(lm_model)

dim(dmg_df)


##--------------Q-Q plot -----------
myqq <- function(x, y, ...) {
  rg <- range(x, y, na.rm = T)
  qqplot(x, y, xlim = rg, ylim = rg, ...)
  abline(0,1,col = 'black')
}

SD_pred <- predict(lm_model,dmg_df)
myqq(dmg_df$SD,SD_pred,xlab = "empirical quantile",
     ylab = "predicted quantile",
     main = "Q-Q plot from the regression based on equal rank")