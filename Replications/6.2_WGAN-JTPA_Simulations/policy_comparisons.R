######################################
# Author: Alice Yuan Qi
# Fan, Qi & Xu (2025), Figure 4
# Population-level comparisons between the 0.25-EWM, equality-minded, and 1-EWM policies using WGAN-JTPA data in Section 6.2
######################################

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(optimization)
library(DescTools)  # for computing Gini
library(ggplot2)

superpop <- read.csv("WGAN-JTPA.csv")

#### Treatment statuses of different targeting rules ####
# These policy parameters correspond to the optimized welfare for alpha = 0.25 and 1 in truth_comparisons.csv (Table 9 in Appendix I.2)
best_ITR_0.25 <- 37.369427482072-3.07728059342574*superpop$edu-0.000775845455541579*superpop$prevearn>0
best_ITR_1 <- 30.7961522032475-2.4468437984237*superpop$edu+0.000673120273274756*superpop$prevearn>0

#### Helper function that computes the true value of the objective function (Gini SWF) given some combination of linear policy parameters, used for policy learning via simulated annealing. Computation of the Gini SWF is based on Equation (5) in Kitagawa & Tetenov (2021). ####
### "D" is equivalent to "A" in paper, and "d_pop" is equivalent to ITR "pi"
obj_linear <- function(coef) {
  if (counter %% 1000 == 0) {
    print(paste0("Progress: ", 100*counter/total, "%"))  # show progress
  }
  counter <<- counter+1
  d_pop <- coef[1]+coef[2]*edu_std+coef[3]*prevearn_std>0
  outcomes <- (superpop$D==d_pop)*superpop$earnings+(superpop$D!=d_pop)*superpop$earnings_cf
  Gini_SWF <- mean(outcomes)*(1-Gini(outcomes))
  return((-1)*Gini_SWF)  # include negative sign for loss minimization
}

# standardize edu & prevearn so that each has mean=0 and var=1
edu_std <- as.numeric(scale(superpop$edu))
prevearn_std <- as.numeric(scale(superpop$prevearn))

#### Learn true linear rule: pi=1 if LES>0 ####
## Specify parameters for simulated annealing
rate <- 0.8
n_limit <- 600
counter <- 0
total <- ceiling(log(0.1/1000,rate))*n_limit  # for showing progress

set.seed(1234)
obj_sa <- optim_sa(fun = obj_linear, start = c(0,0,0), lower = c(-5,-5,-5), upper = c(5,5,5), control = list(r = rate, nlimit = n_limit))
print(paste0("True = ", -obj_sa$function_value))
best_linear <- obj_sa$par
b0 <- best_linear[1]-best_linear[2]/sd(superpop$edu)*mean(superpop$edu)-best_linear[3]/sd(superpop$prevearn)*mean(superpop$prevearn)
b1 <- best_linear[2]/sd(superpop$edu)
b2 <- best_linear[3]/sd(superpop$prevearn)
print(paste0("Linear coefficients: ", c(b0,b1,b2)))

d_linear <- b0+b1*superpop$edu+b2*superpop$prevearn>0
# d_linear <- 28.23935-2.222134*superpop$edu+0.0001081318*superpop$prevearn>0 # using output from the optimization above
print(paste0("mean(d_linear): ", mean(d_linear)))

#### Quantile Difference Plots ####
## CVaR with alp=0.25 VS equality-minded
outcomes_CVaR <- (superpop$D==best_ITR_0.25)*superpop$earnings + (superpop$D!=best_ITR_0.25)*superpop$earnings_cf
outcomes_mean <- (superpop$D==best_ITR_1)*superpop$earnings + (superpop$D!=best_ITR_1)*superpop$earnings_cf
outcomes_EM <- (superpop$D==d_linear)*superpop$earnings + (superpop$D!=d_linear)*superpop$earnings_cf

q_diff_1 <- numeric(2000)
seq_p <- seq(0,1,by=0.0005)[-2001]  # the last quantile is the max (a single value), which is too noisy
for (i in 1:length(seq_p)) {
  q_diff_1[i] <- quantile(outcomes_CVaR, seq_p[i]) - quantile(outcomes_EM, seq_p[i])
}
d1 <- cbind(seq_p,q_diff_1)

## Equality-minded VS mean-optimal
q_diff_2 <- numeric(2000)
for (i in 1:length(seq_p)) {
  q_diff_2[i] <- quantile(outcomes_EM, seq_p[i]) - quantile(outcomes_mean, seq_p[i])
}
d2 <- cbind(seq_p,q_diff_2)

## CVaR with alp=0.25 VS mean-optimal
q_diff_3 <- numeric(2000)
for (i in 1:length(seq_p)) {
  q_diff_3[i] <- quantile(outcomes_CVaR, seq_p[i]) - quantile(outcomes_mean, seq_p[i])
}
d3 <- cbind(seq_p,q_diff_3)

## All in one plot
jpeg(file="q_diff_all.jpeg", width=9, height=5, units="in", res=300)
ggplot(data = d3, aes(x = seq_p)) + 
  geom_hline(yintercept = 0, linewidth=0.2) +
  geom_ribbon(aes(ymin = 0, ymax = ifelse(q_diff_3 > 0, q_diff_3, 0), fill="0.25-EWM VS 1-EWM"), alpha=0.7) +
  geom_ribbon(aes(ymax = 0, ymin = ifelse(q_diff_3 > 0, 0, q_diff_3), fill="0.25-EWM VS 1-EWM"), alpha=0.7) +
  geom_line(aes(y = q_diff_3), color="cyan2") +
  geom_ribbon(aes(ymin = 0, ymax = ifelse(q_diff_1 > 0, q_diff_1, 0), fill="0.25-EWM VS equality-minded"), alpha=0.7) +
  geom_ribbon(aes(ymax = 0, ymin = ifelse(q_diff_1 > 0, 0, q_diff_1), fill="0.25-EWM VS equality-minded"), alpha=0.7) +
  geom_line(aes(y = q_diff_1), color="dodgerblue3") +
  geom_ribbon(aes(ymin = 0, ymax = ifelse(q_diff_2 > 0, q_diff_2, 0), fill="Equality-minded VS 1-EWM"), alpha=0.7) +
  geom_ribbon(aes(ymax = 0, ymin = ifelse(q_diff_2 > 0, 0, q_diff_2), fill="Equality-minded VS 1-EWM"), alpha=0.7) +
  geom_line(aes(y = q_diff_2), color="royalblue4") +
  scale_fill_manual(values = c("cyan2", "dodgerblue3", "royalblue4")) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, size=18)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text = element_text(size=10)) +
  xlab("Probability level") +
  ylab("Difference between quantiles") +
  ylim(c(-1950,280))
dev.off()
