######################################
# Author: Alice Yuan Qi
# Fan, Qi & Xu (2025), Figure 7
# Population-level comparisons between the 0.25-EWM, equality-minded, and 1-EWM policies using data generating process in Appendix I.3; tau is specified as (I.1) in paper
######################################

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(optimization)
library(DescTools)  # for computing Gini
library(ggplot2)

superpop <- read.csv("data_tau_1.csv")
superpop$Y_cf <- (superpop$A==1)*(superpop$Y-superpop$tau)+(superpop$A==0)*(superpop$Y+superpop$tau)

#### Treatment statuses of different targeting rules ####
# These policy parameters correspond to the optimized welfare for alpha = 0.25 and 1 in truth_comparisons.csv (Table 10 in Appendix I.3)
best_ITR_0.25 <- -3.36+4.49*superpop$X1+4.32*superpop$X2>0
best_ITR_1 <- -3.31+4.51*superpop$X1+4.44*superpop$X2-0.04*superpop$X3+0.02*superpop$X4>0

#### Helper function that computes the true value of the objective function (Gini SWF) given some combination of linear policy parameters, used for policy learning via simulated annealing. Computation of the Gini SWF is based on Equation (5) in Kitagawa & Tetenov (2021). ####
### "d_pop" is equivalent to ITR "pi"
obj_linear <- function(coef) {
  if (counter %% 1000 == 0) {
    print(paste0("Progress: ", 100*counter/total, "%"))  # show progress
  }
  counter <<- counter+1
  d_pop <- coef[1]+coef[2]*superpop$X1+coef[3]*superpop$X2+coef[4]*superpop$X3+coef[5]*superpop$X4>0
  outcomes <- (superpop$A==d_pop)*superpop$Y+(superpop$A!=d_pop)*superpop$Y_cf
  Gini_SWF <- mean(outcomes)*(1-Gini(outcomes))
  return((-1)*Gini_SWF)  # include negative sign for loss minimization
}

#### Learn true equality-minded linear rule: pi=1 if LES>0 ####
## Specify parameters for simulated annealing
rate <- 0.8
n_limit <- 600
counter <- 0
total <- ceiling(log(0.1/1000,rate))*n_limit  # for showing progress

set.seed(123)
obj_sa <- optim_sa(fun = obj_linear, start = c(-3.31, 4.51, 4.44, -0.04, 0.02), lower = c(-5,-5,-5,-5,-5), upper = c(5,5,5,5,5), control = list(r = rate, nlimit = n_limit))
print(paste0("True = ", -obj_sa$function_value))
best_linear <- obj_sa$par
print(best_linear)

best_ITR_EM <- best_linear[1]+best_linear[2]*superpop$X1+best_linear[3]*superpop$X2+best_linear[4]*superpop$X3+best_linear[5]*superpop$X4>0
print(paste0("mean(best_ITR_EM): ", mean(best_ITR_EM)))

#### Quantile Difference Plot ####
outcomes_CVaR <- (superpop$A==best_ITR_0.25)*superpop$Y+(superpop$A!=best_ITR_0.25)*superpop$Y_cf
outcomes_mean <- (superpop$A==best_ITR_1)*superpop$Y+(superpop$A!=best_ITR_1)*superpop$Y_cf
outcomes_EM <- (superpop$A==best_ITR_EM)*superpop$Y+(superpop$A!=best_ITR_EM)*superpop$Y_cf

## CVaR with alp=0.25 VS mean-optimal; the equality-minded policy is identical to mean-optimal
q_diff <- numeric(2000)
seq_p <- seq(0,1,by=0.0005)[-2001]
for (i in 1:length(seq_p)) {
  q_diff[i] <- quantile(outcomes_CVaR, seq_p[i]) - quantile(outcomes_mean, seq_p[i])
}
df <- cbind(seq_p,q_diff)

jpeg(file="q_diff_all.jpeg", width=9, height=5, units="in", res=300)
ggplot(data = df, aes(x = seq_p)) + 
  geom_hline(yintercept = 0, linewidth=0.2) +
  geom_ribbon(aes(ymin = 0, ymax = ifelse(q_diff > 0, q_diff, 0), fill="0.25-EWM VS 1-EWM\nand equality-minded"), alpha=0.7) +
  geom_ribbon(aes(ymax = 0, ymin = ifelse(q_diff > 0, 0, q_diff), fill="0.25-EWM VS 1-EWM\nand equality-minded"), alpha=0.7) +
  geom_line(aes(y = q_diff), color="cyan2") +
  scale_fill_manual(values = c("cyan2")) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, size=18)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text = element_text(size=10)) +
  xlab("Probability level") +
  ylab("Difference between quantiles")
dev.off()
