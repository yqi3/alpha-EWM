######################################
# Author: Alice Yuan Qi
# Fan, Qi & Xu (2025), Tables 4 & 9: 6 alphas (alpha=1 is not used in main text)
# Population truths using WGAN-JTPA data in Section 6.2

# For reasonable run time, we recommend running this script with at least 40 CPUs simultaneously, or breaking the for loop into smaller sub-tasks.
######################################

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(optimization)

#### Helper function that computes the true value of the objective function (generalized Lorenz, can be easily converted to AVaR by dividing alpha) given some combination of linear policy parameters and eta, used for policy learning via simulated annealing ####
### "D" is equivalent to "A" in paper, and "d_pop" is equivalent to ITR "pi"
obj_linear <- function(coef) {
  if (counter %% 10000 == 0) {
    print(paste0("Progress: ", 100*counter/total, "%"))  # show progress
  }
  counter <<- counter+1
  d_pop <- coef[1]+coef[2]*edu_std+coef[3]*prevearn_std>0
  outcomes <- (superpop$D==d_pop)*superpop$earnings+(superpop$D!=d_pop)*superpop$earnings_cf
  return((-1)*alp*mean(outcomes[outcomes<=quantile(outcomes,alp)]))  # include negative sign for loss minimization
}

#### Read superpopulation data ####
superpop <- read.csv("WGAN-JTPA.csv")

# standardize edu & prevearn so that each has mean=0 and var=1
edu_std <- as.numeric(scale(superpop$edu))
prevearn_std <- as.numeric(scale(superpop$prevearn))

#### Learn true linear rule: pi=1 if LES>0 ####
## Specify parameters for simulated annealing
rate <- 0.8
n_limit <- 600
total <- ceiling(log(0.1/1000,rate))*n_limit  # for showing progress

## Vectors for storing results
welfare_values <- numeric(0)
best_linear_0.25 <- numeric(3)
best_linear_0.3 <- numeric(3)
best_linear_0.4 <- numeric(3)
best_linear_0.5 <- numeric(3)
best_linear_0.8 <- numeric(3)
best_linear_1 <- numeric(3)

## Run SA for superpopulation
for (alp in c(0.25,0.3,0.4,0.5,0.8,1)) {
  print(paste0("Learning true optimal policy for alpha = ", alp, "..."))
  set.seed(1234)
  counter <- 0
  
  if (alp == 0.4) {
    start <- best_linear_0.3  # good candidate from previous case, use as initial guess
  } else if (alp == 0.5) {
    start <- best_linear_0.3  # good candidate from previous case, use as initial guess
  } else { start <- numeric(3) }
  
  obj_sa <- optim_sa(fun = obj_linear, start = start, lower = c(-5,-5,-5), upper = c(5,5,5), control = list(r = rate, nlimit = n_limit))
  best_linear <- obj_sa$par
  assign(paste0("best_linear_",alp), best_linear)
  b0 <- best_linear[1]-best_linear[2]/sd(superpop$edu)*mean(superpop$edu)-best_linear[3]/sd(superpop$prevearn)*mean(superpop$prevearn)
  b1 <- best_linear[2]/sd(superpop$edu)
  b2 <- best_linear[3]/sd(superpop$prevearn)
  print(paste0("Linear coefficients: ", c(b0,b1,b2)))
  
  d_linear <- b0+b1*superpop$edu+b2*superpop$prevearn>0
  print(paste0("mean(d_linear): ", mean(d_linear)))
  outcomes <- (superpop$D==d_linear)*superpop$earnings+(superpop$D!=d_linear)*superpop$earnings_cf
  curr_col <- c(mean(outcomes[outcomes<=quantile(outcomes,0.25)]), mean(outcomes[outcomes<=quantile(outcomes,0.3)]), mean(outcomes[outcomes<=quantile(outcomes,0.4)]), mean(outcomes[outcomes<=quantile(outcomes,0.5)]), mean(outcomes[outcomes<=quantile(outcomes,0.8)]), mean(outcomes[outcomes<=quantile(outcomes,1)]))
  print(curr_col)
  welfare_values <- c(welfare_values, curr_col)
  print("----------")
}

comparisons_linear <- matrix(welfare_values, ncol = 6, nrow = 6)
colnames(comparisons_linear) <- c("0.25", "0.3", "0.4", "0.5", "0.8", "1")
rownames(comparisons_linear) <- c("0.25", "0.3", "0.4", "0.5", "0.8", "1")
write.csv(comparisons_linear, "truth_comparisons.csv")
