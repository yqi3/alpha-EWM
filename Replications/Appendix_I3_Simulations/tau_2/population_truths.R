######################################
# Author: Alice Yuan Qi
# Fan, Qi & Xu (2025), Tables 11 & 13: 6 alphas (alpha=1 is not used in main text)
# Population truths using data generating process in Appendix I.3; tau is specified as (I.2) in paper

# For reasonable run time, we recommend running this script with at least 40 CPUs simultaneously, or breaking the for loop into smaller sub-tasks.
######################################

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(optimization)
set.seed(1234)

#### Helper function that computes the true value of the objective function (generalized Lorenz, can be easily converted to AVaR by dividing alpha) given some combination of linear policy parameters and eta, used for policy learning via simulated annealing ####
obj_linear <- function(coef) {
  if (counter %% 10000 == 0) {
    print(paste0("Progress: ", 100*counter/total, "%"))  # show progress
  }
  counter <<- counter+1
  d_pop <- coef[1]+coef[2]*superpop$X1+coef[3]*superpop$X2+coef[4]*superpop$X3+coef[5]*superpop$X4>0
  outcomes <- (superpop$A==d_pop)*superpop$Y+(superpop$A!=d_pop)*superpop$Y_cf
  return((-1)*alp*mean(outcomes[outcomes<=quantile(outcomes,alp)]))  # include negative sign for loss minimization
}

#### Read superpopulation data ####
superpop <- read.csv("data_tau_2.csv")
superpop$Y_cf <- (superpop$A==1)*(superpop$Y-superpop$tau) + (superpop$A==0)*(superpop$Y+superpop$tau)  # compute true counterfactuals

#### Learn true linear rule: pi=1 if LES>0 ####
## Specify parameters for simulated annealing
rate <- 0.8
n_limit <- 600
total <- ceiling(log(0.1/1000,rate))*n_limit  # for showing progress

## Vectors for storing results
welfare_values <- numeric(0)
best_linear_0.25 <- numeric(5)
best_linear_0.3 <- numeric(5)
best_linear_0.4 <- numeric(5)
best_linear_0.5 <- numeric(5)
best_linear_0.8 <- numeric(5)
best_linear_1 <- numeric(5)
start <- numeric(5)

## Run SA for superpopulation
for (alp in c(0.25,0.3,0.4,0.5,0.8,1)) {
  print(paste0("Learning true optimal policy for alpha = ", alp, "..."))
  set.seed(1234)
  counter <- 0
  
  # Good initial conditions
  if (alp == 0.25) {
    start <- c(-3.49, -2.60, -2.54, -0.11, -0.23)
  } else if (alp == 0.3) {
    start <- c(-3.36, -2.64, -2.54, -0.22, -0.14)
  } else if (alp == 0.4) {
    start <- c(-3.41, -2.59, -2.53, -0.17, -0.11)
  } else if (alp == 0.5) {
    start <- c(-3.73, 3.38, 3.32, -0.03, -0.24)
  } else if (alp == 0.8) {
    start <- c(-3.81, 3.46, 3.40, -0.07, -0.28)
  } else {
    start <- c(-3.95, 3.96, 3.81, 0.19, -0.49)
  }
  
  obj_sa <- optim_sa(fun = obj_linear, start = start, lower = c(-5,-5,-5,-5,-5), upper = c(5,5,5,5,5), control = list(r = rate, nlimit = n_limit))
  best_linear <- obj_sa$par
  assign(paste0("best_linear_",alp), best_linear)
  print(paste0("Linear coefficients: ", best_linear))
  
  d_linear <- best_linear[1]+best_linear[2]*superpop$X1+best_linear[3]*superpop$X2+best_linear[4]*superpop$X3+best_linear[5]*superpop$X4>0
  print(paste0("mean(d_linear): ", mean(d_linear)))
  outcomes <- (superpop$A==d_linear)*superpop$Y+(superpop$A!=d_linear)*superpop$Y_cf
  curr_col <- c(mean(outcomes[outcomes<=quantile(outcomes,0.25)]), mean(outcomes[outcomes<=quantile(outcomes,0.3)]), mean(outcomes[outcomes<=quantile(outcomes,0.4)]), mean(outcomes[outcomes<=quantile(outcomes,0.5)]), mean(outcomes[outcomes<=quantile(outcomes,0.8)]), mean(outcomes[outcomes<=quantile(outcomes,1)]))
  print(curr_col)
  welfare_values <- c(welfare_values, curr_col)
  print("----------")
}

comparisons_linear <- matrix(welfare_values, ncol = 6, nrow = 6)
colnames(comparisons_linear) <- c("0.25", "0.3", "0.4", "0.5", "0.8", "1")
rownames(comparisons_linear) <- c("0.25", "0.3", "0.4", "0.5", "0.8", "1")
write.csv(comparisons_linear, "truth_comparisons.csv")
