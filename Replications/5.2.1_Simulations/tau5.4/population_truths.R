######################################
# Author: Alice Yuan Qi
# Fan, Qi & Xu (2024), Table 5: 6 alphas
# Population truths using data generating process in Section 5.2.1; tau is specified as (5.4) in paper

# For reasonable run time, we recommend running this script with at least 40 CPUs simultaneously, or breaking the for loop into smaller sub-tasks.
######################################

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(optimization)
set.seed(1234)

#### Helper function that computes the true value of the objective function (generalized Lorenz) given some combination of linear policy parameters and eta, used for policy learning via simulated annealing ####
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
superpop <- read.csv("data_tau5.4.csv")
superpop$Y_cf <- (superpop$A==1)*(superpop$Y-superpop$tau) + (superpop$A==0)*(superpop$Y+superpop$tau)  # compute true counterfactuals

#### Learn true linear rule: pi=1 if LES>0 ####
## Specify parameters for simulated annealing
rate <- 0.85
n_limit <- 800
total <- ceiling(log(0.1/1000,rate))*n_limit  # for showing progress

## Run SA for superpopulation
for (alp in c(0.25,0.3,0.4,0.5,0.8,1)) {
  print(paste0("Learning true optimal policy for alpha = ", alp, "..."))
  counter <- 0
  obj_sa <- optim_sa(fun = obj_linear, start = c(0,0,0,0,0), lower = c(-5,-5,-5,-5,-5), upper = c(5,5,5,5,5), control = list(r = rate, nlimit = n_limit))
  print(paste0("alpha = ", alp, ", true Lorenz = ", -obj_sa$function_value))
}
