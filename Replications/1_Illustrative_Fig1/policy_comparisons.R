######################################
# Author: Alice Yuan Qi
# Fan, Qi & Xu (2025), Figure 1
# Population-level comparisons between the 0.1-EWM, 0.1-quantile-optimal, equality-minded, and 1-EWM policies using data generating process in Section 1
######################################

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(optimization)
library(DescTools)  # for computing Gini
library(ggplot2)

superpop <- read.csv("data.csv")

#### Helper function that computes the true value of the objective function (alpha-AVaR) given some combination of policy parameters and eta, used for policy learning via simulated annealing ####
obj_fn <- function(coef) {
  if (counter %% 1000 == 0) {
    print(paste0("Progress: ", 100*counter/total, "%"))  # show progress
  }
  counter <<- counter+1
  d_pop <- superpop$x<=coef
  outcomes <- (superpop$A==d_pop)*superpop$y+(superpop$A!=d_pop)*superpop$y_cf
  return((-1)*mean(outcomes[outcomes<=quantile(outcomes,alp)]))  # include negative sign for loss minimization
}

#### Learn true cutoff rule: pi=1 if x>threshold (can skip this part by directly reading best_AVaR_policies.csv) ####
## Specify parameters for simulated annealing
rate <- 0.8
n_limit <- 100
total <- ceiling(log(0.1/1000,rate))*n_limit  # for showing progress

## Df for storing results
results <- matrix(NA, ncol = 2, nrow = 2)
colnames(results) <- c("thres", "AVaR")
rownames(results) <- c("0.1", "1")

## Run SA for superpopulation
for (alp in c(0.1,1)) {
  print(paste0("Learning true optimal policy for alpha = ", alp, "..."))
  set.seed(1234)
  counter <- 0
  
  obj_sa <- optim_sa(fun = obj_fn, start = 0.5, lower = 0, upper = 1, control = list(r = rate, nlimit = n_limit))
  print(paste0("alpha = ", alp, ", true AVaR = ", -obj_sa$function_value))
  print(paste0("thres = ", obj_sa$par))
  results[toString(alp),"thres"] <- obj_sa$par
  results[toString(alp),"AVaR"] <- -obj_sa$function_value
  print("----------")
}

write.csv(results, "best_AVaR_policies.csv", row.names = F)

best_rules <- read.csv("best_AVaR_policies.csv")
best_ITR_0.1 <- superpop$x <= best_rules$thres[1]
best_ITR_1 <- superpop$x <= best_rules$thres[2]

#### Helper function that computes the true value of the objective function (Gini SWF) given some policy parameter, used for policy learning via simulated annealing. Computation of the Gini SWF is based on Equation (5) in Kitagawa & Tetenov (2021). ####
### "d_pop" is equivalent to policy "pi"
obj_EM <- function(cutoff) {
  if (counter %% 1000 == 0) {
    print(paste0("Progress: ", 100*counter/total, "%"))  # show progress
  }
  counter <<- counter+1
  d_pop <- superpop$x<=cutoff
  outcomes <- (superpop$A==d_pop)*superpop$y+(superpop$A!=d_pop)*superpop$y_cf
  Gini_SWF <- mean(outcomes)*(1-Gini(outcomes))
  return((-1)*Gini_SWF)  # include negative sign for loss minimization
}

#### Helper function that computes the true value of the objective function (outcome quantile) given some policy parameter, used for policy learning via simulated annealing. Quantile-optimal policies are studied by Want et al. (2018). ####
### "d_pop" is equivalent to policy "pi"
obj_quantile <- function(cutoff) {
  if (counter %% 1000 == 0) {
    print(paste0("Progress: ", 100*counter/total, "%"))  # show progress
  }
  counter <<- counter+1
  d_pop <- superpop$x<=cutoff
  outcomes <- (superpop$A==d_pop)*superpop$y+(superpop$A!=d_pop)*superpop$y_cf
  
  return((-1)*quantile(outcomes, 0.1))  # include negative sign for loss minimization
}

#### Learn true equality-minded cutoff rule: pi=1 if x<=cutoff ####
## Specify parameters for simulated annealing
rate <- 0.8
n_limit <- 50
counter <- 0
total <- ceiling(log(0.1/1000,rate))*n_limit  # for showing progress

set.seed(123)
obj_sa <- optim_sa(fun = obj_EM, start = mean(superpop$x), lower = min(superpop$x), upper = max(superpop$x), control = list(r = rate, nlimit = n_limit))
print(obj_sa$par)

best_ITR_EM <- superpop$x <= obj_sa$par
print(paste0("mean(best_ITR_EM): ", mean(best_ITR_EM)))

#### Learn true quantile-optimal cutoff rule: pi=1 if x<=cutoff ####
## Specify parameters for simulated annealing
set.seed(123)
counter <- 0
obj_sa <- optim_sa(fun = obj_quantile, start = mean(superpop$x), lower = min(superpop$x), upper = max(superpop$x), control = list(r = rate, nlimit = n_limit))
print(obj_sa$par)

best_ITR_quantile <- superpop$x <= obj_sa$par
print(paste0("mean(best_ITR_quantile): ", mean(best_ITR_quantile)))

#### Density Plot ####
outcomes_AVaR <- (superpop$A==best_ITR_0.1)*superpop$y+(superpop$A!=best_ITR_0.1)*superpop$y_cf
outcomes_mean <- (superpop$A==best_ITR_1)*superpop$y+(superpop$A!=best_ITR_1)*superpop$y_cf
outcomes_EM <- (superpop$A==best_ITR_EM)*superpop$y+(superpop$A!=best_ITR_EM)*superpop$y_cf
outcomes_quantile <- (superpop$A==best_ITR_quantile)*superpop$y+(superpop$A!=best_ITR_quantile)*superpop$y_cf

jpeg(file="density_all.jpeg", width=6, height=5, units="in", res=300)
plot(density(outcomes_AVaR), col="blue", lwd=2, xlab="Post-treatment Y", main = "")
lines(density(outcomes_quantile), col="red", lwd=2)
lines(density(outcomes_EM), col="orange", lwd=2)
lines(density(outcomes_mean), col="green", lwd=2)
legend(x="topright", legend=c("0.1-EWM", "0.1-quantile", "1-EWM", "Equality-minded"),
       col=c("blue", "red", "green", "orange"), lty=1, cex=0.8, lwd=3, y.intersp=1)
dev.off()
