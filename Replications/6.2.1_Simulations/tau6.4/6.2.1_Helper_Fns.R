######################################
# Author: Alice Yuan Qi
# Fan, Qi & Xu (2025)
# Helper functions for simulations in 
# Section 6.2.1
######################################

#### Functions that split data for cross-fitting ####
### make.cvgroup() and make.cvgroup.balanced() below are adopted from Kallus (2023)
make.cvgroup <- function(n, K, right = TRUE) {
  split <- runif(n)
  return(as.numeric(cut(split, quantile(split, probs = seq(0, 1, 1/K)), include.lowest = TRUE, right = right)))
}

make.cvgroup.balanced <- function(data, K, form_t) {
  cvgroup <- numeric(nrow(data))
  cvgroup[data[[form_t]]==1] <- make.cvgroup(sum(data[[form_t]]==1), K, right = TRUE)
  cvgroup[data[[form_t]]==0] <- make.cvgroup(sum(data[[form_t]]==0), K, right = FALSE)
  return(cvgroup)
}

#### Function that builds nuisance outcome models with a given data fold and some value q (eta in paper) ####
get_mu <- function(q, dataframe) {
  dataframe$Y <- (dataframe$Y-q)*(dataframe$Y<q)  # pseudo outcome for constructing nuisances
  X0 <- dataframe[dataframe$A==0, c(1:4)]
  Y0 <- dataframe$Y[dataframe$A==0]
  fmu0 <- regression_forest(X0, Y0)
  X1 <- dataframe[dataframe$A==1, c(1:4)]
  Y1 <- dataframe$Y[dataframe$A==1]
  fmu1 <- regression_forest(X1, Y1)
  return(list(fmu0, fmu1))
}

#### Function that computes the doubly robust scores from alpha, sample data, treatment assignments, cross-fitted/true nuisance estimates, and eta ####
### "d" is equivalent to ITR "pi"
compute_scores <- function(alp, A, d, y, mu1, mu0, p, q) {
  aug_mu0 <- (1-A)*((y-q)*(y<q)-mu0)/(1-p)
  aug1 <- d*A*((y-q)*(y<q)-mu1)/p
  aug2 <- -d*(1-A)*((y-q)*(y<q)-mu0)/(1-p)
  return(mu0*(1-d)+mu1*d+aug_mu0+aug1+aug2+alp*q)
}

#### Function that computes the value of the objective function (generalized Lorenz, can be easily converted to AVaR by dividing alpha) given some combination of linear policy parameters and eta, used for policy learning via simulated annealing ####
### We use K = 2 for ease of computation
obj_linear <- function(coef) {
  if (counter %% 50 == 0) {
    print(paste0("Progress: ", 100*counter/total, "%"))  # show progress
  }
  counter <<- counter+1
  d <- coef[1]+coef[2]*datause$X1+coef[3]*datause$X2+coef[4]*datause$X3+coef[5]*datause$X4>0  # ITR assignment
  q <- coef[6]  # eta in paper
  
  mu0 <- numeric(n)  # record CF mu0 outcomes
  mu1 <- numeric(n)  # record CF mu1 outcomes
  mu_g1 <- get_mu(q, datause_g2)  # construct outcome nuisances using the other group
  mu0[group==1] <- predict(mu_g1[[1]], datause_g1[,c(1:4)])$predictions  # predict on current group
  mu1[group==1] <- predict(mu_g1[[2]], datause_g1[,c(1:4)])$predictions
  
  mu_g2 <- get_mu(q, datause_g1)
  mu0[group==2] <- predict(mu_g2[[1]], datause_g2[,c(1:4)])$predictions
  mu1[group==2] <- predict(mu_g2[[2]], datause_g2[,c(1:4)])$predictions
  
  return((-1)*mean(compute_scores(alp, datause$A, d, datause$Y, mu1, mu0, datause$p, q)))  # include negative sign for loss minimization
}
