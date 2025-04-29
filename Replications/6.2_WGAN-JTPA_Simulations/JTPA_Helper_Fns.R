######################################
# Author: Alice Yuan Qi
# Helper functions for the empirical 
# application: the JTPA Study (1997)
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
  dataframe$Y <- (dataframe$earnings-q)*(dataframe$earnings<q)  # pseudo outcome for constructing nuisances
  X0 <- dataframe[dataframe$D==0, c(3,4)]
  Y0 <- dataframe$Y[dataframe$D==0]
  fmu0 <- regression_forest(X0, Y0)  # predict mu using RF
  X1 <- dataframe[dataframe$D==1, c(3,4)]
  Y1 <- dataframe$Y[dataframe$D==1]
  fmu1 <- regression_forest(X1, Y1)
  
  return(list(fmu0, fmu1))
}

#### Function that computes the doubly robust scores from alpha, sample data, treatment assignments, cross-fitted/true nuisance estimates, and eta ####
### "D" is equivalent to "A" in paper, and "d" is equivalent to ITR "pi"
compute_scores <- function(alp, D, d, y, mu1, mu0, p, q) {
  aug_mu0 <- (1-D)*((y-q)*(y<q)-mu0)/(1-p)
  aug1 <- d*D*((y-q)*(y<q)-mu1)/p
  aug2 <- -d*(1-D)*((y-q)*(y<q)-mu0)/(1-p)
  return(mu0*(1-d)+mu1*d+aug_mu0+aug1+aug2+alp*q)
}

#### Function that computes the value of the objective function (generalized Lorenz, can be easily converted to AVaR by dividing alpha) given some combination of linear policy parameters and eta, used for policy learning via simulated annealing ####
### We use K = 2 for ease of computation
obj_linear <- function(coef) {
  if (counter %% 100 == 0) {
    print(paste0("Progress: ", 100*counter/total, "%"))  # show progress
  }
  counter <<- counter+1
  
  d <- coef[1]+coef[2]*edu_std+coef[3]*prevearn_std>0  # ITR assignment
  q <- coef[4]  # eta in paper
  
  mu0 <- numeric(n)  # record CF mu0 outcomes
  mu1 <- numeric(n)  # record CF mu1 outcomes
  mu_g1 <- get_mu(q, jtpa_g2)  # construct outcome nuisances using the other group
  mu0[group==1] <- predict(mu_g1[[1]], jtpa_g1[,c(3,4)])$predictions  # predict on current group
  mu1[group==1] <- predict(mu_g1[[2]], jtpa_g1[,c(3,4)])$predictions
  
  mu_g2 <- get_mu(q, jtpa_g1)
  mu0[group==2] <- predict(mu_g2[[1]], jtpa_g2[,c(3,4)])$predictions
  mu1[group==2] <- predict(mu_g2[[2]], jtpa_g2[,c(3,4)])$predictions
  
  return((-1)*mean(compute_scores(alp, jtpa$D, d, jtpa$earnings, mu1, mu0, jtpa$p, q)))  # include negative sign for loss minimization
}

#### Function that computes the value of the objective function (generalized Lorenz, can be easily converted to AVaR by dividing alpha) given some combination of cubic policy parameters and eta, used for policy learning via simulated annealing ####
### We use K = 2 for ease of computation
obj_cubic <- function(coef) {
  if (counter %% 100 == 0) {
    print(paste0("Progress: ", 100*counter/total, "%"))  # show progress
  }
  counter <<- counter+1
  d <- coef[1]+coef[2]*edu_std+coef[3]*prevearn_std+coef[4]*edu_2_std+coef[5]*edu_3_std>0  # ITR assignment
  q <- coef[6]  # eta in paper
  
  mu0 <- numeric(n)  # record CF mu0 outcomes
  mu1 <- numeric(n)  # record CF mu1 outcomes
  mu_g1 <- get_mu(q, jtpa_g2)  # construct outcome nuisances using the other group
  mu0[group==1] <- predict(mu_g1[[1]], jtpa_g1[,c(3,4)])$predictions  # predict on current group
  mu1[group==1] <- predict(mu_g1[[2]], jtpa_g1[,c(3,4)])$predictions
  
  mu_g2 <- get_mu(q, jtpa_g1)
  mu0[group==2] <- predict(mu_g2[[1]], jtpa_g2[,c(3,4)])$predictions
  mu1[group==2] <- predict(mu_g2[[2]], jtpa_g2[,c(3,4)])$predictions
  
  return((-1)*mean(compute_scores(alp, jtpa$D, d, jtpa$earnings, mu1, mu0, jtpa$p, q)))  # include negative sign for loss minimization
}

#### Function that computes the value of the objective function (generalized Lorenz, can be easily converted to AVaR by dividing alpha) given some combination of (fixed) policy parameters and (variable) eta, used for policy evaluation via simulated annealing ####
obj_given_d <- function(q) {
  ## d is given outside of the function
  if (counter %% 100 == 0) {
    print(paste0("Progress: ", 100*counter/total, "%"))  # show progress
  }
  counter <<- counter+1
  
  mu0 <- numeric(n)  # record CF mu0 outcomes
  mu1 <- numeric(n)  # record CF mu1 outcomes
  mu_g1 <- get_mu(q, jtpa_g2)  # construct outcome nuisances using the other group
  mu0[group==1] <- predict(mu_g1[[1]], jtpa_g1[,c(3,4)])$predictions  # predict on current group
  mu1[group==1] <- predict(mu_g1[[2]], jtpa_g1[,c(3,4)])$predictions
  
  mu_g2 <- get_mu(q, jtpa_g1)
  mu0[group==2] <- predict(mu_g2[[1]], jtpa_g2[,c(3,4)])$predictions
  mu1[group==2] <- predict(mu_g2[[2]], jtpa_g2[,c(3,4)])$predictions
  
  return((-1)*mean(compute_scores(alp, jtpa$D, d, jtpa$earnings, mu1, mu0, jtpa$p, q)))
}

#### Objective function for the numerical delta method, where the policy is either treating no one or everyone ####
obj_numerical_delta_simple <- function(q) {
  if (counter %% 250 == 0) {
    print(paste0("Progress: ", 100*counter/total, "%"))  # show progress
  }
  counter <<- counter+1
  
  mu0 <- numeric(n)  # record CF mu0 outcomes
  mu1 <- numeric(n)  # record CF mu1 outcomes
  mu_g1 <- get_mu(q, jtpa_g2)  # construct outcome nuisances using the other group
  mu0[group==1] <- predict(mu_g1[[1]], jtpa_g1[,c(3,4)])$predictions  # predict on current group
  mu1[group==1] <- predict(mu_g1[[2]], jtpa_g1[,c(3,4)])$predictions
  
  mu_g2 <- get_mu(q, jtpa_g1)
  mu0[group==2] <- predict(mu_g2[[1]], jtpa_g2[,c(3,4)])$predictions
  mu1[group==2] <- predict(mu_g2[[2]], jtpa_g2[,c(3,4)])$predictions
  
  scores <- compute_scores(alp, jtpa$D, d, jtpa$earnings, mu1, mu0, jtpa$p, q)/alp
  G_hat <- n^(-1/2)*sum(rnorm(n)*(scores-mean(scores)))
  
  return((-1)*(mean(scores)+eps*G_hat))  # include negative sign for loss minimization
}

#### Objective function for the numerical delta method, where the policy is linear in edu and prevearn ####
obj_numerical_delta_linear <- function(coef) {
  if (counter %% 250 == 0) {
    print(paste0("Progress: ", 100*counter/total, "%"))  # show progress
  }
  counter <<- counter+1
  
  d <- coef[1]+coef[2]*edu_std+coef[3]*prevearn_std>0  # ITR assignment
  q <- coef[4]  # eta in paper
  
  mu0 <- numeric(n)  # record CF mu0 outcomes
  mu1 <- numeric(n)  # record CF mu1 outcomes
  mu_g1 <- get_mu(q, jtpa_g2)  # construct outcome nuisances using the other group
  mu0[group==1] <- predict(mu_g1[[1]], jtpa_g1[,c(3,4)])$predictions  # predict on current group
  mu1[group==1] <- predict(mu_g1[[2]], jtpa_g1[,c(3,4)])$predictions
  
  mu_g2 <- get_mu(q, jtpa_g1)
  mu0[group==2] <- predict(mu_g2[[1]], jtpa_g2[,c(3,4)])$predictions
  mu1[group==2] <- predict(mu_g2[[2]], jtpa_g2[,c(3,4)])$predictions
  
  scores <- compute_scores(alp, jtpa$D, d, jtpa$earnings, mu1, mu0, jtpa$p, q)/alp
  G_hat <- n^(-1/2)*sum(rnorm(n)*(scores-mean(scores)))
  
  return((-1)*(mean(scores)+eps*G_hat))  # include negative sign for loss minimization
}

#### Objective function for the numerical delta method, where the policy is linear in edu, prevearn, edu^2, and edu^3 ####
obj_numerical_delta_cubic <- function(coef) {
  if (counter %% 250 == 0) {
    print(paste0("Progress: ", 100*counter/total, "%"))  # show progress
  }
  counter <<- counter+1
  
  d <- coef[1]+coef[2]*edu_std+coef[3]*prevearn_std+coef[4]*edu_2_std+coef[5]*edu_3_std>0  # ITR assignment
  q <- coef[6]  # eta in paper
  
  mu0 <- numeric(n)  # record CF mu0 outcomes
  mu1 <- numeric(n)  # record CF mu1 outcomes
  mu_g1 <- get_mu(q, jtpa_g2)  # construct outcome nuisances using the other group
  mu0[group==1] <- predict(mu_g1[[1]], jtpa_g1[,c(3,4)])$predictions  # predict on current group
  mu1[group==1] <- predict(mu_g1[[2]], jtpa_g1[,c(3,4)])$predictions
  
  mu_g2 <- get_mu(q, jtpa_g1)
  mu0[group==2] <- predict(mu_g2[[1]], jtpa_g2[,c(3,4)])$predictions
  mu1[group==2] <- predict(mu_g2[[2]], jtpa_g2[,c(3,4)])$predictions
  
  scores <- compute_scores(alp, jtpa$D, d, jtpa$earnings, mu1, mu0, jtpa$p, q)/alp
  G_hat <- n^(-1/2)*sum(rnorm(n)*(scores-mean(scores)))
  
  return((-1)*(mean(scores)+eps*G_hat))  # include negative sign for loss minimization
}