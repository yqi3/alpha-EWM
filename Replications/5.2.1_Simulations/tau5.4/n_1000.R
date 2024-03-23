######################################
# Author: Alice Yuan Qi
# Fan, Qi & Xu (2024), Table 5: n = 1000, 6 alpha values
# Simulation results using data generating process in Section 5.2.1; tau is specified as (5.4) in paper

# For exact reproducibility of random forests using grf, this script should be run on the following platform
# R version 4.1.1 (2021-08-10)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS release 6.5 (Final)

# For reasonable run time, we recommend running this script with at least 50 CPUs simultaneously, or breaking the for loop into smaller sub-tasks.
######################################

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(grf)
library(optimization)
source("5.2.1_Helper_Fns.R")

#### Read superpopulation data, specify number of MC trials and sample size ####
superpop <- read.csv("data_tau5.4.csv")
n_samples <- 1000
n <- 1000

#### Simulation ####
for (alp in c(0.25,0.3,0.4,0.5,0.8,1)) {
  ## Match population truth and specify parameters for simulated annealing
  if (alp == 0.25) {
    Lor_value_true <- 2.261504
    rate <- 0.6
    n_limit <- 5
  } else if (alp == 0.3) {
    Lor_value_true <- 2.751794
    rate <- 0.6
    n_limit <- 5
  } else if (alp == 0.4) {
    Lor_value_true <- 3.758295
    rate <- 0.6
    n_limit <- 3
  } else if (alp == 0.5) {
    Lor_value_true <- 4.795932
    rate <- 0.6
    n_limit <- 3
  } else if (alp == 0.8) {
    Lor_value_true <- 8.112018
    rate <- 0.6
    n_limit <- 3
  } else {
    Lor_value_true <- 10.62597
    rate <- 0.6
    n_limit <- 3
  }
  
  set.seed(1234)
  est_Lor <- numeric(n_samples)
  frac_trt <- numeric(n_samples)
  SEs <- numeric(n_samples)
  coverage <- numeric(n_samples)
  
  for (idx in 1:n_samples) {
    print(paste0("alpha = ", alp, ", n = ", n, ", sample ", idx))
    
    #### Prep Data ####
    rand <- sample.int(1e6, n, replace = F)
    datause <- superpop[rand,]
    datause$p <- 2/3  # randomized study with propensity score=2/3
    group <- make.cvgroup.balanced(datause, 2, 'A')  # for CF, K = 2
    datause_g1 <- datause[group==1,]
    datause_g2 <- datause[group==2,]
    
    #### Linear Rule: Optimization over pi Using Simulated Annealing ####
    counter <- 0
    total <- ceiling(log(0.1/1000,rate))*n_limit  # for showing SA progress
    print("Running SA for hat(pi)...")
    obj_sa <- optim_sa(fun = obj_linear, start = c(0,0,0,0,0,as.numeric(quantile(datause$Y,alp))), lower = c(-5,-5,-5,-5,-5,min(datause$Y)), upper = c(5,5,5,5,5,max(datause$Y)), control = list(r = rate, nlimit = n_limit))
    
    #### Estimation & Inference ####
    best_linear <- obj_sa$par
    d_linear <- best_linear[1]+best_linear[2]*datause$X1+best_linear[3]*datause$X2+best_linear[4]*datause$X3+best_linear[5]*datause$X4>0
    frac_trt[idx] <- mean(d_linear)
    datause$d <- d_linear
    
    ## Re-evaluation Using Best Parameters
    mu0 <- numeric(n)
    mu1 <- numeric(n)
    mu_g1 <- get_mu(best_linear[6], datause_g2)  # construct outcome nuisances using the other group
    mu0[group==1] <- predict(mu_g1[[1]], datause_g1[,c(1:4)])$predictions  # predict on current group
    mu1[group==1] <- predict(mu_g1[[2]], datause_g1[,c(1:4)])$predictions
    
    mu_g2 <- get_mu(best_linear[6], datause_g1)
    mu0[group==2] <- predict(mu_g2[[1]], datause_g2[,c(1:4)])$predictions
    mu1[group==2] <- predict(mu_g2[[2]], datause_g2[,c(1:4)])$predictions
    
    scores <- compute_scores(alp, datause$A, d_linear, datause$Y, mu1, mu0, datause$p, best_linear[6])
    Lor_linear <- mean(scores)
    est_Lor[idx] <- Lor_linear
    Lor_var <- var(scores)/n
    Lor_se <- sqrt(Lor_var)
    SEs[idx] <- Lor_se
    CI <- c(Lor_linear-1.96*Lor_se, Lor_linear+1.96*Lor_se)
    
    ## Coverage Rates
    if (Lor_value_true >= CI[1] & Lor_value_true <= CI[2]) {
      coverage[idx] <- 1
    }
    else {
      coverage[idx] <- 0
    }
    
    ## Print Results
    print(paste0("Current est Lor: ", est_Lor[idx]))
    print(paste0("Current var: ", Lor_var))
    print(paste0("Mean SE: ", mean(SEs[1:idx])))
    print(paste0("Median SE: ", median(SEs[1:idx])))
    print(paste0("Avg est Lor: ", mean(est_Lor[1:idx])))
    print(paste0("Avg est frac trt: ", mean(frac_trt[1:idx])))
    print(paste0("Bias: ", mean(est_Lor[1:idx])-Lor_value_true))
    print(paste0("Var: ", var(est_Lor[1:idx])))
    print(paste0("SD: ", sd(est_Lor[1:idx])))
    print(paste0("MSE: ", (mean(est_Lor[1:idx])-Lor_value_true)^2+var(est_Lor[1:idx])))
    print(paste0("Coverage: ", mean(coverage[1:idx])))
    print("----------------")
  }
  
  ## Save Results
  assign(paste0("result_n",n,"_alpha",alp), as.data.frame(cbind("n"=n, "alpha"=alp, "avg_frac_trt"=mean(frac_trt[1:idx]), "bias"=mean(est_Lor[1:idx])-Lor_value_true, "variance"=var(est_Lor[1:idx]), "MSE"=(mean(est_Lor[1:idx])-Lor_value_true)^2+var(est_Lor[1:idx]), "coverage"=mean(coverage[1:idx]))))
}

#### Show Final Results ####
result_n1000_alpha0.25
result_n1000_alpha0.3
result_n1000_alpha0.4
result_n1000_alpha0.5
result_n1000_alpha0.8
result_n1000_alpha1
