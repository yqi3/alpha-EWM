######################################
# Author: Alice Yuan Qi
# Fan, Qi & Xu (2024), Table 6: n = 2000, 6 alpha values
# Simulation results using samples drawn from WGAN-JTPA

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
source("JTPA_Helper_Fns.R")

#### Read superpopulation data, specify number of MC trials, sample size, and parameters for simulated annealing ####
superpop <- read.csv("WGAN-JTPA.csv")
n_samples <- 1000
n <- 2000
rate <- 0.6
n_limit <- 5
total <- ceiling(log(0.1/1000,rate))*n_limit  # for showing SA progress

#### Simulations ####
for (alp in c(0.25,0.3,0.4,0.5,0.8,1)) {
  ## Match Population Truth
  if (alp == 0.25) {
    Lor_value_true <- 279.7646
  } else if (alp == 0.3) {
    Lor_value_true <- 572.3733
  } else if (alp == 0.4) {
    Lor_value_true <- 1384.31
  } else if (alp == 0.5) {
    Lor_value_true <- 2433.763
  } else if (alp == 0.8) {
    Lor_value_true <- 7580.301
  } else {
    Lor_value_true <- 14643.49
  }
  
  set.seed(123)
  est_Lor <- numeric(n_samples)
  frac_trt <- numeric(n_samples)
  SEs <- numeric(n_samples)
  coverage <- numeric(n_samples)
  
  for (idx in 1:n_samples) {
    print(paste0("alpha = ", alp, ", n = ", n, ", sample ", idx))
    
    #### Prep Data ####
    rand <- sample.int(1e6, n, replace = F)
    jtpa <- superpop[rand,]
    jtpa$p <- mean(superpop$D)  # true population propensity score
    group <- make.cvgroup.balanced(jtpa, 2, 'D')  # for CF, K = 2
    jtpa_g1 <- jtpa[group==1,]
    jtpa_g2 <- jtpa[group==2,]
    
    #### Linear Rule ####
    edu_std <- as.numeric(scale(jtpa$edu))
    prevearn_std <- as.numeric(scale(jtpa$prevearn))
    
    #### Optimization over pi Using Simulated Annealing ####
    counter <- 0
    print("Running SA for hat(pi)...")
    obj_sa <- optim_sa(fun = obj_linear, start = c(0,0,0,as.numeric(quantile(jtpa$earnings,alp))), lower = c(-5,-5,-5,min(jtpa$earnings)), upper = c(5,5,5,max(jtpa$earnings)), control = list(r = rate, nlimit = n_limit))
    
    #### Estimation & Inference ####
    best_linear <- obj_sa$par
    d_linear <- best_linear[1]+best_linear[2]*edu_std+best_linear[3]*prevearn_std>0
    frac_trt[idx] <- mean(d_linear)
    jtpa$d <- d_linear
    
    ## Re-evaluation Using Best Parameters
    mu0 <- numeric(n)
    mu1 <- numeric(n)
    mu_g1 <- get_mu(best_linear[4], jtpa_g2)  # construct outcome nuisances using the other group
    mu0[group==1] <- predict(mu_g1[[1]], jtpa_g1[,c(2,3)])$predictions  # predict on current group
    mu1[group==1] <- predict(mu_g1[[2]], jtpa_g1[,c(2,3)])$predictions
    
    mu_g2 <- get_mu(best_linear[4], jtpa_g1)
    mu0[group==2] <- predict(mu_g2[[1]], jtpa_g2[,c(2,3)])$predictions
    mu1[group==2] <- predict(mu_g2[[2]], jtpa_g2[,c(2,3)])$predictions
    
    scores <- compute_scores(alp, jtpa$D, d_linear, jtpa$earnings, mu1, mu0, jtpa$p, best_linear[4])
    Lor_linear <- mean(scores)
    est_Lor[idx] <- Lor_linear
    Lor_var <- var(scores)/n
    Lor_se <- sqrt(Lor_var)
    SEs[idx] <- Lor_se
    CI <- c(Lor_linear-1.96*Lor_se, Lor_linear+1.96*Lor_se)
    
    ## Coverage Rates
    if (Lor_value_true >= CI[1] & Lor_value_true <= CI[2]) {
      coverage[idx] <- 1
    } else {
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
result_n2000_alpha0.25
result_n2000_alpha0.3
result_n2000_alpha0.4
result_n2000_alpha0.5
result_n2000_alpha0.8
result_n2000_alpha1
