######################################
# Author: Alice Yuan Qi
# Tables 2 & 3 in Fan, Qi & Xu (2024)
# Estimated generalized Lorenz for different actual alpha's of interest and alpha's for ITR selection, standard errors are reported as well

# For exact reproducibility of random forests using grf, this script should be run on the following platform
# R version 4.1.1 (2021-08-10)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS release 6.5 (Final)

# For reasonable run time, we recommend running this script with at least 20 CPUs simultaneously. In practice (i.e., not in this particular example), however, researchers have the discretion to adjust the parameters used in simulated annealing as long as there is sufficient evidence of convergence, or use other gradient-free optimization methods.
######################################

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(grf)
library(optimization)
set.seed(1)
source("JTPA_Helper_Fns.R")

#### Prep Data ####
path <- "Tables_1&7_ITR_para/"  # path for reading policy parameters
path_result <- "Tables_2&3_comparisons/"  # path for saving results
jtpa <- read.csv("KT_JTPA.csv")
jtpa$p <- 2/3  # randomized study with propensity score=2/3
n <- nrow(jtpa)
group <- make.cvgroup.balanced(jtpa, 2, 'D')  # K = 2

comparisons_linear <- matrix(0, ncol = 6, nrow = 12)
colnames(comparisons_linear) <- c("0.25", "0.3", "0.4", "0.5", "0.8", "1")
rownames(comparisons_linear) <- c("0.25", "0.25 se", "0.3", "0.3 se", "0.4", "0.4 se", "0.5", "0.5 se", "0.8", "0.8 se", "1", "1 se")

comparisons_cubic <- matrix(0, ncol = 6, nrow = 12)
colnames(comparisons_cubic) <- c("0.25", "0.3", "0.4", "0.5", "0.8", "1")
rownames(comparisons_cubic) <- c("0.25", "0.25 se", "0.3", "0.3 se", "0.4", "0.4 se", "0.5", "0.5 se", "0.8", "0.8 se", "1", "1 se")

jtpa_g1 <- jtpa[group==1,]  # for CF
jtpa_g2 <- jtpa[group==2,]  # for CF
rate <- 0.7
n_limit <- 2
total <- ceiling(log(0.1/1000,rate))*n_limit  # for showing progress

#### Run Comparisons ####
for (alp_ITR in c(0.25, 0.3, 0.4, 0.5, 0.8, 1)) {
  print("------------")
  print(paste0("Running alp_ITR = ", alp_ITR))
  linear_para <- read.csv(paste0(path, alp_ITR,"_linear_para.csv"))
  cubic_para <- read.csv(paste0(path, alp_ITR,"_cubic_para.csv"))
  
  ### Linear rule
  print("Running linear rules...")
  b0 <- linear_para$b0
  b1 <- linear_para$b1
  b2 <- linear_para$b2
  d <- b0+b1*jtpa$edu+b2*jtpa$prevearn>0
  
  for (alp_actual in c(0.25, 0.3, 0.4, 0.5, 0.8, 1)) {
    print(paste0("Running alp_actual = ", alp_actual))
    
    if (alp_actual==alp_ITR) {
      print(linear_para$Lor_linear)
      comparisons_linear[toString(alp_actual), toString(alp_ITR)] <- round(linear_para$Lor_linear,3)  # fills diagonal element
      comparisons_linear[paste0(toString(alp_actual), " se"), toString(alp_ITR)] <- round(linear_para$Lor_se,3)  # fills diagonal se
    } else {
      counter <- 0
      alp <- alp_actual
      obj_sa <- optim_sa(fun = obj_given_d, start = as.numeric(quantile(jtpa$earnings,alp_actual)), lower = min(jtpa$earnings), upper = max(jtpa$earnings), control = list(r = rate, nlimit = n_limit))
      q <- obj_sa$par
      mu0 <- numeric(n)  # record CF mu0 outcomes
      mu1 <- numeric(n)  # record CF mu1 outcomes
      mu_g1 <- get_mu(q, jtpa_g2)  # construct outcome nuisances using the other group
      mu0[group==1] <- predict(mu_g1[[1]], jtpa_g1[,c(3,4)])$predictions  # predict on current group
      mu1[group==1] <- predict(mu_g1[[2]], jtpa_g1[,c(3,4)])$predictions
      
      mu_g2 <- get_mu(q, jtpa_g1)
      mu0[group==2] <- predict(mu_g2[[1]], jtpa_g2[,c(3,4)])$predictions
      mu1[group==2] <- predict(mu_g2[[2]], jtpa_g2[,c(3,4)])$predictions
      
      scores <- compute_scores(alp_actual, jtpa$D, d, jtpa$earnings, mu1, mu0, jtpa$p, q)
      Lor_linear <- mean(scores)
      print(paste0("Lor_linear: ", Lor_linear))
      Lor_var <- var(scores)/n
      Lor_se <- sqrt(Lor_var)
      
      comparisons_linear[toString(alp_actual), toString(alp_ITR)] <- round(Lor_linear,3)  # fills estimate
      comparisons_linear[paste0(toString(alp_actual), " se"), toString(alp_ITR)] <- round(Lor_se,3)  # fills se
    }
  }
  
  ### Cubic rule
  print("Running cubic rules...")
  b0 <- cubic_para$b0
  b1 <- cubic_para$b1
  b2 <- cubic_para$b2
  b3 <- cubic_para$b3
  b4 <- cubic_para$b4
  d <- b0+b1*jtpa$edu+b2*jtpa$prevearn+b3*(jtpa$edu)^2+b4*(jtpa$edu)^3>0
  
  for (alp_actual in c(0.25, 0.3, 0.4, 0.5, 0.8, 1)) {
    print(paste0("Running alp_actual = ", alp_actual))
    
    if (alp_actual==alp_ITR) {
      print(cubic_para$Lor_cubic)
      comparisons_cubic[toString(alp_actual), toString(alp_ITR)] <- round(cubic_para$Lor_cubic,3)  # fills diagonal element
      comparisons_cubic[paste0(toString(alp_actual), " se"), toString(alp_ITR)] <- round(cubic_para$Lor_se,3)  # fills diagonal se
    } else {
      counter <- 0
      alp <- alp_actual
      obj_sa <- optim_sa(fun = obj_given_d, start = as.numeric(quantile(jtpa$earnings,alp_actual)), lower = min(jtpa$earnings), upper = max(jtpa$earnings), control = list(r = rate, nlimit = n_limit))
      q <- obj_sa$par
      mu0 <- numeric(n)  # record CF mu0 outcomes
      mu1 <- numeric(n)  # record CF mu1 outcomes
      mu_g1 <- get_mu(q, jtpa_g2)  # construct outcome nuisances using the other group
      mu0[group==1] <- predict(mu_g1[[1]], jtpa_g1[,c(3,4)])$predictions  # predict on current group
      mu1[group==1] <- predict(mu_g1[[2]], jtpa_g1[,c(3,4)])$predictions
      
      mu_g2 <- get_mu(q, jtpa_g1)
      mu0[group==2] <- predict(mu_g2[[1]], jtpa_g2[,c(3,4)])$predictions
      mu1[group==2] <- predict(mu_g2[[2]], jtpa_g2[,c(3,4)])$predictions
      
      scores <- compute_scores(alp_actual, jtpa$D, d, jtpa$earnings, mu1, mu0, jtpa$p, q)
      Lor_cubic <- mean(scores)
      print(paste0("Lor_cubic: ", Lor_cubic))
      Lor_var <- var(scores)/n
      Lor_se <- sqrt(Lor_var)
      
      comparisons_cubic[toString(alp_actual), toString(alp_ITR)] <- round(Lor_cubic,3)  # fills estimate
      comparisons_cubic[paste0(toString(alp_actual), " se"), toString(alp_ITR)] <- round(Lor_se,3)  # fills se
    }
  }
}

write.csv(comparisons_linear, paste0(path_result, "comparisons_linear.csv"))
write.csv(comparisons_cubic, paste0(path_result, "comparisons_cubic.csv"))
