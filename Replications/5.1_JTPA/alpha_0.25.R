######################################
# Author: Alice Yuan Qi
# Table 1 in Fan, Qi & Xu (2024): alpha = 0.25
# Optimized Lorenz function values and 95%-CIs under simple/linear/linear with edu^2 and edu^3 ITR classes

# For exact reproducibility of random forests using grf, this script should be run on the following platform
# R version 4.1.1 (2021-08-10)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS release 6.5 (Final)
######################################

rm(list=ls())
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(optimization)
library(grf)
set.seed(1234)
source("JTPA_Helper_Fns.R")

#### Prep Data ####
path <- "Tables_1&7_ITR_para/"  # path for saving results
path_training_loss <- "training_losses/"  # path for saving loss plots
jtpa <- read.csv("KT_JTPA.csv")
jtpa$p <- 2/3  # randomized study with known propensity score=2/3
n <- nrow(jtpa)
group <- make.cvgroup.balanced(jtpa, 2, 'D')  # K = 2
jtpa_g1 <- jtpa[group==1,]  # for CF
jtpa_g2 <- jtpa[group==2,]  # for CF


#### --- Targeted policy learning based on generalized Lorenz --- ####
alp <- 0.25

#### Simple rule: d=1 or 0 for all ("d" is equivalent to "pi" in paper) ####
print("Evaluating simple ITRs (first column in Table 1)...")
rate <- 0.6
n_limit <- 50
total <- ceiling(log(0.1/1000,rate))*n_limit  # for showing progress

## Treat all
print("If treating all...")
d <- rep(1,n)

# Run SA to obtain q (eta in paper)
counter <- 0
obj_sa <- optim_sa(fun = obj_given_d, start = as.numeric(quantile(jtpa$earnings,alp)), lower = min(jtpa$earnings), upper = max(jtpa$earnings), control = list(r = rate, nlimit = n_limit))
q <- obj_sa$par

mu0 <- numeric(n)  # record CF mu0 outcomes
mu1 <- numeric(n)  # record CF mu1 outcomes
mu_g1 <- get_mu(q, jtpa_g2)  # construct outcome nuisances using the other group
mu0[group==1] <- predict(mu_g1[[1]], jtpa_g1[,c(3,4)])$predictions  # predict on current group
mu1[group==1] <- predict(mu_g1[[2]], jtpa_g1[,c(3,4)])$predictions

mu_g2 <- get_mu(q, jtpa_g1)
mu0[group==2] <- predict(mu_g2[[1]], jtpa_g2[,c(3,4)])$predictions
mu1[group==2] <- predict(mu_g2[[2]], jtpa_g2[,c(3,4)])$predictions

scores <- compute_scores(alp, jtpa$D, d, jtpa$earnings, mu1, mu0, jtpa$p, q)
Lor_value <- mean(scores)
print(paste0("Lor_value: ", Lor_value))
Lor_var <- var(scores)/n
Lor_se <- sqrt(Lor_var)
simple_CI <- c(Lor_value-1.96*Lor_se, Lor_value+1.96*Lor_se)
print("simple_CI: ")
simple_CI

## Treat none
print("If treating none...")
d <- rep(0,n)

# Run SA to obtain q
counter <- 0
obj_sa <- optim_sa(fun = obj_given_d, start = as.numeric(quantile(jtpa$earnings,alp)), lower = min(jtpa$earnings), upper = max(jtpa$earnings), control = list(r = rate, nlimit = n_limit))
q <- obj_sa$par

mu0 <- numeric(n)  # record CF mu0 outcomes
mu1 <- numeric(n)  # record CF mu1 outcomes
mu_g1 <- get_mu(q, jtpa_g2)  # construct outcome nuisances using the other group
mu0[group==1] <- predict(mu_g1[[1]], jtpa_g1[,c(3,4)])$predictions  # predict on current group
mu1[group==1] <- predict(mu_g1[[2]], jtpa_g1[,c(3,4)])$predictions

mu_g2 <- get_mu(q, jtpa_g1)
mu0[group==2] <- predict(mu_g2[[1]], jtpa_g2[,c(3,4)])$predictions
mu1[group==2] <- predict(mu_g2[[2]], jtpa_g2[,c(3,4)])$predictions

scores <- compute_scores(alp, jtpa$D, d, jtpa$earnings, mu1, mu0, jtpa$p, q)
Lor_value_baseline <- mean(scores)
print(paste0("Lor_value_baseline: ", Lor_value_baseline))
Lor_var_baseline <- var(scores)/n
Lor_se_baseline <- sqrt(Lor_var)
simple_CI_baseline <- c(Lor_value_baseline-1.96*Lor_se_baseline, Lor_value_baseline+1.96*Lor_se_baseline)
print("simple_CI_baseline: ")
simple_CI_baseline

#### Linear rule: d=1 if b0+b1*edu_std+b2*prevearn_std>0 ####
# Standardize edu & prevearn so that each has mean=0 and var=1
edu_std <- as.numeric(scale(jtpa$edu))
prevearn_std <- as.numeric(scale(jtpa$prevearn))

print("Learning linear ITRs (second column in Table 1)...")
set.seed(1234)
rate <- 0.75
n_limit <- 100
total <- ceiling(log(0.1/1000,rate))*n_limit  # for showing progress
counter <- 0

# Run simulated annealing for optimization
obj_sa <- optim_sa(fun = obj_linear, start = c(0,0,0,as.numeric(quantile(jtpa$earnings,alp))), lower = c(-5,-5,-5,min(jtpa$earnings)), upper = c(5,5,5,max(jtpa$earnings)), control = list(r = rate, nlimit = n_limit), trace = T)
pdf(file=paste0(path_training_loss, alp, "linear.pdf"), width=8, height=4)
plot(obj_sa)  # can inspect if the loss function shows convergence
dev.off()

## Re-evaluate welfare for the best policy parameters
best_linear <- obj_sa$par
b0 <- best_linear[1]-best_linear[2]/sd(jtpa$edu)*mean(jtpa$edu)-best_linear[3]/sd(jtpa$prevearn)*mean(jtpa$prevearn)
b1 <- best_linear[2]/sd(jtpa$edu)
b2 <- best_linear[3]/sd(jtpa$prevearn)
q <- best_linear[4]
print("Linear coefficients: ")
c(b0,b1,b2,q)

d_linear <- b0+b1*jtpa$edu+b2*jtpa$prevearn>0
print(paste0("mean(d_linear): ", mean(d_linear)))

mu0 <- numeric(n)  # record CF mu0 outcomes
mu1 <- numeric(n)  # record CF mu1 outcomes

mu_g1 <- get_mu(q, jtpa_g2)  # construct outcome nuisances using the other group
mu0[group==1] <- predict(mu_g1[[1]], jtpa_g1[,c(3,4)])$predictions  # predict on current group
mu1[group==1] <- predict(mu_g1[[2]], jtpa_g1[,c(3,4)])$predictions

mu_g2 <- get_mu(q, jtpa_g1)
mu0[group==2] <- predict(mu_g2[[1]], jtpa_g2[,c(3,4)])$predictions
mu1[group==2] <- predict(mu_g2[[2]], jtpa_g2[,c(3,4)])$predictions

scores <- compute_scores(alp, jtpa$D, d_linear, jtpa$earnings, mu1, mu0, jtpa$p, q)
Lor_linear <- mean(scores)
print(paste0("Lor_linear: ", Lor_linear))
Lor_var <- var(scores)/n
Lor_se <- sqrt(Lor_var)
Lor_linear_CI <- c(Lor_linear-1.96*Lor_se, Lor_linear+1.96*Lor_se)
print("Lor_linear_CI:")
Lor_linear_CI

print("Normalized coefficients: ")
if (b2 < 0) {
  print(c(-b0/b2,-b1/b2,-1))
  write.csv(data.frame("b0"=-b0/b2, "b1"=-b1/b2, "b2"=-1, "q"=q, "Lor_linear"=Lor_linear, "Lor_se"=Lor_se, "CI"=paste0("(",toString(round(Lor_linear_CI,3)),")")), paste0(path,alp,"_linear_para.csv"), row.names = F)  # for plotting
} else {
  print(c(b0/b2,b1/b2,1))
  write.csv(data.frame("b0"=b0/b2, "b1"=b1/b2, "b2"=1, "q"=q, "Lor_linear"=Lor_linear, "Lor_se"=Lor_se, "CI"=paste0("(",toString(round(Lor_linear_CI,3)),")")), paste0(path,alp,"_linear_para.csv"), row.names = F)
}

#### Linear with cubic terms: d=1 if b0+b1*edu_std+b2*prevearn_std+b3*edu_2_std+b4*edu_3_std>0 ####
# standardize edu^2 & edu^3 so that each has mean=0 and var=1
edu_2_std <- as.numeric(scale((jtpa$edu)^2))
edu_3_std <- as.numeric(scale((jtpa$edu)^3))

print("Learning linear ITRs with squared and cubic edu (third column in Table 1)...")
set.seed(1234)
rate <- 0.8
n_limit <- 200
total <- ceiling(log(0.1/1000,rate))*n_limit  # for showing progress
counter <- 0

# Run simulated annealing for optimization
obj_sa <- optim_sa(fun = obj_cubic, start = c(0,0,0,0,0,as.numeric(quantile(jtpa$earnings,alp))), lower = c(-5,-5,-5,-5,-5,min(jtpa$earnings)), upper = c(5,5,5,5,5,max(jtpa$earnings)), control = list(r = rate, nlimit = n_limit), trace = T)
pdf(file=paste0(path_training_loss, alp, "cubic.pdf"), width=8, height=4)
plot(obj_sa)  # can inspect if the loss function shows convergence
dev.off()

## Re-evaluate welfare for the best policy parameters
best_cubic <- obj_sa$par
b0 <- best_cubic[1]-best_cubic[2]/sd(jtpa$edu)*mean(jtpa$edu)-best_cubic[3]/sd(jtpa$prevearn)*mean(jtpa$prevearn)-best_cubic[4]/sd((jtpa$edu)^2)*mean((jtpa$edu)^2)-best_cubic[5]/sd((jtpa$edu)^3)*mean((jtpa$edu)^3)
b1 <- best_cubic[2]/sd(jtpa$edu)
b2 <- best_cubic[3]/sd(jtpa$prevearn)
b3 <- best_cubic[4]/sd((jtpa$edu)^2)
b4 <- best_cubic[5]/sd((jtpa$edu)^3)
q <- best_cubic[6]
print("Cubic coefficients: ")
c(b0,b1,b2,b3,b4,q)

d_cubic <- b0+b1*jtpa$edu+b2*jtpa$prevearn+b3*((jtpa$edu)^2)+b4*((jtpa$edu)^3)>0
print(paste0("mean(d_cubic): ", mean(d_cubic)))

mu0 <- numeric(n)  # record CF mu0 outcomes
mu1 <- numeric(n)  # record CF mu1 outcomes

mu_g1 <- get_mu(q, jtpa_g2)  # construct outcome nuisances using the other group
mu0[group==1] <- predict(mu_g1[[1]], jtpa_g1[,c(3,4)])$predictions  # predict on current group
mu1[group==1] <- predict(mu_g1[[2]], jtpa_g1[,c(3,4)])$predictions

mu_g2 <- get_mu(q, jtpa_g1)
mu0[group==2] <- predict(mu_g2[[1]], jtpa_g2[,c(3,4)])$predictions
mu1[group==2] <- predict(mu_g2[[2]], jtpa_g2[,c(3,4)])$predictions

scores <- compute_scores(alp, jtpa$D, d_cubic, jtpa$earnings, mu1, mu0, jtpa$p, q)
Lor_cubic <- mean(scores)
print(paste0("Lor_cubic: ", Lor_cubic))
Lor_var <- var(scores)/n
Lor_se <- sqrt(Lor_var)
Lor_cubic_CI <- c(Lor_cubic-1.96*Lor_se, Lor_cubic+1.96*Lor_se)
print("Lor_cubic_CI:")
Lor_cubic_CI

print("Normalized coefficients: ")
if (b2 < 0) {
  print(c(-b0/b2,-b1/b2,-b3/b2,-b4/b2,-1))
  write.csv(data.frame("b0"=-b0/b2, "b1"=-b1/b2, "b3"=-b3/b2, "b4"=-b4/b2, "b2"=-1, "q"=q, "Lor_cubic"=Lor_cubic, "Lor_se"=Lor_se, "CI"=paste0("(",toString(round(Lor_cubic_CI,3)),")")), paste0(path,alp,"_cubic_para.csv"), row.names = F)  # for plotting
} else {
  print(c(b0/b2,b1/b2,b3/b2,b4/b2,1))
  write.csv(data.frame("b0"=b0/b2, "b1"=b1/b2, "b3"=b3/b2, "b4"=b4/b2, "b2"=1, "q"=q, "Lor_cubic"=Lor_cubic, "Lor_se"=Lor_se, "CI"=paste0("(",toString(round(Lor_cubic_CI,3)),")")), paste0(path,alp,"_cubic_para.csv"), row.names = F)
}

