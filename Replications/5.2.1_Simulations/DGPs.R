######################################
# Author: Alice Yuan Qi
# Data generating processes for Fan, Qi & Xu (2024), Section 5.2.1
# Specifications of Y & tau are adopted from Athey & Wager (2021), Section 5.2, with minor adjustments
######################################

rm(list=ls())
set.seed(1234)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#### tau as (5.3) in paper ####
X <- matrix(rnorm(1e6*4), 1e6, 4)
eps <- rnorm(1e6)
A <- rbernoulli(1e6, p=2/3)
get_tau <- function(x) {
  return((max(x[1],0)+max(x[2],0)-1)/2)
}
tau <- apply(X, 1, get_tau)

d <- cbind(X,tau,A,eps)
get_Y <- function(x_tau_a_e) {
  return(10+max(x_tau_a_e[3]+x_tau_a_e[4],0)+x_tau_a_e[5]*x_tau_a_e[6]+x_tau_a_e[7])
}
Y <- apply(d, 1, get_Y)
d <- cbind(d,Y)
colnames(d)[1:4] <- c("X1", "X2", "X3", "X4")
write.csv(d, "tau5.3/data_tau5.3.csv", row.names = F)

#### tau as (5.4) in paper ####
X <- matrix(rnorm(1e6*4), 1e6, 4)
eps <- rnorm(1e6)
A <- rbernoulli(1e6, p=2/3)
get_tau <- function(x) {
  return(sign(x[1]*x[2])/2)
}
tau <- apply(X, 1, get_tau)

d <- cbind(X,tau,A,eps)
get_Y <- function(x_tau_a_e) {
  return(10+max(x_tau_a_e[3]+x_tau_a_e[4],0)+x_tau_a_e[5]*x_tau_a_e[6]+x_tau_a_e[7])
}
Y <- apply(d, 1, get_Y)
d <- cbind(d,Y)
colnames(d)[1:4] <- c("X1", "X2", "X3", "X4")
write.csv(d, "tau5.4/data_tau5.4.csv", row.names = F)
