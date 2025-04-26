######################################
# Author: Alice Yuan Qi
# Data generating process for Fan, Qi & Xu (2025), Section 1
# Specifications are adopted from Wang et al. (2018), with minor adjustments
######################################

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
set.seed(1234)

#### Prep data ####
n <- 1e6
x <- runif(n)
A <- rbinom(n, size = 1, prob = 0.5)
eps <- rnorm(n)
y <- 20+3*A+x-5*A*x+(1+A+2*A*x)*eps
y_cf <- 20+3*(1-A)+x-5*(1-A)*x+(1+(1-A)+2*(1-A)*x)*eps
superpop <- data.frame(cbind(x,A,y,y_cf))
write.csv(superpop, "data.csv", row.names = F)
