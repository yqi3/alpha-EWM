######################################
# Author: Alice Yuan Qi
# Figures 1 & 2 in Fan, Qi & Xu (2025)
# Optimal treatment regions, 5 alpha values × 2 ITR classes = 10 plots
######################################

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)

path <- "Figs_2&3_plots/"  # for saving plots
path_ITR <- "Tables_1&7_ITR_para/"  # for reading policy parameters
plot_ub <- 20000  # for better visualization
for (alp in c(0.25, 0.3, 0.4, 0.5, 0.8)) {
  linear_para <- read.csv(paste0(path_ITR, alp,"_linear_para.csv"))
  cubic_para <- read.csv(paste0(path_ITR, alp,"_cubic_para.csv"))
  
  ### Linear rule
  b0 <- linear_para$b0
  b1 <- linear_para$b1
  b2 <- linear_para$b2
  slope <- -b1/b2
  intercept <- -b0/b2
  
  plot_data <- read.csv("plot_data.csv")  # helper dataframe for plotting
  jpeg(file=paste0(path, alp, "linear.jpeg"), width=6, height=4, units="in", res=200)
  if (b2 < 0){
    # Line is upper bound
    plot_data$ub <- (slope*plot_data$edu+intercept>=0)*(slope*plot_data$edu+intercept)
    plot_data$ub <- (plot_data$ub>=plot_ub)*plot_ub + (plot_data$ub<plot_ub)*plot_data$ub
    plot_data <- rbind(plot_data, c((plot_ub-intercept)/slope, 0, 0, plot_ub))  # add helper point for plotting trt region
    plot_data <- rbind(plot_data, c((0-intercept)/slope, 0, 0, 0))  # add a helper point for plotting trt region
    print(ggplot(plot_data, aes(edu, prevearn)) +
            geom_point(aes(size=n)) +
            scale_size_continuous(range = c(-1, 8)) +
            geom_ribbon(data=plot_data, aes(ymin = 0, ymax = ub), fill = "yellow", alpha = .2) +
            ylim(0,plot_ub) +
            xlim(7,18) +
            theme_bw() +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            labs(title=paste0("Linear class, \u03B1=", alp)))
  } else {
    # Line is lower bound
    plot_data$lb <- (slope*plot_data$edu+intercept>=0)*(slope*plot_data$edu+intercept)
    plot_data <- rbind(plot_data, c((plot_ub-intercept)/slope, 0, 0, plot_ub))  # add helper point for plotting trt region
    plot_data <- rbind(plot_data, c((0-intercept)/slope, 0, 0, 0))  # add a helper point for plotting trt region
    
    print(ggplot(plot_data, aes(edu, prevearn)) +
            geom_point(aes(size=n)) +
            scale_size_continuous(range = c(-1, 8)) +
            geom_ribbon(data=plot_data, aes(ymin = lb, ymax = plot_ub), fill = "yellow", alpha = .2) +
            ylim(0,plot_ub) +
            xlim(7,18) +
            theme_bw() +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            labs(title=paste0("Linear class, \u03B1=", alp)))
  }
  dev.off()
  
  ### Cubic rule
  b0 <- cubic_para$b0
  b1 <- cubic_para$b1
  b2 <- cubic_para$b2
  b3 <- cubic_para$b3
  b4 <- cubic_para$b4
  
  # Make a helper grid for plotting
  edu_grid <- seq(7,18,0.005)
  helper <- data.frame("edu"=edu_grid, "prevearn"=rep(0,length(edu_grid)), "n"=rep(0,length(edu_grid)))
  plot_data <- read.csv("plot_data.csv")
  plot_data <- rbind(plot_data, helper)
  intercept <- -b0/b2
  coef_1 <- -b1/b2
  coef_2 <- -b3/b2
  coef_3 <- -b4/b2
  jpeg(file=paste0(path, alp, "cubic.jpeg"), width=6, height=4, units="in", res=200)
  
  if (b2 < 0) {
    # Curve is upper bound
    plot_data$ub <- (coef_1*plot_data$edu+coef_2*(plot_data$edu)^2+coef_3*(plot_data$edu)^3+intercept>=0)*(coef_1*plot_data$edu+coef_2*(plot_data$edu)^2+coef_3*(plot_data$edu)^3+intercept)
    plot_data$ub <- (plot_data$ub>=plot_ub)*plot_ub + (plot_data$ub<plot_ub)*plot_data$ub
    
    print(ggplot(plot_data, aes(edu, prevearn)) +
            geom_point(aes(size=n)) +
            scale_size_continuous(range = c(-1, 8)) +
            geom_ribbon(data=plot_data, aes(ymin = 0, ymax = ub), fill = "yellow", alpha = .2) +
            ylim(0, plot_ub) +
            xlim(7,18) +
            theme_bw() +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            labs(title=paste0("Linear class with squared and cubic edu, \u03B1=", alp)))
  } else {
    # Curve is lower bound
    plot_data$lb <- (coef_1*plot_data$edu+coef_2*(plot_data$edu)^2+coef_3*(plot_data$edu)^3+intercept>=0)*(coef_1*plot_data$edu+coef_2*(plot_data$edu)^2+coef_3*(plot_data$edu)^3+intercept)
    
    print(ggplot(plot_data, aes(edu, prevearn)) +
            geom_point(aes(size=n)) +
            scale_size_continuous(range = c(-1, 8)) +
            geom_ribbon(data=plot_data, aes(ymin = lb, ymax = plot_ub), fill = "yellow", alpha = .2) +
            ylim(0, plot_ub) +
            xlim(7,18) +
            theme_bw() +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            labs(title=paste0("Linear class with squared and cubic edu, \u03B1=", alp)))
  }
  dev.off()
}