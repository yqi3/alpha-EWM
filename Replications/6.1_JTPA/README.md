# Replication files for Fan, Qi, and Xu (2025), Section 6.1
## Empirical Application: The JTPA Study
Data files: 
- KT_JTPA.csv: JTPA data organized by Kitagawa & Tenenov (2018). The original data source is Bloom et al. (1997).
- plot_data.csv: Helper data frame that bins observations in KT_JTPA.csv by a grid of (edu, prevearn) values. Used for plotting treatment regions (see plot.R). 

Scripts and results folders:
- alpha_0.25.R, alpha_0.3.R, alpha_0.4.R, alpha_0.5.R, alpha_0.8.R: R scripts for learning optimal policies targeting the worst-affected alpha = 0.25, 0.3, 0.4, 0.5, and 0.8, respectively. Results are saved in folder "Tables_1&5_ITR_para." Figures of training losses are saved in folder "training_losses."
- plot.R: R script for plotting the corresponding treatment regions. Figures are saved in folder "Figs_2&3_plots."
- comparisons_with_se.R: R script for comparing the estimated welfare for different actual alpha's of interest and alpha's for policy selection. Results are saved in folder "Tables_2&3_comparisons."
- JTPA_Helper_Fns.R: Helper functions called by other scripts in this folder. 

Note: We implement random forests using the [`grf`](https://grf-labs.github.io/grf/reference/index.html) package, whose reproducibility depends on both seed and computing platform. Exact reproducibility is therefore not guaranteed. 