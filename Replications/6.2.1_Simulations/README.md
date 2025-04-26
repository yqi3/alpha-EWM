# Replication files for Fan, Qi & Xu (2025), Section 6.2.1
## Two Simulation Studies Based on Data Generating Processes in Athey & Wager (2021)
Script and subfolders:
- DGPs.R: R script for generating superpopulation data according to the DGPs in Section 6.2.1. Data files are saved as "data_tau6.3.csv" and "data_tau6.4.csv" in the corresponding subfolders. These data files are zipped on Github and need to be unzipped for use.
- Within subfolders "tau6.3" and "tau6.4"
  - population_truths.R: R script for computing the population true welfare under six specifications of alpha (alpha=1 is not used in the main text). The truths are already incorporated in the simulation scripts, so there is no need to rerun this file before running the simulations. This script also produces Tables 10 and 11 in Appendix I.2 (truth_comparisons.csv in each subfolder).
  - policy_comparisons.R: R script that performs population-level comparisons between the 0.25-EWM, equality-minded, and 1-EWM policies and replicates Figures 5 and 6 (q_diff_all.jpeg in each subfolder).
  - n_300.R, n_500.R, n_1000.R, n_1500.R: R scripts for Monte Carlo simulations (of four sample sizes) that replicate Tables 4 and 5 in paper.
  - 6.2.1_Helper_Fns.R: Helper functions called by other scripts in this folder.

Note: We implement random forests using the [`grf`](https://grf-labs.github.io/grf/reference/index.html) package, whose reproducibility depends on both seed and computing platform. Exact reproducibility is therefore not guaranteed. 