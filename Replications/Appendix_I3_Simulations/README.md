# Replication files for Fan, Qi, and Xu (2025), Appendix I.3
## Two Simulation Studies Based on Data Generating Processes in Athey & Wager (2021)
Script and subfolders:
- DGPs.R: R script for generating superpopulation data according to the DGPs in Section I.3. Data files are saved as "data_tau_1.csv" and "data_tau_2.csv" in the corresponding subfolders. These data files are zipped on Github and need to be unzipped for use.
- Within subfolders "tau_1" and "tau_2"
  - population_truths.R: R script for computing the population true welfare under six specifications of alpha (alpha=1 is not used in the main text). The truths are already incorporated in the simulation scripts, so there is no need to rerun this file before running the simulations. This script also produces Tables 10 and 11 in Appendix I.3 (truth_comparisons.csv in each subfolder).
  - policy_comparisons.R: R script that performs population-level comparisons between the 0.25-EWM, equality-minded, and 1-EWM policies and replicates Figures 7 and 8 (q_diff_all.jpeg in each subfolder).
  - n_300.R, n_500.R, n_1000.R, n_1500.R: R scripts for Monte Carlo simulations (of four sample sizes) that replicate Tables 12 and 13 in Appendix I.3.
  - I3_Helper_Fns.R: Helper functions called by other scripts in this folder.

Note: We implement random forests using the [`grf`](https://grf-labs.github.io/grf/reference/index.html) package, whose reproducibility depends on both seed and computing platform. Exact reproducibility is therefore not guaranteed. 