# Replication files for Fan, Qi & Xu (2024), Section 5.2.1
## Two Simulation Studies Based on Data Generating Processes in Athey & Wager (2021)
Script and subfolders:
- DGPs.R: R script for generating superpopulation data according to the DGPs in Section 5.2.1. Data files are saved as "data_tau5.3.csv" and "data_tau5.4.csv" in the corresponding subfolders. These data files are zipped on Github and needs to be unzipped for use.
- Within subfolders "tau5.3" and "tau5.4"
  - population_truths.R: R script for computing the population true welfares under six specifications of alpha. The truths are already incorporated in the simulation scripts.
  - n_300.R, n_500.R, n_1000.R, n_1500_part1.R, n_1500_part2.R: R scripts for Monte Carlo simulations (of four sample sizes) that replicate Tables 4 and 5 in paper.
  - 5.2.1_Helper_Fns.R: Helper functions called by other scripts in this folder.

Note: We implement random forests using the [`grf`](https://grf-labs.github.io/grf/reference/index.html) package, whose replicability depends on both seed and computing platform. See the header of each R script.