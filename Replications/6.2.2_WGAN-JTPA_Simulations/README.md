# Replication files for Fan, Qi & Xu (2025), Section 6.2.2
## Simulated JTPA Data using Wasserstein Generative Adversarial Networks (WGANs)
The [wgan package repository](https://github.com/gsbDBI/ds-wgan) contains installation instructions and usage examples.

Data files: 
- KT_JTPA.csv: JTPA data organized by Kitagawa & Tenenov (2018). The original data source is Bloom et al. (1997).
- WGAN-JTPA.csv: Simulated JTPA-based superpopulation data using WGAN (see WGAN.py).

Scripts:
- WGAN.py: Python script for generating superpopulation data based on the original JTPA data. The generated data is saved as "WGAN-JTPA.csv," and summary statistics and comparative plots are saved in subfolder "comparison_plots."
- population_truths.R: R script for computing the population true welfares under six specifications of alpha (alpha=1 is not used in the main text). The truths are already incorporated in the simulation scripts, so there is no need to rerun this file before running the simulations. This script also produces Tables 13 in Appendix I.2 (truth_comparisons.csv).
- policy_comparisons.R: R script that performs population-level comparisons between the 0.25-EWM, equality-minded, and 1-EWM policies and replicates Figure 4 (q_diff_all.jpeg).
- n_2000.R, n_5000.R, n_10000.R: R scripts for Monte Carlo simulations (of three sample sizes) that replicate Table 6 in paper.
- JTPA_Helper_Fns.R: Helper functions called by other scripts in this folder. 

Note: We implement random forests using the [`grf`](https://grf-labs.github.io/grf/reference/index.html) package, whose reproducibility depends on both seed and computing platform. Exact reproducibility is therefore not guaranteed. 