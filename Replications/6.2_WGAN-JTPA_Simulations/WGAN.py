'''
Author: Alice Yuan Qi
This file uses WGAN to generate simulation data that looks like
the data from the JTPA data used by Kitagawa & Tetenov (2018).
The pre-treatment covariates are years of education (edu) and
pre-program earnings (prevearn) and the outcome variable is an
applicant's earnings 30 months after the assignment (earnins).
See Bloom et al. (1997) for details on the training samples
and Athey et al. (2021) for WGAN.
'''

import wgan
import pandas as pd
from copy import copy
import matplotlib.pyplot as plt

# Function for saving the correlation plots after fake data is generated
def save_corr(df_real, df_fake, figsize=5, path=""):
    if "source" in list(df_real.columns): df_real = df_real.drop("source", axis=1)
    if "source" in list(df_fake.columns): df_fake = df_fake.drop("source", axis=1)
    df_real.insert(0, "source", "real"), df_fake.insert(0, "source", "fake")
    common_cols = [c for c in df_real.columns if c in df_fake.columns]
    df_joined = pd.concat([df_real[common_cols], df_fake[common_cols]], axis=0, ignore_index=True)
    df_real, df_fake = df_real.drop("source", axis=1), df_fake.drop("source", axis=1)
    common_cols = [c for c in df_real.columns if c in df_fake.columns]
    fig1 = plt.figure(figsize=(figsize * 2, figsize * 1))
    s1 = [fig1.add_subplot(1, 2, i) for i in range(1, 3)]
    s1[0].set_xlabel("real")
    s1[1].set_xlabel("fake")
    s1[0].matshow(df_real[common_cols].corr())
    s1[1].matshow(df_fake[common_cols].corr())
    fig1.savefig(path+'_corr.png')

# Read sample data
df = pd.read_csv('KT_JTPA.csv')

# First train X | A, then train Y | (X, A)
continuous_vars_0 = ["prevearn"]
continuous_lower_bounds_0 = {"prevearn": 0}
categorical_vars_0 = ["edu"]
context_vars_0 = ["D"]  # "D" is equivalent to "A" in paper

continuous_vars_1 = ["earnings"]
continuous_lower_bounds_1 = {"earnings": 0}
categorical_vars_1 = []
context_vars_1 = ["D", "edu", "prevearn"]

data_wrappers = [wgan.DataWrapper(df, continuous_vars_0, categorical_vars_0, 
                                  context_vars_0, continuous_lower_bounds_0),
                 wgan.DataWrapper(df, continuous_vars_1, categorical_vars_1, 
                                  context_vars_1, continuous_lower_bounds_1)]

specs = [wgan.Specifications(dw, batch_size=4096, max_epochs=1000, critic_lr=1e-3, generator_lr=1e-3,
                             print_every=100) for dw in data_wrappers]

generators = [wgan.Generator(spec) for spec in specs]
critics = [wgan.Critic(spec) for spec in specs]

x, context = data_wrappers[0].preprocess(df)
wgan.train(generators[0], critics[0], x, context, specs[0])

x, context = data_wrappers[1].preprocess(df)
wgan.train(generators[1], critics[1], x, context, specs[1])

# Simulated superpopulation size = 1e6
df_generated = data_wrappers[0].apply_generator(generators[0], df.sample(int(1e6), replace=True))
df_generated = data_wrappers[1].apply_generator(generators[1], df_generated)  # contains generated X & Y

# Apply the generator for Y | (X, A) on Y | (X, 1-A) to obtain the population counterfactuals
df_generated_cf = copy(df_generated)
df_generated_cf["D"] = 1 - df_generated_cf["D"]
df_generated["earnings_cf"] = data_wrappers[1].apply_generator(generators[1], df_generated_cf)["earnings"]  # add counterfactual Ys
df_generated.to_csv("WGAN-JTPA.csv", index=False)

# Compare WGAN-JTPA to original JTPA
wgan.compare_dfs(df, df_generated, 
                 scatterplot=dict(x=["edu", "prevearn"], y=["earnings"], samples=500, smooth=0),
                 table_groupby=["D"],
                 histogram=dict(variables=["earnings", "earnings", "edu", "prevearn"], nrow=2, ncol=2),
                 figsize=5, save=True, path="comparison_plots/")
save_corr(df, df_generated, figsize=5, path="comparison_plots/")
plt.show(block=True)
