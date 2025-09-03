#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  8 11:23:19 2025

@author: plkn
"""

# Imports
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pingouin as pg
import scipy.stats
import numpy as np
import scipy.special


def ez_diffusion(MRT, VRT, Pc, s = 0.1):
    """
    Translates the R function get.vaTer to Python.
    
    Parameters:
    - Pc: Proportion correct (0 < Pc < 1)
    - VRT: Variance of RT (typically correct trials only)
    - MRT: Mean RT (in seconds, correct trials only)
    - s: Scaling parameter for DDM noise (default: 0.1)

    Returns:
    - v: Drift rate
    - a: Boundary separation
    - Ter: Non-decision time
    """
    s2 = s ** 2

    # Basic validation
    if Pc <= 0 or Pc == 0.5 or Pc >= 1:
        raise ValueError("Pc must be between 0 and 1 and not equal to 0.5 for this method to work.")

    # Logit of Pc
    L = scipy.special.logit(Pc)  # same as qlogis in R

    # Compute intermediate term
    x = L * (L * Pc**2 - L * Pc + Pc - 0.5) / VRT

    # Drift rate
    v = np.sign(Pc - 0.5) * s * (x ** 0.25)

    # Boundary separation
    a = s2 * L / v

    # y = -v * a / s²
    y = -v * a / s2

    # Mean decision time
    MDT = (a / (2 * v)) * (1 - np.exp(y)) / (1 + np.exp(y))

    # Non-decision time
    Ter = MRT - MDT

    return a, v, Ter

# Paths
path_in = "/mnt/data_dump/pixelcheck/1_behavioral_data/"
path_out = "/mnt/data_dump/pixelcheck/"

# Load data
fn = os.path.join(path_in, "behavioral_data.csv")
df = pd.read_csv(fn)

# Check for outlier rts
df['rt_flag'] = np.where(df['rt'] > 1200, 1, 0)

# Summarize for subjects
df_outlier = df.groupby(["id"])[["rt_flag"]].sum().reset_index()

# Identify subjects to remove
outlier_id = df_outlier.id[df_outlier.rt_flag > 10].values

# Remove outlier 
df = df[~df['id'].isin(outlier_id)]

# Exclude
df = df[(df.trial_in_block > 0)]

# Transform to seconds
df.rt = df.rt / 1000

# Rename cols
df.rename(columns={'ma_condition': 'agency'}, inplace=True)
df.rename(columns={'reward_condition': 'reward'}, inplace=True)

# Make categorial
df['agency'] = pd.Categorical(
    df['agency'].replace({1: 'neutral', 2: 'self', 3: 'other'}),
    categories=['neutral', 'self', 'other'],
    ordered=True
)
df['reward'] = pd.Categorical(
    df['reward'].replace({0: 'low', 1: 'high'}),
    categories=['low', 'high'],
    ordered=True
)

# Get correct only
df_rt = df[df.accuracy == 1]

# Make rt mean df
df_rt_mean = df_rt.groupby(["id", "agency", "reward"])[["rt"]].mean().reset_index()
df_rt_mean.rename(columns={'rt': 'rt_mean'}, inplace=True)

# Make rt variance df
df_rt_var = df_rt.groupby(["id", "agency", "reward"])[["rt"]].var().reset_index()
df_rt_var.rename(columns={'rt': 'rt_var'}, inplace=True)

# Binarize accuracy
df_acc = df
df_acc['accuracy'] = df_acc['accuracy'].apply(lambda x: 1 if x == 1 else 0)

# Make accuracy df
df_acc = df_acc.groupby(["id", "agency", "reward"])[["accuracy"]].mean().reset_index()

# Make common dataframe
df_grouped = df_rt_mean
df_grouped["rt_var"] = df_rt_var.rt_var
df_grouped["accuracy"] = df_acc.accuracy

# Calculate EZ-diffusion parameters
df_grouped["drift_rate"] = 0
df_grouped["boundary_seperation"] = 0
df_grouped["non_decision_time"] = 0

for idx, row in df_grouped.iterrows():
    
    a, v, t0 = ez_diffusion(row["rt_mean"], row["rt_var"], row["accuracy"])
    df_grouped.at[idx, "drift_rate"] = v
    df_grouped.at[idx, "boundary_seperation"] = a
    df_grouped.at[idx, "non_decision_time"] = t0
    

df_long = df_grouped.melt(id_vars=["id", "agency", "reward"], var_name='variable', value_name='value')

# Plot measures
sns.relplot(
    data=df_long,
    x='reward', y='value', hue="agency",
    col='variable',
    kind='line', col_wrap=3, facet_kws={'sharey': False}
)
plt.show()



# Define your dependent variables
dvs = ['rt_mean', 'accuracy', 'drift_rate', 'boundary_seperation', 'non_decision_time']

anova_results_all = []
posthoc_results_all = []

for dv in dvs:
    
    # rm ANOVA
    anova_res = pg.rm_anova(dv=dv,
                            within=['agency', 'reward'],
                            subject='id',
                            data=df_grouped,
                            detailed=True)
    anova_res['DV'] = dv
    anova_results_all.append(anova_res)

    # Check for significant main effect of agency
    sig_agency = (anova_res.loc[anova_res['Source'] == 'agency', 'p-GG-corr'] < 0.05).any()
    if sig_agency:
        post_agency = pg.pairwise_ttests(dv=dv,
                                         within='agency',
                                         subject='id',
                                         data=df_grouped,
                                         padjust='fdr_bh',
                                         effsize='cohen')
        post_agency['DV'] = dv
        post_agency['Effect'] = 'agency'
        posthoc_results_all.append(post_agency)

    # Check for significant interaction
    sig_interaction = (anova_res.loc[anova_res['Source'] == 'agency * reward', 'p-GG-corr'] < 0.05).any()
    if sig_interaction:
        
        # Simple effects: agency within each reward
        for rw in df_grouped['reward'].unique():
            simple_agency = pg.pairwise_ttests(dv=dv,
                                               within='agency',
                                               subject='id',
                                               data=df_grouped[df_grouped['reward'] == rw],
                                               padjust='fdr_bh',
                                               effsize='cohen')
            simple_agency['DV'] = dv
            simple_agency['Effect'] = f'agency_within_reward_{rw}'
            posthoc_results_all.append(simple_agency)

        # Simple effects: reward within each agency
        for ag in df_grouped['agency'].unique():
            simple_reward = pg.pairwise_ttests(dv=dv,
                                               within='reward',
                                               subject='id',
                                               data=df_grouped[df_grouped['agency'] == ag],
                                               padjust='fdr_bh',
                                               effsize='cohen')
            simple_reward['DV'] = dv
            simple_reward['Effect'] = f'reward_within_agency_{ag}'
            posthoc_results_all.append(simple_reward)

# Combine into dataframes
anova_results_df = pd.concat(anova_results_all, ignore_index=True)
posthoc_results_df = pd.concat(posthoc_results_all, ignore_index=True) if posthoc_results_all else pd.DataFrame()






