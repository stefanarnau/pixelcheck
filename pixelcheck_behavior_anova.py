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
df = df[(df.trial_in_block > 0) & (df.block_nr > 3)]

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


# ANOVA
anova_results_rt = pg.rm_anova(dv='rt_mean', within=['agency', 'reward'], subject='id', data=df_grouped, detailed=True)
anova_results_acc = pg.rm_anova(dv='accuracy', within=['agency', 'reward'], subject='id', data=df_grouped, detailed=True)
anova_results_drift_rate = pg.rm_anova(dv='drift_rate', within=['agency', 'reward'], subject='id', data=df_grouped, detailed=True)
anova_results_boundary_seperation = pg.rm_anova(dv='boundary_seperation', within=['agency', 'reward'], subject='id', data=df_grouped, detailed=True)
anova_results_non_decision_time = pg.rm_anova(dv='non_decision_time', within=['agency', 'reward'], subject='id', data=df_grouped, detailed=True)






