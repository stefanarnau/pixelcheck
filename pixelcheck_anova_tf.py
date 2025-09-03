#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  2 11:56:19 2025

@author: plkn
"""
# Imports
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pingouin as pg

# Paths
path_in_cuelocked = "/mnt/data_dump/pixelcheck/4_tf_results_cuelocked/"
path_in_resplocked = "/mnt/data_dump/pixelcheck/4_tf_results_resplocked/"

# Error monitoring | frontal theta resplocked =====================================

# Load data
fn = os.path.join(path_in_resplocked, "error_monitoring_frontal_theta.csv")
df = pd.read_csv(fn, header=None)

# rename cols and factor levels
df.columns = ['id', 'agency', 'reward', 'value']
df['agency'] = df['agency'].replace({1: 'sff', 2: 'oth'})
df['reward'] = df['reward'].replace({1: 'low', 2: 'high'})
df['id'] = df['id'].astype('category')

# Plot measures
plt_error_monitoring = sns.relplot(
    data=df,
    x='reward',
    y='value',
    hue="agency",
    kind='line',
)
plt_error_monitoring.fig.subplots_adjust(top=0.9)
plt_error_monitoring.fig.suptitle('frontal theta error monitoring')
plt.show()

# Anova
aov_error_monitoring = pg.rm_anova(
    dv="value",
    within=["agency", "reward"],
    subject="id",
    data=df,
    detailed=True,
)

# Flip trials | frontal theta resplocked =====================================

# Load data
fn = os.path.join(path_in_resplocked, "fliptrials_frontal_theta.csv")
df = pd.read_csv(fn, header=None)

# rename cols and factor levels
df.columns = ['id', 'agency', 'reward', 'value']
df['agency'] = df['agency'].replace({1: 'sff', 2: 'oth'})
df['reward'] = df['reward'].replace({1: 'low', 2: 'high'})
df['id'] = df['id'].astype('category')

# Plot measures
plt_flip_trials = sns.relplot(
    data=df,
    x='reward',
    y='value',
    hue="agency",
    kind='line',
)
plt_flip_trials.fig.subplots_adjust(top=0.9)
plt_flip_trials.fig.suptitle('frontal theta flip trials')
plt.show()

# Anova
aov_flip_trials = pg.rm_anova(
    dv="value",
    within=["agency", "reward"],
    subject="id",
    data=df,
    detailed=True,
)

# Motor preparation | centro-lateral beta cuelocked =====================================

# Load data
fn = os.path.join(path_in_cuelocked, "motorprep_beta.csv")
df = pd.read_csv(fn, header=None)

# rename cols and factor levels
df.columns = ['id', 'agency', 'reward', 'value']
df['agency'] = df['agency'].replace({1: 'neu', 2: 'slf', 3: 'oth'})
df['reward'] = df['reward'].replace({0: 'low', 1: 'high'})
df['id'] = df['id'].astype('category')

# Plot measures
plt_motorprep = sns.relplot(
    data=df,
    x='reward',
    y='value',
    hue="agency",
    kind='line',
)
plt_motorprep.fig.subplots_adjust(top=0.9)
plt_motorprep.fig.suptitle('beta motor preparation')
plt.show()

# Anova
aov_motorprep = pg.rm_anova(
    dv="value",
    within=["agency", "reward"],
    subject="id",
    data=df,
    detailed=True,
)

posthoc_agency_motorprep = pg.pairwise_ttests(dv="value",
                                   within='agency',
                                   subject='id',
                                   data=df,
                                   padjust='bh',
                                   effsize='cohen')


# Attention preparation | Occipital alpha cuelocked =====================================

# Load data
fn = os.path.join(path_in_cuelocked, "attentionprep_alpha.csv")
df = pd.read_csv(fn, header=None)

# rename cols and factor levels
df.columns = ['id', 'agency', 'reward', 'value']
df['agency'] = df['agency'].replace({1: 'neu', 2: 'slf', 3: 'oth'})
df['reward'] = df['reward'].replace({0: 'low', 1: 'high'})
df['id'] = df['id'].astype('category')

# Plot measures
plt_attentionprep = sns.relplot(
    data=df,
    x='reward',
    y='value',
    hue="agency",
    kind='line',
)
plt_attentionprep.fig.subplots_adjust(top=0.9)
plt_attentionprep.fig.suptitle('alpha attention preparation')
plt.show()

# Anova
aov_attentionprep = pg.rm_anova(
    dv="value",
    within=["agency", "reward"],
    subject="id",
    data=df,
    detailed=True,
)

posthoc_agency_attentionprep = pg.pairwise_ttests(dv="value",
                                   within='agency',
                                   subject='id',
                                   data=df,
                                   padjust='bh',
                                   effsize='cohen')





