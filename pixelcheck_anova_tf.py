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

# Load data
#fn = os.path.join(path_in_cuelocked, "sensory_alpha.csv")
fn = os.path.join(path_in_cuelocked, "motor_beta.csv")
df = pd.read_csv(fn, header=None)

# rename cols and factor levels
df.columns = ['id', 'agency', 'reward', 'value']
df['agency'] = df['agency'].replace({1: 'neu', 2: 'slf', 3: 'oth'})
df['reward'] = df['reward'].replace({0: 'low', 1: 'high'})
df['id'] = df['id'].astype('category')

# Plot measures
sns.relplot(
    data=df,
    x='reward',
    y='value',
    hue="agency",
    kind='line',
)
plt.show()

# Anova
aov = pg.rm_anova(
    dv="value",
    within=["agency", "reward"],
    subject="id",
    data=df,
    detailed=True,
)