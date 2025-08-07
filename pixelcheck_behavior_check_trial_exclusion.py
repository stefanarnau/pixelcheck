#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 16:27:06 2025

@author: plkn
"""

# Imports
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Paths
path_in = "/mnt/data_dump/pixelcheck/1_behavioral_data/"
path_out = "/mnt/data_dump/pixelcheck/"

# Load data
fn = os.path.join(path_in, "behavioral_data.csv")
df = pd.read_csv(fn)

rows = []

# Iterate subjects
for subject in set(df.id.values):
    
    # Filter
    df_subject = df[df.id == subject]
        
    # Iterate trial exclusion
    for trial_exclusion in range(90):
        
        # Filter
        df_start_excluded = df_subject[df_subject.trial_in_block > trial_exclusion]
    
        # Iterate agency and reward
        for agency in range(1, 4):
            for reward in range(2):
                
                # Filter
                tmp = df_start_excluded[(df_start_excluded.ma_condition == agency) & (df_start_excluded.reward_condition == reward)]
                
                # iterate block exclusion
                for block_exclusion in range(4):
                    
                    #Filter
                    tmp2 = tmp[tmp["block_nr"] > block_exclusion]
                    
                    # iterate exclude post flips
                    for post_flip_exclusion in range(2):
                        
                        if post_flip_exclusion == 1:
                            tmp3 = tmp2[tmp2.last_fb_flipped == 0]
                        else:
                            tmp3 = tmp2
       
                
                        # Init row dict
                        row = {}
                        row["id"] = subject
                        row["post_flip_exclusion"] = post_flip_exclusion
                        row["block_exclusion"] = block_exclusion
                        row["trial_exclusion"] = trial_exclusion
                        row["agency"] = ["neu", "slf", "oth"][agency - 1]
                        row["reward"] = ["lo", "hi"][reward]
                        row["rt"] = tmp2[tmp2.accuracy == 1].rt.values.mean()
                        row["accuracy"] = tmp2[tmp2.accuracy == 1].shape[0] / tmp2.shape[0]
                        rows.append(row)
                    
df_means = pd.DataFrame(rows)


sns.set(style="whitegrid")
g = sns.relplot(
    data=df_means,
    x="trial_exclusion",         
    y="rt",     
    hue="agency",  
    style="reward",
    kind="line",      
    col="block_exclusion",       
    row="post_flip_exclusion",    
    markers=False,
    ci=None
)
g.fig.subplots_adjust(top=0.9)
g.fig.suptitle("RT dependency on trial exclusion criteria")

sns.set(style="whitegrid")
g = sns.relplot(
    data=df_means,
    x="trial_exclusion",         
    y="accuracy",     
    hue="agency",  
    style="reward",
    kind="line",      
    col="block_exclusion",       
    row="post_flip_exclusion",    
    markers=False,
    ci=None
)
g.fig.subplots_adjust(top=0.9)
g.fig.suptitle("Accuracy dependency on trial exclusion criteria")