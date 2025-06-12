#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 16:27:06 2025

@author: plkn
"""

# Imports
import os
import re
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pingouin as pg

# Paths
path_in = "/mnt/data_dump/pixelcheck/pixelcheck_ratings/"
path_out = "/mnt/data_dump/pixelcheck/"

# Subjects
subject_list = list(
    set(
        [f[:3] for f in os.listdir(path_in) if os.path.isfile(os.path.join(path_in, f))]
    )
)

# Loop subjects
df = []
for subject in subject_list:

    # Get subject id number as list
    subject_id = [float(subject[2:])]

    # Read trial log
    fn = os.path.join(path_in, subject + "_pixelCheckLog.txt")
    with open(fn) as file:
        lines_trial_log = file.readlines()

    # Get trial list as array. Cols are trial number, block number, ma condition
    trial_list = [line.split(" ")[:3] for line in lines_trial_log[6:-6]]
    trial_list = np.array([[float(num) for num in nums] for nums in trial_list])

    # Read ratings
    fn = os.path.join(path_in, subject + "_thoughtProbeLog.txt")
    with open(fn) as file:
        lines_thought_probes = file.readlines()

    # Get ratings as list
    ratings_list = [line.split(" ")[:-1] for line in lines_thought_probes[4:]]

    # Iterate ratings
    for rating in ratings_list:

        # get trial number
        rating_trial_nr, rating_values = float(rating[0]), [
            float(x) for x in rating[1:]
        ]
        
        # Get position in block
        position_in_block = np.mod(rating_trial_nr, 120)

        # Get block nr
        block_nr = trial_list[trial_list[:, 0] == rating_trial_nr, 1][0]

        # Get ma_cond
        ma_cond = trial_list[trial_list[:, 0] == rating_trial_nr, 2][0]

        # Compile dict
        df.append(
            {
                "id": float(subject),
                "trial_nr": rating_trial_nr,
                "block_nr": block_nr - 1,
                "position_in_block": position_in_block,
                "ma_cond": ma_cond,
                "focus": rating_values[0],
                "self": rating_values[1],
                "computer": rating_values[2],
                "performance": rating_values[3],
            }
        )

# Create dataframe
df = pd.DataFrame(df)

# Remove block 0 and 1 (practice and fist neutral)
df = df[~df["block_nr"].isin([0, 1])]

# Remove ratings from first 40 trials
df = df[df["position_in_block"] > 40]

# Rename values
df["ma_cond"] = df["ma_cond"].replace({1: "neutral", 2: "self", 3: "other"})

# Average topo df across ids
df = df.groupby(["id", "ma_cond"])["focus", "self", "computer", "performance"].mean().reset_index()


# Long format with satisfaction self versus other
df_long = pd.melt(
    df,
    id_vars=["id", "ma_cond"],
    value_vars=["self", "computer"],
    var_name="satisfaction",
    value_name="rating",
)


# Set the figure size (30x30 cm converted to inches: 1 inch = 2.54 cm)
fig, ax = plt.subplots(figsize=(30 / 2.54, 30 / 2.54))


# Create the lineplot
sns.lineplot(data=df_long, x="ma_cond", y="rating", hue="satisfaction", estimator="mean", ax=ax)

# Set title
ax.set_title("satisfaction with own versus computer performance")

# Optional: Tweak other elements if needed
ax.set_xlabel("feedback condition")
ax.set_ylabel("Rating")

# Show the plot
plt.tight_layout()
plt.show()

# Save the figure with high DPI
fn = os.path.join(path_out, "ratings.png")
fig.savefig(fn, dpi=600) 

# Repeated-measures ANOVA
aov = pg.rm_anova(
    dv="rating",
    within=["ma_cond", "satisfaction"],
    subject="id",
    data=df_long,
    detailed=True,
)

print(aov)
