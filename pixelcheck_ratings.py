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

# Rename values
df["ma_cond"] = df["ma_cond"].replace({1: "neutral", 2: "self", 3: "other"})


# Long format with confidence self versus other
df_long = pd.melt(
    df,
    id_vars=["id", "ma_cond"],
    value_vars=["self", "computer"],
    var_name="confidence",
    value_name="rating",
)

sns.lineplot(data=df_long, x="ma_cond", y="rating", hue="confidence", estimator="mean")
plt.title("Lineplot of dv1 and dv2 by ma_cond")
plt.show()


# Repeated-measures ANOVA
aov = pg.rm_anova(
    dv="rating",
    within=["ma_cond", "confidence"],
    subject="id",
    data=df_long,
    detailed=True,
)

print(aov)
