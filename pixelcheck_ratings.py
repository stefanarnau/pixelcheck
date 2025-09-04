#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 16:27:06 2025

@author: plkn
"""

# Imports
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pingouin as pg

# Paths
path_in = "/mnt/data_dump/pixelcheck/pixelcheck_ratings/"
path_out = "/mnt/data_dump/pixelcheck/"

# Subjects
subject_list = [
    "001",
    "003",
    "005",
    "006",
    "007",
    "009",
    "010",
    "011",
    "012",
    "013",
    "014",
    "015",
    "016",
    "017",
    "018",
    "019",
    "020",
    "021",
    "022",
    "023",
    "024",
    "026",
    "027",
    "028",
    "029",
    "030",
    "031",
    "032",
    "033",
    "034",
    "035",
    "036",
    "037",
    "038",
    "041",
]

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

# Remove block 0
df = df[~df["block_nr"].isin([0])]

# Remove ratings from first 40 trials
#df = df[df["position_in_block"] > 40]

# Rename values
df["ma_cond"] = df["ma_cond"].replace({1: "neutral", 2: "self", 3: "other"})

#
df = (
    df.groupby(["id", "ma_cond"])[["focus", "self", "computer", "performance"]]
    .mean()
    .reset_index()
)

df['SELF'] = df['self']

# Repeated-measures ANOVAs for individual ratings
aov_focus = pg.rm_anova(
    dv="focus",
    within=["ma_cond"],
    subject="id",
    data=df,
    detailed=True,
)
post_focus = pg.pairwise_tests(dv="focus",
                                 within='ma_cond',
                                 subject='id',
                                 data=df,
                                 padjust='bh',
                                 effsize='cohen')

aov_self = pg.rm_anova(
    dv="SELF",
    within=["ma_cond"],
    subject="id",
    data=df,
    detailed=True,
)
post_self = pg.pairwise_tests(dv="SELF",
                                 within='ma_cond',
                                 subject='id',
                                 data=df,
                                 padjust='bh',
                                 effsize='cohen')

aov_computer = pg.rm_anova(
    dv="computer",
    within=["ma_cond"],
    subject="id",
    data=df,
    detailed=True,
)
post_computer = pg.pairwise_tests(dv="computer",
                                 within='ma_cond',
                                 subject='id',
                                 data=df,
                                 padjust='bh',
                                 effsize='cohen')

aov_performance = pg.rm_anova(
    dv="performance",
    within=["ma_cond"],
    subject="id",
    data=df,
    detailed=True,
)
post_performance = pg.pairwise_tests(dv="performance",
                                 within='ma_cond',
                                 subject='id',
                                 data=df,
                                 padjust='bh',
                                 effsize='cohen')

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
sns.lineplot(
    data=df_long, x="ma_cond", y="rating", hue="satisfaction", estimator="mean", ax=ax
)

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

# Long format with satisfaction self versus other
df_long_all = pd.melt(
    df,
    id_vars=["id", "ma_cond"],
    value_vars=["focus", "self", "computer", "performance"],
    var_name="measure",
    value_name="rating",
)


# Plot measures
sns.relplot(
    data=df_long_all,
    x='ma_cond', y='rating',
    col='measure',
    kind='line', col_wrap=3, facet_kws={'sharey': False}
)
plt.show()