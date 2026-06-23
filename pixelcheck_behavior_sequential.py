# Imports
import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.stats.anova import AnovaRM
import numpy as np

# Paths
path_in = "/home/plkn/repos/pixelcheck/"

# Load data
df = pd.read_csv(os.path.join(path_in, "behavioral_data.csv"))

# Accuracy on the previous trial within the same block
df["last_accuracy"] = (
    df.groupby(["id", "block_nr"])["accuracy"]
      .shift(1)
)

# First trial in each block has no previous-trial accuracy
df.loc[df["trial_in_block"] == 1, "last_accuracy"] = np.nan

# Keep only MA conditions 2 and 3
dat = df.loc[df["ma_condition"].isin([2, 3])].copy()

# Remove errors if desired
dat = dat.loc[dat["accuracy"] == 1].copy()

# Convert to categorical
for col in ["last_fb_flipped", "ma_condition", "reward_condition"]:
    dat[col] = dat[col].astype("category")
    
# Ensure trial order within each block
df = df.sort_values(["id", "block_nr", "trial_in_block"]).copy()


# Filter out post error trials
#dat = dat.loc[dat["last_accuracy"] == 1].copy()


# Subject-level cell means
plotdat = (
    dat.groupby(
        ["id",
         "last_fb_flipped",
         "ma_condition",
         "reward_condition"],
        observed=True,
        as_index=False
    )
    .agg(rt=("rt", "mean"),
         n_trials=("rt", "size"))
)

# Optional: check missing cells
cell_counts = (
    plotdat
    .groupby("id")
    .size()
)

print(cell_counts.value_counts())
print("Subjects with incomplete cells:")
print(cell_counts[cell_counts < 8])

# Keep only subjects with all 2 × 2 × 2 cells
complete_ids = cell_counts[cell_counts == 8].index
plotdat_complete = plotdat[plotdat["id"].isin(complete_ids)].copy()

# Repeated-measures ANOVA on subject averages
aov = AnovaRM(
    data=plotdat_complete,
    depvar="rt",
    subject="id",
    within=["last_fb_flipped", "ma_condition", "reward_condition"]
)

res = aov.fit()
print(res)


g = sns.catplot(
    data=plotdat_complete,
    x="last_fb_flipped",
    y="rt",
    hue="ma_condition",
    col="reward_condition",
    kind="point",
    errorbar=("se", 1),
    height=5,
    aspect=1
)

g.set_axis_labels("Previous feedback flipped", "Mean RT (ms)")
plt.show()