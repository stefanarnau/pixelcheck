#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 10 15:23:19 2025

@author: plkn
"""

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
path_in = "/mnt/data_dump/pixelcheck/4_tf_results/"

# Load behavior
fn = os.path.join(path_in, "behavior.csv")
df = pd.read_csv(fn, header=0)

pivot = df.pivot(index=['id'], columns=['agency', 'reward'], values=['rt_mean', 'accuracy']).reset_index()


# Load data
fn = os.path.join(path_in, "error_monitoring_frontal_theta.csv")
df1 = pd.read_csv(fn, header=None)

# Load data
fn = os.path.join(path_in, "fliptrials_frontal_theta.csv")
df2 = pd.read_csv(fn, header=None)

# Load data
fn = os.path.join(path_in, "motorprep_beta.csv")
df3 = pd.read_csv(fn, header=None)

# Load data
fn = os.path.join(path_in, "attentionprep_alpha.csv")
df = pd.read_csv(fn, header=None)


# rename cols and factor levels
df.columns = ["id", "agency", "reward", "value"]
df["agency"] = df["agency"].replace({1: "neu", 2: "slf", 3: "oth"})
df["reward"] = df["reward"].replace({0: "low", 1: "high"})
df["id"] = df["id"].astype("category")





# Load data
fn = os.path.join(path_in, "ratings.csv")
df6 = pd.read_csv(fn, header=None)