#!/usr/bin/env python
# coding: utf-8

# ## Compare Fat Globule Phenotypes
# 
# **Gregory Way, 2019**

# In[1]:


import os
import pandas as pd
import numpy as np
from scipy import stats

import plotnine as gg


# In[2]:


def compare_groups(df, ffa, day, features):
    # Perform an independent t-test for given combination of ffa and day
    day = str(day)
    ffa = int(ffa)
    comp_df = (
        df
        .query("Metadata_FFA == @ffa")
        .query("Metadata_diff_day == @day")
        .reset_index(drop=True)
    )
    a_df = comp_df.query("Metadata_category == 1").loc[:, cp_features].reset_index(drop=True)
    b_df = comp_df.query("Metadata_category == 0").loc[:, cp_features].reset_index(drop=True)
    
    t_stat, p_val = stats.ttest_ind(a=a_df, b=b_df)
    
    t_result_df = pd.DataFrame(t_stat).transpose()
    t_result_df.columns = cp_features
    t_result_df = t_result_df.melt().rename({"value": "t_stat",
                                             "variable": "cp_feature"}, axis="columns")

    p_result_df = pd.DataFrame(p_val).transpose()
    p_result_df.columns = cp_features
    p_result_df = p_result_df.melt().rename({"value": "p_val",
                                             "variable": "cp_feature"}, axis="columns")


    result_df = (
        t_result_df
        .merge(p_result_df, on="cp_feature")
        .assign(FFA=ffa, diff_day=day)
        .sort_values(by="cp_feature")
        .reset_index(drop=True)
    )

    return result_df


# In[3]:


# Load labels
file = os.path.join("data", "category_labels.csv")
label_df = pd.read_csv(file)
label_df.columns = ["Metadata_{}".format(x) for x in label_df.columns]

label_df


# In[4]:


# Load profiles
file = os.path.join("data", "batch1_batch3_combined_normalized_variable_selected.tsv")
df = pd.read_csv(file, sep='\t').drop(["Metadata_Plate",
                                       "Metadata_Assay_Plate_Barcode",
                                       "Metadata_well_position"], axis="columns")

print(df.shape)
df.head()


# In[5]:


aggregate_cols = ["Metadata_Plate_Map_Name", "Metadata_cell_line", "Metadata_patient",
                  "Metadata_FFA", "Metadata_diff_day", "Metadata_Batch"]

agg_df = df.groupby(aggregate_cols).mean().reset_index()
agg_df.loc[:, "Metadata_patient"] = ["m{}".format(x.strip("PAC_")) for x in agg_df.Metadata_patient]

# Merge labels
group_df = label_df.merge(agg_df, left_on="Metadata_IID", right_on="Metadata_patient", how="inner")

print(group_df.shape)
group_df.head()


# In[6]:


bodipy_features = [x for x in group_df.columns if "BODIPY" in x]

print(len(bodipy_features))
bodipy_features


# In[7]:


group_df.Metadata_patient.value_counts()


# In[8]:


group_df.Metadata_cell_line.value_counts()


# In[9]:


# Split visceral and subcutaneous
vc_df = group_df.query("Metadata_cell_line == 'vc'")
sc_df = group_df.query("Metadata_cell_line == 'sc'")


# In[10]:


vc_df.Metadata_diff_day.value_counts()


# ## Calculate all comparisons

# In[11]:


cp_features = [x for x in vc_df.columns if not x.startswith("Metadata_")]


# In[12]:


all_results = []
for ffa in [0, 1]:
    for day in [0, 3, 14]:
        result_df = compare_groups(df=vc_df, ffa=ffa, day=day, features=cp_features)
        all_results.append(result_df)


# In[13]:


all_results_df = pd.concat(all_results).sort_values(by="cp_feature").reset_index(drop=True)
all_results_df = all_results_df.assign(neg_log_10_p=-1 * np.log10(all_results_df.p_val))

print(all_results_df.shape)
all_results_df.head()


# In[14]:


# Output results
file = os.path.join("results", "fat_globule_analysis_results.tsv")
all_results_df.to_csv(file, sep='\t', index=False)


# ## Quick Visualization

# In[15]:


gg.ggplot(all_results_df, gg.aes(x="t_stat", y="neg_log_10_p")) +     gg.facet_grid("FFA~diff_day") +     gg.geom_point()

