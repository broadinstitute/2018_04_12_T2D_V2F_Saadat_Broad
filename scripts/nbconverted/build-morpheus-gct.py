#!/usr/bin/env python
# coding: utf-8

# ## Generate Morpheus Input Data
# 
# **Gregory Way, 2019**
# 
# Use this script to concatenate all of the cell painting data into one `.gct` file for input into morpheus.

# In[1]:


import os
import pandas as pd


# In[2]:


get_ipython().run_line_magic('load_ext', 'rpy2.ipython')


# In[3]:


def load_files(plate_dir_list, file_grep="normalized_variable_selected.csv"):
    # Build full cell painting dataset
    df_list = []
    for plate_dir in plate_dirs:
        if ".DS_Store" in plate_dir:
            continue
        plate_files = os.listdir(plate_dir)
        for plate_file in plate_files:
            if file_grep in plate_file:
                plate_file = os.path.join(plate_dir, plate_file)
                df = pd.read_csv(plate_file)
                print("reading {} with profile count: {}".format(plate_file, df.shape[0]))
                df_list.append(df)
                
    cp_df = pd.concat(df_list).reset_index(drop=True)
    return cp_df


# ## Batch 1

# In[4]:


batch_id = "2019_04_16_Batch1"
backend_dir = os.path.join("..", "..", "backend", batch_id)

plate_dirs = [os.path.join(backend_dir, x) for x in os.listdir(backend_dir)]


# In[5]:


cp_df = load_files(plate_dirs, file_grep="normalized_variable_selected.csv")
cp_df.Metadata_diff_day = cp_df.Metadata_diff_day.astype(str)

print(cp_df.shape)
cp_df.head()


# In[6]:


# Extract out the cell painting features selected in this batch
batch_1_features = [x for x in cp_df.columns if not x.startswith("Metadata_")]


# In[7]:


# Output combined file
file = os.path.join("data", "combined_normalized_variable_selected.tsv")
cp_df.to_csv(file, index=False, sep='\t')


# # Collapse Data for Morpheus Heatmaps

# In[8]:


# Extract out day 15 data
cp_15_df = cp_df.query("Metadata_diff_day in ['15', '15+iso']").reset_index(drop=True)
print(cp_15_df.shape)

cp_non15_df = cp_df.query("Metadata_diff_day not in ['15', '15+iso']").reset_index(drop=True)
print(cp_non15_df.shape)


# In[9]:


replicate_cols = ["Metadata_Plate",
                  "Metadata_cell_line",
                  "Metadata_patient",
                  "Metadata_FFA",
                  "Metadata_diff_day"]

cp_15_collapsed_df = cp_15_df.groupby(replicate_cols).mean().reset_index()
cp_non15_collapsed_df = cp_non15_df.groupby(replicate_cols).mean().reset_index()


# ## Use `write_gct.R` to build the Moxrpheus Input

# In[10]:


get_ipython().run_cell_magic('R', '-i cp_df -i batch_id -i backend_dir -i cp_15_df -i cp_non15_df -i cp_15_collapsed_df -i cp_non15_collapsed_df', '\nlibrary(dplyr)\nlibrary(magrittr)\n\nfile <- file.path("..", "cytominer_scripts", "write_gct.R")\nsource(file)\n\noutput <- file.path("results", "morpheus",\n                    paste0("full_", batch_id, "_morpheus.gct"))\nchannels <- NULL\ncreate_row_annotations <- TRUE\nfeature_regex <- "^Nuclei_|^Cells_|^Cytoplasm_"\n\n# Step 1: Output combined gct file\nwrite_gct(x = cp_df,\n          path = output,\n          channels = channels,\n          create_row_annotations = create_row_annotations,\n          feature_regex = feature_regex)\n\n# Step 2: Output specific plate combinations\n# (replicate collapsed and non replicate collapsed)\noutput <- file.path("results", "morpheus",\n                    paste0("full_", batch_id, "_morpheus_day15.gct"))\n\nwrite_gct(x = cp_15_df,\n          path = output,\n          channels = channels,\n          create_row_annotations = create_row_annotations,\n          feature_regex = feature_regex)\n\noutput <- file.path("results", "morpheus",\n                    paste0("full_", batch_id, "_morpheus_nonday15.gct"))\n\nwrite_gct(x = cp_non15_df,\n          path = output,\n          channels = channels,\n          create_row_annotations = create_row_annotations,\n          feature_regex = feature_regex)\n\noutput <- file.path("results", "morpheus",\n                    paste0("full_", batch_id, "_morpheus_day15_collapsed.gct"))\n\nwrite_gct(x = cp_15_collapsed_df,\n          path = output,\n          channels = channels,\n          create_row_annotations = create_row_annotations,\n          feature_regex = feature_regex)\n\noutput <- file.path("results", "morpheus",\n                    paste0("full_", batch_id, "_morpheus_nonday15_collapsed.gct"))\n\nwrite_gct(x = cp_non15_collapsed_df,\n          path = output,\n          channels = channels,\n          create_row_annotations = create_row_annotations,\n          feature_regex = feature_regex)\n\n# Step 3: Also generate and write individual gct files\nplate_cols <- readr::cols(\n    .default = readr::col_double(),\n    Metadata_Plate = readr::col_character(),\n    Metadata_Well = readr::col_character(),\n    Metadata_Assay_Plate_Barcode = readr::col_character(),\n    Metadata_Plate_Map_Name = readr::col_character(),\n    Metadata_well_position = readr::col_character(),\n    Metadata_cell_line = readr::col_character(),\n    Metadata_patient = readr::col_character(),\n    Metadata_FFA = readr::col_character(),\n    Metadata_diff_day = readr::col_character()\n)\n\nall_plate_dirs <- list.files(backend_dir, full.names = TRUE)\nfor (plate_dir in all_plate_dirs) {\n    plate_file <- list.files(plate_dir, full.names = FALSE, pattern = "normalized_variable_selected")[1]\n    full_plate_file <- file.path(plate_dir, plate_file)\n\n    df <- readr::read_csv(full_plate_file, col_types = plate_cols)\n    \n    output_file <- file.path("results", "morpheus",\n                             paste0(tools::file_path_sans_ext(plate_file),\n                                    "_", batch_id, "_morpheus.gct"))\n    write_gct(x = df,\n              path = output_file,\n              channels = channels,\n              create_row_annotations = create_row_annotations,\n              feature_regex = feature_regex)\n}')


# ## Batch 2

# In[11]:


batch_id = "2019_06_11_Batch2"
file = os.path.join("data", "merged_profiles_{}.tsv.gz".format(batch_id))

cp_df = pd.read_csv(file, sep='\t')
print(cp_df.shape)
cp_df.head(3)


# In[12]:


# Merge replicate cols
replicate_cols = ["Metadata_Plate",
                  "Metadata_cell_line",
                  "Metadata_condition_O2",
                  "Metadata_treatment"]

cp_collapse_df = cp_df.groupby(replicate_cols).mean().reset_index()


# ## Use `write_gct.R` to build the Morpheus Input

# In[13]:


get_ipython().run_cell_magic('R', '-i cp_df -i batch_id -i cp_collapse_df', '\nlibrary(dplyr)\nlibrary(magrittr)\n\nfile <- file.path("..", "cytominer_scripts", "write_gct.R")\nsource(file)\n\noutput <- file.path("results", "morpheus",\n                    paste0("full_", batch_id, "_morpheus.gct"))\noutput_collapsed <- file.path("results", "morpheus",\n                              paste0("collapsed_", batch_id, "_morpheus.gct"))\n\nchannels <- NULL\ncreate_row_annotations <- TRUE\nfeature_regex <- "^Nuclei_|^Cells_|^Cytoplasm_"\n\n# Step 1: Output combined gct file\nwrite_gct(x = cp_df,\n          path = output,\n          channels = channels,\n          create_row_annotations = create_row_annotations,\n          feature_regex = feature_regex)\n\n# Step 2: Output replicate collapsed gct file\nwrite_gct(x = cp_collapse_df,\n          path = output_collapsed,\n          channels = channels,\n          create_row_annotations = create_row_annotations,\n          feature_regex = feature_regex)')


# ## Batch 3

# In[14]:


batch_id = "2019_08_06_Batch3"
backend_dir = os.path.join("..", "..", "backend", batch_id)

plate_dirs = [os.path.join(backend_dir, x) for x in os.listdir(backend_dir)]


# In[15]:


cp_df = load_files(plate_dirs, file_grep="normalized_variable_selected.csv")
cp_df.Metadata_diff_day = cp_df.Metadata_diff_day.astype(str)
cp_df.Metadata_patient = cp_df.Metadata_patient.astype(str)

print(cp_df.shape)
cp_df.head()


# In[16]:


# Extract out the cell painting features selected in this batch
batch_3_features = [x for x in cp_df.columns if not x.startswith("Metadata_")]


# In[17]:


# Output combined file
file = os.path.join("data", "{}_combined_normalized_variable_selected.tsv".format(batch_id))
cp_df.to_csv(file, index=False, sep='\t')


# In[18]:


replicate_cols = ["Metadata_Plate",
                  "Metadata_cell_line",
                  "Metadata_patient",
                  "Metadata_FFA",
                  "Metadata_diff_day"]

cp_collapsed_df = cp_df.groupby(replicate_cols).mean().reset_index()


# ## Use `write_gct.R` to build the Moxrpheus Input

# In[19]:


get_ipython().run_cell_magic('R', '-i cp_df -i batch_id -i cp_df -i cp_collapsed_df', '\nlibrary(dplyr)\nlibrary(magrittr)\n\nfile <- file.path("..", "cytominer_scripts", "write_gct.R")\nsource(file)\n\noutput <- file.path("results", "morpheus",\n                    paste0("full_", batch_id, "_morpheus.gct"))\noutput_collapsed <- file.path("results", "morpheus",\n                              paste0("collapsed_", batch_id, "_morpheus.gct"))\n\nchannels <- NULL\ncreate_row_annotations <- TRUE\nfeature_regex <- "^Nuclei_|^Cells_|^Cytoplasm_"\n\n# Step 1: Output combined gct file\nwrite_gct(x = cp_df,\n          path = output,\n          channels = channels,\n          create_row_annotations = create_row_annotations,\n          feature_regex = feature_regex)\n\n# Step 2: Output replicate collapsed gct file\nwrite_gct(x = cp_collapsed_df,\n          path = output_collapsed,\n          channels = channels,\n          create_row_annotations = create_row_annotations,\n          feature_regex = feature_regex)')


# ## Combine Batch 1 and Batch 3 Profiles

# In[20]:


combined_batch_features = list(set(batch_1_features + batch_3_features))
len(combined_batch_features)


# In[21]:


batch_id = "2019_04_16_Batch1"
backend_dir = os.path.join("..", "..", "backend", batch_id)

plate_dirs = [os.path.join(backend_dir, x) for x in os.listdir(backend_dir)]


# In[22]:


# Build full cell painting dataset
batch_one_df = load_files(plate_dirs, file_grep="normalized.csv").assign(Metadata_Batch="batch_one")
batch_one_df.Metadata_diff_day = batch_one_df.Metadata_diff_day.astype(str)

meta_cols = [x for x in batch_one_df.columns if x.startswith("Metadata")]

batch_one_df = batch_one_df.loc[:, meta_cols + combined_batch_features]

print(batch_one_df.shape)
batch_one_df.head(2)


# In[23]:


batch_id = "2019_08_06_Batch3"
backend_dir = os.path.join("..", "..", "backend", batch_id)

plate_dirs = [os.path.join(backend_dir, x) for x in os.listdir(backend_dir)]


# In[24]:


batch_three_df = load_files(plate_dirs, file_grep="normalized.csv").assign(Metadata_Batch="batch_three")
batch_three_df.Metadata_diff_day = batch_three_df.Metadata_diff_day.astype(str)
batch_three_df.Metadata_patient = batch_three_df.Metadata_patient.astype(str)

meta_cols = [x for x in batch_three_df.columns if x.startswith("Metadata")]

batch_three_df = batch_three_df.loc[:, meta_cols + combined_batch_features]

print(batch_three_df.shape)
batch_three_df.head(2)


# In[25]:


combined_df = pd.concat([batch_one_df, batch_three_df], axis="rows").dropna(axis="columns")
combined_df.Metadata_patient = combined_df.Metadata_patient.astype(str)

print(combined_df.shape)
combined_df.head(2)


# In[26]:


# Recode ER features to BODIPY
combined_df.columns = [x.replace("_ER", "_BODIPY") if "_ER" in x else x for x in combined_df.columns]


# In[27]:


# Output combined file
file = os.path.join("data", "batch1_batch3_combined_normalized_variable_selected.tsv")
combined_df.to_csv(file, index=False, sep='\t')


# ## Setup Data Specifically for Morpheus

# In[28]:


# Load labels
file = os.path.join("data", "category_labels.csv")
label_df = pd.read_csv(file)
label_df.columns = ["Metadata_{}".format(x) for x in label_df.columns]

label_df


# In[29]:


combined_df.loc[:, "Metadata_patient"] = ["m{}".format(x.strip("PAC_")) for x in combined_df.Metadata_patient]

# Merge labels
combined_df = label_df.merge(combined_df, left_on="Metadata_IID", right_on="Metadata_patient", how="inner")

print(combined_df.shape)
combined_df.head()


# In[30]:


combined_df = combined_df.query("Metadata_FFA == 0")


# In[31]:


replicate_cols = [
    "Metadata_Batch",
    "Metadata_Plate",
    "Metadata_cell_line",
    "Metadata_patient",
    "Metadata_FFA",
    "Metadata_diff_day"
]

combined_collapsed_df = combined_df.groupby(replicate_cols).mean().reset_index()

print(combined_collapsed_df.shape)
combined_collapsed_df.head(2)


# ## Use `write_gct.R` to output files for Morpheus

# In[32]:


get_ipython().run_cell_magic('R', '-i combined_df -i combined_collapsed_df', '\nlibrary(dplyr)\nlibrary(magrittr)\n\nfile <- file.path("..", "cytominer_scripts", "write_gct.R")\nsource(file)\n\noutput <- file.path("results", "morpheus",\n                    paste0("full_combined_batch1_batch3_morpheus.gct"))\noutput_collapsed <- file.path("results", "morpheus",\n                              paste0("collapsed_full_combined_batch1_batch3_morpheus.gct"))\n\nchannels <- NULL\ncreate_row_annotations <- TRUE\nfeature_regex <- "^Nuclei_|^Cells_|^Cytoplasm_"\n\n# Step 1: Output combined gct file\nwrite_gct(x = combined_df,\n          path = output,\n          channels = channels,\n          create_row_annotations = create_row_annotations,\n          feature_regex = feature_regex)\n\n# Step 2: Output replicate collapsed gct file\nwrite_gct(x = combined_collapsed_df,\n          path = output_collapsed,\n          channels = channels,\n          create_row_annotations = create_row_annotations,\n          feature_regex = feature_regex)')

