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


# ## Batch 1

# In[3]:


batch_id = "2019_04_16_Batch1"
backend_dir = os.path.join("..", "..", "backend", batch_id)

plate_dirs = [os.path.join(backend_dir, x) for x in os.listdir(backend_dir)]


# In[4]:


# Build full cell painting dataset
df_list = []
for plate_dir in plate_dirs:
    plate_files = os.listdir(plate_dir)
    for plate_file in plate_files:
        if "normalized_variable_selected.csv" in plate_file:
            plate_file = os.path.join(plate_dir, plate_file)
            df = pd.read_csv(plate_file)
            print("reading {} with profile count: {}".format(plate_file, df.shape[0]))
            df_list.append(df)


# In[5]:


# Combine into a single file
cp_df = pd.concat(df_list).reset_index(drop=True)
cp_df.Metadata_diff_day = cp_df.Metadata_diff_day.astype(str)

print(cp_df.shape)
cp_df.head()


# In[6]:


# Output combined file
file = os.path.join("data", "combined_normalized_variable_selected.tsv")
cp_df.to_csv(file, index=False, sep='\t')


# # Collapse Data for Morpheus Heatmaps

# In[7]:


# Extract out day 15 data
cp_15_df = cp_df.query("Metadata_diff_day in ['15', '15+iso']").reset_index(drop=True)
print(cp_15_df.shape)

cp_non15_df = cp_df.query("Metadata_diff_day not in ['15', '15+iso']").reset_index(drop=True)
print(cp_non15_df.shape)


# In[8]:


replicate_cols = ["Metadata_Plate",
                  "Metadata_cell_line",
                  "Metadata_patient",
                  "Metadata_FFA",
                  "Metadata_diff_day"]

cp_15_collapsed_df = cp_15_df.groupby(replicate_cols).mean().reset_index()
cp_non15_collapsed_df = cp_non15_df.groupby(replicate_cols).mean().reset_index()


# ## Use `write_gct.R` to build the Moxrpheus Input

# In[9]:


get_ipython().run_cell_magic('R', '-i cp_df -i batch_id -i backend_dir -i cp_15_df -i cp_non15_df -i cp_15_collapsed_df -i cp_non15_collapsed_df', '\nlibrary(dplyr)\nlibrary(magrittr)\n\nfile <- file.path("..", "cytominer_scripts", "write_gct.R")\nsource(file)\n\noutput <- file.path("results", "morpheus",\n                    paste0("full_", batch_id, "_morpheus.gct"))\nchannels <- NULL\ncreate_row_annotations <- TRUE\nfeature_regex <- "^Nuclei_|^Cells_|^Cytoplasm_"\n\n# Step 1: Output combined gct file\nwrite_gct(x = cp_df,\n          path = output,\n          channels = channels,\n          create_row_annotations = create_row_annotations,\n          feature_regex = feature_regex)\n\n# Step 2: Output specific plate combinations\n# (replicate collapsed and non replicate collapsed)\noutput <- file.path("results", "morpheus",\n                    paste0("full_", batch_id, "_morpheus_day15.gct"))\n\nwrite_gct(x = cp_15_df,\n          path = output,\n          channels = channels,\n          create_row_annotations = create_row_annotations,\n          feature_regex = feature_regex)\n\noutput <- file.path("results", "morpheus",\n                    paste0("full_", batch_id, "_morpheus_nonday15.gct"))\n\nwrite_gct(x = cp_non15_df,\n          path = output,\n          channels = channels,\n          create_row_annotations = create_row_annotations,\n          feature_regex = feature_regex)\n\noutput <- file.path("results", "morpheus",\n                    paste0("full_", batch_id, "_morpheus_day15_collapsed.gct"))\n\nwrite_gct(x = cp_15_collapsed_df,\n          path = output,\n          channels = channels,\n          create_row_annotations = create_row_annotations,\n          feature_regex = feature_regex)\n\noutput <- file.path("results", "morpheus",\n                    paste0("full_", batch_id, "_morpheus_nonday15_collapsed.gct"))\n\nwrite_gct(x = cp_non15_collapsed_df,\n          path = output,\n          channels = channels,\n          create_row_annotations = create_row_annotations,\n          feature_regex = feature_regex)\n\n# Step 3: Also generate and write individual gct files\nplate_cols <- readr::cols(\n    .default = readr::col_double(),\n    Metadata_Plate = readr::col_character(),\n    Metadata_Well = readr::col_character(),\n    Metadata_Assay_Plate_Barcode = readr::col_character(),\n    Metadata_Plate_Map_Name = readr::col_character(),\n    Metadata_well_position = readr::col_character(),\n    Metadata_cell_line = readr::col_character(),\n    Metadata_patient = readr::col_character(),\n    Metadata_FFA = readr::col_character(),\n    Metadata_diff_day = readr::col_character()\n)\n\nall_plate_dirs <- list.files(backend_dir, full.names = TRUE)\nfor (plate_dir in all_plate_dirs) {\n    plate_file <- list.files(plate_dir, full.names = FALSE, pattern = "normalized_variable_selected")[1]\n    full_plate_file <- file.path(plate_dir, plate_file)\n\n    df <- readr::read_csv(full_plate_file, col_types = plate_cols)\n    \n    output_file <- file.path("results", "morpheus",\n                             paste0(tools::file_path_sans_ext(plate_file),\n                                    "_", batch_id, "_morpheus.gct"))\n    write_gct(x = df,\n              path = output_file,\n              channels = channels,\n              create_row_annotations = create_row_annotations,\n              feature_regex = feature_regex)\n}')


# ## Batch 2

# In[10]:


batch_id = "2019_06_11_Batch2"
file = os.path.join("data", "merged_profiles_{}.tsv.gz".format(batch_id))

cp_df = pd.read_csv(file, sep='\t')
print(cp_df.shape)
cp_df.head(3)


# In[11]:


# Merge replicate cols
replicate_cols = ["Metadata_Plate",
                  "Metadata_cell_line",
                  "Metadata_condition_O2",
                  "Metadata_treatment"]

cp_collapse_df = cp_df.groupby(replicate_cols).mean().reset_index()


# ## Use `write_gct.R` to build the Morpheus Input

# In[12]:


get_ipython().run_cell_magic('R', '-i cp_df -i batch_id -i cp_collapse_df', '\nlibrary(dplyr)\nlibrary(magrittr)\n\nfile <- file.path("..", "cytominer_scripts", "write_gct.R")\nsource(file)\n\noutput <- file.path("results", "morpheus",\n                    paste0("full_", batch_id, "_morpheus.gct"))\noutput_collapsed <- file.path("results", "morpheus",\n                              paste0("collapsed_", batch_id, "_morpheus.gct"))\n\nchannels <- NULL\ncreate_row_annotations <- TRUE\nfeature_regex <- "^Nuclei_|^Cells_|^Cytoplasm_"\n\n# Step 1: Output combined gct file\nwrite_gct(x = cp_df,\n          path = output,\n          channels = channels,\n          create_row_annotations = create_row_annotations,\n          feature_regex = feature_regex)\n\n# Step 2: Output replicate collapsed gct file\nwrite_gct(x = cp_collapse_df,\n          path = output_collapsed,\n          channels = channels,\n          create_row_annotations = create_row_annotations,\n          feature_regex = feature_regex)')


# ## Batch 3

# In[13]:


batch_id = "2019_08_06_Batch3"
backend_dir = os.path.join("..", "..", "backend", batch_id)

plate_dirs = [os.path.join(backend_dir, x) for x in os.listdir(backend_dir)]


# In[14]:


# Build full cell painting dataset
df_list = []
for plate_dir in plate_dirs:
    plate_files = os.listdir(plate_dir)
    for plate_file in plate_files:
        if "normalized_variable_selected.csv" in plate_file:
            plate_file = os.path.join(plate_dir, plate_file)
            df = pd.read_csv(plate_file)
            print("reading {} with profile count: {}".format(plate_file, df.shape[0]))
            df_list.append(df)


# In[15]:


# Combine into a single file
cp_df = pd.concat(df_list).reset_index(drop=True)
cp_df.Metadata_diff_day = cp_df.Metadata_diff_day.astype(str)

print(cp_df.shape)
cp_df.head()


# In[16]:


# Output combined file
file = os.path.join("data", "{}_combined_normalized_variable_selected.tsv".format(batch_id))
cp_df.to_csv(file, index=False, sep='\t')


# In[17]:


replicate_cols = ["Metadata_Plate",
                  "Metadata_cell_line",
                  "Metadata_patient",
                  "Metadata_FFA",
                  "Metadata_diff_day"]

cp_collapsed_df = cp_df.groupby(replicate_cols).mean().reset_index()


# ## Use `write_gct.R` to build the Moxrpheus Input

# In[18]:


get_ipython().run_cell_magic('R', '-i cp_df -i batch_id -i cp_df -i cp_collapsed_df', '\nlibrary(dplyr)\nlibrary(magrittr)\n\nfile <- file.path("..", "cytominer_scripts", "write_gct.R")\nsource(file)\n\noutput <- file.path("results", "morpheus",\n                    paste0("full_", batch_id, "_morpheus.gct"))\noutput_collapsed <- file.path("results", "morpheus",\n                              paste0("collapsed_", batch_id, "_morpheus.gct"))\n\nchannels <- NULL\ncreate_row_annotations <- TRUE\nfeature_regex <- "^Nuclei_|^Cells_|^Cytoplasm_"\n\n# Step 1: Output combined gct file\nwrite_gct(x = cp_df,\n          path = output,\n          channels = channels,\n          create_row_annotations = create_row_annotations,\n          feature_regex = feature_regex)\n\n# Step 2: Output replicate collapsed gct file\nwrite_gct(x = cp_collapsed_df,\n          path = output_collapsed,\n          channels = channels,\n          create_row_annotations = create_row_annotations,\n          feature_regex = feature_regex)')

