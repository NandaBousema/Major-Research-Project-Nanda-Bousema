#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  9 08:48:02 2025

@author: nandabousema
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 19 18:42:49 2025

@author: nandabousema
"""

import pandas as pd
import openpyxl as px
import matplotlib.pyplot as plt
import os
import cv2
import math
from skimage.measure import label, regionprops
from skimage.color import label2rgb
import numpy as np
from scipy import ndimage as ndi
from scipy import stats
from skimage.feature import peak_local_max
from skimage.segmentation import watershed
from skimage.segmentation import clear_border
from statistics import mean
import seaborn as sns
from matplotlib.ticker import FormatStrFormatter

#pip install pandas matplotlib openpyxl > dit is voor in de terminal

#%%
"""
Open excel files: seperate one for Ficoll, DMEM and AA together.
"""

file_path = 'file.xlsx'  # Replace with your Excel file's path
df_AA_DMEM = pd.read_excel(file_path, header=None)

# Inspect the data
print(df_AA_DMEM.head())

file_path = 'file.xlsx'  # Replace with your Excel file's path
df_Ficoll = pd.read_excel(file_path, header=None)

# Inspect the data
print(df_Ficoll.head())

#%%
"""
Calculate Blank corrected OD data for DMEM and AA
"""
# Select blank OD data
OD_Blank_AA_DMEM = df_AA_DMEM.iloc[55, 1:4]
OD_Blank_means_AA_DMEM = OD_Blank_AA_DMEM.mean()    

print(OD_Blank_AA_DMEM)
print(OD_Blank_means_AA_DMEM)


# Compute the mean across replicates for each row
mean_blankOD_AA_DMEM = OD_Blank_means_AA_DMEM

print("Mean Blank OD AA DMEM:", mean_blankOD_AA_DMEM)

"""
Calculate Blank corrected OD data for Ficoll
"""
# Select blank OD data
OD_Blank_Ficoll = df_Ficoll.iloc[66, 1:4]
OD_Blank_means_Ficoll = OD_Blank_Ficoll.mean()    

print(OD_Blank_Ficoll)
print(OD_Blank_means_Ficoll)


# Compute the mean across replicates for each row
mean_blankOD_Ficoll = OD_Blank_means_Ficoll

print("Mean Blank OD Ficoll:", mean_blankOD_Ficoll)

#%%
"""
Make dataframe for standard curve AA DMEM
"""
# Concentrations for the standard curve and calculate means
Con_Standard_AA_DMEM = df_AA_DMEM.iloc[26:33, 1:4]     # Rows, Columns
Con_Standard_means_AA_DMEM = Con_Standard_AA_DMEM.mean(axis=1).reset_index(drop=True)    

# Standard OD data and calculate means
OD_Standard_cols_AA_DMEM = df_AA_DMEM.iloc[48:55, 1:4]     # Rows 1–8, Columns 2–4

# Blank Corrected Data
OD_Standard_AA_DMEM = OD_Standard_cols_AA_DMEM - mean_blankOD_AA_DMEM

OD_Standard_means_AA_DMEM = OD_Standard_AA_DMEM.mean(axis=1).reset_index(drop=True)  

standard_df_AA_DMEM = pd.DataFrame({
    'Collagen Concentration (µg/mL)' : Con_Standard_means_AA_DMEM,
    'Absorbance (OD)' : OD_Standard_means_AA_DMEM 
    
})

# Dtype from above is object, error='coerce' will turn invalid entries like text of empry strings into NaN
standard_df_AA_DMEM['Absorbance (OD)'] = pd.to_numeric(standard_df_AA_DMEM['Absorbance (OD)'], errors='coerce')
standard_df_AA_DMEM['Collagen Concentration (µg/mL)'] = pd.to_numeric(standard_df_AA_DMEM['Collagen Concentration (µg/mL)'], errors='coerce')

# Drop rows with NaNs in either column
standard_df_AA_DMEM = standard_df_AA_DMEM.dropna(subset=['Absorbance (OD)', 'Collagen Concentration (µg/mL)'])


# Plot standard curve with best-fit line to predict concentrations later on
# Extract values
x = standard_df_AA_DMEM['Collagen Concentration (µg/mL)'].values
y = standard_df_AA_DMEM['Absorbance (OD)'].values

# Linear regression
slope_AA_DMEM, intercept_AA_DMEM, r_value_AA_DMEM, p_value_AA_DMEM, std_err_AA_DMEM = stats.linregress(x, y)

"""
Make dataframe for standard curve Ficoll
"""
# Concentrations for the standard curve and calculate means
Con_Standard_Ficoll = df_Ficoll.iloc[37:44, 1:4]     # Rows, Columns
Con_Standard_means_Ficoll = Con_Standard_Ficoll.mean(axis=1).reset_index(drop=True)     

# Standard OD data and calculate means
OD_Standard_cols_Ficoll = df_Ficoll.iloc[59:66, 1:4]     # Rows 1–8, Columns 2–4

# Blank Corrected Data
OD_Standard_Ficoll = OD_Standard_cols_Ficoll - mean_blankOD_Ficoll

OD_Standard_means_Ficoll = OD_Standard_Ficoll.mean(axis=1).reset_index(drop=True)    

standard_df_Ficoll = pd.DataFrame({
    'Collagen Concentration (µg/mL)' : Con_Standard_means_Ficoll,
    'Absorbance (OD)' : OD_Standard_means_Ficoll 
    
})

# Dtype from above is object, error='coerce' will turn invalid entries like text of empry strings into NaN
standard_df_Ficoll['Absorbance (OD)'] = pd.to_numeric(standard_df_Ficoll['Absorbance (OD)'], errors='coerce')
standard_df_Ficoll['Collagen Concentration (µg/mL)'] = pd.to_numeric(standard_df_Ficoll['Collagen Concentration (µg/mL)'], errors='coerce')

# Drop rows with NaNs in either column
standard_df_Ficoll = standard_df_Ficoll.dropna(subset=['Absorbance (OD)', 'Collagen Concentration (µg/mL)'])

# Plot standard curve with best-fit line to predict concentrations later on
# Extract values
x = standard_df_Ficoll['Collagen Concentration (µg/mL)'].values
y = standard_df_Ficoll['Absorbance (OD)'].values

# Linear regression
slope_ficoll, intercept_ficoll, r_value_ficoll, p_value_ficoll, std_err_ficoll = stats.linregress(x, y)


#%%
"""
Make dataframe for results CDM samples AA DMEM Ficoll
"""

# Mean OD values CDM samples and blank corrected
OD_CDM_DMEM = {
   '4000' : df_AA_DMEM.iloc[48:51, 4],
   '800' : df_AA_DMEM.iloc[48:51, 5],
   '500' : df_AA_DMEM.iloc[48:51, 6],
   '300' : df_AA_DMEM.iloc[48:51, 7],
   '200' : df_AA_DMEM.iloc[48:51, 8],
   '100' : df_AA_DMEM.iloc[48:51, 9],
   '50' : df_AA_DMEM.iloc[48:51, 10],
   '10' : df_AA_DMEM.iloc[48:51, 11],
   '5' : df_AA_DMEM.iloc[48:51, 12]
}

OD_DMEM_corrected = {
    label: (values - mean_blankOD_AA_DMEM)
    for label, values in OD_CDM_DMEM.items()
}

OD_CDM_AA = {
   '4000' : df_AA_DMEM.iloc[52:55, 4],
   '800' : df_AA_DMEM.iloc[52:55, 5],
   '500' : df_AA_DMEM.iloc[52:55, 6],
   '300' : df_AA_DMEM.iloc[52:55, 7],
   '200' : df_AA_DMEM.iloc[52:55, 8],
   '100' : df_AA_DMEM.iloc[52:55, 9],
   '50' : df_AA_DMEM.iloc[52:55, 10],
   '10' : df_AA_DMEM.iloc[52:55, 11],
   '5' : df_AA_DMEM.iloc[52:55, 12]
}

OD_AA_corrected = {
    label: (values - mean_blankOD_AA_DMEM)
    for label, values in OD_CDM_AA.items()
}


OD_CDM_Ficoll = {
   '4000' : df_Ficoll.iloc[59:62, 4],
   '800' : df_Ficoll.iloc[59:62, 5],
   '500' : df_Ficoll.iloc[59:62, 6],
   '300' : df_Ficoll.iloc[59:62, 7],
   '200' : df_Ficoll.iloc[59:62, 8],
   '100' : df_Ficoll.iloc[59:62, 9],
   '50' : df_Ficoll.iloc[59:62, 10],
   '10' : df_Ficoll.iloc[59:62, 11],
   '5' : df_Ficoll.iloc[59:62, 12]
}

OD_Ficoll_corrected = {
    label: (values - mean_blankOD_Ficoll)
    for label, values in OD_CDM_Ficoll.items()
}


# 1. Build tidy DataFrame from corrected OD dicts
rows = []
for cond_name, corrected_dict in zip(
    ['DMEM', 'AA', 'Ficoll'],
    [OD_DMEM_corrected, OD_AA_corrected, OD_Ficoll_corrected]
):
    for conc_label, values in corrected_dict.items():
        for v in values:
            rows.append({
                'CDM_conc_label': conc_label,
                'Condition': cond_name,
                'OD': v
            })

df_all_data = pd.DataFrame(rows)

# 2. Predict collagen concentrations
def predict_conc(row):
    if row['Condition'] == 'DMEM':
        return (row['OD'] - intercept_AA_DMEM) / slope_AA_DMEM
    elif row['Condition'] == 'AA':
        return (row['OD'] - intercept_AA_DMEM) / slope_AA_DMEM
    elif row['Condition'] == 'Ficoll':
        return (row['OD'] - intercept_ficoll) / slope_ficoll
    else:
        return np.nan

df_all_data['Collagen_µg_mL'] = df_all_data.apply(predict_conc, axis=1)

# 3. Exclude negative concentrations
df_all_data = df_all_data[df_all_data['Collagen_µg_mL'] >= 0]


# 5. Group by concentration and condition to calculate mean, std, CV
df_stats = df_all_data.groupby(
    ['CDM_conc_label', 'Condition'], sort=False
)['Collagen_µg_mL'].agg(['mean','std']).reset_index()

df_stats_clean = df_stats.dropna(subset=['std']).copy()

# Create the CV column
df_stats_clean['CV_%'] = (df_stats_clean['std'] / df_stats_clean['mean']) * 100

print(df_stats_clean)


#%%
"""
Make a barplot
"""

# Rename columns for clarity in plotting
df_plot = df_stats_clean.rename(columns={'mean': 'Mean', 'std': 'Std', 'CDM_conc_label': 'Concentration'})

# Sort concentrations numerically (optional)
def conc_to_numeric(label):
    return float(label.split()[0])  # assumes labels like '5 µg/mL'

df_plot['Conc_numeric'] = df_plot['Concentration'].apply(conc_to_numeric)
df_plot = df_plot.sort_values('Conc_numeric')

# Set style
sns.set(style="whitegrid")
palette = sns.color_palette("flare", n_colors=3)

plt.figure(figsize=(8,6))
ax = sns.barplot(
    data=df_plot,
    x='Concentration',
    y='Mean',
    hue='Condition',
    palette=palette,
    ci=None
)

# Add error bars manually
for container, cond in zip(ax.containers, df_plot['Condition'].unique()):
    for bar, (_, row) in zip(container, df_plot[df_plot['Condition'] == cond].iterrows()):
        bar_height = bar.get_height()
        ax.errorbar(
            bar.get_x() + bar.get_width()/2,
            bar_height,
            yerr=row['Std'],
            color='black',
            capsize=4,
            elinewidth=2.5,
            capthick=2
        )

plt.xticks(fontsize=20, weight='bold')
plt.yticks(fontsize=20, weight='bold')
plt.ylabel("Collagen Concentration (µg/mL)", fontsize=20, weight='bold')
plt.xlabel("CDM Concentration (µg/mL)", fontsize=20, weight='bold')
plt.title("Collagen Content in CDM across Different Conditions", fontsize=20, weight='bold')
legend = plt.legend(title="Condition", title_fontsize=16, fontsize=16, loc="upper left")
plt.setp(legend.get_title(), fontweight='bold')
plt.show()


