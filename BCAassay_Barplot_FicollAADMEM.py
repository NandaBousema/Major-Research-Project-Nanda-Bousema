#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  9 08:53:02 2025

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
Open: one excel file with Ficoll and one for DMEM and AA together
"""


file_path = 'file.xlsx'  
df_AA_DMEM = pd.read_excel(file_path, header=None)

# Inspect the data
print(df_AA_DMEM.head())

file_path = 'file.xlsx' 
df_Ficoll = pd.read_excel(file_path, header=None)

# Inspect the data
print(df_Ficoll.head())

#%%
"""
Calculate Blank corrected OD data for DMEM and AA
"""
# Select blank OD data
OD_Blank_AA_DMEM = df_AA_DMEM.iloc[25:28, 4]
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
OD_Blank_Ficoll = df_Ficoll.iloc[25:28, 4]
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
Con_Standard_AA_DMEM = df_AA_DMEM.iloc[14:22, 1:4].apply(pd.to_numeric, errors="coerce").dropna(how="any")
Con_Standard_means_AA_DMEM = Con_Standard_AA_DMEM.mean(axis=1).reset_index(drop=True)  

# Standard OD data and calculate means
OD_Standard_cols_AA_DMEM = df_AA_DMEM.iloc[25:33, 1:4].apply(pd.to_numeric, errors="coerce").dropna(how="any")  

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
Con_Standard_Ficoll = df_Ficoll.iloc[14:22, 1:4].apply(pd.to_numeric, errors="coerce").dropna(how="any")  
Con_Standard_means_Ficoll = Con_Standard_Ficoll.mean(axis=1).reset_index(drop=True)   

# Standard OD data and calculate means
OD_Standard_cols_Ficoll = df_Ficoll.iloc[25:33, 1:4].apply(pd.to_numeric, errors="coerce").dropna(how="any")    

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

"""
DMEM
"""

OD_CDM_DMEM = {
    '800': pd.to_numeric(df_AA_DMEM.iloc[25:28, 5], errors='coerce').dropna(),
    '500' : pd.to_numeric(df_AA_DMEM.iloc[25:28, 6], errors='coerce').dropna(),
    '300' : pd.to_numeric(df_AA_DMEM.iloc[25:28, 7], errors='coerce').dropna(),
    '200' : pd.to_numeric(df_AA_DMEM.iloc[25:28, 8], errors='coerce').dropna(),
    '100' : pd.to_numeric(df_AA_DMEM.iloc[25:28, 9], errors='coerce').dropna(),
    '50' : pd.to_numeric(df_AA_DMEM.iloc[25:28, 10], errors='coerce').dropna(),
    '25'  : pd.to_numeric(df_AA_DMEM.iloc[25:28, 11], errors='coerce').dropna(),
    '10'  : pd.to_numeric(df_AA_DMEM.iloc[25:28, 12], errors='coerce').dropna(),
    '5'   : pd.Series(pd.to_numeric(df_AA_DMEM.iloc[28, 5:8], errors='coerce')).dropna()
}


OD_DMEM_corrected = {
    label: (values - mean_blankOD_AA_DMEM)
    for label, values in OD_CDM_DMEM.items()
}

Conc_CDM_DMEM = {
    '800': pd.to_numeric(df_AA_DMEM.iloc[14:17, 5], errors='coerce').dropna(),
    '500' : pd.to_numeric(df_AA_DMEM.iloc[14:17, 6], errors='coerce').dropna(),
    '300' : pd.to_numeric(df_AA_DMEM.iloc[14:17, 7], errors='coerce').dropna(),
    '200' : pd.to_numeric(df_AA_DMEM.iloc[14:17, 8], errors='coerce').dropna(),
    '100' : pd.to_numeric(df_AA_DMEM.iloc[14:17, 9], errors='coerce').dropna(),
    '50' : pd.to_numeric(df_AA_DMEM.iloc[14:17, 10], errors='coerce').dropna(),
    '25'  : pd.to_numeric(df_AA_DMEM.iloc[14:17, 11], errors='coerce').dropna(),
    '10'  : pd.to_numeric(df_AA_DMEM.iloc[14:17, 12], errors='coerce').dropna(),
    '5'   : pd.Series(pd.to_numeric(df_AA_DMEM.iloc[17, 5:8], errors='coerce')).dropna()
}

"""
AA
"""

OD_CDM_AA = {
    '800': pd.to_numeric(df_AA_DMEM.iloc[30:33, 4], errors='coerce').dropna(),
    '500' : pd.to_numeric(df_AA_DMEM.iloc[30:33, 5], errors='coerce').dropna(),
    '300' : pd.to_numeric(df_AA_DMEM.iloc[30:33, 6], errors='coerce').dropna(),
    '200' : pd.to_numeric(df_AA_DMEM.iloc[30:33, 7], errors='coerce').dropna(),
    '100' : pd.to_numeric(df_AA_DMEM.iloc[30:33, 8], errors='coerce').dropna(),
    '50' : pd.to_numeric(df_AA_DMEM.iloc[30:33, 9], errors='coerce').dropna(),
    '25'  : pd.to_numeric(df_AA_DMEM.iloc[30:33, 10], errors='coerce').dropna(),
    '10'  : pd.to_numeric(df_AA_DMEM.iloc[30:33, 11], errors='coerce').dropna(),
    '5'   : pd.to_numeric(df_AA_DMEM.iloc[30:33, 12], errors='coerce').dropna()
}


OD_AA_corrected = {
    label: (values - mean_blankOD_AA_DMEM)
    for label, values in OD_CDM_AA.items()
}

Conc_CDM_AA = {
    '800': pd.to_numeric(df_AA_DMEM.iloc[19:22, 4], errors='coerce').dropna(),
    '500' : pd.to_numeric(df_AA_DMEM.iloc[19:22, 5], errors='coerce').dropna(),
    '300' : pd.to_numeric(df_AA_DMEM.iloc[19:22, 6], errors='coerce').dropna(),
    '200' : pd.to_numeric(df_AA_DMEM.iloc[19:22, 7], errors='coerce').dropna(),
    '100' : pd.to_numeric(df_AA_DMEM.iloc[19:22, 8], errors='coerce').dropna(),
    '50' : pd.to_numeric(df_AA_DMEM.iloc[19:22, 9], errors='coerce').dropna(),
    '25'  : pd.to_numeric(df_AA_DMEM.iloc[19:22, 10], errors='coerce').dropna(),
    '10'  : pd.to_numeric(df_AA_DMEM.iloc[19:22, 11], errors='coerce').dropna(),
    '5'   : pd.to_numeric(df_AA_DMEM.iloc[19:22, 12], errors='coerce').dropna()
}

"""
Ficoll
"""

OD_CDM_Ficoll = {
    '800': pd.to_numeric(df_Ficoll.iloc[30:33, 4], errors='coerce').dropna(),
    '500' : pd.to_numeric(df_Ficoll.iloc[30:33, 5], errors='coerce').dropna(),
    '300' : pd.to_numeric(df_Ficoll.iloc[30:33, 6], errors='coerce').dropna(),
    '200' : pd.to_numeric(df_Ficoll.iloc[30:33, 7], errors='coerce').dropna(),
    '100' : pd.to_numeric(df_Ficoll.iloc[30:33, 8], errors='coerce').dropna(),
    '50' : pd.to_numeric(df_Ficoll.iloc[30:33, 9], errors='coerce').dropna(),
    '25'  : pd.to_numeric(df_Ficoll.iloc[30:33, 10], errors='coerce').dropna(),
    '10'  : pd.to_numeric(df_Ficoll.iloc[30:33, 11], errors='coerce').dropna(),
    '5'   : pd.to_numeric(df_Ficoll.iloc[30:33, 12], errors='coerce').dropna()
}

OD_Ficoll_corrected = {
    label: (values - mean_blankOD_Ficoll)
    for label, values in OD_CDM_Ficoll.items()
}

Conc_CDM_Ficoll = {
    '800': pd.to_numeric(df_Ficoll.iloc[19:22, 4], errors='coerce').dropna(),
    '500' : pd.to_numeric(df_Ficoll.iloc[19:22, 5], errors='coerce').dropna(),
    '300' : pd.to_numeric(df_Ficoll.iloc[19:22, 6], errors='coerce').dropna(),
    '200' : pd.to_numeric(df_Ficoll.iloc[19:22, 7], errors='coerce').dropna(),
    '100' : pd.to_numeric(df_Ficoll.iloc[19:22, 8], errors='coerce').dropna(),
    '50' : pd.to_numeric(df_Ficoll.iloc[19:22, 9], errors='coerce').dropna(),
    '25'  : pd.to_numeric(df_Ficoll.iloc[19:22, 10], errors='coerce').dropna(),
    '10'  : pd.to_numeric(df_Ficoll.iloc[19:22, 11], errors='coerce').dropna(),
    '5'   : pd.to_numeric(df_Ficoll.iloc[19:22, 12], errors='coerce').dropna()
}

#%%

# Example: your datasets
datasets = [
    ('AA', OD_AA_corrected, Conc_CDM_AA),
    ('Ficoll', OD_Ficoll_corrected, Conc_CDM_Ficoll),
    ('DMEM', OD_DMEM_corrected, Conc_CDM_DMEM)
]

# Initialize an empty DataFrame
df_all = pd.DataFrame()

for condition, od_dict, conc_dict in datasets:
    for cdm_conc in od_dict.keys():  # e.g., '4000', '800', etc.
        od_values = od_dict[cdm_conc]
        prot_values = conc_dict[cdm_conc]

        # Make a temporary DataFrame
        df_temp = pd.DataFrame({
            'Condition': condition,
            'CDM_Concentration': int(cdm_conc),  # convert to int if you like
            'Protein_Concentration': prot_values.values,
            'OD': od_values.values
        })

        # Append to main DataFrame
        df_all = pd.concat([df_all, df_temp], ignore_index=True)

# Optional: remove negative OD or protein concentrations
df_all = df_all[(df_all['OD'] >= 0) & (df_all['Protein_Concentration'] >= 0)]

# Check
print(df_all)



# Group by condition and concentration, then calculate mean, std, CV
stats_protein = df_all.groupby(['Condition', 'CDM_Concentration'])['Protein_Concentration'].agg(
    mean_protein='mean',
    std_protein='std'
).reset_index()

# Calculate coefficient of variation
stats_protein['CV_protein'] = stats_protein['std_protein'] / stats_protein['mean_protein'] * 100

stats_protein


#%%
"""
Plot CDM sample data in the plot of the standard curve with seaborn
"""

from statannotations.Annotator import Annotator
import matplotlib.pyplot as plt
import seaborn as sns

plt.figure(figsize=(10, 8))
# Barplot
sns.set(style="whitegrid")
ax = sns.barplot(
    data=df_all,
    x='CDM_Concentration',
    y='Protein_Concentration',
    hue='Condition',
    palette='flare',
    ci='sd',
    capsize=0.1,
    errwidth=2.5
)

# Define pairs to compare per concentration
pairs = []
for conc in df_all['CDM_Concentration'].unique():
    conditions = df_all['Condition'].unique()
    for i in range(len(conditions)):
        for j in range(i+1, len(conditions)):
            pairs.append( ((conc, conditions[i]), (conc, conditions[j])) )

# Create the annotator
annot = Annotator(
    ax, pairs,
    data=df_all,
    x='CDM_Concentration',
    y='Protein_Concentration',
    hue='Condition'
)
annot.configure(
    test='t-test_ind', 
    text_format='star', 
    loc='inside', 
    verbose=2,
    hide_non_significant=True,  # <-- this hides all "ns"
    fontsize=25
)

# Compute tests and annotate in one step
annot.apply_and_annotate()

plt.xticks(fontsize=20, weight='bold', color='black')
plt.yticks(fontsize=20, weight='bold')
plt.xlabel("CDM Concentration (µg/mL)", fontsize=20, weight='bold')
plt.ylabel("Total Protein Concentration (µg/mL)", fontsize=20, weight='bold')
plt.title("Total Protein Concentration in CDM across Different Conditions", fontsize=20, weight='bold')
legend = plt.legend(title="Condition", title_fontsize=16, fontsize=16)
plt.setp(legend.get_title(), fontweight='bold')
plt.tight_layout()
plt.show()
























