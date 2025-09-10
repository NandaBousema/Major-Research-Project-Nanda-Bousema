#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  9 08:50:14 2025

@author: nandabousema
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 19 10:05:30 2025

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
Open
"""

file_path = 'file.xlsx'
#pip install xlrd>=2.0.1
df = pd.read_excel(file_path, header=None)

# Inspect the data
print(df.head())

#%%
# Calculate Blank corrected OD data
# Select blank OD data
OD_Blank_cols = df.iloc[43, 1:4]
OD_Blank_means = OD_Blank_cols.mean()    

print(OD_Blank_cols)
print(OD_Blank_means)


# Compute the mean across replicates for each row
mean_blankOD = OD_Blank_means

print("Mean Blank OD:", mean_blankOD)


#%%
"""
Make dataframe for standard curve
"""
# Concentrations for the standard curve and calculate means
Con_Standard_cols = df.iloc[47:54, 1:4]     # Rows, Columns
Con_Standard_means = Con_Standard_cols.mean(axis=1).reset_index(drop=True)    

# Standard OD data and calculate means
OD_Standard_cols = df.iloc[36:43, 1:4]     # Rows 1–8, Columns 2–4

# Blank Corrected Data
OD_Standard = OD_Standard_cols - mean_blankOD

OD_Standard_means = OD_Standard.mean(axis=1).reset_index(drop=True)    

standard_df = pd.DataFrame({
    'Collagen Concentration (µg/mL)' : Con_Standard_means,
    'Absorbance (OD)' : OD_Standard_means 
    
})

# Dtype from above is object, error='coerce' will turn invalid entries like text of empry strings into NaN
standard_df['Absorbance (OD)'] = pd.to_numeric(standard_df['Absorbance (OD)'], errors='coerce')
standard_df['Collagen Concentration (µg/mL)'] = pd.to_numeric(standard_df['Collagen Concentration (µg/mL)'], errors='coerce')

# Drop rows with NaNs in either column
standard_df = standard_df.dropna(subset=['Absorbance (OD)', 'Collagen Concentration (µg/mL)'])

#%%
"""
Create plot standard curve with seaborn
"""

# Choose palette
colors = sns.color_palette("flare", 6)
my_color = colors[0]  # the first color, the first color is 0, second is 1, etc. so last one is 8
sns.palplot(colors)

# Set Seaborn theme
sns.set_theme(style="darkgrid", context = "talk")

# Plot standard curve with best-fit line
# Extract values
x = standard_df['Collagen Concentration (µg/mL)'].values
y = standard_df['Absorbance (OD)'].values

# Linear regression
slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)


# Create scatterplot and fit linear regression line
plt.figure(figsize=(8, 6))
sns.regplot(
    data=standard_df,
    x='Collagen Concentration (µg/mL)',
    y='Absorbance (OD)',
    scatter=True,
    ci=None,          
    color='black',
    marker='X',
    line_kws={'color': 'black', 'linestyle': '--'}
)

# Axis limits and ticks (adjust as needed)
plt.xlim(0, 60)
plt.ylim(-0.01, 0.06)
plt.xticks(np.arange(0, 65, 20))

# Force x-axis tick labels to have 1 decimal place
ax = plt.gca()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# Customization
plt.title('Standard Curve - Hydroxyproline Assay')
plt.xlabel('Collagen Concentration (µg/mL)')
plt.ylabel('Absorbance (OD)')
plt.tight_layout()
plt.show()

#%%
"""
Make dataframe for results CDM samples
"""

# Mean OD values CDM samples and blank corrected
OD_CDM = {
   'CDM 10.000 µg/mL' : df.iloc[36:39, 4].mean(),
   'CDM 2000 µg/mL' : df.iloc[36:39, 5].mean(),
   'CDM 800 µg/mL' : df.iloc[36:39, 6].mean(),
   'CDM 500 µg/mL' : df.iloc[36:39, 7].mean(),
   'CDM 300 µg/mL' : df.iloc[36:39, 8].mean(),
   'CDM 200 µg/mL' : df.iloc[36:39, 9].mean(),
   'CDM 100 µg/mL' : df.iloc[36:39, 10].mean(),
   'CDM 50 µg/mL' : df.iloc[36:39, 11].mean(),
   'CDM 10 µg/mL' : df.iloc[36:39, 12].mean(),
}

OD_CDM_means_corrected = {
    label: (values - mean_blankOD).mean()
    for label, values in OD_CDM.items()
}


CDM_df = pd.DataFrame({
    'Absorbance (OD)' : OD_CDM_means_corrected
})


CDM_OD_values = CDM_df['Absorbance (OD)'].values

# Calculate collagen concentrations based on OD values and standard curve
CDM_predicted_conc = (CDM_OD_values - intercept) / slope
CDM_df['Collagen Concentration (µg/mL)'] = CDM_predicted_conc
CDM_df['Type'] = 'CDM'  

#%%
"""
Plot CDM sample data in the plot of the standard curve with seaborn
"""

# Color palette per group
standard_colors = sns.color_palette("husl", 9)

# Palette for CDM samples (e.g., "Set2")
CDM_colors = sns.color_palette("flare", len(CDM_df))[::-1] # one color per CDM sample

# Seaborn theme
sns.set_theme(style="whitegrid", context="talk")

# Create figure
plt.figure(figsize=(14, 8))

# Plot standard curve with regression line
sns.regplot(
    data=standard_df,
    x='Collagen Concentration (µg/mL)',
    y='Absorbance (OD)',
    scatter=True,
    ci=None,
    color="black", 
    marker='X',
    line_kws={'color': 'black', 'linestyle': '--', 'alpha': 1, 'label': 'Standard Curve'},
    scatter_kws={'zorder': 1}
)

# Plot standard points separately with label for legend
plt.scatter(
    standard_df['Collagen Concentration (µg/mL)'],
    standard_df['Absorbance (OD)'],
    color="black",
    marker='X',
    s=100,
    label='Standards',  # <-- This creates the legend entry for standards
    zorder=5
)


# Plot CDM samples one-by-one to get individual legend entries
for i, (idx, row) in enumerate(CDM_df.iterrows()):
    plt.scatter(
        row['Collagen Concentration (µg/mL)'],
        row['Absorbance (OD)'],
        color=CDM_colors[i],  
        marker='o',
        s=150,
        label=idx,
        zorder=10
    )

plt.legend(title='CDM Samples', bbox_to_anchor=(1.05, 1), loc='upper left') 

# Axis limits and ticks
plt.xlim(0, 60)
plt.ylim(-0.01, 0.06)
plt.xticks(np.arange(0, 65, 20))

# Format tick labels to 1 decimal
ax = plt.gca()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# Labels and title
plt.title('Results Hydroxyproline Assay on CDM from FaDus', fontsize=20, weight='bold')
plt.xlabel('Collagen Concentration (µg/mL)', fontsize=20, weight='bold')
plt.ylabel('Absorbance (OD)', fontsize=20, weight='bold')
plt.legend()
plt.tight_layout()

ax = plt.gca()

# Bold tick labels
for label in ax.get_xticklabels():
    label.set_fontweight('bold')
for label in ax.get_yticklabels():
    label.set_fontweight('bold')
    
plt.show()




























