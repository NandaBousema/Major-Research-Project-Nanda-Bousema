#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  9 08:55:53 2025

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

"""
Sample I used on this day = Ficoll SD dissolved in PBS+ for quite some time
"""

#%%
"""
Open
"""

file_path = '/Users/nandabousema/Library/CloudStorage/OneDrive-Persoonlijk/Internship_masterRMT/CDM/Protein assays/BCA assay/20250818_BCAassay_AA_DMEM.xlsx'  # Replace with your Excel file's path
#pip install xlrd>=2.0.1
df_AADMEM = pd.read_excel(file_path, header=None)

# Inspect the data
print(df_AADMEM.head())

file_path = '/Users/nandabousema/Library/CloudStorage/OneDrive-Persoonlijk/Internship_masterRMT/CDM/Protein assays/BCA assay/20250818_BCAassay_Ficoll.xlsx'  # Replace with your Excel file's path
#pip install xlrd>=2.0.1
df_Ficoll = pd.read_excel(file_path, header=None)

# Inspect the data
print(df_Ficoll.head())

#%%
"""
AA and DMEM
"""
# Calculate Blank corrected OD data
# Select blank OD data
OD_Blank_cols_AD = df_AADMEM.iloc[25:28,4]
OD_Blank_means_AD = OD_Blank_cols_AD.mean()    

print(OD_Blank_cols_AD)
print(OD_Blank_means_AD)


# Compute the mean across replicates for each row
mean_blankOD_AD = OD_Blank_means_AD

print("Mean Blank OD:", mean_blankOD_AD)


#%%
"""
Make dataframe for standard curve
"""
# Concentrations for the standard curve and calculate means
Con_Standard_cols_AD = df_AADMEM.iloc[14:22, 1:4]    
Con_Standard_cols_AD = Con_Standard_cols_AD.apply(pd.to_numeric, errors="coerce")
Con_Standard_means_AD = Con_Standard_cols_AD.mean(axis=1, skipna=True).reset_index(drop=True)     #Merge standard columns into one column with the means of these standards.


# Standard OD data and calculate means
OD_Standard_cols_AD = df_AADMEM.iloc[25:33, 1:4]
OD_Standard_cols_AD = OD_Standard_cols_AD.apply(pd.to_numeric, errors="coerce")

# Blank Corrected Data
OD_Standard_AD = OD_Standard_cols_AD - mean_blankOD_AD

OD_Standard_means_AD = OD_Standard_AD.mean(axis=1).reset_index(drop=True)     #Merge standard columns into one column with the means of these standards.

standard_df_AD = pd.DataFrame({
    'Protein Concentration (µg/mL)' : Con_Standard_means_AD,
    'Absorbance (OD)' : OD_Standard_means_AD 
    
})

# Dtype from above is object, error='coerce' will turn invalid entries like text of empry strings into NaN
standard_df_AD['Absorbance (OD)'] = pd.to_numeric(standard_df_AD['Absorbance (OD)'], errors='coerce')
standard_df_AD['Protein Concentration (µg/mL)'] = pd.to_numeric(standard_df_AD['Protein Concentration (µg/mL)'], errors='coerce')

# Drop rows with NaNs in either column
standard_df_AD = standard_df_AD.dropna(subset=['Absorbance (OD)', 'Protein Concentration (µg/mL)'])

#%%
"""
Create plot standard curve with seaborn
"""

# Choose palette
colors = sns.color_palette("flare", 6)
my_color = colors[0]  # the first color, the first color is 0, second is 1, etc. so last one is 8
#sns.palplot(colors)

# Set Seaborn theme
sns.set_theme(style="darkgrid", context = "talk")

# Plot standard curve with best-fit line
# Extract values
x = standard_df_AD['Protein Concentration (µg/mL)'].values
y = standard_df_AD['Absorbance (OD)'].values

# Linear regression
slope_AD, intercept_AD, r_value_AD, p_value_AD, std_err_AD = stats.linregress(x, y)


# Extended x range (e.g., up to 150)
x_extended = np.linspace(0, 320, 200)
y_extended = slope_AD * x_extended + intercept_AD

# Create scatterplot and fit linear regression line
plt.figure(figsize=(14, 8))
#sns.regplot(
 #   data=standard_df,
  #  x='Protein Concentration (µg/mL)',
   # y='Absorbance (OD)',
    #scatter=True,
    #ci=None,           # Removes confidence interval
    #color='black',
    #marker='X',
    #line_kws={'color': 'black', 'linestyle': '--'}
#)

# Scatter plot (standard curve points)
sns.scatterplot(
    data=standard_df_AD,
    x='Protein Concentration (µg/mL)',
    y='Absorbance (OD)',
    color='black',
    marker='X',
    s=150,
    label='Standard Points'
)

# Plot extended regression line
plt.plot(
    x_extended,
    y_extended,
    linestyle='--',
    color='black',
    label='Best-fit Line'
)

# Axis limits and ticks (adjust as needed)
plt.xlim(0, 30)
plt.ylim(0, 3.5)
plt.xticks(np.arange(0, 30, 5))

# Force x-axis tick labels to have 1 decimal place
ax = plt.gca()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# Customization
plt.title('Standard Curve - BCA Assay')
plt.xlabel('Protein Concentration (µg/mL)')
plt.ylabel('Absorbance (OD)')
plt.tight_layout()
plt.show()

#%%
"""
AA - Make dataframe for results CDM samples BCA
"""

OD_CDM_AA = {
    '800': pd.to_numeric(df_AADMEM.iloc[30:33, 4], errors='coerce').dropna(),
    '500' : pd.to_numeric(df_AADMEM.iloc[30:33, 5], errors='coerce').dropna(),
    '300' : pd.to_numeric(df_AADMEM.iloc[30:33, 6], errors='coerce').dropna(),
    '200' : pd.to_numeric(df_AADMEM.iloc[30:33, 7], errors='coerce').dropna(),
    '100' : pd.to_numeric(df_AADMEM.iloc[30:33, 8], errors='coerce').dropna(),
    '50' : pd.to_numeric(df_AADMEM.iloc[30:33, 9], errors='coerce').dropna(),
    '25'  : pd.to_numeric(df_AADMEM.iloc[30:33, 10], errors='coerce').dropna(),
    '10'  : pd.to_numeric(df_AADMEM.iloc[30:33, 11], errors='coerce').dropna(),
    '5'   : pd.to_numeric(df_AADMEM.iloc[30:33, 12], errors='coerce').dropna()
}


OD_CDM_corrected_AA = {
    label: (values - mean_blankOD_AD)
    for label, values in OD_CDM_AA.items()
}

Conc_CDM_AA = {
    '800': pd.to_numeric(df_AADMEM.iloc[19:22, 4], errors='coerce').dropna(),
    '500' : pd.to_numeric(df_AADMEM.iloc[19:22, 5], errors='coerce').dropna(),
    '300' : pd.to_numeric(df_AADMEM.iloc[19:22, 6], errors='coerce').dropna(),
    '200' : pd.to_numeric(df_AADMEM.iloc[19:22, 7], errors='coerce').dropna(),
    '100' : pd.to_numeric(df_AADMEM.iloc[19:22, 8], errors='coerce').dropna(),
    '50' : pd.to_numeric(df_AADMEM.iloc[19:22, 9], errors='coerce').dropna(),
    '25'  : pd.to_numeric(df_AADMEM.iloc[19:22, 10], errors='coerce').dropna(),
    '10'  : pd.to_numeric(df_AADMEM.iloc[19:22, 11], errors='coerce').dropna(),
    '5'   : pd.to_numeric(df_AADMEM.iloc[19:22, 12], errors='coerce').dropna()
}


labels = list(Conc_CDM_AA.keys())

CDM_df_AA = pd.DataFrame({
    "CDM Concentration (µg/mL)": labels,
    "Protein Concentration (µg/mL)": [Conc_CDM_AA[l].mean() for l in labels],
    "Absorbance (OD)": [OD_CDM_corrected_AA[l].mean() for l in labels]
})

CDM_df_AA_stats = pd.DataFrame({
    "CDM Concentration (µg/mL)": labels,
    "Protein Concentration (µg/mL)": [np.mean(Conc_CDM_AA[l]) for l in labels],
    "Protein Concentration SD": [np.std(Conc_CDM_AA[l], ddof=1) for l in labels],  # sample std
    "Absorbance (OD)": [np.mean(OD_CDM_corrected_AA[l]) for l in labels],
    "Absorbance (OD) SD": [np.std(OD_CDM_corrected_AA[l], ddof=1) for l in labels]
})

CDM_df_AA_stats["Protein Concentration (µg/mL)"] = pd.to_numeric(CDM_df_AA_stats["Protein Concentration (µg/mL)"], errors='coerce')
CDM_df_AA_stats["Absorbance (OD)"] = pd.to_numeric(CDM_df_AA_stats["Absorbance (OD)"], errors='coerce')

print(CDM_df_AA_stats)


# If you want to check the whole dataframe run the following:
# Show all rows
pd.set_option('display.max_rows', None)

# Show all columns
pd.set_option('display.max_columns', None)

# Show full column width (if long text)
pd.set_option('display.max_colwidth', None)

# Now print the dataframe
print(CDM_df_AA_stats)

#%%
"""
Plot CDM sample data in the plot of the standard curve with seaborn
"""

# Reverse color palette: darkest for highest concentration
CDM_colors = sns.color_palette("flare", len(CDM_df_AA_stats))[::-1]

sns.set_theme(style="whitegrid", context="talk")
plt.figure(figsize=(10, 4))

# Plot standard curve points
sns.scatterplot(
    data=standard_df_AD,
    x='Protein Concentration (µg/mL)',
    y='Absorbance (OD)',
    color='black',
    marker='X',
    s=150
)

# Plot regression line
plt.plot(
    x_extended,
    y_extended,
    linestyle='--',
    color='black',
    label='Standard Curve'
)

# Plot CDM points with error bars and proper colors
for i, row in enumerate(CDM_df_AA_stats.itertuples()):
    plt.errorbar(
        x=row._2,           # x = mean protein concentration
        y=row._4,             # y = mean OD
        yerr=row._5,           # vertical error bar
        fmt='o',                   # marker style
        color=CDM_colors[i],       # color from reversed palette
        markersize=10,
        capsize=5,
        label=f"{int(row._1)} µg/mL"  # label = actual CDM concentration
    )

# Axis limits and ticks
plt.xlim(-2, 22)
plt.ylim(0, 3.5)
plt.xticks(np.arange(-2, 22, 2))
ax = plt.gca()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# Labels and title
plt.title('FaDu AA CDM Total Protein Concentration', fontsize=20, weight='bold')
plt.xlabel('Protein Concentration (µg/mL)', fontsize=20, weight='bold')
plt.ylabel('Absorbance (OD)', fontsize=20, weight='bold')

# Bold tick labels
for label in ax.get_xticklabels():
    label.set_fontweight('bold')
for label in ax.get_yticklabels():
    label.set_fontweight('bold')

plt.show()

#%%
"""
DMEM - Make dataframe for results CDM samples BCA
"""

OD_CDM_DMEM = {
    '800': pd.to_numeric(df_AADMEM.iloc[25:28, 5], errors='coerce').dropna(),
    '500' : pd.to_numeric(df_AADMEM.iloc[25:28, 6], errors='coerce').dropna(),
    '300' : pd.to_numeric(df_AADMEM.iloc[25:28, 7], errors='coerce').dropna(),
    '200' : pd.to_numeric(df_AADMEM.iloc[25:28, 8], errors='coerce').dropna(),
    '100' : pd.to_numeric(df_AADMEM.iloc[25:28, 9], errors='coerce').dropna(),
    '50' : pd.to_numeric(df_AADMEM.iloc[25:28, 10], errors='coerce').dropna(),
    '25'  : pd.to_numeric(df_AADMEM.iloc[25:28, 11], errors='coerce').dropna(),
    '10'  : pd.to_numeric(df_AADMEM.iloc[25:28, 12], errors='coerce').dropna(),
    '5'   : pd.Series(pd.to_numeric(df_AADMEM.iloc[28, 5:8], errors='coerce')).dropna()
}


OD_DMEM_corrected = {
    label: (values - mean_blankOD_AD)
    for label, values in OD_CDM_DMEM.items()
}

Conc_CDM_DMEM = {
    '800': pd.to_numeric(df_AADMEM.iloc[14:17, 5], errors='coerce').dropna(),
    '500' : pd.to_numeric(df_AADMEM.iloc[14:17, 6], errors='coerce').dropna(),
    '300' : pd.to_numeric(df_AADMEM.iloc[14:17, 7], errors='coerce').dropna(),
    '200' : pd.to_numeric(df_AADMEM.iloc[14:17, 8], errors='coerce').dropna(),
    '100' : pd.to_numeric(df_AADMEM.iloc[14:17, 9], errors='coerce').dropna(),
    '50' : pd.to_numeric(df_AADMEM.iloc[14:17, 10], errors='coerce').dropna(),
    '25'  : pd.to_numeric(df_AADMEM.iloc[14:17, 11], errors='coerce').dropna(),
    '10'  : pd.to_numeric(df_AADMEM.iloc[14:17, 12], errors='coerce').dropna(),
    '5'   : pd.Series(pd.to_numeric(df_AADMEM.iloc[17, 5:8], errors='coerce')).dropna()
}

labels = list(Conc_CDM_DMEM.keys())

CDM_df_DMEM = pd.DataFrame({
    "CDM Concentration (µg/mL)": labels,
    "Protein Concentration (µg/mL)": [Conc_CDM_DMEM[l].mean() for l in labels],
    "Absorbance (OD)": [OD_DMEM_corrected[l].mean() for l in labels]
})

CDM_df_DMEM_stats = pd.DataFrame({
    "CDM Concentration (µg/mL)": labels,
    "Protein Concentration (µg/mL)": [np.mean(Conc_CDM_DMEM[l]) for l in labels],
    "Protein Concentration SD": [np.std(Conc_CDM_DMEM[l], ddof=1) for l in labels],  # sample std
    "Absorbance (OD)": [np.mean(OD_DMEM_corrected[l]) for l in labels],
    "Absorbance (OD) SD": [np.std(OD_DMEM_corrected[l], ddof=1) for l in labels]
})

CDM_df_DMEM_stats["Protein Concentration (µg/mL)"] = pd.to_numeric(CDM_df_DMEM_stats["Protein Concentration (µg/mL)"], errors='coerce')
CDM_df_DMEM_stats["Absorbance (OD)"] = pd.to_numeric(CDM_df_DMEM_stats["Absorbance (OD)"], errors='coerce')

print(CDM_df_DMEM_stats)


# If you want to check the whole dataframe run the following:
# Show all rows
pd.set_option('display.max_rows', None)

# Show all columns
pd.set_option('display.max_columns', None)

# Show full column width (if long text)
pd.set_option('display.max_colwidth', None)

# Now print the dataframe
print(CDM_df_DMEM_stats)

#%%
"""
Plot CDM sample data in the plot of the standard curve with seaborn
"""
# Reverse color palette: darkest for highest concentration
CDM_colors = sns.color_palette("flare", len(CDM_df_DMEM_stats))[::-1]

sns.set_theme(style="whitegrid", context="talk")
plt.figure(figsize=(10, 4))

# Plot standard curve points
sns.scatterplot(
    data=standard_df_AD,
    x='Protein Concentration (µg/mL)',
    y='Absorbance (OD)',
    color='black',
    marker='X',
    s=150
)

# Plot regression line
plt.plot(
    x_extended,
    y_extended,
    linestyle='--',
    color='black',
    label='Standard Curve'
)

# Plot CDM points with error bars and proper colors
for i, row in enumerate(CDM_df_DMEM_stats.itertuples()):
    plt.errorbar(
        x=row._2,           # x = mean protein concentration
        y=row._4,             # y = mean OD
        yerr=row._5,           # vertical error bar
        fmt='o',                   # marker style
        color=CDM_colors[i],       # color from reversed palette
        markersize=10,
        capsize=5,
        label=f"{int(row._1)} µg/mL"  # label = actual CDM concentration
    )

# Axis limits and ticks
plt.xlim(-2, 22)
plt.ylim(0, 3.5)
plt.xticks(np.arange(-2, 22, 2))
ax = plt.gca()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# Labels and title
plt.title('FaDu DMEM CDM Total Protein Concentration', fontsize=20, weight='bold')
plt.xlabel('Protein Concentration (µg/mL)', fontsize=20, weight='bold')
plt.ylabel('Absorbance (OD)', fontsize=20, weight='bold')


# Bold tick labels
for label in ax.get_xticklabels():
    label.set_fontweight('bold')
for label in ax.get_yticklabels():
    label.set_fontweight('bold')

plt.show()



#%%
"""
Ficoll 
"""
# Calculate Blank corrected OD data
# Select blank OD data
OD_Blank_cols_F = df_Ficoll.iloc[25:28,4]
OD_Blank_means_F = OD_Blank_cols_F.mean()    

print(OD_Blank_cols_F)
print(OD_Blank_means_F)


# Compute the mean across replicates for each row
mean_blankOD_F = OD_Blank_means_F

print("Mean Blank OD:", mean_blankOD_F)


#%%
"""
Make dataframe for standard curve
"""
# Concentrations for the standard curve and calculate means
Con_Standard_cols_F = df_Ficoll.iloc[14:22, 1:4]    
Con_Standard_cols_F = Con_Standard_cols_F.apply(pd.to_numeric, errors="coerce")
Con_Standard_means_F = Con_Standard_cols_F.mean(axis=1, skipna=True).reset_index(drop=True)     #Merge standard columns into one column with the means of these standards.


# Standard OD data and calculate means
OD_Standard_cols_F = df_Ficoll.iloc[25:33, 1:4]
OD_Standard_cols_F = OD_Standard_cols_F.apply(pd.to_numeric, errors="coerce")

# Blank Corrected Data
OD_Standard_F = OD_Standard_cols_F - mean_blankOD_F

OD_Standard_means_F = OD_Standard_F.mean(axis=1).reset_index(drop=True)     #Merge standard columns into one column with the means of these standards.

standard_df_F = pd.DataFrame({
    'Protein Concentration (µg/mL)' : Con_Standard_means_F,
    'Absorbance (OD)' : OD_Standard_means_F 
    
})

# Dtype from above is object, error='coerce' will turn invalid entries like text of empry strings into NaN
standard_df_F['Absorbance (OD)'] = pd.to_numeric(standard_df_F['Absorbance (OD)'], errors='coerce')
standard_df_F['Protein Concentration (µg/mL)'] = pd.to_numeric(standard_df_F['Protein Concentration (µg/mL)'], errors='coerce')

# Drop rows with NaNs in either column
standard_df_F = standard_df_F.dropna(subset=['Absorbance (OD)', 'Protein Concentration (µg/mL)'])

#%%
"""
Create plot standard curve with seaborn
"""

# Choose palette
colors = sns.color_palette("flare", 6)
my_color = colors[0]  # the first color, the first color is 0, second is 1, etc. so last one is 8
#sns.palplot(colors)

# Set Seaborn theme
sns.set_theme(style="darkgrid", context = "talk")

# Plot standard curve with best-fit line
# Extract values
x = standard_df_F['Protein Concentration (µg/mL)'].values
y = standard_df_F['Absorbance (OD)'].values

# Linear regression
slope_F, intercept_F, r_value_F, p_value_F, std_err_F = stats.linregress(x, y)


# Extended x range (e.g., up to 150)
x_extended = np.linspace(0, 320, 200)
y_extended = slope_F * x_extended + intercept_F

# Create scatterplot and fit linear regression line
plt.figure(figsize=(14, 8))
#sns.regplot(
 #   data=standard_df,
  #  x='Protein Concentration (µg/mL)',
   # y='Absorbance (OD)',
    #scatter=True,
    #ci=None,           # Removes confidence interval
    #color='black',
    #marker='X',
    #line_kws={'color': 'black', 'linestyle': '--'}
#)

# Scatter plot (standard curve points)
sns.scatterplot(
    data=standard_df_F,
    x='Protein Concentration (µg/mL)',
    y='Absorbance (OD)',
    color='black',
    marker='X',
    s=150,
    label='Standard Points'
)

# Plot extended regression line
plt.plot(
    x_extended,
    y_extended,
    linestyle='--',
    color='black',
    label='Best-fit Line'
)

# Axis limits and ticks (adjust as needed)
plt.xlim(0, 30)
plt.ylim(0, 3.5)
plt.xticks(np.arange(0, 30, 5))

# Force x-axis tick labels to have 1 decimal place
ax = plt.gca()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# Customization
plt.title('Standard Curve - BCA Assay')
plt.xlabel('Protein Concentration (µg/mL)')
plt.ylabel('Absorbance (OD)')
plt.tight_layout()
plt.show()



#%%
"""
Ficoll - Make dataframe for results CDM samples BCA
"""

OD_CDM_F = {
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

OD_F_corrected = {
    label: (values - mean_blankOD_F)
    for label, values in OD_CDM_F.items()
}

Conc_CDM_F = {
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

labels = list(Conc_CDM_F.keys())

CDM_df_F = pd.DataFrame({
    "CDM Concentration (µg/mL)": labels,
    "Protein Concentration (µg/mL)": [Conc_CDM_F[l].mean() for l in labels],
    "Absorbance (OD)": [OD_F_corrected[l].mean() for l in labels]
})

CDM_df_F_stats = pd.DataFrame({
    "CDM Concentration (µg/mL)": labels,
    "Protein Concentration (µg/mL)": [np.mean(Conc_CDM_F[l]) for l in labels],
    "Protein Concentration SD": [np.std(Conc_CDM_F[l], ddof=1) for l in labels],  # sample std
    "Absorbance (OD)": [np.mean(OD_F_corrected[l]) for l in labels],
    "Absorbance (OD) SD": [np.std(OD_F_corrected[l], ddof=1) for l in labels]
})

CDM_df_F_stats["Protein Concentration (µg/mL)"] = pd.to_numeric(CDM_df_F_stats["Protein Concentration (µg/mL)"], errors='coerce')
CDM_df_F_stats["Absorbance (OD)"] = pd.to_numeric(CDM_df_F_stats["Absorbance (OD)"], errors='coerce')

print(CDM_df_F_stats)


# If you want to check the whole dataframe run the following:
# Show all rows
pd.set_option('display.max_rows', None)

# Show all columns
pd.set_option('display.max_columns', None)

# Show full column width (if long text)
pd.set_option('display.max_colwidth', None)

# Now print the dataframe
print(CDM_df_F_stats)

#%%
"""
Plot CDM sample data in the plot of the standard curve with seaborn - without legend
"""

# Reverse color palette: darkest for highest concentration
CDM_colors = sns.color_palette("flare", len(CDM_df_F_stats))[::-1]

sns.set_theme(style="whitegrid", context="talk")
plt.figure(figsize=(10, 4))

# Plot standard curve points
sns.scatterplot(
    data=standard_df_F,
    x='Protein Concentration (µg/mL)',
    y='Absorbance (OD)',
    color='black',
    marker='X',
    s=150
)

# Plot regression line
plt.plot(
    x_extended,
    y_extended,
    linestyle='--',
    color='black',
    label='Standard Curve'
)

# Plot CDM points with error bars and proper colors
for i, row in enumerate(CDM_df_F_stats.itertuples()):
    plt.errorbar(
        x=row._2,           # x = mean protein concentration
        y=row._4,             # y = mean OD
        yerr=row._5,           # vertical error bar
        fmt='o',                   # marker style
        color=CDM_colors[i],       # color from reversed palette
        markersize=10,
        capsize=5,
        label=f"{int(row._1)} µg/mL"  # label = actual CDM concentration
    )

# Axis limits and ticks
plt.xlim(-2, 22)
plt.ylim(0, 3.5)
plt.xticks(np.arange(-2, 22, 2))
ax = plt.gca()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# Labels and title
plt.title('FaDu Ficoll CDM Total Protein Concentration', fontsize=20, weight='bold')
plt.xlabel('Protein Concentration (µg/mL)', fontsize=20, weight='bold')
plt.ylabel('Absorbance (OD)', fontsize=20, weight='bold')

# Bold tick labels
for label in ax.get_xticklabels():
    label.set_fontweight('bold')
for label in ax.get_yticklabels():
    label.set_fontweight('bold')

plt.show()

#%%
"""
Plot CDM sample data in the plot of the standard curve with seaborn - with legend
"""

# Reverse color palette: darkest for highest concentration
CDM_colors = sns.color_palette("flare", len(CDM_df_F_stats))[::-1]

sns.set_theme(style="whitegrid", context="talk")
plt.figure(figsize=(10, 7))

# Plot standard curve points
sns.scatterplot(
    data=standard_df_F,
    x='Protein Concentration (µg/mL)',
    y='Absorbance (OD)',
    color='black',
    marker='X',
    s=150,
    label='Standards'
)

# Plot regression line
plt.plot(
    x_extended,
    y_extended,
    linestyle='--',
    color='black',
    label='Standard Curve'
)

# Plot CDM points with error bars and proper colors
for i, row in enumerate(CDM_df_F_stats.itertuples()):
    plt.errorbar(
        x=row._2,           # x = mean protein concentration
        y=row._4,             # y = mean OD
        yerr=row._5,           # vertical error bar
        fmt='o',                   # marker style
        color=CDM_colors[i],       # color from reversed palette
        markersize=10,
        capsize=5,
        label=f"{int(row._1)} µg/mL"  # label = actual CDM concentration
    )

# Axis limits and ticks
plt.xlim(-2, 22)
plt.ylim(0, 3.5)
plt.xticks(np.arange(-2, 22, 2))
ax = plt.gca()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# Labels and title
plt.title('FaDu Ficoll CDM Total Protein Concentration', fontsize=20, weight='bold')
plt.xlabel('Protein Concentration (µg/mL)', fontsize=20, weight='bold')
plt.ylabel('Absorbance (OD)', fontsize=20, weight='bold')

# Legend
leg = plt.legend(title='CDM Concentration (µg/mL)', loc='upper left', prop={'size': 16, 'weight': 'bold'})
leg.get_title().set_fontsize(16)  # title font size
leg.get_title().set_fontweight('bold')  # title bold
plt.tight_layout()

# Bold tick labels
for label in ax.get_xticklabels():
    label.set_fontweight('bold')
for label in ax.get_yticklabels():
    label.set_fontweight('bold')

plt.show()





























