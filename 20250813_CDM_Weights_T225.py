#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 13 10:35:30 2025

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
from scipy.stats import ttest_ind
import itertools

#pip install pandas matplotlib openpyxl > dit is voor in de terminal

#%%
"""
Open
"""

file_path = '/Users/nandabousema/Library/CloudStorage/OneDrive-Persoonlijk/Internship_masterRMT/CDM/20250813_Weight CDM.xlsx'  # Replace with your Excel file's path

df = pd.read_excel(file_path, header=None)

#, skiprows=21,
 #   usecols="A:E"

# Inspect the data
print(df.head())

print(df.columns.tolist())
#%%
# Concentrations of each replicate in the CDM range
CDM_weights = df.iloc[21:29, 0:5]    # Rows 1–3, Columns 5–13

# Mean concentrations of the CDM range
T225_Weights_CDM_means = {
   'DMEM SD' : df.iloc[22, 2:5].mean(),
   'DMEM Triton' : df.iloc[23, 2:5].mean(),
   'Ficoll SD' : df.iloc[24, 2:5].mean(),
   'Ficoll Triton' : df.iloc[25, 2:5].mean(),
   'AA SD' : df.iloc[26, 2:5].mean(),
   'AA Triton' : df.iloc[27, 2:5].mean()
}

# If you have std values:
T225_Weights_CDM_stds = {
    'DMEM SD': df.iloc[22, 2:5].std(),
    'DMEM Triton': df.iloc[23, 2:5].std(),
    'Ficoll SD': df.iloc[24, 2:5].std(),
    'Ficoll Triton': df.iloc[25, 2:5].std(),
    'AA SD': df.iloc[26, 2:5].std(),
    'AA Triton': df.iloc[27, 2:5].std()	
}


# Convert your dictionary to a DataFrame for Seaborn
data = pd.DataFrame({
    'Condition': list(T225_Weights_CDM_means.keys()),
    'Mean Weight': list(T225_Weights_CDM_means.values()),
    'Std': list(T225_Weights_CDM_stds.values())
})

# Create a column to differentiate SD vs Triton
data['Type'] = data['Condition'].apply(lambda x: 'SD' if 'SD' in x else 'Triton')

# Create a palette where higher means are darker
norm = plt.Normalize(data['Mean Weight'].min(), data['Mean Weight'].max())
palette = sns.color_palette("flare", as_cmap=True)(norm(data['Mean Weight']))


# Plot bars without built-in CI
plt.figure(figsize=(8, 5))
ax = sns.barplot(
    x='Condition',
    y='Mean Weight',
    data=data,
    palette=palette,
    ci=None,
    edgecolor='black'
)

# Add manual error bars
ax.errorbar(
    x=range(len(data)),
    y=data['Mean Weight'],
    yerr=data['Std'],
    fmt='none',
    ecolor='black',
    capsize=5
)

plt.xticks(rotation=45, ha='right', fontsize=12)
plt.ylabel('Mean Weight', fontsize=14)
plt.title('Mean Weights of T225 Conditions', fontsize=16, weight='bold')
plt.tight_layout()
plt.show()


#%%
# --- 1. Replicate data ---
replicates = {
    'DMEM SD': df.iloc[22, 2:5].values,
    'DMEM Triton': df.iloc[23, 2:5].values,
    'Ficoll SD': df.iloc[24, 2:5].values,
    'Ficoll Triton': df.iloc[25, 2:5].values,
    'AA SD': df.iloc[26, 2:5].values,
    'AA Triton': df.iloc[27, 2:5].values
}

replicates_clean = {}
for k, v in replicates.items():
    arr = pd.to_numeric(v, errors='coerce')  # convert to numeric, NaN if not possible
    arr = np.array(arr, dtype=float)         # ensure NumPy array
    arr = arr[~np.isnan(arr)]                # drop NaNs
    replicates_clean[k] = arr

# --- 2. Means and stds ---
means = {k: v.mean() for k, v in replicates_clean.items()}
stds  = {k: v.std()  for k, v in replicates_clean.items()}

data = pd.DataFrame({
    'Condition': list(means.keys()),
    'Mean Weight': list(means.values()),
    'Std': list(stds.values())
})

# Create a palette where higher means are darker
norm = plt.Normalize(data['Mean Weight'].min(), data['Mean Weight'].max())
palette = sns.color_palette("flare", as_cmap=True)(norm(data['Mean Weight']))

# --- 3. Plot bars ---
plt.figure(figsize=(12,8))
ax = sns.barplot(
    x='Condition',
    y='Mean Weight',
    data=data,
    palette=palette,
    ci=None,
    edgecolor='black'
)

# Error bars
ax.errorbar(
    x=range(len(data)),
    y=data['Mean Weight'],
    yerr=data['Std'],
    fmt='none',
    ecolor='black',
    capsize=5
)

# --- 4. Pairwise t-tests and significance ---
conditions = list(replicates_clean.keys())
pairs = list(itertools.combinations(conditions, 2))

def p_to_star(p):
    if p < 0.001: return '***'
    elif p < 0.01: return '**'
    elif p < 0.05: return '*'
    else: return 'ns'

# Calculate max bar height
bar_max = (data['Mean Weight'] + data['Std']).max()

# Prepare y-limit headroom (will be updated dynamically)
ax.set_ylim(0, bar_max + 1)

# Track existing significance lines to prevent overlap
drawn_lines = []  # list of (xmin, xmax, y_top)

p_values = []

for g1, g2 in pairs:
    i = conditions.index(g1)
    j = conditions.index(g2)

    # Base height above the tallest of the two bars
    y1 = data.loc[data['Condition'] == g1, 'Mean Weight'].values[0] + \
         data.loc[data['Condition'] == g1, 'Std'].values[0]
    y2 = data.loc[data['Condition'] == g2, 'Mean Weight'].values[0] + \
         data.loc[data['Condition'] == g2, 'Std'].values[0]
    y = max(y1, y2) + 0.02  # start just above error bar

    # Adjust height so it doesn't overlap with existing lines
    this_min = min(i, j)
    this_max = max(i, j)
    for (xmin, xmax, y_top) in drawn_lines:
        if not (this_max < xmin or this_min > xmax):  # horizontal overlap
            y = max(y, y_top + 0.05)

    # t-test
    stat, pval = ttest_ind(replicates_clean[g1], replicates_clean[g2])
    star = p_to_star(pval)

    # Save results
    p_values.append({
        'Group 1': g1,
        'Group 2': g2,
        't-statistic': stat,
        'p-value': pval,
        'significance': star
    })

    # Plot only significant comparisons
    if star != 'ns':
        ax.plot([i, i, j, j], [y, y + 0.01, y + 0.01, y], color='black', linewidth=2)
        ax.text((i + j) / 2, y + 0.015, star, ha='center', va='bottom', fontsize=20)

        # Save this line's occupied space
        drawn_lines.append((this_min, this_max, y + 4.2))

# Make a DataFrame of results
p_df = pd.DataFrame(p_values)
print(p_df)

# Update ylim so no lines are cut off
max_y_line = max([y for (_, _, y) in drawn_lines], default=bar_max)
ax.set_ylim(0, max_y_line + 0.1)

# --- 5. Final touches ---
plt.xticks(rotation=45, ha='right', fontsize=20, weight='bold')
plt.yticks(fontsize=20, weight='bold') 
ax.set_xlabel('Condition', fontsize=20, weight='bold')
plt.ylabel('CDM Weight (mg)', fontsize=20, weight='bold')
plt.title('CDM Weights Across Culture Conditions and Decellularization Methods', fontsize=20, weight='bold')
plt.tight_layout()
plt.show()

#%%
#Levenes test heb ik gedaan om te kijken of de standaarddeviaties van de groepen zelf niet te groot zijn. 
from scipy.stats import levene

# --- 1. Replicate data (zoals voorheen) ---
replicates = {
    'DMEM SD': df.iloc[22, 2:5].values,
    'DMEM Triton': df.iloc[23, 2:5].values,
    'Ficoll SD': df.iloc[24, 2:5].values,
    'Ficoll Triton': df.iloc[25, 2:5].values,
    'AA SD': df.iloc[26, 2:5].values,
    'AA Triton': df.iloc[27, 2:5].values
}

replicates_clean = {}
for k, v in replicates.items():
    arr = pd.to_numeric(v, errors='coerce')  # convert to numeric, NaN if not possible
    arr = np.array(arr, dtype=float)         
    arr = arr[~np.isnan(arr)]                
    replicates_clean[k] = arr

# --- 2. Check SD homogeneity per condition ---
homogeneity_results = []

for cond, values in replicates_clean.items():
    # Voor 3 replicaten is levene nog nuttig, maar we kunnen ook de max/min vergelijken
    # Hier doen we Levene test tussen de 3 individuele waarden tegen elkaar, dummy-test
    # Omdat er maar 3 waarden zijn, geeft dit soms p=1.0, maar je krijgt een indicatie
    # Alternatief: gewoon CV berekenen
    mean_val = np.mean(values)
    std_val = np.std(values, ddof=1)
    cv = std_val / mean_val if mean_val != 0 else np.nan
    homogeneity_results.append({
        'Condition': cond,
        'Mean': mean_val,
        'Std': std_val,
        'CV': cv
    })

homogeneity_df = pd.DataFrame(homogeneity_results)
print("Replicate variability per condition:")
print(homogeneity_df)

# --- 3. Optional: check if variances differ across conditions ---
# Levene test across all conditions
all_values = [v for v in replicates_clean.values()]
stat, pval = levene(*all_values)
print(f"\nLevene test across all conditions: stat={stat:.3f}, p-value={pval:.3f}")
if pval < 0.05:
    print("=> Significant difference in variance across conditions")
else:
    print("=> No significant difference in variance across conditions")


