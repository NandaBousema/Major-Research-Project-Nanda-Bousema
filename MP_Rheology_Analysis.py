#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 15:14:27 2024

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
from skimage.feature import peak_local_max
from skimage.segmentation import watershed
from skimage.segmentation import clear_border
from statistics import mean
import seaborn as sns
from scipy.interpolate import CubicSpline


#pip install pandas matplotlib openpyxl > dit is voor in de terminal
#%%
"""
Open
"""

file_path = 'file.xls'  # Replace with your Excel file's path
#pip install xlrd>=2.0.1
df = pd.read_excel(file_path)

# Inspect the data
print(df.head())

#%%
"""
Select data from excel file

Temperature data
"""

# Select temperature data
start_row = 116 - 1
end_row = 281 - 1 
temp_rows = range(start_row, end_row + 1, 5)

temperature_sweep = df.iloc[temp_rows, 0] 

# Rename Temperature Indexes 
new_index = [f"{i}" for i in range(4, 4 + len(temperature_sweep))]  # Creëer nieuwe indexen
temperature_sweep.index = new_index  # Hernoem de index van temperatures

print(temperature_sweep)

"""
Storage modulus data
"""

# Select the storage modulus data
start_row = 119 - 1
end_row = 284 - 1 
storage_modulus_rows = range(start_row, end_row + 1, 5)

storage_modulus = df.iloc[storage_modulus_rows, 0] 
#storage_modulus = [abs(value) + 1e-1 if value <= 0 else value for value in storage_modulus]

print(storage_modulus)

"""
Loss modulus data
"""

# Select loss modulus data
start_row = 119 - 1
end_row = 284 - 1 
loss_modulus_rows = range(start_row, end_row + 1, 5)

loss_modulus = df.iloc[loss_modulus_rows, 1] 
#loss_modulus = [abs(value) + 1e-1 if value <= 0 else value for value in loss_modulus]

print(loss_modulus)




#%%
"""
Plot both loss modulus and storage modulus data in one graph with intersections
"""
sns.color_palette("flare")
full_palette = sns.color_palette("flare", 6)


# Zorg ervoor dat temperatuur, loss en stor numeriek zijn
loss_modulus = np.array(loss_modulus).astype(float)
storage_modulus = np.array(storage_modulus).astype(float)
temperature_sweep.index = np.array (temperature_sweep.index).astype(float)

# Interpolation for finding intersections
loss_modulus_interp = CubicSpline(temperature_sweep.index, loss_modulus, extrapolate=True)
storage_modulus_interp = CubicSpline(temperature_sweep.index, storage_modulus, extrapolate=True)

temperature_sweep_dense = np.linspace(min(temperature_sweep.index), max(temperature_sweep.index), 10000)
y_diff_dense = loss_modulus_interp(temperature_sweep_dense) - storage_modulus_interp(temperature_sweep_dense)

sign_changes = np.where(np.diff(np.sign(y_diff_dense)))[0]
intersect_x = temperature_sweep_dense[sign_changes]
intersect_y = loss_modulus_interp(intersect_x)

# ----- Plot with black background -----
fig, ax = plt.subplots(figsize=(8, 6))

# Black background
fig.patch.set_facecolor('white')
ax.set_facecolor('white')

# Plot lines
ax.plot(temperature_sweep.index, storage_modulus, label="Storage Modulus G'", color=full_palette[0], linewidth=3)
ax.plot(temperature_sweep.index, loss_modulus, label='Loss Modulus G"', color=full_palette[5], linewidth=3)

# Intersections
ax.scatter(intersect_x[0], intersect_y[0], color='red', edgecolor='red', zorder=5, label="Intersection")

# Grid, ticks, labels, legend
ax.grid(True, color='gray')
ax.set_xticks(np.arange(0, 38, 4))

# Axis labels and title
ax.set_title("Rheology 90 kDa Gelatin MPs - 43% Ethanol", fontsize=20, weight='bold', color='black')
ax.set_xlabel('Temperature Sweep (°C)', fontsize=20, weight='bold', color='black')
ax.set_ylabel('Pascal (Pa)', fontsize=20, weight='bold', color='black')

# Tick label styling
ax.tick_params(axis='x', colors='black', labelsize=20)
ax.tick_params(axis='y', colors='black', labelsize=20)
for label in ax.get_xticklabels():
    label.set_fontsize(20)
    label.set_fontweight('bold')
    label.set_color('black')
for label in ax.get_yticklabels():
    label.set_fontsize(20)
    label.set_fontweight('bold')
    label.set_color('black')

# Legend styling
legend = ax.legend(facecolor='white', edgecolor='black', fontsize=20)
for text in legend.get_texts():
    text.set_color('black')

# White spines (borders)
for spine in ax.spines.values():
    spine.set_color('black')

plt.tight_layout()
plt.show()


#%%

"""
Make a graph with log values
"""

# Put the data in 1 graph
plt.figure(figsize=(10, 10))
plt.plot(temperature_sweep.index, storage_modulus, label="Storage Modulus $G'$", linestyle='-')
plt.plot(temperature_sweep.index, loss_modulus, label='Loss Modulus G"', linestyle='--')

plt.xticks(ticks=np.arange(4, 38, 1))  # Ticks van 4 t/m 37 in stappen van 1


# Transform y-axis into log-scale 
plt.yscale('log')

# Limits for the y-axis 
plt.ylim(0.001, 1000)

# Graph details
plt.title('Rheology 16-12-2024 43% ethanol 90kDa MPs - 40mm')
plt.xlabel('Temperature Sweep (°C)')
plt.ylabel('Log Modulus')
plt.legend()
plt.grid(True)
plt.show()

"""
Calculate intersections
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


# Zorg ervoor dat temperatuur, loss en stor numeriek zijn
loss_modulus = np.array(loss_modulus).astype(float)
storage_modulus = np.array(storage_modulus).astype(float)
temperature_sweep.index = np.array (temperature_sweep.index).astype(float)


# Gebruik Cubic Spline voor niet-lineaire interpolatie
loss_modulus_interp = CubicSpline(temperature_sweep.index, loss_modulus, extrapolate=True)
storage_modulus_interp = CubicSpline(temperature_sweep.index, storage_modulus, extrapolate=True)

# Bereken de verschillen
temperature_sweep_dense = np.linspace(min(temperature_sweep.index), max(temperature_sweep.index), 10000)  # Hoge resolutie
y_diff_dense = loss_modulus_interp(temperature_sweep_dense) - storage_modulus_interp(temperature_sweep_dense)

# Zoek de snijpunten
sign_changes = np.where(np.diff(np.sign(y_diff_dense)))[0]
intersect_x = temperature_sweep_dense[sign_changes]
intersect_y = loss_modulus_interp(intersect_x)  # Of stor_interp(intersect_x), beide zijn gelijk bij snijpunten

# Zet de intersectie-waarden om naar logaritmisch
intersect_y_log = np.log10(intersect_y)  # Of np.log(intersect_y) afhankelijk van je voorkeur

# Vermenigvuldig de negatieve logaritmen met -1 om ze positief te maken
#corrected_intersect_y_log = np.where(intersect_y_log < 0, intersect_y_log * -1, intersect_y_log)

#print("Originele logaritmes:", intersect_y_log)
#print("Aangepaste logaritmes:", corrected_intersect_y_log)


# Print de snijpunten
for i in range(len(intersect_x)):
    print(f"Intersect {i+1}: ({intersect_x[i]}, {intersect_y_log[i]})")

# Plot
plt.figure(figsize=(10, 7))
plt.plot(temperature_sweep.index, storage_modulus, label="Stor", color='blue')
plt.plot(temperature_sweep.index, loss_modulus, label='Loss', color='orange')

# Markeer de snijpunten
plt.scatter(intersect_x, intersect_y_log, color='black', zorder=5, label="Intersections")

# Stel x-as limieten in om 4 tot 37 te tonen
plt.xticks(ticks=np.arange(4, 38, 1))

# Logaritmische schaal
plt.yscale('log')
plt.ylim(0.001, 1000)

# Grafiekdetails
plt.title("Rheology 16-12-2024 43% ethanol 90kDa MPs - 40mm")
plt.xlabel('Temperature Sweep (°C)')
plt.ylabel('Log Modulus')
plt.legend()
plt.grid(True)
plt.show()



















