#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  5 19:41:52 2025

@author: nandabousema
"""

import os
import cv2
import matplotlib.pyplot as plt
import math
from skimage.measure import label, regionprops
from skimage.color import label2rgb
import numpy as np
from scipy import ndimage as ndi
from skimage.feature import peak_local_max
from skimage.segmentation import watershed
from skimage.segmentation import clear_border
from statistics import mean
import pandas as pd
import seaborn as sns


#%%
"""
Open

"""

list_folders = ["folder where data is stored"]



def open_files(folder):
    image_files = []
    subtitles = []
    for file in os.listdir(folder):
        if file.endswith('.jpg') or file.endswith('.tif'):
           image_files.append(os.path.join(folder, file))
           subtitles.append(os.path.splitext(os.path.basename(file))[0])
    return image_files, subtitles

image_files_20250605, subtitles_20250605 = open_files(list_folders[0])




"""
Plot
"""

# Read the image
#img = cv2.imread(image_files_25012024[0], cv2.IMREAD_GRAYSCALE)

# Display the image
#plt.imshow(img, cmap='gray')
#plt.show()

#This only shows one image and not all of them


"""
Plot all from a file
"""

def plot_images_grid(image_files, title):
    # Calculate the grid size: square root of the number of images
    grid_size = math.ceil(math.sqrt(len(image_files)))
    
    # Create a new figure with a title
    fig = plt.figure(figsize=(20, 20))
    fig.suptitle(title, fontsize=16)
    
    # Add each image to the plot
    for i in range(len(image_files)):
        # Read the image
        img = cv2.imread(image_files[i], cv2.IMREAD_GRAYSCALE)
        image_title = os.path.splitext(os.path.basename(image_files[i]))[0]
        ax = fig.add_subplot(grid_size, grid_size, i + 1)
        ax.imshow(img, cmap='gray')
        ax.axis('on')  # Show axes
        ax.set_title(image_title, fontsize=10)
        

    
    plt.tight_layout()
    plt.show()
    

plot_images_grid(image_files_20250605, 'MP sizes Type A gelatin')

#%%
"""
Binarize
"""


def binarize_images(image_files):
    binarized_images = []
    for image_file in image_files:
        
        # Read the image
        img = cv2.imread(image_file, cv2.IMREAD_GRAYSCALE)
        
        # Perform Gaussian blurring
        img = cv2.GaussianBlur(img, (15, 15), 0)

        # Apply binary thresholding
        ret, binarized = cv2.threshold(img, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
        binarized = cv2.bitwise_not(binarized)
        binarized_images.append(binarized) # invert the binarization, to get the MPs instead of the background
    return binarized_images


bin_20250605 = binarize_images(image_files_20250605)

#In the code above the microparticles are white (225) and the background is black (0). so the image converted into binary (0 or 225). 


"""
Plot
"""

def plot_bin_grid(images, title, subtitles):
    plt.style.use('default')
    # Calculate the grid size: square root of the number of images
    grid_size = math.ceil(math.sqrt(len(images)))
    
    # Create a new figure with a title
    fig = plt.figure(figsize=(20, 20))
    fig.suptitle(title, fontsize=16)
    
    # Add each image to the plot
    for i, image in enumerate(images):
        ax = fig.add_subplot(grid_size, grid_size, i + 1)
        ax.imshow(image)
        ax.axis('on')  # Show axes
        ax.set_title(subtitles[i], fontsize=10)
    
    plt.tight_layout()
    plt.show()

plot_bin_grid(bin_20250605, 'Masks', subtitles_20250605)



#%%
"""
Label
"""

#'code' to exclude the watershed function:
def label_unique_MP(binarized_files):
    labels_list = []
    #labels_overlay = []
    for file in binarized_files:
        # Distance transform
        #distance = ndi.distance_transform_edt(file)

        # Perform the watershed
        #coords = peak_local_max(distance, min_distance=20, labels = file)
        #mask = np.zeros(distance.shape, dtype=bool)
        #mask[tuple(coords.T)] = True
        #markers, _ = ndi.label(mask)
        #labels = watershed(-distance, markers, mask=file)
        labels = label(file)
        
        # Remove labels on the border
        labels = clear_border(labels)

        # Visualize the results
        #image_label_overlay = label2rgb(labels, image=file)
        
        labels_list.append(labels)
        #labels_overlay.append(image_label_overlay)
    
    return labels_list



labelled_20250605 = label_unique_MP(bin_20250605)


#Here a list of label images where each unique microparticle gets a distinct label (1,2,3,...) like assigining a number to each particle. 

"""
Plot
"""


plot_bin_grid(labelled_20250605, 'Labels', subtitles_20250605)



#%%
"""
Measure minor and major axes
NB: 20X with Olympus microscope of the confocal room, has a scale of 4424 pixels / mm = 4.424 pixels / µm
"""


def find_axes(labelled_images, image_files, scale):

    file_nb = 0
    results = []
    
    # ANALYSE
    for labelled in labelled_images:
        
        major_axes = []
        minor_axes = []
        
        for region in regionprops(labelled):
            # Append the length of major and minor axes of each particle to the lists
            major_axes.append(region.major_axis_length/scale)
            minor_axes.append(region.minor_axis_length/scale)
            
    
        # ASSOCIATE TO ETHANOL
        file_name = image_files[file_nb]
        # List of numbers to check for the percentage of ethanol
        numbers = ["30.0", "31.0", "36.0", "46.0", "45.5", "45.0", "44.5", "44.0", "43.5", "43.0", "42.5", "42.0", "41.5", "41.0"]

        for number in numbers:
             if number in file_name:
                 ethanol = number
                 break  # Exit the loop once a number is found
        file_nb += 1
    
        # FINAL RESULT
        results.append([ethanol, major_axes, minor_axes])
    
    
    return results

# scale = 4424 pixels / mm ie scale = 4.424 pixels / µm

results_20250605 = find_axes(labelled_20250605, image_files_20250605, scale = 4.424)


#%%
"""
Plot
"""

# 1. Create dataframes

def dataframes(results):
    # df1 for single values
    df1 = pd.DataFrame({
        'ethanol': [ethanol for ethanol, major_axes, minor_axes, in results for _ in range(len(major_axes))],
        'major_axes': [value for _, major_axes, _,  in results for value in major_axes],
        'minor_axes': [value for _, _, minor_axes, in results for value in minor_axes]
    })

    major_axis_av = df1.groupby('ethanol')['major_axes'].mean().reset_index()
    minor_axis_av = df1.groupby('ethanol')['minor_axes'].mean().reset_index()
    
    max_major_axis = df1.groupby('ethanol')['major_axes'].max().reset_index()
    max_minor_axis = df1.groupby('ethanol')['minor_axes'].max().reset_index()

    # df2 for averages
    df2 = pd.DataFrame({
        'ethanol': major_axis_av['ethanol'],
        'major_axis_av': major_axis_av['major_axes'],
        'minor_axis_av': minor_axis_av['minor_axes'],
        'max_major_axis': max_major_axis['major_axes'],
        'max_minor_axis': max_minor_axis['minor_axes']
    })

    return df1, df2


df1_20250605, df2_20250605 = dataframes(results_20250605)



# 2. Plot dataframes

# The codes behind the # for df1 and df2 are for when you use multiple dataframes. So, e.g. for when you want to plot multiple MP ranges in one plot. 
#df1 = [df1_25012024, df1_08022024, df1_12032024, df1_13032024]
#df1_all = pd.concat(df1)
df1_all = df1_20250605

#df2 = [df2_25012024, df2_08022024, df2_12032024, df2_13032024]
#df2_all = pd.concat(df2, ignore_index=True)
df2_all = df2_20250605 

# Ensuring the ethanol column is numeric
df1_all["ethanol"] = df1_all["ethanol"].astype(float)
df2_all["ethanol"] = df2_all["ethanol"].astype(float)

# Defining the order of ethanol categories
ethanol_order = [30.0, 31.0, 36.0, 43.0, 44.0, 45.0, 45.5, 46.0]

# Converting ethanol column to categorical with order
df1_all['ethanol'] = pd.Categorical(df1_all['ethanol'], categories=ethanol_order, ordered=True)
df2_all['ethanol'] = pd.Categorical(df2_all['ethanol'], categories=ethanol_order, ordered=True)


plt.style.use('ggplot')
fig, axis = plt.subplots(figsize=(12, 6)) # width=12 inches, height=6 inches


sns.violinplot(data=df1_all, x="ethanol", y="major_axes", color = 'lightgray', ax=axis, order=ethanol_order)
sns.stripplot(x='ethanol', y='major_axes', data=df1_all, palette='gist_rainbow', ax=axis, order=ethanol_order)
sns.stripplot(data=df2_all, x='ethanol', y='major_axis_av',  color='red', ax=axis, order=ethanol_order)

# Set the y-axis limits
axis.set_ylim([df1_all['major_axes'].min() - 30, df1_all['major_axes'].max() + 30])


# Change the size and make the axes titles bold
axis.set_xlabel('% Ethanol', fontsize=16, weight='bold') 
axis.set_ylabel('Major Axes (µm)', fontsize=16, weight='bold')  
plt.title("MP Sizes Type A Gelatin with Different Ethanol Percentages")


# Customizing x-axis tick positions and labels 
axis.set_xticks(range(len(ethanol_order)))
axis.set_xticklabels(ethanol_order, fontsize=12, weight='bold')
plt.show()

#%%
# 2. Plot dataframes


# Set base Seaborn theme (will be overridden by manual color settings)
sns.set_theme(style="whitegrid", context="talk")

# Create the full HUSL palette and select specific colors
full_palette = sns.color_palette("husl", 10)
selected_colors = [full_palette[i] for i in [8, 4, 7]]

# The codes behind the # for df1 and df2 are for when you use multiple dataframes. So, e.g. for when you want to plot multiple MP ranges in one plot. 
#df1 = [df1_25012024, df1_08022024, df1_12032024, df1_13032024]
#df1_all = pd.concat(df1)
df1_all = df1_20250605

#df2 = [df2_25012024, df2_08022024, df2_12032024, df2_13032024]
#df2_all = pd.concat(df2, ignore_index=True)
df2_all = df2_20250605 

# Ensuring the ethanol column is numeric
df1_all["ethanol"] = df1_all["ethanol"].astype(float)
df2_all["ethanol"] = df2_all["ethanol"].astype(float)

ethanol_order = [30.0, 31.0, 36.0, 43.0, 44.0, 45.0, 45.5, 46.0]

df1_all['ethanol'] = pd.Categorical(df1_all['ethanol'], categories=ethanol_order, ordered=True)
df2_all['ethanol'] = pd.Categorical(df2_all['ethanol'], categories=ethanol_order, ordered=True)

# Create the figure and axis
fig, axis = plt.subplots(figsize=(12, 6))

# Set black background for figure and axes
fig.patch.set_facecolor('white')
axis.set_facecolor('white')

# Plot data
sns.violinplot(
    data=df1_all,
    x="ethanol",
    y="major_axes",
    color='lightgray',
    ax=axis,
    order=ethanol_order
)
sns.stripplot(
    x='ethanol',
    y='major_axes',
    data=df1_all,
    palette=selected_colors,
    ax=axis,
    order=ethanol_order
)
sns.stripplot(
    data=df2_all,
    x='ethanol',
    y='major_axis_av',
    color='red',
    ax=axis,
    order=ethanol_order
)

# Adjust y-axis limits
axis.set_ylim([
    df1_all['major_axes'].min() - 30,
    df1_all['major_axes'].max() + 30
])

# Set grid color for visibility on black
#axis.grid(True, color='gray')

# Format y-axis ticks to show one decimal place
axis.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f'{y:.1f}'))

# Axis labels and title
axis.set_xlabel('% Ethanol', fontsize=20, weight='bold', color='black')
axis.set_ylabel('Major Axes (µm)', fontsize=20, weight='bold', color='black')
plt.title("40 kDa Gelatin - MP Size vs % Ethanol", fontsize=20, weight='bold', color='black')

# X-axis tick labels
axis.set_xticks(range(len(ethanol_order)))
axis.set_xticklabels(ethanol_order, fontsize=20, weight='bold', color='black')

# Bold white tick labels
for label in axis.get_xticklabels():
    label.set_fontweight('bold')
    label.set_color('black')

for label in axis.get_yticklabels():
    label.set_fontsize(20)
    label.set_fontweight('bold')
    label.set_color('black')

# Set white spines (borders)
for spine in axis.spines.values():
    spine.set_color('black')

plt.show()


#%% 
"""
Statistics - test simple ANOVA 1 - not valid if statistical assumptions (independence, normality, homoscedasticity) not verified
Significance level < 0.05
"""

import scipy.stats as stats

# Compute the p-value
ethanol_groups = [group["major_axes"].values for name, group in df1_all.groupby("ethanol")]
f_val, p_val = stats.f_oneway(*ethanol_groups)

print(f"ANOVA results: F={f_val}, p={p_val}")


#%%
"""
Statistics - whole
Significance level < 0.05
"""

from scipy.stats import norm as normaldensityfunction
from scipy.stats import f_oneway, ttest_ind, levene, kstest, tukey_hsd, mannwhitneyu, kruskal
import scikit_posthocs as sp




# ================= Check the statistical assumptions ==================  
# The statistical assumptions are independence, normality, homoscedasticity.

def check_assumptions(data_for_statistics):
    # INDEPENDENCE
    text0 = "The independence has to be checked by the user"
    
    # NORMALITY
    normality = []
    for x in data_for_statistics: # we compute the normality on each group studied
        res = kstest(x, normaldensityfunction.cdf)
        normality.append(res[1]) # extract the p-value
    check_normality = all(x > 0.05 for x in normality)
    if check_normality is True:
        text1 = "Kolmogorov-Smirnov test indicates that the Normality is respected for all the samples" 
    else:
        text1 = "Kolmogorov-Smirnov test indicates that the Normality is NOT respected for all the samples"
        
    # HOMOSCEDASTICITY
    res, p_var = levene(*data_for_statistics)# usually, samples are small. We thus use Levene's test. Note the use of * in the call that unpacks the list, which is the same as calling with multiple argument
    if p_var > 0.05:
        text2 = "Levene's test indicates equal variances between populations"
        check_homoscedasticity = True
    else:
        text2 = "Levene's test indicates that the variances DIFFER between populations"
        check_homoscedasticity = False
    
    return text0, text1, text2, check_normality, check_homoscedasticity






# ================= STATISTICS computing ==================  

def statistics(data_for_statistics, check_normality, check_homoscedasticity):
    
    labels_stat = []
    # ======== 2 groups =========
    
    if len(data_for_statistics) == 2: # t-test between 2 independent groups
        x_list = [0,1] # list of the columns on which add the stats annotations
        
        # PARAMETRIC test
        if check_normality is True and check_homoscedasticity is True:
            statistics_name = "Student's t-test"
            statistic, p_value = ttest_ind(data_for_statistics[0], data_for_statistics[1])
            if 0.01 < p_value < 0.05:
                label_stat = "*"
            elif 0.001 < p_value < 0.01:
                label_stat = "**"
            elif p_value < 0.001:
                label_stat = "***"
            else:
                label_stat = "ns (p>0.05)"

                
        # NON PARAMETRIC test
        else:
            statistics_name = "Mann-Whitney U test"
            statistic, p_value = mannwhitneyu(data_for_statistics[0], data_for_statistics[1]) # to cross-check the statistic: https://planetcalc.com/7858/
            if 0.01 < p_value < 0.05:
                label_stat = "*"
            elif 0.001 < p_value < 0.01:
                label_stat = "**"
            elif p_value < 0.001:
                label_stat = "***"
            else:
                label_stat = "ns (p>0.05)"
                
        post_hoc_name = "N/A"
        labels_stat.append(label_stat)
        
        
    # ========= More than 2 groups =========  
       
    elif len(data_for_statistics) > 2:
           
        # PARAMETRIC test
           if check_normality is True and check_homoscedasticity is True:
               statistics_name = "One-way ANOVA"
               statistic, pvalue = f_oneway(*data_for_statistics)
               label_stat = []
               x_list = []
               p_value = []
               
               # post-hoc test
               if pvalue < 0.05: 
                   tukey = tukey_hsd(*data_for_statistics)
                   post_hoc_name = "Tukey"
                   for i in range (len(tukey.pvalue)):
                       p_list_i = tukey.pvalue[i] # tukey needs normality and homoscedasticity
                       # In case of unequal sample sizes, the tukey test uses the Tukey-Kramer method
                       for j in range (len(p_list_i)):
                           pval = p_list_i[j] # p_value for the comparison between groups i and j.
                           p_value.append(f"{pval:.1E}")
                           if 0.01 < pval < 0.05:
                               label_stat = "*"

                           elif 0.001 < pval < 0.01:
                               label_stat = "*"

                           elif pval < 0.001:
                               label_stat = "*"
 
                           else:
                               label_stat = "ns"
                           labels_stat.append(label_stat)
       
           # NON PARAMETRIC test
           else:
               statistics_name = "Kruskal-Wallis"
               statistic, pvalue = kruskal(*data_for_statistics) # pvalue from the kruskal-wallis test
               label_stat = []
               p_value = [] # list of all the p_values from the dunn test
               k = 0
               if pvalue < 0.05: # post-hoc test
                   dunn = sp.posthoc_dunn(data_for_statistics, p_adjust = 'bonferroni') # as Kruskal-Wallis test above, we use a multiple comparison procedure on ranks: Dunn's test. To avoid multiple-testing bias, we perform a bonferroni correction.
                   post_hoc_name = "Dunn"
                   for i in range (len(dunn)):
                       for j in range (i, len(dunn)): # as repeats in the dataframe, following the diagonal
                           pval = dunn.iat[i,j] # p_value for the comparison between groups i and j; dunn dataframe, ie use .iat
                           p_value.append(f"{pval:.1E}")
                           if 0.01 < pval < 0.05:
                               label_stat = "*"

                           elif 0.001 < pval < 0.01:
                               label_stat = "**"
                               #y, h, col = data_parameter['Data'].max() + 30 + k, 10, 'k' # determine the height of the annotations

                           elif pval < 0.001:
                               label_stat = "***"

                           else:
                               label_stat = "ns"
                           labels_stat.append(label_stat)
                           
               else:
                   label_stat = "ns (p>0.05)"
                   p_value = "ns (p>0.05)"
                   post_hoc_name = "N/A"
       
    # Number of samples, for indication under the plot
    number_samples = []
    for i in data_for_statistics:
        n = len(i)
        number_samples.append(n)
        
    text3 = "The statistics are calculated via " + str(statistics_name) + " and the post-hoc test via " + str(post_hoc_name)
        
    return statistics_name, text3, p_value, post_hoc_name, labels_stat, number_samples

# Prepare data for statistics
data_for_statistics = [group["major_axes"].values for name, group in df1_all.groupby("ethanol")] # groups by the ethanol column, and gives a list of numpy arrays, where each array contains the "major_axes" values for a specific "ethanol" group
# Goes from lowest ethanol concentration to highest

# Check assumptions
text0, text1, text2, check_normality, check_homoscedasticity = check_assumptions(data_for_statistics)
print(text0)
print(text1)
print(text2)

# Compute statistics
statistics_name, text3, p_value, post_hoc_name, labels_stat, number_samples = statistics(data_for_statistics, check_normality, check_homoscedasticity)
print(text3)
print(p_value)
print(labels_stat),
print(number_samples)




#%%
"""
Plot with STATISTICS

1. Preparation
"""


def get_max_major_axis(df):
    # Group by 'ethanol' and get the max of 'max_major_axis' for each group
    grouped = df.groupby('ethanol')['max_major_axis'].max()
    
    # Convert the Series to a DataFrame and reset the index
    maxima = pd.DataFrame(grouped).reset_index()
    
    return maxima

maxima = get_max_major_axis(df2_all)


df2_all_major_axis_av = df2_all.groupby('ethanol')['major_axis_av'].mean().reset_index()
#df2_all_major_axis_av['ethanol'] = df2_all_major_axis_av['ethanol'] - df1_all['ethanol'].min() # for plotting the scatterplot, as sns does with indexes ie 1 for group 1 of ethanol, 2 for group 2 of ethanol, etc



# RELIRE FROM HERE

def display_statistics(df1_all, maxima, p_value, labels_stat):

    ethanol_groups = np.sort(df1_all['ethanol'].unique()) # to put it from lowest concentration to highest, same order as before in data_for_statistics
    maxima_groups = np.array(maxima['max_major_axis'])

    # Initialize lists to store the results
    mean_ethanol = []
    max_ethanol = []
    p_values = []
    labels = []
    groups = []

    # Initialize a counter for the current index in the p_value and labels_stat lists
    counter = 0

    # Loop over each pair of groups
    for i in range(len(ethanol_groups)):
        for j in range(i, len(ethanol_groups)):
            # Calculate the mean ethanol concentration between the two groups
            mean_ethanol.append(np.mean([ethanol_groups[i], ethanol_groups[j]]))
            max_ethanol.append(np.max([maxima_groups[i], maxima_groups[j]])) # for the y of the line plotted between 2 groups
            
            # Get the p-value and label for this pair of groups
            p_values.append(p_value[counter])
            labels.append(labels_stat[counter])
            groups.append([ethanol_groups[i], ethanol_groups[j]])
            
            # Increment the counter
            counter += 1

    # Create the new DataFrame
    df = pd.DataFrame({
        'Mean Ethanol Concentration': mean_ethanol,
        'p-value': p_values,
        'Statistical Label': labels,
        'Group': groups,
        'Max y': max_ethanol
    })
    
    # Filter the DataFrame to keep only the rows where the p-value is not '1.0E+00'
    df_selected = df[df['p-value'] != '1.0E+00']

    return df, df_selected


df, df_selected = display_statistics(df1_all, maxima, p_value, labels_stat)





"""
Plot with STATISTICS

2. Plotting
"""


def plot_ethanol_concentration(df1_all, df2_all_major_axis_av, df_selected, my_color, my_label_size, style, title, text):
    y_start = df_selected['Max y'].max()

    plt.style.use(style)
    fig, axis = plt.subplots()
    #axis.set_facecolor(background)

    x_count = 0
    count = 10
    y_line = y_start - 10
    for i, row in df_selected.iterrows():
        x_start = row['Group'][0] - df1_all['ethanol'].min()
        x_end = row['Group'][1] - df1_all['ethanol'].min()
        x_count += 1
        y_line += count

        plt.plot([x_start, x_end], [y_line+count, y_line+count], lw=1.5, c= my_color)

        x = row['Mean Ethanol Concentration'] - df1_all['ethanol'].min()
        y = y_line + count + 1
        statistical_label = row['Statistical Label']
        plt.text(x, y, statistical_label, ha='center', va='bottom', color = my_color, fontsize=12, fontweight='bold')

    sns.violinplot(data=df1_all, x="ethanol", y="major_axes", color = 'lavender', ax=axis)
    sns.stripplot(x='ethanol', y='major_axes', data=df1_all, palette='gist_rainbow', ax=axis)
    sns.scatterplot(data=df2_all_major_axis_av, x='ethanol', y='major_axis_av',  color='crimson', s=70, ax=axis, zorder = 3)

    # Cosmetics
    axis.set_xlabel('% Ethanol concentration', weight='bold', size = my_label_size, labelpad = 10)
    axis.set_ylabel('Major axis (µm)', weight='bold', size = my_label_size, labelpad = 10)
    [axis.spines[side].set_linewidth(2) for side in ['left','right','top','bottom']]
    axis.tick_params(axis='both', which='major', labelsize=my_label_size, width=2)

    plt.ylim(-20, y_line + count + 20)
    
    #axis.text(3, 1, 
             #text,
             #size=12, ha="center", transform=axis.transAxes)
    plt.title(title, fontsize=my_label_size + 5, color=my_color, pad = 40, weight='bold')

    plt.show()



plot_ethanol_concentration(df1_all, df2_all_major_axis_av, df_selected, my_color = "black", my_label_size = 20, style = 'default', title = "Size of the MPs vs Ethanol Concentration", text = text3)


#%%
"""
Plot with Statistics 
In dark background
"""

plot_ethanol_concentration(df1_all, df2_all_major_axis_av, df_selected, my_color = "white", my_label_size = 20, style = 'dark_background', title = "Size of the MPs vs Ethanol Concentration", text = text3)







