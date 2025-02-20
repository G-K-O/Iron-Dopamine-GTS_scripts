#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 13:42:13 2024

@author: gkotsoulias
"""
from scipy.stats import f_oneway
import os
import numpy as np
import seaborn
from scipy.io import savemat
import nibabel as nib
import matplotlib.pyplot as plt
from scipy import stats
import scipy.io
import statsmodels.api as sm
import pylab as py

"""
Load .mat file with info
"""
mat = scipy.io.loadmat('/----/---/---.mat')


"""
Assign variables and region
"""
Controls  = mat['Controls']
GTS = mat['GTS']
name = 'variable_name'


"""
Remove outliers based on percentiles, if needed.
The outliers are given an impossible value here,
in comparison to the metrics e.g. 1000.
"""
q75,q25 = np.percentile(Controls,[75,25])
intr_qr = q75-q25
max_1 = q75+(1.5*intr_qr)
min_1 = q25-(1.5*intr_qr)
Controls[Controls < min_1] = 1000
Controls[Controls > max_1] = 1000 

q75,q25 = np.percentile(GTS,[75,25])
intr_qr = q75-q25
max_1 = q75+(1.5*intr_qr)
min_1 = q25-(1.5*intr_qr)
GTS[GTS < min_1] = 1000
GTS[GTS > max_1] = 1000 


# Colors
BG_WHITE = "#fbf9f4"
BG_WHITE_2 = "#ffffff"
GREY_LIGHT = "#b4aea9"
GREY50 = "#7F7F7F"
BLUE_DARK = "#1B2838"
BLUE = "#2a475e"
BLACK = "#282724"
GREY_DARK = "#747473"
RED_DARK = "#850e00"

# Colors taken from Dark2 palette in RColorBrewer R library
COLOR_SCALE = ["#850e00", "#3776ab"]

# Horizontal positions for the violins. 
POSITIONS = [0, 1]

# Create jittered version of "x" (which is only 0, 1, and 2)
import scipy.stats as st
jitter = 0.04
# Use always the data that are not outliers
x_data = [np.array([i] * len(d)) for i, d in enumerate([GTS[GTS<1000], Controls[Controls<1000]])]
x_jittered = [x + st.t(df=6, scale=jitter).rvs(len(x)) for x in x_data]

fig, ax = plt.subplots(figsize= (5.7,8.27), dpi = 300)

# Some layout stuff ----------------------------------------------
# Background color
fig.patch.set_facecolor(BG_WHITE_2)
ax.set_facecolor(BG_WHITE_2)

# Horizontal lines that are used as scale reference
for h in HLINES:
    ax.axhline(h, color=GREY50, ls=(0, (5, 5)), alpha=0.8, zorder=0)

# Add violins ----------------------------------------------------
violins = ax.violinplot(
    [GTS[GTS<1000], Controls[Controls<1000]], 
    positions=POSITIONS,
    widths=0.4,
    bw_method="scott",
    showmeans=True, 
    showmedians=False,
    showextrema=False
)

# Customize violins (remove fill, customize line, etc.)
for pc in violins["bodies"]:
    pc.set_facecolor("none")
    pc.set_edgecolor(BLACK)
    pc.set_linewidth(1.4)
    pc.set_alpha(1)
    

# Add boxplots ---------------------------------------------------
medianprops = dict(
    linewidth=0.001, 
    color=GREY_DARK,
    solid_capstyle="butt"
)
boxprops = dict(
    linewidth=2, 
    color=GREY_DARK
)

ax.boxplot(
    [GTS[GTS<1000], Controls[Controls<1000]],
    positions=POSITIONS, 
    showfliers = False, # Do not show the outliers beyond the caps.
    showcaps = False,   # Do not show the caps
    medianprops = medianprops,
    whiskerprops = boxprops,
    boxprops = boxprops
)

# Add jittered dots ----------------------------------------------
for x, y, color in zip(x_jittered, [GTS[GTS<1000], Controls[Controls<1000]], COLOR_SCALE):
    ax.scatter(x, y, s = 100, color=color, alpha=0.7)
    
    
imname = ('/---/---/'+name+'.png')
plt.savefig(imname)   
 
#"""
#Calculate p-vals and welch p-vals
#"two-sided" tests for PET, no hypothesis
#"less" for 7T as we hypothesize decreased iron
#"""
#pval_welch = stats.ttest_ind(GTS[GTS<1000],Controls[Controls<1000],alternative='two-sided', equal_var = False) #alternative='less'  
#mean_GTS = GTS[GTS<1000].mean()
#std_GTS = GTS[GTS<1000].std()
#mean_Controls = Controls[Controls<1000].mean()
#std_Controls = Controls[Controls<1000].std()
#sm.stats.fdrcorrection([-,-,-,-,-,-,-], alpha=0.05, method='indep', is_sorted=False)
   

