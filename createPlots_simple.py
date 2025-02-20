#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 13:16:38 2023

@author: gkotsoulias
"""

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
mat = scipy.io.loadmat('/---/---/---.mat')

"""
Assign variables and region
"""
#HERE ADJUST DEPENDED ON WHAT IS BEING PLOTTED. QSM, DECOMPOSE, R2*, BPnd
Region = 'Striatum'
Controls  = mat['Controls_BPnd'] 
GTS = mat['GTS_BPnd']


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

"""
Create plot and save figure
"""
my_plot = seaborn.violinplot(data=[GTS[GTS<1000], Controls[Controls<1000]], 
                   bw=0.65,width=0.5, inner= 'box',orient= 'v', palette="Set1", linewidth=1.5, saturation= 0.5)
my_plot.set(ylabel = "BPnd", title =Region)
seaborn.set(rc={'figure.figsize':(5.7,8.27)})
imname = '/---/'+Region+'---_violinPLots.png'
plt.savefig(imname) 


"""
Calculate p-vals and welch p-vals
"two-sided" tests for PET, no hypothesis
"less" for 7T as we hypothesize decreased iron
"""
pval_welch = stats.ttest_ind(GTS[GTS<1000],Controls[Controls<1000],alternative='two-sided', equal_var = False) #alternative='less'  
mean_GTS = GTS[GTS<1000].mean()
std_GTS = GTS[GTS<1000].std()
mean_Controls = Controls[Controls<1000].mean()
std_Controls = Controls[Controls<1000].std()

#sm.stats.fdrcorrection([-,-,-,-,-,-,-], alpha=0.05, method='indep', is_sorted=False)


