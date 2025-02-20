#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 17:05:12 2024

@author: gkotsoulias
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 14:59:11 2023

@author: gkotsoulias
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 12:04:38 2023

@author: gkotsoulias
"""

from scipy.io import savemat
import mat73
import scipy.io
import ants
import seaborn as sb
from matplotlib import pyplot as plt
import pandas as pd
import mat73
import scipy.stats as scs
import numpy as np


Patients_7T_codes =  [4,6,8,23,12,11,14,7,15,16,18,20,21,17,19,22,25] #Mp013, excluded movement, Mp009 no PET


PUTS_7T = [33,26,20,22,22,24,28,27,1000,1000,18,23,28,1000,26,27,1000,21,13,23,28,27,37,28,33]
ATQ_MOTOR_FREQ_7T = np.round([1.357142857, 0.615384615, 1.214285714,0.928571429, 2.846153846, 0.428571429, 0.461538462, 2.285714286, 1000, 1000, 1.928571429,0.571428571,0.357142857,1000,1.5,1.142857143,0.857142857,1.785714286,1.285714286,0.785714286,0.714285714,1.428571429,0.714285714,1.285714286,2.571428571],2)
ATQ_MOTOR_INT_7T = np.round([1.428571429,0.46153846,0.5,0.714285714,1.769230769,0.285714286,0.642857143,1.214285714,1000,1000,1.142857143,0.428571429,0.714285714,1000,0.928571429,0.785714286,0.642857143,0.928571429,0.571428571,0.5,1.214285714,0.785714286,1000,0.928571429,1.714285714],2)
ATQ_VOCAL_FREQ_7T = np.round([0.636363636,0.307692308,0.076923077,1.416666667,1.384615385,0.307692308,0,0.461538462,1000,1000,0.583333333,0.307692308,0.333333333,1000,1.230769231,0.230769231,0.15384615,1.307692308,1.153846154,0.615384615,0.307692308,0.307692308,0,1.384615385,2.384615385],2)
ATQ_VOCAL_INT_7T = np.round([0.454545455,0.307692308,0.076923077,1,0.923076923,0.076923077,0,0.307692308,1000,1000,0.461538462,0.230769231,0.666666667,1000,1,0.307692308,0.153846154,0.538461538,0.615384615,0.153846154,0.846153846,0.615384615,0,0.923076923,1.461538462],2)

#Put the metrics into order of PET
PUTS_3T  = np.zeros((17))
ATQ_MOTOR_FREQ_3T  = np.zeros((17))
ATQ_MOTOR_INT_3T  = np.zeros((17))
ATQ_VOCAL_FREQ_3T  = np.zeros((17))
ATQ_VOCAL_INT_3T  = np.zeros((17))

for i in range(17):
    PUTS_3T[i] = PUTS_7T[Patients_7T_codes[i]-1]
    ATQ_MOTOR_FREQ_3T[i] = ATQ_MOTOR_FREQ_7T[Patients_7T_codes[i]-1]
    ATQ_MOTOR_INT_3T[i] = ATQ_MOTOR_INT_7T[Patients_7T_codes[i]-1]
    ATQ_VOCAL_FREQ_3T[i] = ATQ_VOCAL_FREQ_7T[Patients_7T_codes[i]-1]
    ATQ_VOCAL_INT_3T[i] = ATQ_VOCAL_INT_7T[Patients_7T_codes[i]-1]
    
#IMPORTANT: Manually appended also by GTSPp017, that had no 7T measurement!!! 
"""
Create all needed YGTSS vectors------------------------------------------------
"""

GTS_BPnd_ATQ_PUTS  = np.zeros((12,18))
mat = scipy.io.loadmat('/data/pt_02518/7T_Processing/Stats_Results/PET/Caudate/PET_BPnd_Caudate_stats.mat')
GTS_BPnd_ATQ_PUTS[0,:] = mat['GTS_BPnd']
mat = scipy.io.loadmat('/data/pt_02518/7T_Processing/Stats_Results/PET/Pallidum/PET_BPnd_Pallidum_stats.mat')
GTS_BPnd_ATQ_PUTS[1,:] = mat['GTS_BPnd']  
mat = scipy.io.loadmat('/data/pt_02518/7T_Processing/Stats_Results/PET/Putamen/PET_BPnd_Putamen_stats.mat')
GTS_BPnd_ATQ_PUTS[2,:] = mat['GTS_BPnd']    
mat = scipy.io.loadmat('/data/pt_02518/7T_Processing/Stats_Results/PET/Striatum/PET_BPnd_Striatum_stats.mat')
GTS_BPnd_ATQ_PUTS[3,:] = mat['GTS_BPnd']        
mat = scipy.io.loadmat('/data/pt_02518/7T_Processing/Stats_Results/PET/Thalamus/PET_BPnd_Thalamus_stats.mat')
GTS_BPnd_ATQ_PUTS[4,:] = mat['GTS_BPnd']           
GTS_BPnd_ATQ_PUTS[5,:] = PUTS_3T           
GTS_BPnd_ATQ_PUTS[6,:] =  ATQ_MOTOR_FREQ_3T
GTS_BPnd_ATQ_PUTS[7,:] =  ATQ_MOTOR_INT_3T
GTS_BPnd_ATQ_PUTS[8,:] =  (ATQ_MOTOR_FREQ_3T + ATQ_MOTOR_INT_3T)/2
GTS_BPnd_ATQ_PUTS[9,:] =  ATQ_VOCAL_FREQ_3T
GTS_BPnd_ATQ_PUTS[10,:] =   ATQ_VOCAL_INT_3T     
GTS_BPnd_ATQ_PUTS[11,:] =   (ATQ_VOCAL_FREQ_3T+ATQ_VOCAL_INT_3T)/2

final_data =  np.delete(GTS_BPnd_ATQ_PUTS, [], axis=1) 

def calculate_pvalues(df):
    dfcols = pd.DataFrame(columns=df.columns)
    pvalues = dfcols.transpose().join(dfcols, how='outer')
    for r in df.columns:
        for c in df.columns:
            tmp = df[df[r].notnull() & df[c].notnull()]
            pvalues[r][c] = round(pearsonr(tmp[r], tmp[c])[1], 4)
    return pvalues

# Print heatmap
final_data_T = final_data.transpose()
frame_patients = pd.DataFrame(data=final_data_T,index=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15], columns=['Caudate','Pallidum','Putamen','STR','Thalamus','PUTS','ATQ_MOTOR_FREQ','ATQ_MOTOR_INT','ATQ_MOTOR_MEAN','ATQ_VOCAL_FREQ','ATQ_VOCAL_INT','ATQ_VOCAL_MEAN'])  
corr = frame_patients.corr(method = 'pearson')
corr
pvals = calculate_pvalues(frame_patients)
pvals
plt.figure(figsize=(20,20), dpi =300)
sb.heatmap(corr,annot=True, square = 'True',cmap="vlag")
imname = '/---/---/ATQ_PUTS_vs_BPnd_pearsonCoeff.png'
plt.savefig(imname)      

plt.figure(figsize=(20,20), dpi =300)
sb.heatmap(corr,annot=pvals, square = 'True',cmap="rocket_r",vmin=0, vmax=1)
imname = '/---/---/ATQ_PUTS_vs_BPnd_pvals.png'
plt.savefig(imname)   


