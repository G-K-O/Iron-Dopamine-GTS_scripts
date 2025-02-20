#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 14:30:38 2024

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
from scipy.stats import pearsonr
import pandas as pd


"""
Create all needed YGTSS vectors------------------------------------------------
"""

PUTS_7T = [33,26,20,22,22,24,28,27,1000,1000,18,23,28,1000,26,27,1000,21,13,23,28,27,37,28,33]

ATQ_MOTOR_FREQ = np.round([1.357142857, 0.615384615, 1.214285714,0.928571429, 2.846153846, 0.428571429, 0.461538462, 2.285714286, 1000, 1000, 1.928571429,0.571428571,0.357142857,1000,1.5,1.142857143,0.857142857,1.785714286,1.285714286,0.785714286,0.714285714,1.428571429,0.714285714,1.285714286,2.571428571],2)
ATQ_MOTOR_INT = np.round([1.428571429,0.46153846,0.5,0.714285714,1.769230769,0.285714286,0.642857143,1.214285714,1000,1000,1.142857143,0.428571429,0.714285714,1000,0.928571429,0.785714286,0.642857143,0.928571429,0.571428571,0.5,1.214285714,0.785714286,1000,0.928571429,1.714285714],2)

ATQ_VOCAL_FREQ = np.round([0.636363636,0.307692308,0.076923077,1.416666667,1.384615385,0.307692308,0,0.461538462,1000,1000,0.583333333,0.307692308,0.333333333,1000,1.230769231,0.230769231,0.15384615,1.307692308,1.153846154,0.615384615,0.307692308,0.307692308,0,1.384615385,2.384615385],2)
ATQ_VOCAL_INT = np.round([0.454545455,0.307692308,0.076923077,1,0.923076923,0.076923077,0,0.307692308,1000,1000,0.461538462,0.230769231,0.666666667,1000,1,0.307692308,0.153846154,0.538461538,0.615384615,0.153846154,0.846153846,0.615384615,0,0.923076923,1.461538462],2)
 
GTS_QSMvsMetrics  = np.zeros((15,25)) 

mat = scipy.io.loadmat('/data/pt_02518/7T_Processing/Stats_Results/T2star_7T/Caudate_T2star_stats_7T.mat')
GTS_QSMvsMetrics[0,:] = mat['GTS']
mat = scipy.io.loadmat('/data/pt_02518/7T_Processing/Stats_Results/T2star_7T/Pallidum_T2star_stats_7T.mat')
GTS_QSMvsMetrics[1,:] = mat['GTS']  
mat = scipy.io.loadmat('/data/pt_02518/7T_Processing/Stats_Results/T2star_7T/Putamen_T2star_stats_7T.mat')
GTS_QSMvsMetrics[2,:] = mat['GTS']    
mat = scipy.io.loadmat('/data/pt_02518/7T_Processing/Stats_Results/T2star_7T/RN_T2star_stats_7T.mat')
GTS_QSMvsMetrics[3,:] = mat['GTS']
mat = scipy.io.loadmat('/data/pt_02518/7T_Processing/Stats_Results/T2star_7T/SN_T2star_stats_7T.mat')
GTS_QSMvsMetrics[4,:] = mat['GTS']  
mat = scipy.io.loadmat('/data/pt_02518/7T_Processing/Stats_Results/T2star_7T/Striatum_T2star_stats_7T.mat')
GTS_QSMvsMetrics[5,:] = mat['GTS']        
mat = scipy.io.loadmat('/data/pt_02518/7T_Processing/Stats_Results/T2star_7T/Thalamus_T2star_stats_7T.mat')
GTS_QSMvsMetrics[6,:] = mat['GTS']
mat = scipy.io.loadmat('/data/pt_02518/7T_Processing/Stats_Results/T2star_7T/StN_T2star_stats_7T.mat')
GTS_QSMvsMetrics[7,:] = mat['GTS']
GTS_QSMvsMetrics[8,:] =   ATQ_MOTOR_FREQ         
GTS_QSMvsMetrics[9,:] =   ATQ_MOTOR_INT
GTS_QSMvsMetrics[10,:] =   (ATQ_MOTOR_FREQ+ATQ_MOTOR_INT)/2
GTS_QSMvsMetrics[11,:] =   ATQ_VOCAL_FREQ         
GTS_QSMvsMetrics[12,:] =   ATQ_VOCAL_INT
GTS_QSMvsMetrics[13,:] =   (ATQ_VOCAL_FREQ+ATQ_VOCAL_INT)/2
GTS_QSMvsMetrics[14,:] =  PUTS_7T

final_data =  np.delete(GTS_QSMvsMetrics, [6,8,9, 12,13, 20], axis=1)  
#GTSMp007,GTSMp013, GTSMp021 excluded due to bad quality,the rest,no ATQ-PUTS data


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
frame_patients = pd.DataFrame(data=final_data_T,index=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19], columns=['Caudate','Pallidum','Putamen','RN','SN','STR','Thalamus','StN','ATQ_MOTOR_FREQ','ATQ_MOTOR_INT','ATQ_MOTOR_MEAN','ATQ_VOCAL_FREQ','ATQ_VOCAL_INT','ATQ_VOCAL_MEAN','PUTS'])  
corr = frame_patients.corr(method = 'pearson')
corr
pvals = calculate_pvalues(frame_patients)
pvals
plt.figure(figsize=(20,20), dpi =300)
sb.heatmap(corr,annot=True, square = 'True',cmap="vlag")
imname = '/---/---/ATQ_PUTS_vs_PCS_pearsonsCoeff.png'
plt.savefig(imname)      
    


plt.figure(figsize=(20,20), dpi =300)
sb.heatmap(corr,annot=pvals, square = 'True',cmap="rocket_r",vmin=0, vmax=1)
imname = '/---/---/ATQ_PUTS_vs_T2star_pvals.png'
plt.savefig(imname)   


