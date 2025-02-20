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
from scipy.stats import pearsonr
import pandas as pd


"""
Create all needed YGTSS vectors------------------------------------------------
"""

YGTSS_Motor_7T =           [11,18,12,13,15,15,17,16,7,1000,19,9,8 ,14,19,21,17,15,14,10,18,18,12,14,18]
YGTSS_Vocal_7T =           [10,13, 0,14,14,12, 0,8, 0,1000,10,7,9 ,10,18,11,13, 9,16, 6,13,22, 0,14,14]
YGTSS_SubjImpairment_7T =  [10,20,10,10,20, 0,10,20,0,1000,20,0,30,20,40,10,10,20,10,30,30,30,20,10,30]

YGTSS_Full_7T = []
for i in range(len(YGTSS_Motor_7T)):
    YGTSS_Full_7T.append(YGTSS_Motor_7T[i] + YGTSS_Vocal_7T[i])

YGTSS_Total_7T = []
for i in range(len(YGTSS_Motor_7T)):
    YGTSS_Total_7T.append(YGTSS_Motor_7T[i] + YGTSS_Vocal_7T[i]+ YGTSS_SubjImpairment_7T[i])

mat = scipy.io.loadmat('/---/---/Caudate_QSM_stats_7T_refVentricles.mat') #/---/---/ to be changed with directory
GTS_QSMvsMetrics[0,:] = mat['GTS']
mat = scipy.io.loadmat('/---/---/Pallidum_QSM_stats_7T_refVentricles.mat')
GTS_QSMvsMetrics[1,:] = mat['GTS']  
mat = scipy.io.loadmat('/---/---/Putamen_QSM_stats_7T_refVentricles.mat')
GTS_QSMvsMetrics[2,:] = mat['GTS']    
mat = scipy.io.loadmat('/---/---/RN_QSM_stats_7T_refVentricles.mat')
GTS_QSMvsMetrics[3,:] = mat['GTS']
mat = scipy.io.loadmat('/---/---/SN_QSM_stats_7T_refVentricles.mat')
GTS_QSMvsMetrics[4,:] = mat['GTS']  
mat = scipy.io.loadmat('/---/---/Striatum_QSM_stats_7T_refVentricles.mat')
GTS_QSMvsMetrics[5,:] = mat['GTS']        
mat = scipy.io.loadmat('/---/---/Thalamus_QSM_stats_7T_refVentricles.mat')
GTS_QSMvsMetrics[6,:] = mat['GTS']
mat = scipy.io.loadmat('/---/---/StN_QSM_stats_7T_refVentricles.mat')
GTS_QSMvsMetrics[7,:] = mat['GTS']


GTS_QSMvsMetrics[8,:] = YGTSS_Motor_7T        
GTS_QSMvsMetrics[9,:] = YGTSS_Vocal_7T
GTS_QSMvsMetrics[10,:] = YGTSS_Full_7T        
GTS_QSMvsMetrics[11,:] = YGTSS_Total_7T 
final_data =  np.delete(GTS_QSMvsMetrics, [ 6, 9, 12, 20], axis=1)  #GTSMp007,GTSMp013, GTSMp021 excluded due to bad quality,GTSMp010 bad quality and no YGTSS

def calculate_pvalues(df):
    dfcols = pd.DataFrame(columns=df.columns)
    pvalues = dfcols.transpose().join(dfcols, how='outer')
    for r in df.columns:
        for c in df.columns:
            tmp = df[df[r].notnull() & df[c].notnull()]
            pvalues[r][c] = round(scs.pearsonr(tmp[r], tmp[c])[1], 4)
    return pvalues

# Print heatmap
final_data_T = final_data.transpose()
frame_patients = pd.DataFrame(data=final_data_T,index=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21], columns=['Caudate','Pallidum','Putamen','RN','SN','STR','Thalamus','StN','YGTSS Motor','YGTSS Vocal','YGTSS Full','YGTSS Total'])  
corr = frame_patients.corr(method = 'pearson')
corr
pvals = calculate_pvalues(frame_patients)
pvals
plt.figure(figsize=(10,10), dpi = 100)
sb.heatmap(corr,annot=True, square = 'True',cmap="vlag")   
    
plt.figure(figsize=(10,10), dpi =100)
sb.heatmap(corr,annot=pvals, square = 'True',cmap="rocket_r",vmin=0, vmax=1)
 

      

data = pd.DataFrame(data={'YGTSS Total': final_data[0,:], 'iron': final_data[5,:] })  
sb.set_style("darkgrid")
g=sb.lmplot(x="BPnd", y="YGTSS Motor",palette='Set1',height=5, data=data)
sb.set_context("notebook", font_scale=1.5, rc={"lines.linewidth": 1.5})
g = (g.set_axis_labels("YGTSS Motor","CAUDATE BPnd").set(xlim=(8,22), ylim=(1.48, 2.73)))
imname = '/---/---/CAUDATE_BPnd_Vs_BPnd_correlationGraph.png'
plt.savefig(imname)         

data1 = pd.DataFrame(data={'x': final_data[3,:], 'y': final_data[5,:] })  #,'Weights': Weights
p = sb.regplot(x="x", y="y",data=data1)

#calculate slope and intercept of regression equation
slope, intercept, r, g, sterr = scipy.stats.linregress(x=p.get_lines()[0].get_xdata(),
                                                       y=p.get_lines()[0].get_ydata())
print(intercept, slope)

    
