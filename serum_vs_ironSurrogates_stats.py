#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script is for serum vs iron
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

#VALUES = 1000 are just like this for easy manual exclusion
#Patients_age = [33,29,34,25,28,24,56,25,30,44,28,21,37,24,27,39,38,37,21,20,20,23,25,40,26]
#Patients_gender = [1,1,1,1,1,1,2,1,2,1,1,1,1,1,2,1,2,1,1,1,1,1,2,1,1]
#gts_ferittin_7T = np.array([67.5,81.3,61,157,76.7,1106,69.1,1000,25.3,142,116,75.3,1000,
#                             124,1000,194,61.8,312,44.7,85.2,27.2,65.6,38.4,215,167])
#gts_iron_7T = np.array([18.9,15.5,18.2,22.7,19.1,15.8,10.2,1000,13.1,13.9,17.4,16.3,1000,
#                        18.5,1000,18.6,16.9,16.9,24.3,16.2,16.8,21.2,12.6,12.7,17.2])
#gts_Transferrin_7T = np.array([2.4,2.3,2.3,2.1,2.5,2.1,2.3,1000,2.3,2.5,2.1,2.3,1000,
#                                2.7,1000,2.2,2.8,2.3,2.3,2.5,3.1,2.3,2.6,2.2,2.2])
#gts_TransfSat_7T = np.array([31.3,26.8,31.5,43,30.4,29.9,17.7,1000,22.7,22.1,33,28.2,1000,
#                              27.3,1000,33.6,24,29.2,42,25.8,21.6,36.7,19.3,23,31.1])
#gts_SolTransfRec_7T = np.array([2.63,3.56,2.68,2.47,2.5,3.11,4.45,1000,3.92,2.37,2.63,2.58,
#                                 1000,2.48,1000,2.76,2.29,2.89,2.8,2.69,3.47,2.66,3.97,2.4,2.11])

GTS_QSMvsMetrics_patients  = np.zeros((13,25))

mat1 = scipy.io.loadmat('/---/---/Caudate_QSM_stats_7T_refVentricles.mat')
GTS_QSMvsMetrics_patients[0,:] = mat1['GTS']
mat2 = scipy.io.loadmat('/---/---/Pallidum_QSM_stats_7T_refVentricles.mat')
GTS_QSMvsMetrics_patients[1,:] = mat2['GTS']  
mat3 = scipy.io.loadmat('/---/---/Putamen_QSM_stats_7T_refVentricles.mat')
GTS_QSMvsMetrics_patients[2,:] = mat3['GTS']    
mat4 = scipy.io.loadmat('/---/---/RN_QSM_stats_7T_refVentricles.mat')
GTS_QSMvsMetrics_patients[3,:] = mat4['GTS']
mat5 = scipy.io.loadmat('/---/---/SN_QSM_stats_7T_refVentricles.mat')
GTS_QSMvsMetrics_patients[4,:] = mat5['GTS']  
mat8 = scipy.io.loadmat('/---/---/Striatum_QSM_stats_7T_refVentricles.mat')
GTS_QSMvsMetrics_patients[5,:] = (mat8['GTS'])  
mat6 = scipy.io.loadmat('/---/---/Thalamus_QSM_stats_7T_refVentricles.mat')
GTS_QSMvsMetrics_patients[6,:] = mat6['GTS']
mat7 = scipy.io.loadmat('/---/---/StN_QSM_stats_7T_refVentricles.mat')
GTS_QSMvsMetrics_patients[7,:] = mat7['GTS']
GTS_QSMvsMetrics_patients[8,:] =      gts_ferittin_7T      
GTS_QSMvsMetrics_patients[9,:] =   gts_iron_7T
GTS_QSMvsMetrics_patients[10,:] =   gts_Transferrin_7T
GTS_QSMvsMetrics_patients[11,:] =    gts_TransfSat_7T        
GTS_QSMvsMetrics_patients[12,:] =   gts_SolTransfRec_7T

#replace 1116 (GTSMp006) with the mean ferittin and kept it in as there was error with the measurement in serum
final_data_patients =  np.delete(GTS_QSMvsMetrics_patients, [], axis=1) ##IMPORTANT: manually exclude all that have blood measurements issues and iron surrogates bad quality

def calculate_pvalues(df):
    dfcols = pd.DataFrame(columns=df.columns)
    pvalues = dfcols.transpose().join(dfcols, how='outer')
    for r in df.columns:
        for c in df.columns:
            tmp = df[df[r].notnull() & df[c].notnull()]
            pvalues[r][c] = round(pearsonr(tmp[r], tmp[c])[1], 4)
    return pvalues

# Print heatmap
final_data_T = final_data_patients.transpose()
frame_patients = pd.DataFrame(data=final_data_T,index=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19], columns=['Caudate','Pallidum','Putamen','RN','SN','STR','Thalamus','StN','Ferittin','Iron','Transferin','TransfSat','SolTransfRec'])  
corr = frame_patients.corr(method = 'pearson')
corr
pvals = calculate_pvalues(frame_patients)
pvals
plt.figure(figsize=(20,20), dpi =300)
sb.heatmap(corr,annot=True, square = 'True',cmap="vlag")
imname = '/---/---/Blood_vs_QSM_patients_pearson.png'
plt.savefig(imname)      
    
plt.figure(figsize=(20,20), dpi =300)
sb.heatmap(corr,annot=pvals, square = 'True',cmap="rocket_r",vmin=0, vmax=1)
imname = '/---/---/Blood_vs_QSM_patients_pvals.png'
plt.savefig(imname)   





'''
CONTROLS---------------------------------------------------------------------------------------------------------------------------
'''

Controls_age = np.array([28,36,30,31,30,31,36,33,26,25,34,40,35,36,38,30,41,25,33,38,34,24,34,43,38,27,37,31,30,28,31,49,27,22,23,28,20,23,29,28])
Controls_gender = np.array([1,2,1,1,2,1,1,1,1,1,2,2,1,1,1,1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,2,2,1]) #1: male, 2:female

ctrl_ferittin_7T = np.array([291,18.2,17.9,275,18.1,132,153,79.7,144,116,50.7,90.2,19,14,128,16.9,61.,41.7,92,105,229,
                             38.6,125,228,68.1,31,436 ,116,59.6,205,70.6,174,19.2,63.8,64.5,7,127,14.1,30.2 ,111])
ctrl_iron_7T = np.array([18.7,19.5,10.8,14.4,17,13.1,32.6,15.4,15.7,13.3,14.6,15.8,8.,13.6,18.4,17.6,26.2,20.5, 22.7,42.6,
                         16.3,13.5,37.1,16.3,23.2,17.3,10,12.2,26.5,20.1,17.5,12.3,20.6,15.8,15,5.3,35,14.4,25.6,26.2])
ctrl_Transferrin_7T = np.array([2.6,2.4,3.1,2.2,2.5,1.8,2.5,2.5,2.5,2.9,2.1,2.9,3.1,2.7,2.9,2.8,2.6,2.7,2.3,2.8,
                                2.6,2.4,2.6,2.7,2.9,2.6,1.9,2.6,2.3,2.5,2.3,2.4,2.9,2.6,2.2,3.8,2.8,2.8,2.3,2.4])
ctrl_TransfSat_7T = np.array([28.6,32.3,13.9,26.1,27.1,29,51.9,24.5,25,18.3,27.7,21.7,11.2,20,25.3,25,40.1,30.2,39.3,60.6,25,
                              22.4,56.8,24,31.8,26.5,20.9,18.7,45.9,32,30.3,20.4,28.3,24.2,27.1,5.6,49.8,20.5,44.3,43.4])
ctrl_SolTransfRec_7T = np.array([2.77,2.52,4.48,2.28,2.3,2.52,2.76,2.48,2.86,2.42,2.2,2.01,3.47,2.46,2.6,3.01,2.44,3.81,2.36,2.93,
                                 3.17,2.97,2.67,2.41,2.76,2.98,3.07,3.51,2.2,2.97,2.63,2.31,3.07,3.33,2.56,3.97,2.24,4.57,2.59,2.2])

GTS_QSMvsMetrics_ctrl  = np.zeros((13,40))
mat = scipy.io.loadmat('/---/---/Caudate_DECOMPOSE_stats_7T_corrRange.mat')
GTS_QSMvsMetrics_ctrl[0,:] = mat['Controls']
mat = scipy.io.loadmat('/---/---/Pallidum_DECOMPOSE_stats_7T_corrRange.mat')
GTS_QSMvsMetrics_ctrl[1,:] = mat['Controls']  
mat = scipy.io.loadmat('/---/---/Putamen_DECOMPOSE_stats_7T_corrRange.mat')
GTS_QSMvsMetrics_ctrl[2,:] = mat['Controls']    
mat = scipy.io.loadmat('/---/---/RN_DECOMPOSE_stats_7T_corrRange.mat')
GTS_QSMvsMetrics_ctrl[3,:] = mat['Controls']
mat = scipy.io.loadmat('/---/---/SN_DECOMPOSE_stats_7T_corrRange.mat')
GTS_QSMvsMetrics_ctrl[4,:] = mat['Controls']  
mat = scipy.io.loadmat('/---/---/Striatum_DECOMPOSE_stats_7T_corrRange.mat')
GTS_QSMvsMetrics_ctrl[5,:] = mat['Controls']        
mat = scipy.io.loadmat('/---/---/Thalamus_DECOMPOSE_stats_7T_corrRange.mat')
GTS_QSMvsMetrics_ctrl[6,:] = mat['Controls']   
mat = scipy.io.loadmat('/---/---/StN_DECOMPOSE_stats_7T_corrRange.mat')
GTS_QSMvsMetrics_ctrl[7,:] = mat['Controls']   
GTS_QSMvsMetrics_ctrl[8,:] =      ctrl_ferittin_7T      
GTS_QSMvsMetrics_ctrl[9,:] =   ctrl_iron_7T
GTS_QSMvsMetrics_ctrl[10,:] =   ctrl_Transferrin_7T
GTS_QSMvsMetrics_ctrl[11,:] =    ctrl_TransfSat_7T        
GTS_QSMvsMetrics_ctrl[12,:] =  ctrl_SolTransfRec_7T

final_data_ctrl =  np.delete(GTS_QSMvsMetrics_ctrl, [], axis=1) ##IMPORTANT: manually exclude all that have blood measurements issues and iron surrogates bad quality


def calculate_pvalues(df):
    dfcols = pd.DataFrame(columns=df.columns)
    pvalues = dfcols.transpose().join(dfcols, how='outer')
    for r in df.columns:
        for c in df.columns:
            tmp = df[df[r].notnull() & df[c].notnull()]
            pvalues[r][c] = round(scs.pearsonr(tmp[r], tmp[c])[1], 4)
    return pvalues

# Print heatmap
final_data_T = final_data_ctrl.transpose()
frame_patients = pd.DataFrame(data=final_data_T,index=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36], columns=['Caudate','Pallidum','Putamen','RN','SN','STR','Thalamus','StN','Ferittin','Iron','Transferin','TransfSat','SolTransfRec'])  
corr = frame_patients.corr(method = 'pearson')
corr
pvals = calculate_pvalues(frame_patients)
pvals
plt.figure(figsize=(20,20), dpi =300)
sb.heatmap(corr,annot=True, square = 'True',cmap="vlag")
imname = '/data/hu_gkotsoulias/Desktop/Blood_vs_T2star_CTRL.png'
plt.savefig(imname)      
    
plt.figure(figsize=(20,20), dpi =300)
sb.heatmap(corr,annot=pvals, square = 'True',cmap="rocket_r",vmin=0, vmax=1)
imname = '/data/hu_gkotsoulias/Desktop/Blood_vs_T2star_CTRL_pvals.png'
plt.savefig(imname)   

'''
MERGING STATS
'''

## Merge all data for final stats
final_data_QSM_all = np.concatenate((final_data_ctrl,final_data_patients),axis=1)


##Load all data if analysis is done again
#mat = scipy.io.loadmat("/---/---/finalAnalysis_data_all_forHeatmaps.mat")
#final_data_QSM_all = mat['final_data_QSM_all']



final_data_T = final_data_PCS_all.transpose()
frame_patients = pd.DataFrame(data=final_data_T,index=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54],columns=['Caudate','Pallidum','Putamen','RN','SN','STR','Thalamus','StN','Ferittin','Iron','Transferin','TransfSat','SolTransfRec'])  
corr = frame_patients.corr(method = 'pearson')
corr
pvals = calculate_pvalues(frame_patients)
pvals
plt.figure(figsize=(20,20), dpi =100)
sb.heatmap(corr,annot=True, square = 'True',cmap="vlag")
imname = '/---/---/Blood_vs_PCS_allParticipants_pearson.png'
plt.savefig(imname)      
    
plt.figure(figsize=(20,20), dpi =100)
sb.heatmap(corr,annot=pvals, square = 'True',cmap="rocket_r",vmin=0, vmax=1)
imname = '/---/---/Blood_vs_PCS_allParticipants_pvals.png'
plt.savefig(imname)   


Weights = [2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
Weights_colour = []
for i in range(np.size(Weights)):
    if (Weights[i] == 1):
        Weights_colour.append("maroon")
    else:
        Weights_colour.append("cornflowerblue")

data = pd.DataFrame(data={'PUT QSM': final_data_T[:,2], 'FERRITIN': final_data_T[:,8] })  #,'Weights': Weights
data['color']= Weights_colour

plt.figure(figsize=(5,5), dpi =500)
sb.set_style("darkgrid")
sb.set_context("paper", font_scale=1.5, rc={"lines.linewidth": 2})
g = sb.regplot(x="FERRITIN", y="PUT PCS", data=data,scatter_kws={"s": 60,"color": "black",'facecolors':data['color']})#,hue="Weights"
g = (g.set(xlim=(-50, 480), ylim=(0.030, 0.09)))#-49, 270
imname = '/---/---/QSM_Putamen_Vs_Ferritin_correlationGraph.png'
plt.savefig(imname)   

data1 = pd.DataFrame(data={'x':  final_data_T[:,2], 'y': final_data_T[:,8] })  #,'Weights': Weights
p = sb.regplot(x="y", y="x",data=data1)

#calculate slope and intercept of regression equation
slope, intercept, r, g, sterr = scipy.stats.linregress(x=p.get_lines()[0].get_xdata(),
                                                       y=p.get_lines()[0].get_ydata())
print(intercept, slope)

