#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 19:15:04 2023

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

"""
Some subjects had to be excluded due to lack of quality data either on 7T or
PET. Refer to manuscript for full details.
"""
Controls_7T_codes =  [13,25,9,17,19,35,22,34,40,33,29,36,18,15,38,39,37]
Patients_7T_codes =  [4,6,8,12,11,14,7,15,16,18,20,21,17,19,22,25]


#Here load the contrast QSM, PCS or T2* for Caudate, Putamen, Pallidum
mat = scipy.io.loadmat('/---/---/Caudate_QSM_stats_7T_refVentricles.mat') 
Controls  = mat['Controls']
GTS = mat['GTS']
Controls_Caud_QSM  = np.zeros((1,17))
GTS_Caud_QSM   = np.zeros((1,16))

for i in range(17):
    print('GTSMc'+str(Controls_7T_codes[i])+', index in matrix:'+str(Controls_7T_codes[i]-1)+', value:'+str(Controls[0,Controls_7T_codes[i]-1]))
    Controls_Caud_QSM[0,i] = Controls[0,Controls_7T_codes[i]-1]
print('----------------------------------------------------------')
for i in range(16):
    print('GTSMp'+str(Patients_7T_codes[i])+', index in matrix:'+str(Patients_7T_codes[i]-1)+', value:'+str(GTS[0,Patients_7T_codes[i]-1]))
    GTS_Caud_QSM[0,i] = GTS[0,Patients_7T_codes[i]-1]
    
    
    
mat = scipy.io.loadmat('/---/---/Pallidum_QSM_stats_7T_refVentricles.mat')
Controls  = mat['Controls']
GTS = mat['GTS']
Controls_Pall_QSM   = np.zeros((1,17))
GTS_Pall_QSM  = np.zeros((1,16))

for i in range(17):
    print('GTSMc'+str(Controls_7T_codes[i])+', index in matrix:'+str(Controls_7T_codes[i]-1)+', value:'+str(Controls[0,Controls_7T_codes[i]-1]))
    Controls_Pall_QSM[0,i] = Controls[0,Controls_7T_codes[i]-1]
print('----------------------------------------------------------')
for i in range(16):
    print('GTSMp'+str(Patients_7T_codes[i])+', index in matrix:'+str(Patients_7T_codes[i]-1)+', value:'+str(GTS[0,Patients_7T_codes[i]-1]))
    GTS_Pall_QSM[0,i] = GTS[0,Patients_7T_codes[i]-1]
    

mat = scipy.io.loadmat('/---/---/Putamen_QSM_stats_7T_refVentricles.mat')
Controls  = mat['Controls']
GTS = mat['GTS']
Controls_Puta_QSM   = np.zeros((1,17))
GTS_Puta_QSM  = np.zeros((1,16))

for i in range(17):
    print('GTSMc'+str(Controls_7T_codes[i])+', index in matrix:'+str(Controls_7T_codes[i]-1)+', value:'+str(Controls[0,Controls_7T_codes[i]-1]))
    Controls_Puta_QSM[0,i] = Controls[0,Controls_7T_codes[i]-1]
print('----------------------------------------------------------')
for i in range(16):
    print('GTSMp'+str(Patients_7T_codes[i])+', index in matrix:'+str(Patients_7T_codes[i]-1)+', value:'+str(GTS[0,Patients_7T_codes[i]-1]))
    GTS_Puta_QSM[0,i] = GTS[0,Patients_7T_codes[i]-1]   
    
   
Controls_Puta_Pall_Caud = np.zeros((3,51))
Controls_Puta_Pall_Caud[0,0:17] = Controls_Caud_QSM
Controls_Puta_Pall_Caud[0,17:34] = Controls_Pall_QSM
Controls_Puta_Pall_Caud[0,34:51] = Controls_Puta_QSM

GTS_Puta_Pall_Caud = np.zeros((3,48))
GTS_Puta_Pall_Caud[0,0:16] = GTS_Caud_QSM
GTS_Puta_Pall_Caud[0,16:32] = GTS_Pall_QSM
GTS_Puta_Pall_Caud[0,32:48] = GTS_Puta_QSM



"""
PET----------------------------------------------------------------------------
"""
mat = scipy.io.loadmat('/data/pt_02518/7T_Processing/Stats_Results/PET/Caudate/PET_BPnd_Caudate_stats.mat')
Controls_BPnd_Caud  = mat['Controls_BPnd']
GTS_BPnd_Caud = mat['GTS_BPnd']

mat = scipy.io.loadmat('/data/pt_02518/7T_Processing/Stats_Results/PET/Pallidum/PET_BPnd_Pallidum_stats.mat')
Controls_BPnd_Pall  = mat['Controls_BPnd']
GTS_BPnd_Pall = mat['GTS_BPnd']

mat = scipy.io.loadmat('/data/pt_02518/7T_Processing/Stats_Results/PET/Putamen/PET_BPnd_Putamen_stats.mat')
Controls_BPnd_Puta  = mat['Controls_BPnd']
GTS_BPnd_Puta = mat['GTS_BPnd']

Controls_Puta_Pall_Caud[1,0:17] = Controls_BPnd_Caud
Controls_Puta_Pall_Caud[1,17:34] = Controls_BPnd_Pall
Controls_Puta_Pall_Caud[1,34:51] = Controls_BPnd_Puta

GTS_Puta_Pall_Caud[1,0:16] = GTS_BPnd_Caud
GTS_Puta_Pall_Caud[1,16:32] = GTS_BPnd_Pall
GTS_Puta_Pall_Caud[1,32:48] = GTS_BPnd_Puta

ctrl_QSM =  np.delete(Controls_Puta_Pall_Caud, [], axis=1) 
gts_QSM =  np.delete(GTS_Puta_Pall_Caud, [], axis=1)
dataset_QSM = np.concatenate((ctrl_QSM,gts_QSM),axis=1)
Weights = np.transpose(np.concatenate((np.ones(np.shape(ctrl_QSM)[1])*2,np.ones(np.shape(gts_QSM)[1])),axis=0))

data_BPnd_vs_QSM = pd.DataFrame(data={'BPnd':dataset_QSM[1,0:np.shape(Weights)[0]], 'QSM': dataset_QSM[0,0:np.shape(Weights)[0]], 'Weights': Weights[0:np.shape(Weights)[0]] })  
   
sb.set_style("darkgrid")
sb.lmplot(x="BPnd", y="QSM",hue="Weights",
             palette='Set1',
             height=10, data=data_BPnd_vs_QSM)
sb.set_context("notebook", font_scale=1.5, rc={"lines.linewidth": 1.5})
imname = '/---/---/QSM_Vs_BPnd_correlationGraph.png'
plt.savefig(imname) 


"""
Correlation metrics------------------------------------------------------------
"""
Controls_corrcoef = scipy.stats.pearsonr(dataset_QSM[1,0:47],dataset_QSM[0,0:47])
GTS_corrcoef = scipy.stats.pearsonr(dataset_QSM[1,47:88],dataset_QSM[0,47:88])

Controls_spearmancoef =  scipy.stats.spearmanr(dataset_QSM[1,0:47],dataset_QSM[0,0:47])
GTS_spearmancoef = scipy.stats.spearmanr(dataset_QSM[1,47:88],dataset_QSM[0,47:88])







    
