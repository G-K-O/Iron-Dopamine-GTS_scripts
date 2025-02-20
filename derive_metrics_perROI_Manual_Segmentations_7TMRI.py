#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  3 17:33:49 2023

@author: gkotsoulias
"""

import os
import numpy as np
import seaborn
from scipy.io import savemat
import nibabel as nib
import matplotlib.pyplot as plt
from scipy import stats


"""
Controls-----------------------------------------------------------------------
"""
Controls_age = [28,36,30,31,30,31,36,33,26,25,34,40,35,36,38,30,41,25,33,38,34,24,34,43,38,27,37,31,30,28,31,49,27,22,23,28,20,23,29,28]
Controls_age_mean = np.average(Controls_age)
Controls_std_mean = np.std(Controls_age)
Controls_gender = [1,2,1,1,2,1,1,1,1,1,2,2,1,1,1,1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,2,2,1] #1: male, 2:female
Controls_dropouts = ['GTSMc001','GTSMc013','GTSMc026','GTSMc030','GTSMc015']

full_SN_QSM_Controls = np.zeros([1,40])

j=0
for i in range(1,41):
    if i<10:
        subject = 'GTSMc00'+str(i)
    else:
        subject = 'GTSMc0'+str(i)
    dir1 = '/data/pt_02518/7T_Processing/Controls_Cohort/'+subject+'/QSM2STDSpace/'
    if os.path.exists(dir1): 
       if subject in Controls_dropouts:
           print('DROPOUT: ', subject)
           full_SN_QSM_Controls[0,j] = 1000 #Give an impossible value to the ones we dont include
           j=j+1
       else:
           print('ADDED IN STATS: ', subject)
           ## Load volumes------------------------------------------------------------
           QSM2_dir = os.path.join(dir1,'QSM_iLSQR_mean_5echoes.nii.gz') #CHANGE TYPE OF CONTRAST HERE BETWEEN QSM, PCS, DCS, T2*
           ROIs_dir = os.path.join(dir1,'ManualSegmentations/ROI_SN_Manual.nii.gz') #CHANGE the ROI here, RN, SN , StN  
           in_vol1 = nib.load(QSM2_dir)
           QSM = in_vol1.get_fdata()
       
           in_vol2 = nib.load(ROIs_dir)      
           ROIs = in_vol2.get_fdata()
        
           full_SN_QSM_Controls[0,j] = QSM[ROIs!= 0].mean()*100
           j=j+1
      

"""
Patients-----------------------------------------------------------------------
"""
Patients_age = [33,29,34,25,28,24,56,25,30,44,28,21,37,24,27,39,38,37,21,20,20,23,25,40,26]
Patients_age_mean = np.average(Patients_age)
Patients_std_mean = np.std(Patients_age)
Patients_gender = [1,1,1,1,1,1,2,1,2,1,1,1,1,1,2,1,2,1,1,1,1,1,2,1,1] #1: male, 2:female 
Patient_dropouts = ['GTSMp010','GTSMp013','GTSMp021', 'GTSMp007']


full_SN_QSM_GTS = np.zeros([1,25])

j=0
for i in range(1,26):
    if i<10:
        subject = 'GTSMp00'+str(i)
    else:
        subject = 'GTSMp0'+str(i)
    dir1 = '/data/pt_02518/7T_Processing/GTS_Cohort/'+subject+'/QSM2STDSpace/'
    if os.path.exists(dir1): 
       if subject in Patient_dropouts:
           print('DROPOUT: ', subject)
           full_SN_QSM_GTS[0,j] = 1000 #Give an impossible value to the ones we dont include
           j=j+1
       else:
           print('ADDED IN STATS: ', subject)
           ## Load volumes------------------------------------------------------------
           QSM2_dir = os.path.join(dir1,'QSM_iLSQR_mean_5echoes.nii.gz') #CHANGE TYPE OF CONTRAST HERE BETWEEN QSM, PCS, DCS, T2*
           ROIs_dir = os.path.join(dir1,'ManualSegmentations/ROI_SN_Manual.nii.gz') #CHANGE the ROI here, RN, SN , StN   
           in_vol1 = nib.load(QSM2_dir)
           QSM = in_vol1.get_fdata()
       
           in_vol2 = nib.load(ROIs_dir)      
           ROIs = in_vol2.get_fdata()
        
           full_SN_QSM_GTS[0,j] = QSM[ROIs!= 0].mean()*100
           j=j+1
    
  
dic_mat = {"GTS": GTS, "Controls": Controls}
savemat("/---/---/***ROI***_***CONTRAST_NAME***_7TMRI.mat", dic_mat)





#

