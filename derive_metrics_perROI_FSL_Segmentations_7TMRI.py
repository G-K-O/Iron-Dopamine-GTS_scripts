#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 23 16:23:15 2023

@author: gkotsoulias
"""


import os
import numpy as np
import seaborn
from scipy.io import savemat
import nibabel as nib
import matplotlib.pyplot as plt
from scipy import stats
import ants
from scipy import ndimage
"""
Controls-----------------------------------------------------------------------
"""
Controls_age = [28,36,30,31,30,31,36,33,26,25,34,40,35,36,38,30,41,25,33,38,34,24,34,43,38,27,37,31,30,28,31,49,27,22,23,28,20,23,29,28]
Controls_age_mean = np.average(Controls_age)
Controls_std_mean = np.std(Controls_age)
Controls_gender = [1,2,1,1,2,1,1,1,1,1,2,2,1,1,1,1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,2,2,1] #1: male, 2:female
Controls_dropouts = ['GTSMc001','GTSMc013','GTSMc026','GTSMc030','GTSMc015'] 

Controls = np.zeros([1,40])

j=0
for i in range(1,41):
    if i<10:
        subject = 'GTSMc00'+str(i)
    else:
        subject = 'GTSMc0'+str(i)
    dir1 = '/---/---/'+subject+'/'
    if os.path.exists(dir1): 
       if subject in Controls_dropouts:
           print('DROPOUT: ', subject)
           full_Thalamus_T2star_Controls[0,j] = 1000 Give an impossible value to the ones we dont include
           j=j+1
       else:
           print('ADDED IN STATS: ', subject)

           ## Load volumes------------------------------------------------------------
           T2star_dir = os.path.join(dir1,'T2star/6echoes_ARLO_R2starmap.nii.gz') #LOAD CONTRAST
           ROIs_R_dir = os.path.join(dir1,'FSL_SubcorticalSeg/FSL_Seg_-R_Thal_corr.nii.gz')  #LOAD SEGMENTATION RIGHT HEMISPHERE
           ROIs_L_dir = os.path.join(dir1,'FSL_SubcorticalSeg/FSL_Seg_-L_Thal_corr.nii.gz') #LOAD SEGMENTATION LEFT HEMISPHERE

           contrast = ants.image_read(T2star_dir).numpy()
           ROIs_L_erod = ants.image_read(ROIs_L_dir).numpy()
           ROIs_R_erod = ants.image_read(ROIs_R_dir).numpy()
           ROIs_erod = ROIs_L_erod+ROIs_R_erod
           
           ROIs_erod_mask = np.zeros([266,266,176]) #hardcoded to resolution for bin mask
           ROIs_erod_mask[ROIs_erod!= 0]=1         
           ROIs_erod_mask = ndimage.binary_erosion(ROIs_erod_mask,iterations=2) #indicative erosion, adjust
          
           Controls[0,j] = contrast[ROIs_erod_mask!= 0].mean()
           j=j+1
        
        
"""
Patients-----------------------------------------------------------------------
"""
Patients_age = [33,29,34,25,28,24,56,25,30,44,28,21,37,24,27,39,38,37,21,20,20,23,25,40,26]
Patients_age_mean = np.average(Patients_age)
Patients_std_mean = np.std(Patients_age)
Patients_gender = [1,1,1,1,1,1,2,1,2,1,1,1,1,1,2,1,2,1,1,1,1,1,2,1,1] #1: male, 2:female 
Patient_dropouts = ['GTSMp010','GTSMp013','GTSMp021', 'GTSMp007'] 


GTS = np.zeros([1,25])

j=0
for i in range(1,26):
    if i<10:
        subject = 'GTSMp00'+str(i)
    else:
        subject = 'GTSMp0'+str(i)
    dir1 = '/---/---/'+subject+'/'
    if os.path.exists(dir1): 
       if subject in Patient_dropouts:
           print('DROPOUT: ', subject)
           full_Thalamus_T2star_GTS[0,j] = 1000 #Give an impossible value to the ones we dont include
           j=j+1
       else:
           print('ADDED IN STATS: ', subject)

           ## Load volumes------------------------------------------------------------
           T2star_dir = os.path.join(dir1,'T2star/6echoes_ARLO_R2starmap.nii.gz') #LOAD CONTRAST
           ROIs_R_dir = os.path.join(dir1,'FSL_SubcorticalSeg/FSL_Seg_-R_Thal_corr.nii.gz')  #LOAD SEGMENTATION RIGHT HEMISPHERE
           ROIs_L_dir = os.path.join(dir1,'FSL_SubcorticalSeg/FSL_Seg_-L_Thal_corr.nii.gz') #LOAD SEGMENTATION LEFT HEMISPHERE

           contrast = ants.image_read(T2star_dir).numpy()
           ROIs_L_erod = ants.image_read(ROIs_L_dir).numpy()
           ROIs_R_erod = ants.image_read(ROIs_R_dir).numpy()
           ROIs_erod = ROIs_L_erod+ROIs_R_erod
           
           ROIs_erod_mask = np.zeros([266,266,176]) #hardcoded to resolution for bin mask
           ROIs_erod_mask[ROIs_erod!= 0]=1         
           ROIs_erod_mask = ndimage.binary_erosion(ROIs_erod_mask,iterations=2) #indicative erosion, adjust
          
           GTS[0,j] = contrast[ROIs_erod_mask!= 0].mean()
           j=j+1
        
dic_mat = {"GTS": GTS, "Controls": Controls}
savemat("/---/---/***ROI***_***CONTRAST_NAME***_7TMRI.mat", dic_mat)





