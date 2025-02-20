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
Controls_7T_codes =  ['GTSMc013','GTSMc025','GTSMc009','GTSMc017','GTSMc019','GTSMc035','GTSMc022','GTSMc034','GTSMc040','GTSMc033','GTSMc029','GTSMc036','GTSMc018','GTSMc015','GTSMc038','GTSMc039','GTSMc037']
Controls_dropouts = [] 
Controls = np.zeros([1,len(Controls_7T_codes)])

for i in range(len(Controls_7T_codes)):
    print('Processing :', Controls_7T_codes[i]) 
    ## Load volumes-----------------------------------------------------------
    dir1 = '/---/---/'+str(Controls_7T_codes[i])+'/'
    BPnd_dir = dir1+'/T12QSMSTDSpace_Reg/BPnd_reg_to7T.nii.gz'
    R1_dir = dir1+'/T12QSMSTDSpace_Reg/R1_reg_to7T.nii.gz'  
    ROIs_R_dir = dir1+'/FSL_SubcorticalSeg/FSL_Seg_-R_Puta_corr.nii.gz'    
    ROIs_L_dir = dir1+'/FSL_SubcorticalSeg/FSL_Seg_-L_Puta_corr.nii.gz'

    BPnd = ants.image_read(BPnd_dir).numpy()
    R1 = ants.image_read(R1_dir).numpy()
    ROIs_L_erod = ants.image_read(ROIs_L_dir).numpy()
    ROIs_R_erod = ants.image_read(ROIs_R_dir).numpy()
    ROIs_erod = ROIs_L_erod+ROIs_R_erod
           
    ROIs_erod_mask = np.zeros([266,266,176])
    ROIs_erod_mask[ROIs_erod!= 0]=1         
    ROIs_erod_mask = ndimage.binary_erosion(ROIs_erod_mask,iterations=3)
        
    full_Putamen_BPnd_Controls[0,i] = BPnd[ROIs_erod_mask!= 0].mean()

######
#Important note here: Manually append to the controls the mean values of the ROIs for 
#controls GTSPc006, GTSPc0014, GTSPc0019, that had no 7T , thus needed to be separately 
#processed and added (in order not to repeat the whole code and change the processing 
#style). The final number of controls included for PET stats is 19: 1 exclusion for 
#external reasons before measurents, despite data were ok.
#####
    
    
"""
Patients-----------------------------------------------------------------------
"""
Patients_7T_codes = ['GTSMp004','GTSMp006','GTSMp008','GTSMp023','GTSMp012','GTSMp011','GTSMp014','GTSMp007','GTSMp015','GTSMp016','GTSMp018','GTSMp020','GTSMp021','GTSMp017','GTSMp019','GTSMp022','GTSMp025']
Patients_dropouts = ['GTSMp013','GTSMp009'] ##13 EXCLUDED DUE TO BAD QUALITY, 9 PET SYSTEM FAILURE, hence, 18/20 
GTS = np.zeros([1,len(Patients_7T_codes)])

for i in range(len(Patients_7T_codes)):
    print('Processing :', Patients_7T_codes[i]) 
    ## Load volumes-----------------------------------------------------------
    dir1 = '/---/---/'+str(Patients_7T_codes[i])+'/'
    BPnd_dir = dir1+'/T12QSMSTDSpace_Reg/BPnd_reg_to7T.nii.gz'
    R1_dir = dir1+'/T12QSMSTDSpace_Reg/R1_reg_to7T.nii.gz'  
    ROIs_R_dir = dir1+'/FSL_SubcorticalSeg/FSL_Seg_-R_Puta_corr.nii.gz'    
    ROIs_L_dir = dir1+'/FSL_SubcorticalSeg/FSL_Seg_-L_Puta_corr.nii.gz'

    BPnd = ants.image_read(BPnd_dir).numpy()
    R1 = ants.image_read(R1_dir).numpy()
    ROIs_L_erod = ants.image_read(ROIs_L_dir).numpy()
    ROIs_R_erod = ants.image_read(ROIs_R_dir).numpy()
    ROIs_erod = ROIs_L_erod+ROIs_R_erod
           
    ROIs_erod_mask = np.zeros([266,266,176])
    ROIs_erod_mask[ROIs_erod!= 0]=1         
    ROIs_erod_mask = ndimage.binary_erosion(ROIs_erod_mask,iterations=2)
    
    full_Putamen_BPnd_GTS[0,i] = BPnd[ROIs_erod_mask!= 0].mean()

dic_mat = {"GTS_BPnd": GTS, "Controls_BPnd": Controls}
savemat("/data/pt_02518/7T_Processing/PET_BPnd_Putamen_stats.mat", dic_mat)








