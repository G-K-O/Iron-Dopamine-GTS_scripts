#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  3 17:33:49 2023

@author: gkotsoulias
"""

import os
import ants
import numpy as np
from scipy.io import savemat
import nibabel as nib
import matplotlib.pyplot as plt
import seaborn as sb
import ants
import SimpleITK as sitk

Controls = np.zeros([4,41])
Patients = np.zeros([4,26])

"""
Calculate mean of ventricles and full brain for each subject---------------------------------------
"""

j=0
for i in range(1,41):
    if i<10:
        subject = 'GTSMc00'+str(i)
    else:
        subject = 'GTSMc0'+str(i)
        
    dir1 = '/---/'+subject+'/---/'
    dir2 = '/---/'+subject+'/---/'#full brain mask dir
    dir3 = '/---/'+subject+'/---/ManualSegmentations/'
    
    if os.path.exists(dir1): 
       print(dir1)      
       print('ADDED')      

        ## Load volumes------------------------------------------------------------
       QSM1_dir = os.path.join(dir1,'QSM_iLSQR_mean_5echoes.nii.gz')
       ROIs_dir = os.path.join(dir2,'mask.nii.gz')
       ROIs2_dir = os.path.join(dir3,'ROI_Ventricles_Manual.nii.gz')
        
       in_vol1 = nib.load(QSM1_dir)
       QSM = in_vol1.get_fdata()
       
       in_vol2 = nib.load(ROIs_dir)      
       ROIs = in_vol2.get_fdata()
       
       in_vol3 = nib.load(ROIs2_dir)
       ROIs2 = in_vol3.get_fdata()
       

       Controls[0,j] = np.mean(QSM[ROIs!= 0])
       Controls[1,j] = QSM[ROIs!= 0].std()          
       Controls[2,j] = np.mean(QSM[ROIs2!= 0])
       Controls[3,j] = QSM[ROIs2!= 0].std()         
       j=j+1
    else:
       print(dir1)      
       print('EXCLUDED')   
    

j=0;
for i in range(1,26):
    if i<10:
        subject = 'GTSMp00'+str(i)
    else:
        subject = 'GTSMp0'+str(i)
        
    dir1 = '/data/pt_02518/7T_Processing/GTS_Cohort/'+subject+'/QSM2STDSpace/'
    dir2 = '/data/pt_02518/7T_Processing/GTS_Cohort/'+subject+'/QSM2STDSpace/'#full brain mask dir
    dir3 = '/data/pt_02518/7T_Processing/GTS_Cohort/'+subject+'/QSM2STDSpace/ManualSegmentations/'
    
    if os.path.exists(dir1): 
       print(dir1)      
       print('ADDED')      

        ## Load volumes------------------------------------------------------------
       QSM1_dir = os.path.join(dir1,'QSM_iLSQR_mean_5echoes.nii.gz')
       ROIs_dir = os.path.join(dir2,'mask.nii.gz')
       ROIs2_dir = os.path.join(dir3,'ROI_Ventricles_Manual.nii.gz')
        
       in_vol1 = nib.load(QSM1_dir)
       QSM = in_vol1.get_fdata()
       
       in_vol2 = nib.load(ROIs_dir)      
       ROIs = in_vol2.get_fdata()
       
       in_vol3 = nib.load(ROIs2_dir)
       ROIs2 = in_vol3.get_fdata()
       
      
       Patients[0,j] = np.mean(QSM[ROIs!= 0])
       Patients[1,j] = QSM[ROIs!= 0].std()          
       Patients[2,j] = np.mean(QSM[ROIs2!= 0])
       Patients[3,j] = QSM[ROIs2!= 0].std()       
       j=j+1
    else:
       print(dir1)      
       print('EXCLUDED')
       
"""
Calculate and save referenced images with custom header, to avoid issues with orientation
Here, only the reference to ventricles
"""
for i in range(1,41):
    if i<10:
       subject = 'GTSMc00'+str(i)
    else:
        subject = 'GTSMc0'+str(i)
    
    input_dir =  '/data/pt_02518/7T_Processing/Controls_Cohort/'+subject+'/QSM2STDSpace/'
    output_dir = '/data/pt_02518/7T_Processing/Controls_Cohort/'+subject+'/QSM2STDSpace/'
    
    filename_in = 'QSM_iLSQR_mean_5echoes.nii.gz'
    filename_out1= 'QSM_iLSQR_mean_5echoes_REFtoVentricles.nii.gz'
        
    path = os.path.join(input_dir, filename_in)
    out_path1 = os.path.join(output_dir, filename_out1)
        
    in_vol = nib.load(path)
    in_vol_data = in_vol.get_fdata()
    h = in_vol.header
       
        
    h['qform_code'] = 1
    h['sform_code'] = 1
    h.set_sform([[1, 0, 0, 0],
                     [0, 1, 0, 0],
                     [0, 0, 1, 0]])
    h['descrip'] = 'Correctly readable -sform based: CHECKED'
    h['pixdim'] = [1, 0.8, 0.8, 0.8, 0 ,0, 0, 0]
    N = np.shape(in_vol_data)
    h['dim'] = [3, N[0] ,N[1]  ,N[2] ,1 , 1, 0, 0]
    h['quatern_b'] = 0.0
    h['quatern_c'] = 0.0
    h['quatern_d'] = 0.0
    h['qoffset_x'] = 0.0
    h['qoffset_y'] = 0.0
    h['qoffset_z'] = 0.0
    print('----------------------------------------------------')
    print(subject+' normalized by '+str(np.abs(Controls[2,i-1])))
    
    out_vol1 = nib.nifti1.Nifti1Image(in_vol_data-Controls[2,i-1],np.eye(4) , header=h)
    out_vol1.to_filename(out_path1)            
    print('----------------------------------------------------')



for i in range(1,26):
    if i<10:
       subject = 'GTSMp00'+str(i)
    else:
        subject = 'GTSMp0'+str(i)
    
    input_dir =  '/data/pt_02518/7T_Processing/GTS_Cohort/'+subject+'/QSM2STDSpace/'
    output_dir = '/data/pt_02518/7T_Processing/GTS_Cohort/'+subject+'/QSM2STDSpace/'
    
    filename_in = 'QSM_iLSQR_mean_5echoes.nii.gz'
    filename_out1= 'QSM_iLSQR_mean_5echoes_REFtoVentricles.nii.gz'
        
    path = os.path.join(input_dir, filename_in)
    out_path1 = os.path.join(output_dir, filename_out1)
        
    in_vol = nib.load(path)
    in_vol_data = in_vol.get_fdata()
    h = in_vol.header
        
    h['qform_code'] = 1
    h['sform_code'] = 1
    h.set_sform([[1, 0, 0, 0],
                     [0, 1, 0, 0],
                     [0, 0, 1, 0]])
    h['descrip'] = 'Correctly readable -sform based: CHECKED'
    h['pixdim'] = [1, 0.8, 0.8, 0.8, 0 ,0, 0, 0]
    N = np.shape(in_vol_data)
    h['dim'] = [3, N[0] ,N[1]  ,N[2] ,1 , 1, 0, 0]
    h['quatern_b'] = 0.0
    h['quatern_c'] = 0.0
    h['quatern_d'] = 0.0
    h['qoffset_x'] = 0.0
    h['qoffset_y'] = 0.0
    h['qoffset_z'] = 0.0
    print('----------------------------------------------------')
    print(subject+' normalized by '+str(np.abs(Patients[2,i-1])))
    
    out_vol1 = nib.nifti1.Nifti1Image(in_vol_data-Patients[2,i-1],np.eye(4) , header=h)
    out_vol1.to_filename(out_path1)     
    print('----------------------------------------------------')

    
