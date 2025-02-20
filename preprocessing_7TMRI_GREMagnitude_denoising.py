#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 12:54:12 2023

@author: gkotsoulias
"""

import os
import matplotlib.pyplot as plt
import ants
import nibabel as nib
import numpy as np

for i in range(14,22):
    if i<10:
        subject = 'GTSMp00'+str(i)
    else:
        subject = 'GTSMp0'+str(i)
    
    input_dir1 = '/data/pt_02518/7T_Processing/GTS_Cohort/'+subject+'/QSM2STDSpace/Magnitude_Phase/' 
    output_dir = '/data/pt_02518/7T_Processing/GTS_Cohort/'+subject+'/QSM2STDSpace/Magnitude_Phase/'
    
    print(input_dir1)
    print(output_dir)
    
    GRE_Magn_1 = ants.image_read(os.path.join(input_dir1,'Magn_Echo_0001.nii.gz'))
    GRE_Magn_2 = ants.image_read(os.path.join(input_dir1,'Magn_Echo_0002.nii.gz'))
    GRE_Magn_3 = ants.image_read(os.path.join(input_dir1,'Magn_Echo_0003.nii.gz'))
    GRE_Magn_4 = ants.image_read(os.path.join(input_dir1,'Magn_Echo_0004.nii.gz'))
    GRE_Magn_5 = ants.image_read(os.path.join(input_dir1,'Magn_Echo_0005.nii.gz'))
    GRE_Magn_6 = ants.image_read(os.path.join(input_dir1,'Magn_Echo_0006.nii.gz'))
    GRE_Magn_7 = ants.image_read(os.path.join(input_dir1,'Magn_Echo_0007.nii.gz'))
    GRE_Magn_8 = ants.image_read(os.path.join(input_dir1,'Magn_Echo_0008.nii.gz'))
    GRE_Magn_9 = ants.image_read(os.path.join(input_dir1,'Magn_Echo_0009.nii.gz'))
    
    in_vol = nib.load(os.path.join(input_dir1,'Magn_Echo_0001.nii.gz'))
    in_vol_data = in_vol.get_fdata()
    h = in_vol.header
    
    print('STEP 2: Corrections...............................')
    GRE_Magn_biasCorr_1 = ants.n4_bias_field_correction(GRE_Magn_1, shrink_factor=3, convergence={'iters': [50, 50, 50, 50],'tol': 1e-08})                                              
    GRE_Magn_biasCorr_denoised_1 = ants.denoise_image(GRE_Magn_biasCorr_1)

    GRE_Magn_biasCorr_2 = ants.n4_bias_field_correction(GRE_Magn_2, shrink_factor=3, convergence={'iters': [50, 50, 50, 50],'tol': 1e-08})                                              
    GRE_Magn_biasCorr_denoised_2 = ants.denoise_image(GRE_Magn_biasCorr_2)
    
    GRE_Magn_biasCorr_3 = ants.n4_bias_field_correction(GRE_Magn_3, shrink_factor=3, convergence={'iters': [50, 50, 50, 50],'tol': 1e-08})                                              
    GRE_Magn_biasCorr_denoised_3 = ants.denoise_image(GRE_Magn_biasCorr_3)   
    
    GRE_Magn_biasCorr_4 = ants.n4_bias_field_correction(GRE_Magn_4, shrink_factor=3, convergence={'iters': [50, 50, 50, 50],'tol': 1e-08})                                              
    GRE_Magn_biasCorr_denoised_4 = ants.denoise_image(GRE_Magn_biasCorr_4)  
    
    GRE_Magn_biasCorr_5 = ants.n4_bias_field_correction(GRE_Magn_5, shrink_factor=3, convergence={'iters': [50, 50, 50, 50],'tol': 1e-08})                                              
    GRE_Magn_biasCorr_denoised_5 = ants.denoise_image(GRE_Magn_biasCorr_5)
    
    GRE_Magn_biasCorr_6 = ants.n4_bias_field_correction(GRE_Magn_6, shrink_factor=3, convergence={'iters': [50, 50, 50, 50],'tol': 1e-08})                                              
    GRE_Magn_biasCorr_denoised_6 = ants.denoise_image(GRE_Magn_biasCorr_6)   
    
    GRE_Magn_biasCorr_7 = ants.n4_bias_field_correction(GRE_Magn_7, shrink_factor=3, convergence={'iters': [50, 50, 50, 50],'tol': 1e-08})                                              
    GRE_Magn_biasCorr_denoised_7 = ants.denoise_image(GRE_Magn_biasCorr_7)  
    
    GRE_Magn_biasCorr_8 = ants.n4_bias_field_correction(GRE_Magn_8, shrink_factor=3, convergence={'iters': [50, 50, 50, 50],'tol': 1e-08})                                              
    GRE_Magn_biasCorr_denoised_8 = ants.denoise_image(GRE_Magn_biasCorr_8)
    
    GRE_Magn_biasCorr_9 = ants.n4_bias_field_correction(GRE_Magn_9, shrink_factor=3, convergence={'iters': [50, 50, 50, 50],'tol': 1e-08})                                              
    GRE_Magn_biasCorr_denoised_9 = ants.denoise_image(GRE_Magn_biasCorr_9)   

   
    h['qform_code'] = 1
    h['sform_code'] = 1
    #sform is the basis of the file reading
    h.set_sform([[1, 0, 0, 0],
                 [0, 1, 0, 0],
                 [0, 0, 1, 0]])

    h['descrip'] = 'Correctly readable -sform based: CHECKED'
    #Correcting pixel dimensions
    h['pixdim'] = [1, 0.8, 0.8, 0.8, 0 ,0, 0, 0]
    #Orientation will be fixed and equires flip and z axis 90 deg rotation
    
    h['quatern_b'] = 0.0
    h['quatern_c'] = 0.0
    h['quatern_d'] = 0.0
    h['qoffset_x'] = 0.0
    h['qoffset_y'] = 0.0
    h['qoffset_z'] = 0.0
    
    out_vol = nib.nifti1.Nifti1Image(GRE_Magn_biasCorr_denoised_1.numpy(),np.eye(4) , header=h)
    out_vol.to_filename(output_dir+'Magn_Echo_0001_biasCorr_Denoised.nii.gz')  
    
    out_vol = nib.nifti1.Nifti1Image(GRE_Magn_biasCorr_denoised_2.numpy(),np.eye(4) , header=h)
    out_vol.to_filename(output_dir+'Magn_Echo_0002_biasCorr_Denoised.nii.gz')   
    
    out_vol = nib.nifti1.Nifti1Image(GRE_Magn_biasCorr_denoised_3.numpy(),np.eye(4) , header=h)
    out_vol.to_filename(output_dir+'Magn_Echo_0003_biasCorr_Denoised.nii.gz')  
    
    out_vol = nib.nifti1.Nifti1Image(GRE_Magn_biasCorr_denoised_4.numpy(),np.eye(4) , header=h)
    out_vol.to_filename(output_dir+'Magn_Echo_0004_biasCorr_Denoised.nii.gz')   
    
    out_vol = nib.nifti1.Nifti1Image(GRE_Magn_biasCorr_denoised_5.numpy(),np.eye(4) , header=h)
    out_vol.to_filename(output_dir+'Magn_Echo_0005_biasCorr_Denoised.nii.gz')  
    
    out_vol = nib.nifti1.Nifti1Image(GRE_Magn_biasCorr_denoised_6.numpy(),np.eye(4) , header=h)
    out_vol.to_filename(output_dir+'Magn_Echo_0006_biasCorr_Denoised.nii.gz')   
    
    out_vol = nib.nifti1.Nifti1Image(GRE_Magn_biasCorr_denoised_7.numpy(),np.eye(4) , header=h)
    out_vol.to_filename(output_dir+'Magn_Echo_0007_biasCorr_Denoised.nii.gz')  
    
    out_vol = nib.nifti1.Nifti1Image(GRE_Magn_biasCorr_denoised_8.numpy(),np.eye(4) , header=h)
    out_vol.to_filename(output_dir+'Magn_Echo_0008_biasCorr_Denoised.nii.gz')   
    
    out_vol = nib.nifti1.Nifti1Image(GRE_Magn_biasCorr_denoised_9.numpy(),np.eye(4) , header=h)
    out_vol.to_filename(output_dir+'Magn_Echo_0009_biasCorr_Denoised.nii.gz') 


    