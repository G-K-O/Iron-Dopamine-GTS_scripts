#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 14:18:41 2023

@author: gkotsoulias
"""

import os
import matplotlib.pyplot as plt
import ants
import nibabel as nib
import SimpleITK as sitk
import numpy as np


for i in range(1,41):
    if i<10:
       subject = 'GTSMc00'+str(i)
    else:
        subject = 'GTSMc0'+str(i)
    
    input_dir = '/---/'+subject+'/QSM/'
    output_dir = '/---/'+subject+'/QSM2STDSpace/'
    print('Getting data from: ', input_dir)
    print('Saving data to: ',output_dir)
    
    im_name = ['QSM_NDI_mean_5echoes.nii.nii',
               'QSM_NDI_echo_3_Chimap.nii.gz',
               'QSM_NDI_echo_4_Chimap.nii.gz',
               'mask.nii.gz',
               'QSM_iLSQR_echo_4_Chimap.nii',
               'QSM_iLSQR_echo_3_Chimap.nii',
               'QSM_iLSQR_mean_5echoes.nii.nii']

    out_name = ['QSM_NDI_mean_5echoes.nii.gz',
                'QSM_NDI_echo_3.nii.gz',
                'QSM_NDI_echo_4.nii.gz',
                'mask.nii.gz',
                'QSM_iLSQR_echo_4.nii.gz',
                'QSM_iLSQR_echo_3.nii.gz',
                'QSM_iLSQR_mean_5echoes.nii.gz']

    for i in range(7):
        filename_in = im_name[i]
        filename_out = out_name[i]
        
        path = os.path.join(input_dir, filename_in)
        out_path = os.path.join(output_dir, filename_out)
        print(path)
        print(out_path)
        
        in_vol = nib.load(path)
        in_vol_data = in_vol.get_fdata()
        h = in_vol.header
        #print(h)
        
        
        ## Making essential changes to the header
        #qform and sform scanner - multiple softwares require the transforms to 
        #be the same.
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
        N = np.shape(in_vol_data)
        h['dim'] = [3, N[1] ,N[0]  ,N[2] ,1 , 1, 0, 0]
        
        h['quatern_b'] = 0.0
        h['quatern_c'] = 0.0
        h['quatern_d'] = 0.0
        h['qoffset_x'] = 0.0
        h['qoffset_y'] = 0.0
        h['qoffset_z'] = 0.0
        #print(h)
        out_vol = nib.nifti1.Nifti1Image(np.fliplr(np.rot90(in_vol_data, k=1, axes=(0, 1))),np.eye(4) , header=h)
        out_vol.to_filename(out_path)  
        
    print('-------------------------------------------------------------------------------------------------')

        

























