#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 25 16:09:40 2023

@author: gkotsoulias
"""


import mat73
import scipy.io
import ants
import seaborn as sb
from matplotlib import pyplot as plt
import pandas as pd
import mat73
import scipy.stats as scs
import os
import nibabel as nib
import numpy as np

"""
Manually connect the PET ,measurement codes with the 7T MRI codes for controls and patients
"""
controls_PET_codes = ['GTSPc001','GTSPc002','GTSPc003','GTSPc004','GTSPc005','GTSPc007','GTSPc008','GTSPc009','GTSPc010','GTSPc011','GTSPc012','GTSPc013','GTSPc015','GTSPc016','GTSPc017','GTSPc018','GTSPc020']
controls_7T_codes =  ['GTSMc013','GTSMc025','GTSMc009','GTSMc017','GTSMc019','GTSMc035','GTSMc022','GTSMc034','GTSMc040','GTSMc033','GTSMc029','GTSMc036','GTSMc018','GTSMc015','GTSMc038','GTSMc039','GTSMc037']

patients_PET_codes = ['GTSPp001','GTSPp002','GTSPp003','GTSPp005','GTSPp006','GTSPp007','GTSPp008','GTSPp009','GTSPp010','GTSPp011','GTSPp012','GTSPp013','GTSPp014','GTSPp015','GTSPp016','GTSPp018','GTSPp019','GTSPp020']
patients_7T_codes =  ['GTSMp004','GTSMp006','GTSMp008','GTSMp012','GTSMp011','GTSMp013','GTSMp009','GTSMp014','GTSMp007','GTSMp015','GTSMp016','GTSMp018','GTSMp020','GTSMp021','GTSMp017','GTSMp019','GTSMp022','GTSMp025']


for i in range(len(controls_PET_codes)):
    print(controls_PET_codes[i] ,' is pairing with ', controls_7T_codes[i]) 
    input_dir1 = '/---/---/'+controls_PET_codes[i]+'/T1_PET2QSMSTDSpace_Reg/'
    input_dir2 = '/---/---/'+controls_7T_codes[i]+'/T12QSMSTDSpace_Reg/'
    output_dir = '/---/---/'+controls_7T_codes[i]+'/T12QSMSTDSpace_Reg/'

    UNI_7T_dir = os.path.join(input_dir2,'MP2RAGE_UNI_Masked_nonLinReg.nii.gz')
    UNI_3T_dir = os.path.join(input_dir1,'MP2RAGE_UNI_reg.nii.gz')
    BPnd_dir =   os.path.join(input_dir1,'BPnd_reg.nii.gz')
    R1_dir =     os.path.join(input_dir1,'R1_reg.nii.gz')
    
    print(UNI_7T_dir)
    print(UNI_3T_dir)
    print(BPnd_dir)
    print(R1_dir) 
    
    print('STEP 1: Loading...................................')    
    UNI_7T  = ants.image_read(UNI_7T_dir)
    UNI_3T = ants.image_read(UNI_3T_dir)    
    BPnd = ants.image_read(BPnd_dir)
    R1 = ants.image_read(R1_dir)
    
    in_vol = nib.load(UNI_7T_dir)
    in_vol_data = in_vol.get_fdata()
    h = in_vol.header
    
    print('STEP 2: Registration (Init affine).................')
    tx_filename = os.path.join(output_dir,'init_affine_Transform_PET.mat')
    txfile = ants.affine_initializer(UNI_7T, 
                                     UNI_3T, 
                                     search_factor=20, 
                                     radian_fraction=0.3, 
                                     use_principal_axis=False, 
                                     local_search_iterations=20, 
                                     mask=None, 
                                     txfn=tx_filename)
    tx = ants.read_transform(txfile, dimension=2)
        
    ## Apply transform to moving image---------------------------------------------
    init_trasf_UNI_3T =  tx.apply_to_image(UNI_3T, reference=UNI_7T, interpolation='linear')

    print('Post-initial affine comparison----------------------------------------')
    ants.plot(UNI_7T, init_trasf_UNI_3T, axis= 2,overlay_alpha = 0.4, nslices=56)
    
    
    print('STEP 3: Registration (Non-linear).................')
    ## Non-linear registration of images-------------------------------------------
    reg = ants.registration(fixed=UNI_7T, 
                                moving=init_trasf_UNI_3T, 
                                type_of_transform = 'SyN',
                                reg_iterations = [60,60,40] )
        
    UNI_3T_reg = reg['warpedmovout']
    print('Final registration comparisons-----------------------------------------')
    ants.plot(UNI_7T, UNI_3T_reg, axis= 2,overlay_alpha = 0.4, nslices=56)
    
    ## Apply transform to other images (masks from segmentation)-------------------
    fullTrasform = reg[ 'fwdtransforms'] 
    fullTrasform.append(txfile)
    BPnd_reg = ants.apply_transforms( fixed = UNI_7T,moving =  BPnd, transformlist = fullTrasform, 
                                      interpolator  = 'linear', whichtoinvert = [False,False,False])
    R1_reg = ants.apply_transforms( fixed = UNI_7T,moving =  R1, transformlist = fullTrasform, 
                                      interpolator  = 'linear', whichtoinvert = [False,False,False])     
    
    h['qform_code'] = 1
    h['sform_code'] = 1
    #sform is the basis of the file reading
    h.set_sform([[1, 0, 0, 0],
                 [0, 1, 0, 0],
                 [0, 0, 1, 0]])

    h['descrip'] = '9 echoes T2*'
    #Correcting pixel dimensions
    h['pixdim'] = [4, 0.8, 0.8, 0.8, 9 ,0, 0, 0]
    #Orientation will be fixed and equires flip and z axis 90 deg rotation
    
    h['quatern_b'] = 0.0
    h['quatern_c'] = 0.0
    h['quatern_d'] = 0.0
    h['qoffset_x'] = 0.0
    h['qoffset_y'] = 0.0
    h['qoffset_z'] = 0.0
    

    out_dir_1 = os.path.join(output_dir,'MP2RAGE_UNI_reg_3T.nii.gz')
    out_dir_2 = os.path.join(output_dir,'BPnd_reg_to7T.nii.gz')
    out_dir_3 = os.path.join(output_dir,'R1_reg_to7T.nii.gz')
      
    out_vol1= nib.nifti1.Nifti1Image(UNI_3T_reg.numpy(), np.eye(4) , header=h)
    out_vol1.to_filename(out_dir_1)  
    
    out_vol2= nib.nifti1.Nifti1Image(BPnd_reg.numpy(), np.eye(4) , header=h)
    out_vol2.to_filename(out_dir_2)  
    
    out_vol3= nib.nifti1.Nifti1Image(R1_reg.numpy(), np.eye(4) , header=h)
    out_vol3.to_filename(out_dir_3)  





for i in range(len(patients_PET_codes)):
    print(patients_PET_codes[i] ,' is pairing with ', patients_7T_codes[i])
    input_dir1 = '/---/---/'+patients_PET_codes[i]+'/T1_PET2QSMSTDSpace_Reg/'
    input_dir2 = '/---/---/'+patients_7T_codes[i]+'/T12QSMSTDSpace_Reg/'
    output_dir = '/---/---/'+patients_7T_codes[i]+'/T12QSMSTDSpace_Reg/'

    UNI_7T_dir = os.path.join(input_dir2,'MP2RAGE_UNI_Masked_nonLinReg.nii.gz')
    UNI_3T_dir = os.path.join(input_dir1,'MP2RAGE_UNI_reg.nii.gz')
    BPnd_dir =   os.path.join(input_dir1,'BPnd_reg.nii.gz')
    R1_dir =     os.path.join(input_dir1,'R1_reg.nii.gz')
    
    print(UNI_7T_dir)
    print(UNI_3T_dir)
    print(BPnd_dir)
    print(R1_dir) 
    
    print('STEP 1: Loading...................................')    
    UNI_7T  = ants.image_read(UNI_7T_dir)
    UNI_3T = ants.image_read(UNI_3T_dir)    
    BPnd = ants.image_read(BPnd_dir)
    R1 = ants.image_read(R1_dir)
    
    in_vol = nib.load(UNI_7T_dir)
    in_vol_data = in_vol.get_fdata()
    h = in_vol.header
    
    print('STEP 2: Registration (Init affine).................')
    tx_filename = os.path.join(output_dir,'init_affine_Transform_PET.mat')
    txfile = ants.affine_initializer(UNI_7T, 
                                     UNI_3T, 
                                     search_factor=20, 
                                     radian_fraction=0.3, 
                                     use_principal_axis=False, 
                                     local_search_iterations=20, 
                                     mask=None, 
                                     txfn=tx_filename)
    tx = ants.read_transform(txfile, dimension=2)
        
    ## Apply transform to moving image---------------------------------------------
    init_trasf_UNI_3T =  tx.apply_to_image(UNI_3T, reference=UNI_7T, interpolation='linear')

    print('Post-initial affine comparison----------------------------------------')
    ants.plot(UNI_7T, init_trasf_UNI_3T, axis= 2,overlay_alpha = 0.4, nslices=56)
    
    
    print('STEP 3: Registration (Non-linear).................')
    ## Non-linear registration of images-------------------------------------------
    reg = ants.registration(fixed=UNI_7T, 
                                moving=init_trasf_UNI_3T, 
                                type_of_transform = 'SyN',
                                reg_iterations = [60,60,40] )
        
    UNI_3T_reg = reg['warpedmovout']
    print('Final registration comparisons-----------------------------------------')
    ants.plot(UNI_7T, UNI_3T_reg, axis= 2,overlay_alpha = 0.4, nslices=56)
    
    ## Apply transform to other images (masks from segmentation)-------------------
    fullTrasform = reg[ 'fwdtransforms'] 
    fullTrasform.append(txfile)
    BPnd_reg = ants.apply_transforms( fixed = UNI_7T,moving =  BPnd, transformlist = fullTrasform, 
                                      interpolator  = 'linear', whichtoinvert = [False,False,False])
    R1_reg = ants.apply_transforms( fixed = UNI_7T,moving =  R1, transformlist = fullTrasform, 
                                      interpolator  = 'linear', whichtoinvert = [False,False,False])     
    
    h['qform_code'] = 1
    h['sform_code'] = 1
    #sform is the basis of the file reading
    h.set_sform([[1, 0, 0, 0],
                 [0, 1, 0, 0],
                 [0, 0, 1, 0]])

    h['descrip'] = '9 echoes T2*'
    #Correcting pixel dimensions
    h['pixdim'] = [4, 0.8, 0.8, 0.8, 9 ,0, 0, 0]
    #Orientation will be fixed and equires flip and z axis 90 deg rotation
    
    h['quatern_b'] = 0.0
    h['quatern_c'] = 0.0
    h['quatern_d'] = 0.0
    h['qoffset_x'] = 0.0
    h['qoffset_y'] = 0.0
    h['qoffset_z'] = 0.0
    

    out_dir_1 = os.path.join(output_dir,'MP2RAGE_UNI_reg_3T.nii.gz')
    out_dir_2 = os.path.join(output_dir,'BPnd_reg_to7T.nii.gz')
    out_dir_3 = os.path.join(output_dir,'R1_reg_to7T.nii.gz')
      
    out_vol1= nib.nifti1.Nifti1Image(UNI_3T_reg.numpy(), np.eye(4) , header=h)
    out_vol1.to_filename(out_dir_1)  
    
    out_vol2= nib.nifti1.Nifti1Image(BPnd_reg.numpy(), np.eye(4) , header=h)
    out_vol2.to_filename(out_dir_2)  
    
    out_vol3= nib.nifti1.Nifti1Image(R1_reg.numpy(), np.eye(4) , header=h)
    out_vol3.to_filename(out_dir_3)  

