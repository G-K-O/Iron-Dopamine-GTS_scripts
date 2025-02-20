#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 16:35:02 2023

@author: gkotsoulias
"""

import os
import matplotlib.pyplot as plt
import ants
import nibabel as nib
import numpy as np

for i in range(1,26):
    if i<10:
        subject = 'GTSMp00'+str(i)
    else:
        subject = 'GTSMp0'+str(i)
    
    input_dir1 = '/---/'+subject+'/QSM2STDSpace/Magnitude_Phase/' 
    input_dir2 = '/---/'+subject+'/QSM2STDSpace/' 
    input_dir3 = '/---/'+subject+'/T1/'
    output_dir = '/---/'+subject+'/T12QSMSTDSpace_Reg/'
    

    GRE_Magn_dir = os.path.join(input_dir1,'Magn_Echo_0001.nii.gz')
    GRE_mask_dir = os.path.join(input_dir2,'mask.nii.gz')
    MP2RAGE_UNI_dir =  os.path.join(input_dir3,'T1_T1_MP2RAGE_UNI_Masked.nii.gz')
    MP2RAGE_INV2_dir = os.path.join(input_dir3,'T1_T1_MP2RAGE_INV2_Masked.nii.gz')

    print('Inputs---------------------------------------------------------------------')
    print(GRE_Magn_dir)
    print(GRE_mask_dir)
    print(MP2RAGE_UNI_dir)
    print(MP2RAGE_INV2_dir)
    print('--------------------------------------------------------------------------')

    print('STEP 0: All images are already registered rigidly in the scanner anatomical space using ITKSNAP and headers are checked............................')
    print('STEP 1: Loading...................................')
    
    GRE_Magn_QSM  = ants.image_read(GRE_Magn_dir)
    GRE_mask  = ants.image_read(GRE_mask_dir)    
    MP2RAGE_UNI = ants.image_read(MP2RAGE_UNI_dir)
    MP2RAGE_INV2 = ants.image_read(MP2RAGE_INV2_dir)
    
    in_vol = nib.load(GRE_Magn_dir)
    in_vol_data = in_vol.get_fdata()
    h = in_vol.header
    
    print('STEP 2: Corrections...............................')
    GRE_Magn_biasCorr = ants.n4_bias_field_correction(GRE_Magn_QSM, shrink_factor=3, convergence={'iters': [50, 50, 50, 50],'tol': 1e-08})                                              
    GRE_Magn_biasCorr_denoised = ants.denoise_image(GRE_Magn_biasCorr)
    plt.figure()
    plt.imshow(GRE_Magn_biasCorr_denoised[:,100,:])
    
    MP2RAGE_UNI_biasCorr = ants.n4_bias_field_correction(MP2RAGE_UNI, shrink_factor=3, convergence={'iters': [50, 50, 50, 50],'tol': 1e-07})                                              
    MP2RAGE_UNI_biasCorr_denoised = ants.denoise_image(MP2RAGE_UNI_biasCorr)
    plt.figure()
    plt.imshow(MP2RAGE_UNI_biasCorr_denoised[:,100,:])
    
    MP2RAGE_INV2_biasCorr = ants.n4_bias_field_correction(MP2RAGE_INV2, shrink_factor=3, convergence={'iters': [50, 50, 50, 50],'tol': 1e-07})                                              
    MP2RAGE_INV2_biasCorr_denoised = ants.denoise_image(MP2RAGE_INV2_biasCorr)
    plt.figure()
    plt.imshow(MP2RAGE_INV2_biasCorr_denoised[:,100,:])
    
    
    
    print('STEP 3: Registration (Init affine).................')
    tx_filename = os.path.join(output_dir,'init_affine_Transform.mat')
    txfile = ants.affine_initializer(GRE_Magn_biasCorr_denoised, 
                                     MP2RAGE_INV2_biasCorr_denoised, 
                                     search_factor=20, 
                                     radian_fraction=0.3, 
                                     use_principal_axis=False, 
                                     local_search_iterations=20, 
                                     mask=None, 
                                     txfn=tx_filename)
    tx = ants.read_transform(txfile, dimension=2)
        
    ## Apply transform to moving image---------------------------------------------
    init_trasf_INV2 = tx.apply_to_image(MP2RAGE_INV2_biasCorr_denoised, reference=GRE_Magn_biasCorr_denoised, interpolation='linear')
    init_trasf_UNI =  tx.apply_to_image(MP2RAGE_UNI_biasCorr_denoised, reference=GRE_Magn_biasCorr_denoised, interpolation='linear')

    print('Post-initial affine comparison----------------------------------------')
    ants.plot( GRE_Magn_biasCorr_denoised, init_trasf_INV2, axis= 2,overlay_alpha = 0.4, nslices=56)
    
    
    print('STEP 4: Registration (Non-linear).................')
    ## Non-linear registration of images-------------------------------------------
    reg = ants.registration(fixed=GRE_Magn_biasCorr_denoised, 
                                moving=init_trasf_INV2, 
                                type_of_transform = 'SyN',
                                reg_iterations = [100,100,50] )
        
    MP2RAGE_INV2_reg = reg['warpedmovout']
    print('Final registration comparisons-----------------------------------------')
    ants.plot( GRE_Magn_biasCorr_denoised, MP2RAGE_INV2_reg, axis= 2,overlay_alpha = 0.4, nslices=56)
    
    ## Apply transform to other images (masks from segmentation)-------------------
    fullTrasform = reg[ 'fwdtransforms'] 
    fullTrasform.append(txfile)
    MP2RAGE_UNI_biasCorr_denoised_reg = ants.apply_transforms( fixed = GRE_Magn_biasCorr_denoised,
                                                              moving = MP2RAGE_UNI_biasCorr_denoised, transformlist = fullTrasform, 
                                                              interpolator  = 'linear', whichtoinvert = [False,False,False])
        
    h['qform_code'] = 1
    h['sform_code'] = 1
    #sform is the basis of the file reading
    h.set_sform([[1, 0, 0, 0],
                 [0, 1, 0, 0],
                 [0, 0, 1, 0]])

    h['descrip'] = 'Correctly readable -sform based: CHECKED'
    #Correcting pixel dimensions
    h['pixdim'] = [1, 0.8, 0.8, 0.8, 0 ,0, 0, 0]
    
    h['quatern_b'] = 0.0
    h['quatern_c'] = 0.0
    h['quatern_d'] = 0.0
    h['qoffset_x'] = 0.0
    h['qoffset_y'] = 0.0
    h['qoffset_z'] = 0.0
    
    out_vol = nib.nifti1.Nifti1Image(MP2RAGE_UNI_biasCorr_denoised_reg.numpy(),np.eye(4) , header=h)
    out_vol.to_filename(output_dir+'MP2RAGE_UNI_Masked_nonLinReg.nii.gz')  
    
    out_vol = nib.nifti1.Nifti1Image(MP2RAGE_INV2_reg.numpy(),np.eye(4) , header=h)
    out_vol.to_filename(output_dir+'MP2RAGE_INV2_Masked_nonLinReg.nii.gz')   
    
    out_vol = nib.nifti1.Nifti1Image(init_trasf_UNI.numpy(),np.eye(4) , header=h)
    out_vol.to_filename(output_dir+'MP2RAGE_UNI_Masked_AffineReg.nii.gz')  
    
    out_vol = nib.nifti1.Nifti1Image(init_trasf_INV2.numpy(),np.eye(4) , header=h)
    out_vol.to_filename(output_dir+'MP2RAGE_INV2_Masked_AffineReg.nii.gz')      


    
    
    
    












