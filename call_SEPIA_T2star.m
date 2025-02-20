function  call_SEPIA_T2star(infolder,outfolder,mask_folder, algo)
        %% 9 echoes 
        % Input/Output filenames
        % algo = 'Trapezoidal' or 'ARLO'
        input = struct();
        input(1).name = '' ;
        input(2).name =  [infolder 'GRE_Magn_9echoes.nii.gz']  ;
        input(3).name = '' ;
        input(4).name = '/---/---/R2star_SepiaHeader.mat' ;
        output_basename = [outfolder '9echoes_'];
        mask_filename = [mask_folder 'mask.nii.gz'] ;

        % General algorithm parameters
        algorParam = struct();
        algorParam.general.isBET = 0 ;
        algorParam.general.isInvert = 0 ;
        algorParam.general.isRefineBrainMask = 0 ;
        % R2* algorithm parameters
        algorParam.r2s.method =  algo;
        algorParam.r2s.s0mode = '1st echo' ;
        sepiaIO(input,output_basename,mask_filename,algorParam);  


        %% 7 echoes 
        % Input/Output filenames
        input = struct();
        input(1).name = '' ;
        input(2).name =  [infolder 'GRE_Magn_7echoes.nii.gz']  ;
        input(3).name = '' ;
        input(4).name = '/data/p_02518/7T_Measurements/R2star_SepiaHeader_7echoes.mat' ;
        output_basename = [outfolder '7echoes_'];
        mask_filename = [mask_folder 'mask.nii.gz'] ;

        % General algorithm parameters
        algorParam = struct();
        algorParam.general.isBET = 0 ;
        algorParam.general.isInvert = 0 ;
        algorParam.general.isRefineBrainMask = 0 ;
        % R2* algorithm parameters
        algorParam.r2s.method = algo ;
        algorParam.r2s.s0mode = '1st echo' ;
        sepiaIO(input,output_basename,mask_filename,algorParam); 


        %% 6 echoes 
        % Input/Output filenames
        input = struct();
        input(1).name = '' ;
        input(2).name =  [infolder 'GRE_Magn_6echoes.nii.gz']  ;
        input(3).name = '' ;
        input(4).name = '/data/p_02518/7T_Measurements/R2star_SepiaHeader_6echoes.mat' ;
        output_basename = [outfolder '6echoes_Trap'];
        mask_filename = [mask_folder 'mask.nii.gz'] ;

        % General algorithm parameters
        algorParam = struct();
        algorParam.general.isBET = 0 ;
        algorParam.general.isInvert = 0 ;
        algorParam.general.isRefineBrainMask = 0 ;
        % R2* algorithm parameters
        algorParam.r2s.method = algo;
        algorParam.r2s.s0mode = '1st echo' ;
        sepiaIO(input,output_basename,mask_filename,algorParam); 

end