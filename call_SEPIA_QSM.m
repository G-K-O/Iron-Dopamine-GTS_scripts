function  call_SEPIA_QSM(outfolder,headerFolder,echonum)
% Calling Non-linear Dipole Inversion
    sepia_addpath;
    % Input/Output filenames
    input(1).name = [outfolder 'tmp_LocalPhase_echo_' num2str(echonum) '.nii'] ;
    input(2).name = [outfolder 'tmp_Magn_echo_' num2str(echonum) '.nii'] ;
    input(3).name = '' ;
    input(4).name = headerFolder ;
    output_basename = [outfolder 'QSM_iLSQR_echo_' num2str(echonum)];
    mask_filename = [outfolder 'mask.nii.gz'] ;

    % General algorithm parameters
    algorParam = struct();
    algorParam.general.isBET = 0 ;
    algorParam.general.isInvert = 0 ;
    algorParam.general.isRefineBrainMask = 0 ;
    % QSM algorithm parameters
    algorParam.qsm.reference_tissue = 'Brain mask' ;
    algorParam.qsm.method = 'iLSQR STI Suite' ;
    algorParam.qsm.tol = 1 ;
    algorParam.qsm.maxiter = 200 ;
    algorParam.qsm.stepSize = 1 ;

    sepiaIO(input,output_basename,mask_filename,algorParam);
end