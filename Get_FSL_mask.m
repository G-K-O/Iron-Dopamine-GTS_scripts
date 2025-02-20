function [mask] = Get_FSL_mask(outfolder)
    command_FSL_mask = ['/afs/cbs.mpg.de/software/fsl/6.0.3/debian-bullseye-amd64/bin/bet ' outfolder 'tmp_Magn_echo_1.nii ' outfolder ' -R -f 0.25 -g 0 -m'];
    command_rm_files = ['rm -f ' outfolder '.nii.gz'];
    command_rn_files = ['mv ' outfolder '_mask.nii.gz '  outfolder 'mask.nii.gz'];
    system(command_FSL_mask); 
    system(command_rm_files); 
    system(command_rn_files);
    mask = single(niftiread([outfolder 'mask.nii.gz']));
end