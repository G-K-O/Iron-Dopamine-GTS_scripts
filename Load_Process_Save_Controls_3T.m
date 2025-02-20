cd /---/---/In-vivo_data_Analysis
addpath(genpath('/---/---/sepia-master/'))  
addpath(genpath('/---/---/STISuite_V3.0/'))  
addpath('/---/---/matScripts')
%[infolder_list,outfolder_list] = Get_IO_Directory_Lists_3T_PETMR('Controls');

voxelsize = [0.8,0.8,0.8]; 
N = [256, 256,192,2];
TE_ms = [18.5, 25];

for i = 20:length(infolder_list(:,1))
    infolder = infolder_list(i,:); 
    outfolder = outfolder_list(i,:); 
    disp(['IN-Folder: ' infolder])
    disp(['OUT-Folder: ' outfolder])

    magnCombined_allEchoes = zeros(N);
    phaseCombined_allEchoes = zeros(N);
    tissuePhase_allEchoes = zeros(N);
    QSM_NDI_allEchoes = zeros(N);
    QSM_iLSQR_allEchoes = zeros(N);
    slice = 1;
    indecho = 1;
    for indecho = 1:N(4)
        disp(['Echo:' num2str(indecho-1, '%04d')])
        clear mag phase datastore 
        for slice = 1:N(3)
            mag(:,:,:, slice) = niftiread([infolder '/uncombined_channels_echo_' num2str(indecho-1,'%04d') '_partition_' num2str(slice-1, '%04d') '_magnitude.nii']);
            phase(:,:,:,slice) = niftiread([infolder '/uncombined_channels_echo_' num2str(indecho-1,'%04d') '_partition_' num2str(slice-1, '%04d')  '_phase.nii']); 
        end        
        data = permute(mag.*exp(1j*phase), [1,3,4,2]);
        [magnCombined_allEchoes(:,:,:,indecho),phaseCombined_allEchoes(:,:,:, indecho)] = Combine_Coils_JC(data, voxelsize); 
    end

   %% Fix the issue of echo spatial drift
   [optimizer,metric] = imregconfig("monomodal")
    optimizer.MaximumIterations = 300;
    % Find translation between echoes based on magnitude
    tform = imregtform(magnCombined_allEchoes(:,:,:,2),magnCombined_allEchoes(:,:,:,1),"translation",optimizer,metric);
   
    % Apply transform to magn_echo2 and phase_echo2
    mag2_reg = imwarp(magnCombined_allEchoes(:,:,:,2),tform,"OutputView",imref3d(size(magnCombined_allEchoes(:,:,:,1))));
    phase2_reg = imwarp(phaseCombined_allEchoes(:,:,:,2),tform,"OutputView",imref3d(size(magnCombined_allEchoes(:,:,:,1))));
    %imshowpair(magnCombined_allEchoes(:,:,100,1),mag2_reg(:,:,100),"Scaling","joint")
   
    % Fix the matrices before saving
    magnCombined_allEchoes(:,:,:,2) = mag2_reg;
    phaseCombined_allEchoes(:,:,:,2) = phase2_reg;

    save([outfolder '/MagnCombined_allEchoes.mat'],'magnCombined_allEchoes');
    save([outfolder '/PhaseCombined_allEchoes.mat'],'phaseCombined_allEchoes');
    niftiwrite(magnCombined_allEchoes(:,:,:,1),[outfolder 'tmp_Magn_echo_1.nii']);
    niftiwrite(magnCombined_allEchoes(:,:,:,2),[outfolder 'tmp_Magn_echo_2.nii']);

    [mask] = Get_FSL_mask(outfolder);

    indecho = 1;
    for indecho = 1:2
        disp(['Processing for local phase of echo:' num2str(indecho-1, '%04d')])
        tissuePhase_allEchoes(:,:,:,indecho) = V_SHARP(single(phaseCombined_allEchoes(:,:,:,indecho)),mask, 'smvsize',20, 'voxelsize', voxelsize);
    end
    save([outfolder '/LocalPhase_allEchoes.mat'],'tissuePhase_allEchoes');
    niftiwrite(tissuePhase_allEchoes(:,:,:,1),[outfolder 'tmp_LocalPhase_echo_1.nii']);
    niftiwrite(tissuePhase_allEchoes(:,:,:,2),[outfolder 'tmp_LocalPhase_echo_2.nii']);

    indecho = 1;
    for indecho = 1:2
        call_SEPIA_QSM(outfolder,['/data/u_gkotsoulias_software/matScripts/In-vivo_data_Analysis/3T_echo_' num2str(indecho) '_NMR134_SEPIA_HEADER.mat' ],indecho);
        QSM_iLSQR_allEchoes(:,:,:,indecho) = QSM_iLSQR(single(tissuePhase_allEchoes(:,:,:,indecho)), mask,'TE',TE_ms(indecho),'B0',3,'H',[0 0 1],'padsize',[12 12 12],'voxelsize',voxelsize);
        QSM_NDI_allEchoes(:,:,:,indecho) = niftiread([outfolder 'QSM_NDI_echo_' num2str(indecho) '_Chimap.nii.gz']);
    end
    save([outfolder '/QSM_NDI_allEchoes.mat'],'QSM_NDI_allEchoes');
    save([outfolder '/QSM_iLSQR_allEchoes.mat'],'QSM_iLSQR_allEchoes');
    niftiwrite(mean(QSM_NDI_allEchoes, 4),[outfolder 'QSM_NDI_mean.nii.gz']);
    niftiwrite(mean(QSM_iLSQR_allEchoes, 4),[outfolder 'QSM_iLSQR_mean.nii.gz']);
    niftiwrite(QSM_iLSQR_allEchoes(:,:,:,1),[outfolder 'QSM_iLSQR_echo_1_Chimap.nii']);
    niftiwrite(QSM_iLSQR_allEchoes(:,:,:,2),[outfolder 'QSM_iLSQR_echo_2_Chimap.nii']);

    %system(['rm -f ' outfolder 'tmp_*']);
    system(['rm -f ' outfolder 'run*']);  
    system(['rm -f ' outfolder 'sepia*']);
    system(['rm -f ' outfolder '*referenceregion.nii.gz']);
end







