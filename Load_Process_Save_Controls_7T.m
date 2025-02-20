cd /---/---/In-vivo_data_Analysis
addpath(genpath('/---/---/sepia-master/'))  
addpath(genpath('/---/---/STISuite_V3.0/'))  
addpath('/---/---/matlabScripts')
%[infolder_list,outfolder_list] = Get_IO_Directory_Lists('Controls');

voxelsize = [0.8,0.8,0.8]; 
N = [266, 266,176,9];
TE_ms = [5,9.1,13.2,17.3,21.4,25.5, 29.6, 33.7, 37.8];


for i = 39:length(infolder_list(:,1))
    infolder = infolder_list(i,:); 
    outfolder = outfolder_list(i,:); 
    disp(['IN-Folder: ' infolder])
    disp(['OUT-Folder: ' outfolder])

    magnCombined_allEchoes = zeros(N);
    phaseCombined_allEchoes = zeros(N);
    tissuePhase_5Echoes = zeros(N(1),N(2),N(3),5);
    QSM_NDI_5Echoes = zeros(N(1),N(2),N(3),5);
    QSM_iLSQR_5Echoes = zeros(N(1),N(2),N(3),5);
    slice = 1;
    indecho = 1;
    tic
    for indecho = 1:N(4)
        disp(['Echo:' num2str(indecho-1, '%04d')])
        clear mag phase datastore 
        for slice = 1:N(3)
            mag(:,:,:, slice) = niftiread([infolder '/uncombined_partition_' num2str(slice-1,'%04d') '_echo_' num2str(indecho-1, '%04d') '_magnitude.nii']);
            phase(:,:,:,slice) = niftiread([infolder '/uncombined_partition_' num2str(slice-1,'%04d') '_echo_' num2str(indecho-1, '%04d') '_phase.nii']); 
        end

        data = permute(mag.*exp(1j*phase), [1,3,4,2]);
        [magnCombined_allEchoes(:,:,:,indecho),phaseCombined_allEchoes(:,:,:, indecho)] = Combine_Coils_JC(data, voxelsize); 
        %USE phaseCombined_allEchoes(:,:,:, indecho) = MRPhaseUnwrap(permute(phase,[1,3,4,2]),'voxelsize',voxelsize,'padsize',padsize); 
        %IN CASE THE UNCOMBINED CHANNELS WERE NOT RETRIEVED
    end
    save([outfolder '/MagnCombined_allEchoes.mat'],'magnCombined_allEchoes');
    save([outfolder '/PhaseCombined_allEchoes.mat'],'phaseCombined_allEchoes');
    niftiwrite(magnCombined_allEchoes(:,:,:,1),[outfolder 'tmp_Magn_echo_1.nii']);
    niftiwrite(magnCombined_allEchoes(:,:,:,2),[outfolder 'tmp_Magn_echo_2.nii']);
    niftiwrite(magnCombined_allEchoes(:,:,:,3),[outfolder 'tmp_Magn_echo_3.nii']);
    niftiwrite(magnCombined_allEchoes(:,:,:,4),[outfolder 'tmp_Magn_echo_4.nii']);
    niftiwrite(magnCombined_allEchoes(:,:,:,5),[outfolder 'tmp_Magn_echo_5.nii']);

    [mask] = Get_FSL_mask(outfolder);

    indecho = 1;
    for indecho = 1:5
        disp(['Processing for local phase of echo:' num2str(indecho-1, '%04d')])
        tissuePhase_5Echoes(:,:,:,indecho) = V_SHARP(single(phaseCombined_allEchoes(:,:,:, indecho)),mask, 'smvsize',20, 'voxelsize', voxelsize);
    end
    save([outfolder '/LocalPhase_5Echoes.mat'],'tissuePhase_5Echoes');
    niftiwrite(tissuePhase_5Echoes(:,:,:,1),[outfolder 'tmp_LocalPhase_echo_1.nii']);
    niftiwrite(tissuePhase_5Echoes(:,:,:,2),[outfolder 'tmp_LocalPhase_echo_2.nii']);
    niftiwrite(tissuePhase_5Echoes(:,:,:,3),[outfolder 'tmp_LocalPhase_echo_3.nii']);
    niftiwrite(tissuePhase_5Echoes(:,:,:,4),[outfolder 'tmp_LocalPhase_echo_4.nii']);
    niftiwrite(tissuePhase_5Echoes(:,:,:,5),[outfolder 'tmp_LocalPhase_echo_5.nii']);

    indecho = 1;
    for indecho = 1:5 
        call_SEPIA_QSM(outfolder,['/data/u_gkotsoulias_software/matScripts/In-vivo_data_Analysis/7T_echo_' num2str(indecho) '_NMR134_SEPIA_HEADER.mat' ],indecho);
        QSM_iLSQR_5Echoes(:,:,:,indecho) = QSM_iLSQR(single(tissuePhase_5Echoes(:,:,:,indecho)), mask,'TE',TE_ms(indecho),'B0',7,'H',[0 0 1],'padsize',[12 12 12],'voxelsize',voxelsize);
        QSM_NDI_5Echoes(:,:,:,indecho) = niftiread([outfolder 'QSM_NDI_echo_' num2str(indecho) '_Chimap.nii.gz']);
    end
    save([outfolder '/QSM_NDI_5Echoes.mat'],'QSM_NDI_5Echoes');
    save([outfolder '/QSM_iLSQR_5Echoes.mat'],'QSM_iLSQR_5Echoes');
    niftiwrite(mean(QSM_NDI_5Echoes, 4),[outfolder 'QSM_NDI_mean_5echoes.nii.gz']);
    niftiwrite(mean(QSM_iLSQR_5Echoes, 4),[outfolder 'QSM_iLSQR_mean_5echoes.nii.gz']);
    niftiwrite(QSM_iLSQR_5Echoes(:,:,:,3),[outfolder 'QSM_iLSQR_echo_3_Chimap.nii']);
    niftiwrite(QSM_iLSQR_5Echoes(:,:,:,4),[outfolder 'QSM_iLSQR_echo_4_Chimap.nii']);

    system(['rm -f ' outfolder 'tmp_*']);
    system(['rm -f ' outfolder 'run*']);  
    system(['rm -f ' outfolder 'sepia*']);
    system(['rm -f ' outfolder '*referenceregion.nii.gz']);
    system(['rm -f ' outfolder 'QSM_NDI_echo_1_Chimap.nii.gz']); 
    system(['rm -f ' outfolder 'QSM_NDI_echo_2_Chimap.nii.gz']); 
    system(['rm -f ' outfolder 'QSM_NDI_echo_5_Chimap.nii.gz']); 
    toc
end


