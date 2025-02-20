function [magniCombined, phaseCombined,W] = Combine_Coils_JC(img,SpatialRes)
%% read 4D dataset [x y z coils] and combine multicoil, return mag, phase of each echo
%% this multiple coils
% Input
% img        -  4D array, [Nx, Ny, Nz, Ncoil]
% SpatialRes -  resolution, e.g. [1 1 2]
% Output
% magniCombined        -  magnitude image
% phaseCombined        -  phase image
% Chunlei Liu, PhD, 05/26/2014
% J Chen Jun 10 2021 Modified

%% parameters
Imagesize=size(img);
ncoil = 1; 
if length(Imagesize) > 3
    ncoil = Imagesize(4);
end


if ncoil > 1
    %% combinie coils
    % calcualte the magnitude and phase
    Magni0=abs(img);
    Phase0=angle(img);
    
    % calcuate the weights for combining the phase
    % the aim is to derive a weighting function that is smooth and without
    % anatimical information
    magftd=single(zeros(Imagesize));
    phaseuwp=single(zeros(Imagesize));
    
    fprintf('Reading Coil#/%d ',ncoil);
    for icoil=1:ncoil
        fprintf('%d..',icoil);
%         magftd(:,:,:,icoil)=BGRemoval_v2(Magni0(:,:,:,icoil),...
%             ones(Imagesize(1:3)),5,[0 0 6],SpatialRes);
%         phaseuwp(:,:,:,icoil)=MRPhaseUnwrap_v1((Phase0(:,:,:,icoil)),SpatialRes,[0 0 6]);
        magftd(:,:,:,icoil)=BGRemoval_v2(Magni0(:,:,:,icoil),...
            ones(Imagesize(1:3)),10,[12 12 12],SpatialRes); %JC
        phaseuwp(:,:,:,icoil)=MRPhaseUnwrap_v1((Phase0(:,:,:,icoil)),SpatialRes,[12 12 12]);
    end
    fprintf('\n');
    Magfiltered=Magni0-magftd;
    clear magftd
    W=Magfiltered.^2;
    clear Magfiltered
    if length(size(W)) >3
        W=W./repmat(sum(W,4),[1 1 1 Imagesize(4) 1]);
        % phase combination
        phaseCombined=squeeze(sum(W.*phaseuwp,4));
        magniCombined=squeeze(sqrt(sum(abs(img).^2,4)));
    else
        % phase combination
        phaseCombined=phaseuwp;
        magniCombined=abs(img);
    end
else
    phaseCombined = angle(squeeze(img));
    magniCombined = abs(squeeze(img));
end


end


function [PhaseFiltered,UpdatedMask]=BGRemoval_v2(phi0,Mask,R_Filter,padsize,SpatialRes)
                                   
% Wei Li, Duke Uiversity, May 5, 2012
% Wei Li, Updated March. 28, 2011

%% zeropadding to square matrix
% tic
phi0=padarray(phi0,padsize);
Mask=padarray(Mask,padsize);
SS=size(Mask);
%% Iterative Filtering
if mod(SS(3),2)==1
    phi0(1,1,SS(3)+1)=0;
    Mask(1,1,SS(3)+1)=0;
    flagSS3=1;
else 
    flagSS3=0;
end

if SpatialRes(3)==SpatialRes(1)
    IsotropicRes=1;
else
    IsotropicRes=0;
end

switch IsotropicRes
    case 1
%         disp('isotropic resolution')
        PhaseFiltered=zeros(size(phi0));
        for i=1:R_Filter
            [LocalPhase,~]=myRemoval(phi0,1,Mask,i);
            PhaseFiltered=LocalPhase;
            % disp(i)
        end
    case 0 
%         disp('anisotropic resolution')
        NP=size(phi0);
        Z2x=SpatialRes(3)/SpatialRes(1);
        [yy,xx,zz]=meshgrid(1:NP(2),1:NP(1),(1:NP(3))*Z2x-Z2x/2);
        Zres1=2*ceil((NP(3)*Z2x-Z2x/2)/2);
        [yy1,xx1,zz1]=meshgrid(1:NP(2),1:NP(1),1:Zres1);
        PhaseUWPUpsampled=interp3(yy,xx,zz,phi0,yy1,xx1,zz1);
        PhaseUWPUpsampled(isnan(PhaseUWPUpsampled))=0;
        MaskUpsampled=interp3(yy,xx,zz,single(Mask),yy1,xx1,zz1);
        MaskUpsampled(isnan(MaskUpsampled))=0;
        MaskUpsampled=MaskUpsampled>0.5;
        PhiFiltered=zeros(size(MaskUpsampled));
        for i=1:R_Filter
            [LocalPhase,ValidPoint]=myRemoval(PhaseUWPUpsampled,1,MaskUpsampled,i);
            Index=abs(ValidPoint-1)<1e-6;
            PhiFiltered(Index)=LocalPhase(Index);
            if i==2
                UpdatedMask=Index;
            end
%             disp(i)
        end
        PhiFiltered=PhiFiltered.*UpdatedMask;
        ResidualPhaseUpsampled=PhaseUWPUpsampled-PhiFiltered;
        ResidualPhase=interp3(yy1,xx1,zz1,ResidualPhaseUpsampled,yy,xx,zz);
        ResidualPhase(isnan(ResidualPhase))=0;
        PhaseFiltered=(phi0-ResidualPhase); 

        UpdatedMask=interp3(yy1,xx1,zz1,single(UpdatedMask),yy,xx,zz);
        UpdatedMask(isnan(UpdatedMask))=0;
        UpdatedMask=(UpdatedMask>0.5);
end

UpdatedMask=0;
%%
PhaseFiltered=PhaseFiltered(padsize(1)+1:end-padsize(1),padsize(2)+1:end-padsize(2),padsize(3)+1:end-padsize(3));
UpdatedMask=UpdatedMask(padsize(1)+1:end-padsize(1),padsize(2)+1:end-padsize(2),padsize(3)+1:end-padsize(3));
% PhaseFiltered=PhaseFiltered.*UpdatedMask;

if flagSS3==1
    %UpdatedMask(:,:,end)=[];
    PhaseFiltered(:,:,end)=[];
end

% toc
return
end

%% Filtering function
function [LocalPhase,ValidPoint]=myRemoval(cleanphase,control,Mask,n)

kernal3D=ball3D(n);
kernal3D=kernal3D/sum(kernal3D(:));
NP=size(cleanphase);  
F0=zeros(NP);
SS=size(F0);
F0(SS(1)/2+1-n:SS(1)/2+1+n,SS(2)/2+1-n:SS(2)/2+1+n,SS(3)/2+1-n:SS(3)/2+1+n)=kernal3D;%--------------------------------------------------------------------
%F0(round(SS(1))/2+1-n:round(SS(1)/2)+1+n,round(SS(2)/2)+1-n:round(SS(2)/2)+1+n,round(SS(3)/2)+1-n:round(SS(3)/2)+1+n)=kernal3D;
F0=fftnc(F0*sqrt(length(F0(:))));

switch control
    case 1
        LocalPhase=cleanphase-ifftnc(F0.*fftnc(cleanphase));
        ValidPoint=ifftnc(F0.*fftnc(Mask));
    case 2
        cleanphase=cleanphase.*Mask;
        LocalPhase=fftnc(cleanphase);
        F0=1-F0;
        F0(NP(1)/2+1,NP(2)/2+1,NP(3)/2+1)=1e9;
        F0=1./F0;
        F0(NP(1)/2+1,NP(2)/2+1,NP(3)/2+1)=0;
        LocalPhase=LocalPhase.*F0;
        LocalPhase=ifftnc(LocalPhase);
        ValidPoint=1;
end
end

function phi=MRPhaseUnwrap_v1(theta,SpatialRes,padsize)
% Wei Li, Duke University
% Last Modified on July 9 2012.
theta=padarray(theta,padsize);
%theta = theta(:,1:234,:);%----------------------------------------------------------------------------------------------------------------------------
SS=size(theta);
if mod(SS(3),2)==1
    theta(1,1,SS(3)+1)=0;
    flagSS3=1;
else 
    flagSS3=0;
end
NP=size(theta); % Number of pixel
[yy,xx,zz]=meshgrid(1:NP(2),1:NP(1),1:NP(3));
FOV=SpatialRes.*NP;
xx=(xx-NP(1)/2-1)/FOV(1);
yy=(yy-NP(2)/2-1)/FOV(2);
zz=(zz-NP(3)/2-1)/FOV(3);
k2=(xx).^2+(yy).^2+(zz).^2;
LPTheta=cos(theta).*ifftnc(k2.*fftnc(sin(theta)));
LPTheta=LPTheta-sin(theta).*ifftnc(k2.*fftnc(cos(theta)));
k2(NP(1)/2+1,NP(2)/2+1,NP(3)/2+1)=1e6;%-------------------------------------------------------------------------------------------------------------------
%k2(round(NP(1))/2+1,round(NP(2))/2+1.5,round(NP(3))/2+1)=1e6;
phi=ifftnc(fftnc(LPTheta)./k2);
phi=real(phi);
phi=phi(padsize(1)+1:end-padsize(1),padsize(2)+1:end-padsize(2),padsize(3)+1:end-padsize(3));

if flagSS3==1
    phi(:,:,end)=[];
end
end

function [index1,index2] =ball3D(r)
% Wei Li, PhD
% Brain Imaging And Analysis Center, Duke Uiversity.

[xx,yy,zz]=meshgrid(-r:r,-r:r,-r:r);
v2=xx.^2+yy.^2++zz.^2;
index1=(v2<=(r.^2+r/3)) & (v2>=0.1);
index2=(v2<=(r.^2+r/3)) & (v2>=0.1);
end