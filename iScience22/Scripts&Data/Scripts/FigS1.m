%% Induced power spctrum of Local field potentials for spatial stimuli, relative to Each frequency means LEFT & RIGHT Hemisphere

clc; clear all; close all;
     
default_dir     = pwd;
params.Fs       = 500; % sampling frequency
params.fpass    = [0.5 100]; % band of frequencies to be kept
params. tapers  = [9 11]; % taper parameters
movingwin       = [.5 0.1];
 
Task       = {'Spatial','Feature'};
Training   = {'PRETRAINING','POSTTRAINING'};
Area       = {'Mid-Dorsal','Posterior-Dorsal','Posterior-Ventral'};
Hemisphere ={'LEFT_HEMILFP','RIGHT_HEMILFP'};

  
fig = figure; hold on
set(fig, 'Position', [30 40 850 650]);

%% MD LEFT

TEM       = importdata(['ALLTRIALSLFP/ELV_',Task{1},'_',Training{1},'_',Area{1},'_LFP.mat']);
S_Induced = zeros(66,51);

for tri=1:size(TEM,1)
    TEMP          = TEM(tri,:);
    TEMPBaseline_ = TEM(tri,1:500);
    
    [STRI,t,f]              = mtspecgramc(TEMP,movingwin,params);
    [STRI_Baseline,t_b,f_b] = mtspecgramc(TEMPBaseline_,movingwin,params);
    for ii=1:size(STRI,2)
        SInducedDb(:,ii) = STRI(:,ii)./ mean(STRI_Baseline(:,ii));
    end
    S_Induced=[S_Induced+SInducedDb];
    clear STRI STRI_Baseline TEMP TEMPBaseline_ SInducedDb;
end
SPre_Induced = 10*log10(S_Induced./size(TEM,1));    
T            = linspace(-2,5,length(t));
h            = subplot(2,3,1); ScriptforPowS(T(10:end),f,SPre_Induced(10:66,:),'PRE-Induced (MD)',h);
clear TEM S_Induced;




TEM       = importdata(['ALLTRIALSLFP/',Hemisphere{1},'_',Task{1},'_',Area{1},'_LFP.mat']);
S_Induced = zeros(66,51);

for tri=1:size(TEM,1)
    TEMP          = TEM(tri,:);
    TEMPBaseline_ = TEM(tri,1:500);
    
    [STRI,t,f]              = mtspecgramc(TEMP,movingwin,params);
    [STRI_Baseline,t_b,f_b] = mtspecgramc(TEMPBaseline_,movingwin,params);
    for ii=1:size(STRI,2)
        SInducedDb(:,ii) = STRI(:,ii)./ mean(STRI_Baseline(:,ii));
    end
    S_Induced=[S_Induced+SInducedDb];
    clear STRI STRI_Baseline TEMP TEMPBaseline_ SInducedDb;
end
SPostL_Induced = 10*log10(S_Induced./size(TEM,1));  
h              = subplot(2,3,2); ScriptforPowS(T(10:end),f,SPostL_Induced(10:66,:),'POST-L-Induced (MD)',h);
clear TEM S_Induced;



%%%% Difference
Diff = SPostL_Induced(10:end,:)-SPre_Induced(10:end,:);
h    = subplot(2,3,3); ScriptforPowDiff(T(10:end),f,Diff,'Difference',h);
clear Diff SPostL_Induced;

%%% Right
h    = subplot(2,3,4); ScriptforPowS(T(10:end),f,SPre_Induced(10:66,:),'PRE-Induced (MD)',h);


TEM       = importdata(['ALLTRIALSLFP/',Hemisphere{2},'_',Task{1},'_',Area{1},'_LFP.mat']);
S_Induced = zeros(66,51);

for tri=1:size(TEM,1)
    TEMP          = TEM(tri,:);
    TEMPBaseline_ = TEM(tri,1:500);
    
    [STRI,t,f]              = mtspecgramc(TEMP,movingwin,params);
    [STRI_Baseline,t_b,f_b] = mtspecgramc(TEMPBaseline_,movingwin,params);
    for ii=1:size(STRI,2)
        SInducedDb(:,ii) = STRI(:,ii)./ mean(STRI_Baseline(:,ii));
    end
    S_Induced=[S_Induced+SInducedDb];
    clear STRI STRI_Baseline TEMP TEMPBaseline_ SInducedDb;
end
SPostR_Induced = 10*log10(S_Induced./size(TEM,1));  
h              = subplot(2,3,5); ScriptforPowS(T(10:end),f,SPostR_Induced(10:66,:),'POST-R-Induced (MD)',h);
clear TEM S_Induced;

%%%% Difference
Diff = SPostR_Induced(10:end,:)-SPre_Induced(10:end,:);
h    = subplot(2,3,6); ScriptforPowDiff(T(10:end),f,Diff,'Difference',h);
clear Diff SPostR_Induced SPre_Induced;

set(gcf,'renderer','Painters');set(gcf,'name',['Induced_Spatial_Left_RightHemiSphere']);
% print(fig,['FigS1.eps'],'-depsc');
% saveas(fig,['FigS1.fig'],'fig' );
% close all;

% 

