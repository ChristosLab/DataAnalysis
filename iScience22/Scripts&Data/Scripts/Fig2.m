%% Induced power spctrum of Local field potentials for spatial stimuli, relative to Each frequency means
default_dir     = pwd;
params.Fs       = 500; % sampling frequency
params.fpass    = [0.5 100]; % band of frequencies to be kept
params. tapers  = [9 11]; % taper parameters
movingwin       = [.5 0.1];
 
Task     = {'Spatial','Feature'};
Training = {'PRETRAINING','POSTTRAINING'};
Area     = {'Mid-Dorsal','Posterior-Dorsal','Posterior-Ventral'};

   
fig = figure; hold on
set(fig, 'Position', [30 40 850 650]);
%% Posterior-Dorsal Pre % Induced (Total Power)
               
TEM       = importdata(['ALLTRIALSLFP/',Task{1},'_',Training{1},'_',Area{2},'_LFP.mat']);
S_Induced = zeros(66,51);
for tri=1:size(TEM,1)
    TEMP          = TEM(tri,:);
    TEMPBaseline_ = TEM(tri,1:500);
    
    [STRI,t,f]               = mtspecgramc(TEMP,movingwin,params);
    [STRI_Beaseline,t_b,f_b] = mtspecgramc(TEMPBaseline_,movingwin,params);
    for ii=1:size(STRI,2)
        SInducedDb(:,ii)  =   STRI(:,ii)./ mean(STRI_Beaseline(:,ii));
    end
    S_Induced=[S_Induced+SInducedDb];
    clear STRI STRI_Beaseline TEMP TEMPBaseline_ SInducedDb;
end
SPre_Induced = 10*log10(S_Induced./size(TEM,1));    
T            = linspace(-2,5,length(t));
h            = subplot(3,3,1); ScriptforPowS(T(10:end),f,SPre_Induced(10:66,:),'PRE-Induced (PD)',h);
clear TEM S_Induced;

% %%% Post-Training              
TEM       = importdata(['ALLTRIALSLFP/',Task{1},'_',Training{2},'_',Area{2},'_LFP.mat']);
S_Induced = zeros(66,51);
for tri=1:size(TEM,1)
    TEMP          = TEM(tri,:);
    TEMPBaseline_ = TEM(tri,1:500);
    
    [STRI,t,f]               = mtspecgramc(TEMP,movingwin,params);
    [STRI_Beaseline,t_b,f_b] = mtspecgramc(TEMPBaseline_,movingwin,params);
    for ii=1:size(STRI,2)
        SInducedDb(:,ii)  =   STRI(:,ii)./ mean(STRI_Beaseline(:,ii));
    end
    S_Induced=[S_Induced+SInducedDb];
    clear STRI STRI_Beaseline TEMP TEMPBaseline_ SInducedDb;
end
SPost_Induced = 10*log10(S_Induced./size(TEM,1));  
h             = subplot(3,3,2); ScriptforPowS(T(10:end),f,SPost_Induced(10:66,:),'POST-Induced (PD)',h);
clear TEM S_Induced;
 
%%%% Difference
Diff = SPost_Induced(10:66,:)-SPre_Induced(10:66,:);
h    = subplot(3,3,3); ScriptforPowDiff(T(10:end),f,Diff,'Difference',h);
clear Diff SPost_Induced SPre_Induced;

%% Mid-Dorsal Pre               %%% Induced (Total Power)
               
TEM       = importdata(['ALLTRIALSLFP/',Task{1},'_',Training{1},'_',Area{1},'_LFP.mat']);
S_Induced = zeros(66,51);
for tri=1:size(TEM,1)
    TEMP          = TEM(tri,:);
    TEMPBaseline_ = TEM(tri,1:500);
    
    [STRI,t,f]               = mtspecgramc(TEMP,movingwin,params);
    [STRI_Beaseline,t_b,f_b] = mtspecgramc(TEMPBaseline_,movingwin,params);
    for ii=1:size(STRI,2)
        SInducedDb(:,ii)  =   STRI(:,ii)./ mean(STRI_Beaseline(:,ii));
    end
    S_Induced=[S_Induced+SInducedDb];
    clear STRI STRI_Beaseline TEMP TEMPBaseline_ SInducedDb;
end
SPre_Induced = 10*log10(S_Induced./size(TEM,1));  
h            = subplot(3,3,4); ScriptforPowS(T(10:end),f,SPre_Induced(10:66,:),'PRE-Induced (MD)',h);
clear TEM S_Induced;

%%% Post-Training              
TEM       = importdata(['ALLTRIALSLFP/',Task{1},'_',Training{2},'_',Area{1},'_LFP.mat']);
S_Induced = zeros(66,51);
for tri=1:size(TEM,1)
    TEMP          = TEM(tri,:);
    TEMPBaseline_ = TEM(tri,1:500);
    
    [STRI,t,f]               = mtspecgramc(TEMP,movingwin,params);
    [STRI_Beaseline,t_b,f_b] = mtspecgramc(TEMPBaseline_,movingwin,params);
    for ii=1:size(STRI,2)
        SInducedDb(:,ii)  =   STRI(:,ii)./ mean(STRI_Beaseline(:,ii));
    end
    S_Induced=[S_Induced+SInducedDb];
    clear STRI STRI_Beaseline TEMP TEMPBaseline_ SInducedDb;
end
SPost_Induced = 10*log10(S_Induced./size(TEM,1));  
h             = subplot(3,3,5); ScriptforPowS(T(10:end),f,SPost_Induced(10:66,:),'POST-Induced (MD)',h);
clear TEM S_Induced;

% %%%% Difference
Diff = SPost_Induced(10:66,:)-SPre_Induced(10:66,:);
h    = subplot(3,3,6); ScriptforPowDiff(T(10:end),f,Diff,'Difference',h);
clear Diff SPost_Induced SPre_Induced;

% %% Posterior-Ventral Pre               
               
TEM       = importdata(['ALLTRIALSLFP/',Task{1},'_',Training{1},'_',Area{3},'_LFP.mat']);
S_Induced = zeros(66,51);
for tri=1:size(TEM,1)
    TEMP          = TEM(tri,:);
    TEMPBaseline_ = TEM(tri,1:500);
    
    [STRI,t,f]               = mtspecgramc(TEMP,movingwin,params);
    [STRI_Beaseline,t_b,f_b] = mtspecgramc(TEMPBaseline_,movingwin,params);
    for ii=1:size(STRI,2)
        SInducedDb(:,ii)  =   STRI(:,ii)./ mean(STRI_Beaseline(:,ii));
    end
    S_Induced=[S_Induced+SInducedDb];
    clear STRI STRI_Beaseline TEMP TEMPBaseline_ SInducedDb;
end
SPre_Induced = 10*log10(S_Induced./size(TEM,1));  
h            = subplot(3,3,7); ScriptforPowS(T(10:end),f,SPre_Induced(10:66,:),'PRE-Induced (PV)',h);
clear TEM S_Induced;
% 
% %%% Post-Training              
TEM       = importdata(['ALLTRIALSLFP/',Task{1},'_',Training{2},'_',Area{3},'_LFP.mat']);
S_Induced = zeros(66,51);
for tri=1:size(TEM,1)
    TEMP          = TEM(tri,:);
    TEMPBaseline_ = TEM(tri,1:500);
    
    [STRI,t,f]           = mtspecgramc(TEMP,movingwin,params);
    [STRI_Beaseline,t_b,f_b] = mtspecgramc(TEMPBaseline_,movingwin,params);
    for ii=1:size(STRI,2)
        SInducedDb(:,ii)  =   STRI(:,ii)./ mean(STRI_Beaseline(:,ii));
    end
    S_Induced=[S_Induced+SInducedDb];
    clear STRI STRI_Beaseline TEMP TEMPBaseline_ SInducedDb;
end
SPost_Induced = 10*log10(S_Induced./size(TEM,1));  
h             = subplot(3,3,8); ScriptforPowS(T(10:end),f,SPost_Induced(10:66,:),'POST-Induced (PV)',h);
clear TEM S_Induced;

%%%% Difference
Diff = SPost_Induced(10:66,:)-SPre_Induced(10:66,:);
h    = subplot(3,3,9); ScriptforPowDiff(T(10:end),f,Diff,'Difference',h);
clear Diff SPost_Induced SPre_Induced;

% set(gcf,'renderer','Painters');set(gcf,'name',['SpatialTask_Induced_ITI']);
% print(fig,['Fig2.eps'],'-depsc');
% saveas(fig,['Fig2.fig'],'fig' );
% close all;

