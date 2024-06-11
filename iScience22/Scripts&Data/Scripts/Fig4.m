%% Evoked power spctrum of Local field potentials for spatial stimuli, relative to Each frequency means% Evoked Powerspectrum % db correct power (minus base power)
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
%% Posterior-Dorsal Pre               
TEM          = importdata(['ALLTRIALSLFP/',Task{1},'_',Training{1},'_',Area{2},'_LFP.mat']);
ERP          = mean(TEM,1);
ERPBaseline_ = mean(TEM(:,1:500));
[EvokedPow,t,f]            = mtspecgramc(ERP,movingwin,params);
[EvokedPow_Baseline_,tb,f] = mtspecgramc(ERPBaseline_,movingwin,params);

for ii=1:size(EvokedPow,2)
    SPre_Evoked(:,ii) =  10*log10(EvokedPow(:,ii)./ mean(EvokedPow_Baseline_(:,ii)));
end
T    = linspace(-2,5,length(t));
SPre = smooth2a(SPre_Evoked,2,2);
h    = subplot(3,3,1); ScriptforPowS(T(10:end),f,SPre(10:end,:),'PRE-Evoked (PD)',h);
clear TEM ERP EvokedPow SPre_Evoked EvokedPow_Baseline_;

%% Posterior-Dorsal Post              
TEM          = importdata(['ALLTRIALSLFP/',Task{1},'_',Training{2},'_',Area{2},'_LFP.mat']);
ERP          = mean(TEM,1);
ERPBaseline_ = mean(TEM(:,1:500));
[EvokedPow,t,f]            = mtspecgramc(ERP,movingwin,params);
[EvokedPow_Baseline_,tb,f] = mtspecgramc(ERPBaseline_,movingwin,params);

for ii=1:size(EvokedPow,2)
    SPost_Evoked(:,ii)  =   10*log10(EvokedPow(:,ii)./ mean(EvokedPow_Baseline_(:,ii)));
end
SPost = smooth2a(SPost_Evoked,2,2);
h     = subplot(3,3,2); ScriptforPowS(T(10:end),f,SPost(10:end,:),'POST-Evoked (PD)',h);

Diff = SPost-SPre;
h    = subplot(3,3,3); ScriptforPowDiff(T(10:end),f,Diff(10:end,:),'Difference',h);
clear Diff SPost SPre TEM ERP EvokedPow SPost_Evoked EvokedPow_Baseline_;



%% Mid-Dorsal Pre               
TEM          = importdata(['ALLTRIALSLFP/',Task{1},'_',Training{1},'_',Area{1},'_LFP.mat']);
ERP          = mean(TEM,1);
ERPBaseline_ = mean(TEM(:,1:500));
[EvokedPow,t,f]            = mtspecgramc(ERP,movingwin,params);
[EvokedPow_Baseline_,tb,f] = mtspecgramc(ERPBaseline_,movingwin,params);

for ii=1:size(EvokedPow,2)
    SPre_Evoked(:,ii)  =   10*log10(EvokedPow(:,ii)./ mean(EvokedPow_Baseline_(:,ii)));
end
SPre = smooth2a(SPre_Evoked,2,2);
h    = subplot(3,3,4); ScriptforPowS(T(10:end),f,SPre(10:end,:),'PRE-Evoked (MD)',h);
clear TEM ERP EvokedPow SPre_Evoked EvokedPow_Baseline_;

%% Mid-Dorsal Post              
TEM          = importdata(['ALLTRIALSLFP/',Task{1},'_',Training{2},'_',Area{1},'_LFP.mat']);
ERP          = mean(TEM,1);
ERPBaseline_ = mean(TEM(:,1:500));
[EvokedPow,t,f]            = mtspecgramc(ERP,movingwin,params);
[EvokedPow_Baseline_,tb,f] = mtspecgramc(ERPBaseline_,movingwin,params);

for ii=1:size(EvokedPow,2)
    SPost_Evoked(:,ii)  =   10*log10(EvokedPow(:,ii)./ mean(EvokedPow_Baseline_(:,ii)));
end
SPost = smooth2a(SPost_Evoked,2,2);
h     = subplot(3,3,5); ScriptforPowS(T(10:end),f,SPost(10:end,:),'POST-Evoked (MD)',h);

Diff = SPost-SPre;
h    = subplot(3,3,6); ScriptforPowDiff(T(10:end),f,Diff(10:end,:),'Difference',h);
clear Diff SPost SPre SPost_Evoked TEM ERP EvokedPow EvokedPow_Baseline_;



%% Posterior-Ventral Pre               
TEM          = importdata(['ALLTRIALSLFP/',Task{1},'_',Training{1},'_',Area{3},'_LFP.mat']);
ERP          = mean(TEM,1);
ERPBaseline_ = mean(TEM(:,1:500));
[EvokedPow,t,f]            = mtspecgramc(ERP,movingwin,params);
[EvokedPow_Baseline_,tb,f] = mtspecgramc(ERPBaseline_,movingwin,params);

for ii=1:size(EvokedPow,2)
    SPre_Evoked(:,ii)  =   10*log10(EvokedPow(:,ii)./ mean(EvokedPow_Baseline_(:,ii)));
end
SPre = smooth2a(SPre_Evoked,2,2);
h    = subplot(3,3,7); ScriptforPowS(T(10:end),f,SPre(10:end,:),'PRE-Evoked (PV)',h);
clear TEM ERP EvokedPow SPre_Evoked EvokedPow_Baseline_;

%% Posterior-Ventral Post              
TEM          = importdata(['ALLTRIALSLFP/',Task{1},'_',Training{2},'_',Area{3},'_LFP.mat']);
ERP          = mean(TEM,1);
ERPBaseline_ = mean(TEM(:,1:500));
[EvokedPow,t,f]            = mtspecgramc(ERP,movingwin,params);
[EvokedPow_Baseline_,tb,f] = mtspecgramc(ERPBaseline_,movingwin,params);

for ii=1:size(EvokedPow,2)
    SPost_Evoked(:,ii)  =  10*log10(EvokedPow(:,ii)./ mean(EvokedPow_Baseline_(:,ii)));
end
SPost= smooth2a(SPost_Evoked,2,2);
h    = subplot(3,3,8); ScriptforPowS(T(10:end),f,SPost(10:end,:),'POST-Evoked (PV)',h);

Diff = SPost-SPre;
h    = subplot(3,3,9); ScriptforPowDiff(T(10:end),f,Diff(10:end,:),'Difference',h);
clear Diff SPost SPre SPost_Evoked TEM ERP EvokedPow EvokedPow_Baseline_;



set(gcf,'renderer','Painters');set(gcf,'name',['SpatialTask_Evoked']);
% print(fig,['Fig4.eps'],'-depsc');
% saveas(fig,['Fig4.fig'],'fig' );
% close all;
