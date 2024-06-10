%% This script for subdivisions spectrogram informative, noninformative LFP sites from Pretraining & Posttraining
clc; clear; close all;
default_dir=pwd;

Task        = {'Spatial','Feature'};
Training    = {'PRETRAINING','POSTTRAINING'};

%%% Areas first AD, Second PD, Third MD, Fourth AV, Fifth PV
Area     = {'Anterior-Dorsal','Posterior-Dorsal','Mid-Dorsal','Anterior-Ventral','Posterior-Ventral'};
AREA     = {{'9','46'},{'8A'},{'8B','9/46'},'45',{'47/12'}};

timeStep    = .02; startTime=-0.748; endTime=4.7520;
timeVector  = timeStep*(startTime/timeStep:endTime/timeStep);
freqsVector = importdata(['Data/Frequency.mat']);


Info= {'Informativetoall','Partial_Informative','NonInformative'}
tas =1;   


%% Post-training

fig = figure; hold on
set(fig, 'Position', [30 40 850 650]);  
%%%% PD
Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/',Task{tas},'_',Training{2},'_',Area{2},'_',Info{1},'_LFP.mat']);
if ~isempty(Pow)    
    h   = subplot(3,3,1); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Task{tas}(1),Training{2}(1:4),Info{1},'-',num2str(Pow.Trials)],h); clear Pow;
end

Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/',Task{tas},'_',Training{2},'_',Area{2},'_',Info{2},'_LFP.mat']);
if ~isempty(Pow)
    h   = subplot(3,3,2); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Task{tas}(1),Training{2}(1:4),Info{2},'-',num2str(Pow.Trials)],h); clear Pow;
end

Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/',Task{tas},'_',Training{2},'_',Area{2},'_',Info{3},'_LFP.mat']);
if ~isempty(Pow)
    h   = subplot(3,3,3); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Task{tas}(1),Training{2}(1:4),Info{3},'-',num2str(Pow.Trials)],h); clear Pow;
end

%%%% MD
Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/',Task{tas},'_',Training{2},'_',Area{3},'_',Info{1},'_LFP.mat']);
if ~isempty(Pow)
    h   = subplot(3,3,4); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Task{tas}(1),Training{2}(1:4),Info{1},'-',num2str(Pow.Trials)],h); clear Pow;
end

Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/',Task{tas},'_',Training{2},'_',Area{3},'_',Info{2},'_LFP.mat']);
if ~isempty(Pow)
    h   = subplot(3,3,5); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Task{tas}(1),Training{2}(1:4),Info{2},'-',num2str(Pow.Trials)],h); clear Pow;
end

Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/',Task{tas},'_',Training{2},'_',Area{3},'_',Info{3},'_LFP.mat']);
if ~isempty(Pow)
    h   = subplot(3,3,6); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Task{tas}(1),Training{2}(1:4),Info{3},'-',num2str(Pow.Trials)],h); clear Pow;
end

%%%%% PV
Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/',Task{tas},'_',Training{2},'_',Area{5},'_',Info{1},'_LFP.mat']);
if ~isempty(Pow)
    h   = subplot(3,3,7); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Task{tas}(1),Training{2}(1:4),Info{1},'-',num2str(Pow.Trials)],h); clear Pow;
end

Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/',Task{tas},'_',Training{2},'_',Area{5},'_',Info{2},'_LFP.mat']);
if ~isempty(Pow)
    h   = subplot(3,3,8); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Task{tas}(1),Training{2}(1:4),Info{2},'-',num2str(Pow.Trials)],h); clear Pow;
end

Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/',Task{tas},'_',Training{2},'_',Area{5},'_',Info{3},'_LFP.mat']);
if ~isempty(Pow)
    h   = subplot(3,3,9); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Task{tas}(1),Training{2}(1:4),Info{3},'-',num2str(Pow.Trials)],h); clear Pow;
end
 
%% Pre-training

fig = figure; hold on
set(fig, 'Position', [30 40 850 650]);  


%%%% PD
Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/',Task{tas},'_',Training{1},'_',Area{2},'_',Info{1},'_LFP.mat']);
if ~isempty(Pow)
    h   = subplot(2,3,1); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Task{tas}(1),Training{1}(1:4),Info{1},'-',num2str(Pow.Trials)],h); clear Pow;
end

Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/',Task{tas},'_',Training{1},'_',Area{2},'_',Info{2},'_LFP.mat']);
if ~isempty(Pow)
    h   = subplot(2,3,2); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Task{tas}(1),Training{1}(1:4),Info{2},'-',num2str(Pow.Trials)],h); clear Pow;
end

Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/',Task{tas},'_',Training{1},'_',Area{2},'_',Info{3},'_LFP.mat']);
if ~isempty(Pow)
    h   = subplot(2,3,3); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Task{tas}(1),Training{1}(1:4),Info{3},'-',num2str(Pow.Trials)],h); clear Pow;
end

%%%% MD
Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/',Task{tas},'_',Training{1},'_',Area{3},'_',Info{1},'_LFP.mat']);
if ~isempty(Pow)
    h   = subplot(2,3,4); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Task{tas}(1),Training{1}(1:4),Info{1},'-',num2str(Pow.Trials)],h); clear Pow;
end

Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/',Task{tas},'_',Training{1},'_',Area{3},'_',Info{2},'_LFP.mat']);
if ~isempty(Pow)
    h   = subplot(2,3,5); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Task{tas}(1),Training{1}(1:4),Info{2},'-',num2str(Pow.Trials)],h); clear Pow;
end

Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/',Task{tas},'_',Training{1},'_',Area{3},'_',Info{3},'_LFP.mat']);
if ~isempty(Pow)
    h   = subplot(2,3,6); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Task{tas}(1),Training{1}(1:4),Info{3},'-',num2str(Pow.Trials)],h); clear Pow;
end




   
