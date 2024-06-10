%% This script for Spatial EPOCH Pretraining-Posttraining (AREA Informative-NonInformative) spectrogram
clc; clear; close all;
default_dir=pwd;

Task        = {'Spatial','Feature'};
Training    = {'PRETRAINING','POSTTRAINING'};
Area        = {'Posterior-Dorsal','Mid-Dorsal','Posterior-Ventral'};
timeStep    = .02; startTime=-0.748; endTime=4.7520;
timeVector  = timeStep*(startTime/timeStep:endTime/timeStep);
freqsVector = importdata(['Data/Frequency.mat']);

Epoch = {'Cue_','Delay1_','Sample_','Delay2_'};

tas=1;     

fig = figure; hold on;
set(fig, 'Position', [30 40 850 650]);  

Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/',Task{tas},'_',Training{2},'_AllRegion_',Epoch{1},'LFP.mat']);
h   = subplot(2,4,1); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Task{tas}(1),Training{2}(1:4),Epoch{1},'-',num2str(Pow.Trials)],h); clear Pow;

Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/',Task{tas},'_',Training{2},'_AllRegion_',Epoch{2},'LFP.mat']);
h   = subplot(2,4,2); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Task{tas}(1),Training{2}(1:4),Epoch{2},'-',num2str(Pow.Trials)],h); clear Pow;

Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/',Task{tas},'_',Training{2},'_AllRegion_',Epoch{3},'LFP.mat']);
h   = subplot(2,4,3); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Task{tas}(1),Training{2}(1:4),Epoch{3},'-',num2str(Pow.Trials)],h); clear Pow;

Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/',Task{tas},'_',Training{2},'_AllRegion_',Epoch{4},'LFP.mat']);
h   = subplot(2,4,4); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Task{tas}(1),Training{2}(1:4),Epoch{4},'-',num2str(Pow.Trials)],h); clear Pow;



Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/',Task{tas},'_',Training{1},'_AllRegion_',Epoch{1},'LFP.mat']);
h   = subplot(2,4,5); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Task{tas}(1),Training{1}(1:4),Epoch{1},'-',num2str(Pow.Trials)],h); clear Pow;

Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/',Task{tas},'_',Training{1},'_AllRegion_',Epoch{2},'LFP.mat']);
h   = subplot(2,4,6); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Task{tas}(1),Training{1}(1:4),Epoch{2},'-',num2str(Pow.Trials)],h); clear Pow;

Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/',Task{tas},'_',Training{1},'_AllRegion_',Epoch{3},'LFP.mat']);
h   = subplot(2,4,7); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Task{tas}(1),Training{1}(1:4),Epoch{3},'-',num2str(Pow.Trials)],h); clear Pow;

Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/',Task{tas},'_',Training{1},'_AllRegion_',Epoch{4},'LFP.mat']);
h   = subplot(2,4,8); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Task{tas}(1),Training{1}(1:4),Epoch{4},'-',num2str(Pow.Trials)],h); clear Pow;
