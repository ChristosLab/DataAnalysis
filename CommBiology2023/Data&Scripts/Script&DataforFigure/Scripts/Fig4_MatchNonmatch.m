%% This script for spectrogram informative, noninformative LFP sites from Match-NonMatch
clc; clear all; close all;
default_dir=pwd;

Task        = {'Spatial','Feature'};
Training    = {'PRETRAINING','POSTTRAINING'};
Area        = {'Posterior-Dorsal','Mid-Dorsal','Posterior-Ventral'};
timeStep    = .02; startTime=-0.748; endTime=4.7520;
timeVector  = timeStep*(startTime/timeStep:endTime/timeStep);
freqsVector = importdata(['Data/Frequency.mat']);


Info= {'Informativetoall','Partial_Informative','NonInformative'}
tas =1;   

fig = figure; hold on
set(fig, 'Position', [30 40 850 650]);  

Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/',Task{tas},'_',Training{2},'_AllRegion_',Info{1},'_LFP_Match.mat']);
h   = subplot(2,3,1); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Task{tas}(1),Training{2}(1:4),Info{1},'-',num2str(Pow.Trials)],h); clear Pow;

Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/',Task{tas},'_',Training{2},'_AllRegion_',Info{2},'_LFP_Match.mat']);
h   = subplot(2,3,2); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Task{tas}(1),Training{2}(1:4),Info{2},'-',num2str(Pow.Trials)],h); clear Pow;

Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/',Task{tas},'_',Training{2},'_AllRegion_',Info{3},'_LFP_Match.mat']);
h   = subplot(2,3,3); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Task{tas}(1),Training{2}(1:4),Info{3},'-',num2str(Pow.Trials)],h); clear Pow;
    
Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/',Task{tas},'_',Training{2},'_AllRegion_',Info{1},'_LFP_NonMatch.mat']);
h   = subplot(2,3,4); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Task{tas}(1),Training{2}(1:4),Info{1},'-',num2str(Pow.Trials)],h); clear Pow;

Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/',Task{tas},'_',Training{2},'_AllRegion_',Info{2},'_LFP_NonMatch.mat']);
h   = subplot(2,3,5); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Task{tas}(1),Training{2}(1:4),Info{2},'-',num2str(Pow.Trials)],h); clear Pow;

Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/',Task{tas},'_',Training{2},'_AllRegion_',Info{3},'_LFP_NonMatch.mat']);
h   = subplot(2,3,6); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Task{tas}(1),Training{2}(1:4),Info{3},'-',num2str(Pow.Trials)],h); clear Pow;
    
