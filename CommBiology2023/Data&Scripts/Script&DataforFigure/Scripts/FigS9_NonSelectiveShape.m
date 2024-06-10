%% This script for spectrogram informative, noninformative LFP sites in nonslective Feature
clc; clear all; close all;
default_dir=pwd;

Training    = {'PRETRAINING','POSTTRAINING'};
timeStep    = .02; startTime=-0.748; endTime=4.7520;
timeVector  = timeStep*(startTime/timeStep:endTime/timeStep);
freqsVector = importdata(['Data/Frequency.mat']);

Info= {'Informative','Partial_Informative','NonInformative'}

fig = figure; hold on
set(fig, 'Position', [30 40 850 650]);  


Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/Feature_',Training{2},'_GP_Selective_Feature_LFP.mat']);
h   = subplot(2,3,1); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Training{2}(1:4),Info{1},'-',num2str(Pow.Trials)],h); clear Pow;
Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/Feature_',Training{2},'_GP_Parselective_Feature_LFP.mat']);
h   = subplot(2,3,2); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Training{2}(1:4),Info{2},'-',num2str(Pow.Trials)],h); clear Pow;
Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/Feature_',Training{2},'_GP_Nonselective_Feature_LFP.mat']);
h   = subplot(2,3,3); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Training{2}(1:4),Info{3},'-',num2str(Pow.Trials)],h); clear Pow;


Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/Feature_',Training{1},'_GP_Selective_Feature_LFP.mat']);
h   = subplot(2,3,4); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Training{1}(1:4),Info{1},'-',num2str(Pow.Trials)],h); clear Pow;
Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/Feature_',Training{1},'_GP_Parselective_Feature_LFP.mat']);
h   = subplot(2,3,5); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Training{1}(1:4),Info{2},'-',num2str(Pow.Trials)],h); clear Pow;
Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/Feature_',Training{1},'_GP_Nonselective_Feature_LFP.mat']);
h   = subplot(2,3,6); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Training{1}(1:4),Info{3},'-',num2str(Pow.Trials)],h); clear Pow;



