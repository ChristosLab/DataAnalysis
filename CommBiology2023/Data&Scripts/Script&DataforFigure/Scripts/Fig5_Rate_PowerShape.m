%% This script for Shape PSTH Pretraining-Posttraining (AllRegionInformative-NonInformative) & spectrogram
clc; clear; close all;
default_dir=pwd;

Task     = {'Spatial','Feature'};
Training = {'PRETRAINING','POSTTRAINING'};
%%% Areas first AD, Second PD, Third MD, Fourth AV, Fifth PV
Area     = {'Anterior-Dorsal','Posterior-Dorsal','Mid-Dorsal','Anterior-Ventral','Posterior-Ventral'};
AREA     = {{'9','46'},{'8A'},{'8B','9/46'},'45',{'47/12'}};

%%%% parameters for psth
bin_width = 0.05;  % 50 milliseconds bin
dur       = [-1, 5]; % duration starts from fixation
step      = 0.02;  %  20 millisecond step 
t = dur(1):step:dur(2)-step; %


tas=2;

fig = figure; hold on
set(fig, 'Position', [30 40 850 650]);

train = 1;
PRESEL = importdata(['Data/PSTH/AllRegion/',Task{tas},'_',Training{train},'_AllRegion_Informative_psth_Upsampled_BESTclass_2FR.mat']);
PRENONSEL  = importdata(['Data/PSTH/AllRegion/',Task{tas},'_',Training{train},'_AllRegion_NonInformative_psth_Upsampled_BESTclass_2FR.mat']);

train = 2;
POSSEL = importdata(['Data/PSTH/AllRegion/',Task{tas},'_',Training{train},'_AllRegion_Informative_psth_Upsampled_BESTclass_2FR.mat']);
POSNONSEL  = importdata(['Data/PSTH/AllRegion/',Task{tas},'_',Training{train},'_AllRegion_NonInformative_psth_Upsampled_BESTclass_2FR.mat']);

NEURONDATA = {PRESEL,PRENONSEL,POSSEL,POSNONSEL};
colormap   = {'b','b', 'r','r'};
linetypes  = {'-',':','-',':'};
Label      = {[Task{tas}(1),'-',Training{1}(1:3),'SEL'],[Task{tas}(1),'-',Training{1}(1:3),'NONSEL'],...
              [Task{tas}(1),'-',Training{2}(1:3),'SEL'],[Task{tas}(1),'-',Training{2}(1:3),'NONSEL']};

subplot(1,3,1); hold on; psth_plot(NEURONDATA, t,'y',linetypes,colormap,Label);hold off; ylim([1.5 11.5]);
clear NEURONDATA PRESEL POSSEL PRENONSEL POSNONSEL;

%%% Match
 train = 1;
PRESEL_M = importdata(['Data/PSTH/AllRegion/',Task{tas},'_',Training{train},'_AllRegion_Informative_psth_UpsampledMatch_BESTclass_2FR.mat']);
PRESEL_NM = importdata(['Data/PSTH/AllRegion/',Task{tas},'_',Training{train},'_AllRegion_Informative_psth_UpsampledNonMatch_BESTclass_2FR.mat']);
 
train = 2;
POSSEL_M = importdata(['Data/PSTH/AllRegion/',Task{tas},'_',Training{train},'_AllRegion_Informative_psth_UpsampledMatch_BESTclass_2FR.mat']);
POSSEL_NM = importdata(['Data/PSTH/AllRegion/',Task{tas},'_',Training{train},'_AllRegion_Informative_psth_UpsampledNonMatch_BESTclass_2FR.mat']);
 
 
NEURONDATA = {PRESEL_M,PRESEL_NM, POSSEL_M,POSSEL_NM};
colormap   = {'g','c','g','c'};
linetypes  = {'-','-',':',':'};
Label      = {[Task{tas}(1),'-',Training{1}(1:3),'M'],[Task{tas}(1),'-',Training{1}(1:3),'NM'],...
              [Task{tas}(1),'-',Training{2}(1:3),'M'],[Task{tas}(1),'-',Training{2}(1:3),'NM']};
subplot(1,3,2); hold on; psth_plot(NEURONDATA, t,'y',linetypes,colormap,Label);hold off;  ylim([3.5 12]); 
clear NEURONDATA PRESEL_M PRESEL_NM POSSEL_M POSSEL_NM;


%%%%% Epoch
train = 1;
PRECUESEL       = importdata(['Data/PSTH/AllRegion/',Task{tas},'_',Training{train},'_AllRegion_Cue_psth_Upsampled_BESTclass_2FR.mat']);
PREDEL1SEL      = importdata(['Data/PSTH/AllRegion/',Task{tas},'_',Training{train},'_AllRegion_Delay1_psth_Upsampled_BESTclass_2FR.mat']);
PRESAMPLE_MSEL  = importdata(['Data/PSTH/AllRegion/',Task{tas},'_',Training{train},'_AllRegion_sample_M_psth_Upsampled_BESTclass_2FR.mat']);
PRESAMPLE_NMSEL = importdata(['Data/PSTH/AllRegion/',Task{tas},'_',Training{train},'_AllRegion_sample_NM_psth_Upsampled_BESTclass_2FR.mat']);
PREDEL2_MSEL    = importdata(['Data/PSTH/AllRegion/',Task{tas},'_',Training{train},'_AllRegion_Delay2M_psth_Upsampled_BESTclass_2FR.mat']);
PREDEL2_NMSEL   = importdata(['Data/PSTH/AllRegion/',Task{tas},'_',Training{train},'_AllRegion_Delay2NM_psth_Upsampled_BESTclass_2FR.mat']);

train = 2;
POSCUESEL       = importdata(['Data/PSTH/AllRegion/',Task{tas},'_',Training{train},'_AllRegion_Cue_psth_Upsampled_BESTclass_2FR.mat']);
POSDEL1SEL      = importdata(['Data/PSTH/AllRegion/',Task{tas},'_',Training{train},'_AllRegion_Delay1_psth_Upsampled_BESTclass_2FR.mat']);
POSSAMPLE_MSEL  = importdata(['Data/PSTH/AllRegion/',Task{tas},'_',Training{train},'_AllRegion_sample_M_psth_Upsampled_BESTclass_2FR.mat']);
POSSAMPLE_NMSEL = importdata(['Data/PSTH/AllRegion/',Task{tas},'_',Training{train},'_AllRegion_sample_NM_psth_Upsampled_BESTclass_2FR.mat']);
POSDEL2_MSEL    = importdata(['Data/PSTH/AllRegion/',Task{tas},'_',Training{train},'_AllRegion_Delay2M_psth_Upsampled_BESTclass_2FR.mat']);
POSDEL2_NMSEL   = importdata(['Data/PSTH/AllRegion/',Task{tas},'_',Training{train},'_AllRegion_Delay2NM_psth_Upsampled_BESTclass_2FR.mat']);


NEURONDATA = {PRECUESEL,PREDEL1SEL,PRESAMPLE_MSEL,PRESAMPLE_NMSEL,PREDEL2_MSEL,PREDEL2_NMSEL,POSCUESEL,POSDEL1SEL,POSSAMPLE_MSEL,POSSAMPLE_NMSEL,POSDEL2_MSEL,POSDEL2_NMSEL};
colormap   = {'r','b','g','c','m','k','r','b','g','c','m','k',};
linetypes  = {'-','-','-','-','-','-',':',':',':',':',':',':'};
Label      = {[Task{tas}(1),'-',Training{1}(1:3),'C'], [Task{tas}(1),'-',Training{1}(1:3),'D1'],...
              [Task{tas}(1),'-',Training{1}(1:3),'S_M'], [Task{tas}(1),'-',Training{1}(1:3),'S_NM'],...
              [Task{tas}(1),'-',Training{1}(1:3),'D2_M'],[Task{tas}(1),'-',Training{1}(1:3),'D2_NM'],...
              [Task{tas}(1),'-',Training{2}(1:3),'C'], [Task{tas}(1),'-',Training{2}(1:3),'D1'],...
              [Task{tas}(1),'-',Training{2}(1:3),'S_M'], [Task{tas}(1),'-',Training{2}(1:3),'S_NM'],...
              [Task{tas}(1),'-',Training{2}(1:3),'D2_M'],[Task{tas}(1),'-',Training{2}(1:3),'D2_NM']}

subplot(1,3,3); hold on; psth_plot(NEURONDATA, t,'y',linetypes,colormap,Label); hold off; ylim([2 17]); 
clear NEURONDATA PRECUESEL PREDEL1SEL PRESAMPLE_MSEL PRESAMPLE_NMSEL PREDEL2_MSEL PREDEL2_NMSEL POSCUESEL POSDEL1SEL POSSAMPLE_MSEL POSSAMPLE_NMSEL POSDEL2_MSEL POSDEL2_NMSEL;


%% Spectrogram


timeStep    = .02; startTime=-0.748; endTime=4.7520;
timeVector  = timeStep*(startTime/timeStep:endTime/timeStep);
freqsVector = importdata(['Data/Frequency.mat']);

Info= {'Informativetoall','Partial_Informative','NonInformative'}
tas =2;   

fig = figure; hold on
set(fig, 'Position', [30 40 850 650]);  

Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/',Task{tas},'_',Training{2},'_AllRegion_',Info{1},'_LFP.mat']);
h   = subplot(2,3,1); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Task{tas}(1),Training{2}(1:4),Info{1},'-',num2str(Pow.Trials)],h); clear Pow;
Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/',Task{tas},'_',Training{2},'_AllRegion_',Info{2},'_LFP.mat']);
h   = subplot(2,3,2); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Task{tas}(1),Training{2}(1:4),Info{2},'-',num2str(Pow.Trials)],h); clear Pow;
Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/',Task{tas},'_',Training{2},'_AllRegion_',Info{3},'_LFP.mat']);
h   = subplot(2,3,3); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Task{tas}(1),Training{2}(1:4),Info{3},'-',num2str(Pow.Trials)],h); clear Pow;


Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/',Task{tas},'_',Training{1},'_AllRegion_',Info{1},'_LFP.mat']);
h   = subplot(2,3,4); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Task{tas}(1),Training{1}(1:4),Info{1},'-',num2str(Pow.Trials)],h); clear Pow;
Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/',Task{tas},'_',Training{1},'_AllRegion_',Info{2},'_LFP.mat']);
h   = subplot(2,3,5); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Task{tas}(1),Training{1}(1:4),Info{2},'-',num2str(Pow.Trials)],h); clear Pow;
Pow = importdata(['Data/PowerSpectrum/TrialbyTrial/',Task{tas},'_',Training{1},'_AllRegion_',Info{3},'_LFP.mat']);
h   = subplot(2,3,6); hold on;  Plot_powspectrum(timeVector,freqsVector,Pow.Induced,[Task{tas}(1),Training{1}(1:4),Info{3},'-',num2str(Pow.Trials)],h); clear Pow;


