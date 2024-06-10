%% This script for Shape PSTH Pretraining-Posttraining (AREA Informative-NonInformative)
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


tas=1;   %% Spatial Task

fig = figure; hold on
set(fig, 'Position', [30 40 850 650]);


%% 'Posterior-Dorsal'

ar=2; train = 1;
PRESEL    = importdata(['Data/PSTH/',Area{ar},'/',Task{tas},'_',Training{train},'_',Area{ar},'_Informative_psth_Upsampled_BESTclass_2FR.mat']);
PRENONSEL = importdata(['Data/PSTH/',Area{ar},'/',Task{tas},'_',Training{train},'_',Area{ar},'_NonInformative_psth_Upsampled_BESTclass_2FR.mat']);

train = 2;
POSSEL    = importdata(['Data/PSTH/',Area{ar},'/',Task{tas},'_',Training{train},'_',Area{ar},'_Informative_psth_Upsampled_BESTclass_2FR.mat']);
POSNONSEL = importdata(['Data/PSTH/',Area{ar},'/',Task{tas},'_',Training{train},'_',Area{ar},'_NonInformative_psth_Upsampled_BESTclass_2FR.mat']);

NEURONDATA = {PRESEL,PRENONSEL,POSSEL,POSNONSEL};
colormap   = {'b','b','r','r'};
linetypes  = {'-',':','-',':'};
Label      = {[Task{tas}(1),'-',Training{1}(1:3),'SEL'],[Task{tas}(1),'-',Training{1}(1:3),'NONSEL'],...                  
              [Task{tas}(1),'-',Training{2}(1:3),'SEL'],[Task{tas}(1),'-',Training{2}(1:3),'NONSEL']};

subplot(1,3,1); hold on; psth_plot(NEURONDATA, t,'y',linetypes,colormap,Label); ylim([1 13]); hold off;
clear NEURONDATA PRESEL POSSEL PRENONSEL POSNONSEL ar;

%% 'Mid-Dorsal'

ar=3; train = 1;
PRESEL    = importdata(['Data/PSTH/',Area{ar},'/',Task{tas},'_',Training{train},'_',Area{ar},'_Informative_psth_Upsampled_BESTclass_2FR.mat']);
PRENONSEL = importdata(['Data/PSTH/',Area{ar},'/',Task{tas},'_',Training{train},'_',Area{ar},'_NonInformative_psth_Upsampled_BESTclass_2FR.mat']);

train = 2;
POSSEL    = importdata(['Data/PSTH/',Area{ar},'/',Task{tas},'_',Training{train},'_',Area{ar},'_Informative_psth_Upsampled_BESTclass_2FR.mat']);
POSNONSEL = importdata(['Data/PSTH/',Area{ar},'/',Task{tas},'_',Training{train},'_',Area{ar},'_NonInformative_psth_Upsampled_BESTclass_2FR.mat']);

NEURONDATA = {PRESEL,PRENONSEL,POSSEL,POSNONSEL};
colormap   = {'b','b','r','r'};
linetypes  = {'-',':','-',':'};
Label      = {[Task{tas}(1),'-',Training{1}(1:3),'SEL'],[Task{tas}(1),'-',Training{1}(1:3),'NONSEL'],...                  
              [Task{tas}(1),'-',Training{2}(1:3),'SEL'],[Task{tas}(1),'-',Training{2}(1:3),'NONSEL']};

subplot(1,3,2); hold on; psth_plot(NEURONDATA, t,'y',linetypes,colormap,Label); ylim([2.5 11]); hold off;
clear NEURONDATA PRESEL POSSEL PRENONSEL POSNONSEL ar;



%% 'Posterior-Ventral'

ar=5; train = 1;
PRESEL    = importdata(['Data/PSTH/',Area{ar},'/',Task{tas},'_',Training{train},'_',Area{ar},'_Informative_psth_Upsampled_BESTclass_2FR.mat']);
PRENONSEL = importdata(['Data/PSTH/',Area{ar},'/',Task{tas},'_',Training{train},'_',Area{ar},'_NonInformative_psth_Upsampled_BESTclass_2FR.mat']);

train = 2;
POSSEL    = importdata(['Data/PSTH/',Area{ar},'/',Task{tas},'_',Training{train},'_',Area{ar},'_Informative_psth_Upsampled_BESTclass_2FR.mat']);
POSNONSEL = importdata(['Data/PSTH/',Area{ar},'/',Task{tas},'_',Training{train},'_',Area{ar},'_NonInformative_psth_Upsampled_BESTclass_2FR.mat']);

NEURONDATA = {PRESEL,PRENONSEL,POSSEL,POSNONSEL};
colormap   = {'b','b','r','r'};
linetypes  = {'-',':','-',':'};
Label      = {[Task{tas}(1),'-',Training{1}(1:3),'SEL'],[Task{tas}(1),'-',Training{1}(1:3),'NONSEL'],...                  
              [Task{tas}(1),'-',Training{2}(1:3),'SEL'],[Task{tas}(1),'-',Training{2}(1:3),'NONSEL']};

subplot(1,3,3); hold on; psth_plot(NEURONDATA, t,'y',linetypes,colormap,Label); ylim([2 12]); hold off;
clear NEURONDATA PRESEL POSSEL PRENONSEL POSNONSEL ar;


