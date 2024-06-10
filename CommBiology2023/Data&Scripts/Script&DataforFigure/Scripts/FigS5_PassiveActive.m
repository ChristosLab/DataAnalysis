%% This script for After Training Passive and Active

clc; clear all; close all;
default_dir=pwd;


timeStep    = .02; startTime=-0.748; endTime=4.7520;
timeVector  = timeStep*(startTime/timeStep:endTime/timeStep);
Fix  = [1:38]; Cue = [39:63]; Delay1 = [64:138]; Sample = [139:163]; Delay2 = [164:238]; Choice = [239:276];

freqsVector = importdata(['Data/Frequency.mat']);
Gamma       = [23:51];

Pow = importdata('Data/PowerSpectrum/SessionBySessionAfterTraining/ELV_Spatial_Passive_AfterTraining_Selective.mat');




for ff=1:length(Pow)  
    SInducedDb       = Pow{ff};   
    SPA_sel_F(ff,:)  = mean(mean(SInducedDb(Fix,Gamma),2));
    SPA_sel_C(ff,:)  = mean(mean(SInducedDb(Cue,Gamma),2)); 
    SPA_sel_D1(ff,:) = mean(mean(SInducedDb(Delay1,Gamma),2)); 
    SPA_sel_S(ff,:)  = mean(mean(SInducedDb(Sample,Gamma),2)); 
    SPA_sel_D2(ff,:) = mean(mean(SInducedDb(Delay2,Gamma),2)); 
    clear SInducedDb; 
end
clear Pow;



Pow = importdata('Data/PowerSpectrum/SessionBySessionAfterTraining/ELV_Spatial_Active_AfterTraining_Selective.mat');


for ff=1:length(Pow)  
    SInducedDb         = Pow{ff}; 
    SPAAC_sel_F(ff,:) = mean(mean(SInducedDb(Fix,Gamma),2));
    SPAAC_sel_C(ff,:) = mean(mean(SInducedDb(Cue,Gamma),2)); 
    SPAAC_sel_D1(ff,:) = mean(mean(SInducedDb(Delay1,Gamma),2));
    
    SPAAC_sel_S(ff,:) = mean(mean(SInducedDb(Sample,Gamma),2)); 
    SPAAC_sel_D2(ff,:) = mean(mean(SInducedDb(Delay2,Gamma),2));
    clear SInducedDb; 
end 
 
 

DATA  = [SPA_sel_F; SPAAC_sel_F; SPA_sel_C; SPAAC_sel_C; SPA_sel_D1; SPAAC_sel_D1;SPA_sel_S; SPAAC_sel_S; SPA_sel_D2; SPAAC_sel_D2];
GROUP = [repmat({'Pas_F'},[length(SPA_sel_F),1]); repmat({'Act_F'},[length(SPAAC_sel_F),1]);... 
         repmat({'Pas_C'},[length(SPA_sel_C),1]); repmat({'Act_C'},[length(SPAAC_sel_C),1]);...
         repmat({'Pas_D'},[length(SPA_sel_D1),1]); repmat({'Act_D'},[length(SPAAC_sel_D1),1]);...
         repmat({'Pas_S'},[length(SPA_sel_S),1]); repmat({'Act_S'},[length(SPAAC_sel_S),1]);...
         repmat({'pas_D2'},[length(SPA_sel_D2),1]); repmat({'Act_D2'},[length(SPAAC_sel_D2),1])];

boxplot (DATA,GROUP,'Colors','mrmrmrmrmr');
title('AFTERTRAINING')
ylabel('Gammapower-Baseline (dB)')
xlabel('Epoch')
ylim([-6 6])
