%% This script for Individual Monkey Gamma power (session)

clc; clear; close all;
default_dir=pwd;

timeStep    = .02; startTime=-0.748; endTime=4.7520;
timeVector  = timeStep*(startTime/timeStep:endTime/timeStep);
Fix  = [1:38]; Cue = [39:63]; Delay1 = [64:138]; Sample = [139:163]; Delay2 = [164:238]; Choice = [239:276];

freqsVector = importdata(['Data//Frequency.mat']);
Gamma       = [23:51];
 
fig = figure; hold on
set(fig, 'Position', [30 40 850 650]); 
flist = importdata('Data/DataForStatics/AllRegion/SessionbySession/ADR_Spatial_PRETRAINING_SELECTIVE_LFP.mat');  
for ff=1:length(flist)
    TEMP = flist{ff};
    for ii=1:size(TEMP,2)
        SInducedDb(:,ii) = 10*log10(TEMP(:,ii));
    end
    SPA_sel_F(ff,:) = mean(mean(SInducedDb(Fix,Gamma),2));
    SPA_sel_C(ff,:) = mean(mean(SInducedDb(Cue,Gamma),2)); 
    SPA_sel_D1(ff,:)= mean(mean(SInducedDb(Delay1,Gamma),2)); 
    
    SPA_sel_S(ff,:) = mean(mean(SInducedDb(Sample,Gamma),2)); 
    SPA_sel_D2(ff,:)= mean(mean(SInducedDb(Delay2,Gamma),2)); 
    clear SInducedDb TEMP; 
end
clear flist;
 
flist = importdata('Data/DataForStatics/AllRegion/SessionbySession/ADR_Spatial_POSTTRAINING_SELECTIVE_LFP.mat');  
for ff=1:length(flist)
    TEMP = flist{ff};
    for ii=1:size(TEMP,2)
        SInducedDb(:,ii) = 10*log10(TEMP(:,ii));
    end
    SPAPos_sel_F(ff,:) = mean(mean(SInducedDb(Fix,Gamma),2));
    SPAPos_sel_C(ff,:) = mean(mean(SInducedDb(Cue,Gamma),2)); 
    SPAPos_sel_D1(ff,:)= mean(mean(SInducedDb(Delay1,Gamma),2));
    
    SPAPos_sel_S(ff,:) = mean(mean(SInducedDb(Sample,Gamma),2)); 
    SPAPos_sel_D2(ff,:)= mean(mean(SInducedDb(Delay2,Gamma),2));
    clear SInducedDb TEMP; 
end
clear flist;


DATA  = [SPA_sel_F; SPAPos_sel_F; SPA_sel_C; SPAPos_sel_C; SPA_sel_D1; SPAPos_sel_D1; SPA_sel_S; SPAPos_sel_S; SPA_sel_D2; SPAPos_sel_D2];
GROUP = [repmat({'PRE_F'},[length(SPA_sel_F),1]); repmat({'POS_F'},[length(SPAPos_sel_F),1]);... 
         repmat({'PRE_C'},[length(SPA_sel_C),1]); repmat({'POS_C'},[length(SPAPos_sel_C),1]);...
         repmat({'PRE_D1'},[length(SPA_sel_D1),1]); repmat({'POS_D1'},[length(SPAPos_sel_D1),1]);...
         repmat({'PRE_S'},[length(SPA_sel_S),1]); repmat({'POS_S'},[length(SPAPos_sel_S),1]);...
         repmat({'PRE_D2'},[length(SPA_sel_D2),1]); repmat({'POS_D2'},[length(SPAPos_sel_D2),1])];


subplot(1,3,1); boxplot (DATA,GROUP,'Colors','brbrbrbrbr');
title('ADR')
ylabel('Gammapower-Baseline (dB)')
xlabel('Epoch')
% ylim([-2 1])

% TF = isoutlier(A)

clear DATA GROUP SPA_sel_F SPAPos_sel_F SPA_sel_C SPAPos_sel_C SPA_sel_D1 SPAPos_sel_D1 SPA_sel_S SPAPos_sel_S SPA_sel_D2 SPAPos_sel_D2;



flist = importdata('Data/DataForStatics/AllRegion/SessionbySession/ELV_Spatial_PRETRAINING_SELECTIVE_LFP.mat');  
for ff=1:length(flist)
    TEMP = flist{ff};
    for ii=1:size(TEMP,2)
        SInducedDb(:,ii) = 10*log10(TEMP(:,ii));
    end
    SPA_sel_F(ff,:) = mean(mean(SInducedDb(Fix,Gamma),2));
    SPA_sel_C(ff,:) = mean(mean(SInducedDb(Cue,Gamma),2)); 
    SPA_sel_D1(ff,:)= mean(mean(SInducedDb(Delay1,Gamma),2)); 
    
    SPA_sel_S(ff,:) = mean(mean(SInducedDb(Sample,Gamma),2)); 
    SPA_sel_D2(ff,:)= mean(mean(SInducedDb(Delay2,Gamma),2)); 
    clear SInducedDb TEMP; 
end
clear flist;

flist = importdata('Data/DataForStatics/AllRegion/SessionbySession/ELV_Spatial_POSTTRAINING_SELECTIVE_LFP.mat');  
for ff=1:length(flist)
    TEMP = flist{ff};
    for ii=1:size(TEMP,2)
        SInducedDb(:,ii) = 10*log10(TEMP(:,ii));
    end
    SPAPos_sel_F(ff,:) = mean(mean(SInducedDb(Fix,Gamma),2));
    SPAPos_sel_C(ff,:) = mean(mean(SInducedDb(Cue,Gamma),2)); 
    SPAPos_sel_D1(ff,:)= mean(mean(SInducedDb(Delay1,Gamma),2));
    
    SPAPos_sel_S(ff,:) = mean(mean(SInducedDb(Sample,Gamma),2)); 
    SPAPos_sel_D2(ff,:)= mean(mean(SInducedDb(Delay2,Gamma),2));
    clear SInducedDb TEMP; 
end
clear flist;

DATA  = [SPA_sel_F; SPAPos_sel_F; SPA_sel_C; SPAPos_sel_C; SPA_sel_D1; SPAPos_sel_D1; SPA_sel_S; SPAPos_sel_S; SPA_sel_D2; SPAPos_sel_D2];
GROUP = [repmat({'PRE_F'},[length(SPA_sel_F),1]); repmat({'POS_F'},[length(SPAPos_sel_F),1]);... 
         repmat({'PRE_C'},[length(SPA_sel_C),1]); repmat({'POS_C'},[length(SPAPos_sel_C),1]);...
         repmat({'PRE_D1'},[length(SPA_sel_D1),1]); repmat({'POS_D1'},[length(SPAPos_sel_D1),1]);...
         repmat({'PRE_S'},[length(SPA_sel_S),1]); repmat({'POS_S'},[length(SPAPos_sel_S),1]);...
         repmat({'PRE_D2'},[length(SPA_sel_D2),1]); repmat({'POS_D2'},[length(SPAPos_sel_D2),1])];

subplot(1,3,2); boxplot (DATA,GROUP,'Colors','brbrbrbrbr');
title('ELV')
ylabel('Gammapower-Baseline (dB)')
xlabel('Epoch')
% ylim([-7 7])

clear DATA GROUP SPA_sel_F SPAPos_sel_F SPA_sel_C SPAPos_sel_C SPA_sel_D1 SPAPos_sel_D1 SPA_sel_S SPAPos_sel_S SPA_sel_D2 SPAPos_sel_D2;




flist = importdata('Data/DataForStatics/AllRegion/SessionbySession/NIN_Spatial_PRETRAINING_SELECTIVE_LFP.mat');  
for ff=1:length(flist)
    TEMP = flist{ff};
    for ii=1:size(TEMP,2)
        SInducedDb(:,ii) = 10*log10(TEMP(:,ii));
    end
    SPA_sel_F(ff,:) = mean(mean(SInducedDb(Fix,Gamma),2));
    SPA_sel_C(ff,:) = mean(mean(SInducedDb(Cue,Gamma),2)); 
    SPA_sel_D1(ff,:)= mean(mean(SInducedDb(Delay1,Gamma),2)); 
    
    SPA_sel_S(ff,:) = mean(mean(SInducedDb(Sample,Gamma),2)); 
    SPA_sel_D2(ff,:)= mean(mean(SInducedDb(Delay2,Gamma),2)); 
    clear SInducedDb TEMP; 
end
clear flist;

flist = importdata('Data/DataForStatics/AllRegion/SessionbySession/NIN_Spatial_POSTTRAINING_SELECTIVE_LFP.mat');  
for ff=1:length(flist)
    TEMP = flist{ff};
    for ii=1:size(TEMP,2)
        SInducedDb(:,ii) = 10*log10(TEMP(:,ii));
    end
    SPAPos_sel_F(ff,:) = mean(mean(SInducedDb(Fix,Gamma),2));
    SPAPos_sel_C(ff,:) = mean(mean(SInducedDb(Cue,Gamma),2)); 
    SPAPos_sel_D1(ff,:)= mean(mean(SInducedDb(Delay1,Gamma),2));
    
    SPAPos_sel_S(ff,:) = mean(mean(SInducedDb(Sample,Gamma),2)); 
    SPAPos_sel_D2(ff,:)= mean(mean(SInducedDb(Delay2,Gamma),2));
    clear SInducedDb TEMP; 
end
clear flist;

DATA  = [SPA_sel_F; SPAPos_sel_F; SPA_sel_C; SPAPos_sel_C; SPA_sel_D1; SPAPos_sel_D1; SPA_sel_S; SPAPos_sel_S; SPA_sel_D2; SPAPos_sel_D2];
GROUP = [repmat({'PRE_F'},[length(SPA_sel_F),1]); repmat({'POS_F'},[length(SPAPos_sel_F),1]);... 
         repmat({'PRE_C'},[length(SPA_sel_C),1]); repmat({'POS_C'},[length(SPAPos_sel_C),1]);...
         repmat({'PRE_D1'},[length(SPA_sel_D1),1]); repmat({'POS_D1'},[length(SPAPos_sel_D1),1]);...
         repmat({'PRE_S'},[length(SPA_sel_S),1]); repmat({'POS_S'},[length(SPAPos_sel_S),1]);...
         repmat({'PRE_D2'},[length(SPA_sel_D2),1]); repmat({'POS_D2'},[length(SPAPos_sel_D2),1])];

subplot(1,3,3); hold on
boxplot (DATA,GROUP,'Colors','brbrbrbrbr');
title('NIN')
ylabel('Gammapower-Baseline (dB)')
xlabel('Epoch')
% ylim([-7 7])
clear DATA GROUP SPA_sel_F SPAPos_sel_F SPA_sel_C SPAPos_sel_C SPA_sel_D1 SPAPos_sel_D1 SPA_sel_S SPAPos_sel_S SPA_sel_D2 SPAPos_sel_D2;

