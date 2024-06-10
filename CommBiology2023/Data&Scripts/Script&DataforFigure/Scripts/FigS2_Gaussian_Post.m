%% This script for Gaussian plot an example:  
clc; clear; close all;
default_dir=pwd;

timeStep    = .02; startTime=-0.748; endTime=4.7520;
timeVector  = timeStep*(startTime/timeStep:endTime/timeStep);
Fix  = [1:38]; Cue = [39:63]; Delay1 = [64:138]; Sample = [139:163]; Delay2 = [164:238]; Choice = [239:276];

freqsVector = importdata(['Data/Frequency.mat']);
Gamma       = [23:51];

%% Allregion

SPA_sel = importdata(['Data/DataForStatics/AllRegion/SessionbySession/Spatial_POSTTRAINING_SELECTIVE_LFP.mat']);
for ff=1:length(SPA_sel)
    TEMP = SPA_sel{ff};
    for ii=1:size(TEMP,2)
        SInducedDb(:,ii) = 10*log10(TEMP(:,ii));
    end
    SPA_sel_F(ff,:) = mean(mean(SInducedDb(Fix,Gamma),2)); 
    SPA_sel_Cu(ff,:) = mean(mean(SInducedDb(Cue,Gamma),2)); 
    SPA_sel_D1(ff,:) = mean(mean(SInducedDb(Delay1,Gamma),2));
    SPA_sel_S(ff,:) = mean(mean(SInducedDb(Sample,Gamma),2)); 
    SPA_sel_D2(ff,:) = mean(mean(SInducedDb(Delay2,Gamma),2));     
    clear TEMP SInducedDb;
end
clear SPA_sel;

SPA_par    = importdata(['Data/DataForStatics/AllRegion/SessionbySession/Spatial_POSTTRAINING_PARTIAL_SELECTIVE_LFP.mat']);
for ff=1:length(SPA_par)
    TEMP = SPA_par{ff};
    for ii=1:size(TEMP,2)
        SInducedDb(:,ii) = 10*log10(TEMP(:,ii));
    end
      
    SPA_par_F(ff,:) = mean(mean(SInducedDb(Fix,Gamma),2)); 
    SPA_par_Cu(ff,:) = mean(mean(SInducedDb(Cue,Gamma),2)); 
    SPA_par_D1(ff,:) = mean(mean(SInducedDb(Delay1,Gamma),2)); 
    SPA_par_S(ff,:) = mean(mean(SInducedDb(Sample,Gamma),2)); 
    SPA_par_D2(ff,:) = mean(mean(SInducedDb(Delay2,Gamma),2)); 
    clear TEMP SInducedDb;
end
clear SPA_par;

SPA_nonsel = importdata(['Data/DataForStatics/AllRegion/SessionbySession/Spatial_POSTTRAINING_NON_SELECTIVE_LFP.mat']);
for ff=1:length(SPA_nonsel)
    TEMP = SPA_nonsel{ff};
    for ii=1:size(TEMP,2)
        SInducedDb(:,ii) = 10*log10(TEMP(:,ii));
    end
     
    SPA_nonsel_F(ff,:) = mean(mean(SInducedDb(Fix,Gamma),2)); 
    SPA_nonsel_Cu(ff,:) = mean(mean(SInducedDb(Cue,Gamma),2)); 

    SPA_nonsel_D1(ff,:) = mean(mean(SInducedDb(Delay1,Gamma),2)); 
    SPA_nonsel_S(ff,:) = mean(mean(SInducedDb(Sample,Gamma),2)); 
    SPA_nonsel_D2(ff,:) = mean(mean(SInducedDb(Delay2,Gamma),2)); 
    clear TEMP SInducedDb;
end
clear SPA_nonsel;

%%

fig = figure; hold on
set(fig, 'Position', [30 40 850 650]); 
subplot(1,5,1);hold on;
h1 = histfit(SPA_sel_F);
h2 = histfit(SPA_par_F);
h3 = histfit(SPA_nonsel_F);
legend('Sel','','Par','','Nonsel')
title('Fix');


subplot(1,5,2);hold on;
h1 = histfit(SPA_sel_Cu);
h2 = histfit(SPA_par_Cu);
h3 = histfit(SPA_nonsel_Cu);
legend('Sel','','Par','','Nonsel')
title('Cue');




subplot(1,5,3);hold on;
h1 = histfit(SPA_sel_D1);
h2 = histfit(SPA_par_D1);
h3 = histfit(SPA_nonsel_D1);
legend('Sel','','Par','','Nonsel')
title('Delay1');



subplot(1,5,4);hold on;
h1 = histfit(SPA_sel_S);
h2 = histfit(SPA_par_S);
h3 = histfit(SPA_nonsel_S);
legend('Sel','','Par','','Nonsel')
title('Sample');




subplot(1,5,5);hold on;
h1 = histfit(SPA_sel_D2);
h2 = histfit(SPA_par_D2);
h3 = histfit(SPA_nonsel_D2);
legend('Sel','','Par','','Nonsel')
title('Delay2');



