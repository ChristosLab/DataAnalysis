%% For MSNG only

clear all; close all; clc; 

BrainArea ='PPC';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NS Neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NS_BS='NS';
[~,~,raw_NS] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\NS_BS_Database_20231109a.xlsx',['MSNG_',BrainArea,'_',NS_BS]);
if strcmp(raw_NS{1,1}, 'Filename')
    raw_NS=raw_NS(2:end, :);
end
%%
%%% Change Parameter %%%
% Area = 'PFC'
TI = 0.05; % sliding time interval(sec) for ROC
TW = 0.5; % Time Window (sec) for ROC
% TW = 3; % Time Window (sec) for ROC
startt = -1;
% startt = 0.5+ TW/2;
endt= 5;
% endt= 0.5+ TW/2;
numer=2; % use neurons with more than this number of errors >numer
%%%%%%%%%%%%%%%%%%%%%%%%
    
time = startt-TW:TI:endt-TW;
results_allHigh_NS(1:length(raw_NS),1:length(time))=NaN;
results_allLow_NS(1:length(raw_NS),1:length(time))=NaN;

%%
for n = 1:size(raw_NS)
    filename = raw_NS{n,1};%correct trial name
    [ROCarea_all_High, ntrsHigh, err_ntrsHigh,ROCarea_all_Low, ntrsLow, err_ntrsLow] = FS_Neuron_Data_ROCarray_onefileErrtria_MNM(filename,TI,TW,startt,endt,numer);
    
    results_allHigh_NS(n,1:length(time)) =ROCarea_all_High;
   
    allntrs_High(n,1:2) = [ntrsHigh err_ntrsHigh];
  
end
%%
%%% Population analysis %%%

popresults1_NS = nanmean(results_allHigh_NS); % Calculate average ROC value in each bin
sem1_NS = nanstd(results_allHigh_NS)./sqrt(mean(sum(~isnan(results_allHigh_NS))));

nBest_NS = mean(sum(~isnan(results_allHigh_NS)));

delay_roc_NS  = mean(results_allHigh_NS(:, 32:91), 2);
delay_roc_NS  = delay_roc_NS(find(~isnan(delay_roc_NS)));
sample_roc_NS = mean(results_allHigh_NS(:, 92:101), 2);
sample_roc_NS  = sample_roc_NS(find(~isnan(sample_roc_NS)));

%%
clearvars -except results_allHigh_NS popresults1_NS sem1_NS nBest_NS BrainArea delay_roc_NS sample_roc_NS TW

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BS Neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NS_BS='BS';
[~,~,raw_BS] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\NS_BS_Database_20231109a.xlsx',['MSNG_',BrainArea,'_',NS_BS]);
if strcmp(raw_BS{1,1}, 'Filename')
    raw_BS=raw_BS(2:end, :);
end
%%
%%% Change Parameter %%%
% Area = 'PFC'
TI = 0.05; % sliding time interval(sec) for ROC
TW = 0.5; % Time Window (sec) for ROC
startt = -1;
endt=5;
numer=2; % use neurons with more than this number of errors >numer
%%%%%%%%%%%%%%%%%%%%%%%%
    
time = startt-TW:TI:endt-TW;
results_allHigh_BS(1:length(raw_BS),1:length(time))=NaN;
results_allLow_BS(1:length(raw_BS),1:length(time))=NaN;

%%
for n = 1:size(raw_BS)
    filename = raw_BS{n,1};%correct trial name
    [ROCarea_all_High, ntrsHigh, err_ntrsHigh,ROCarea_all_Low, ntrsLow, err_ntrsLow] = FS_Neuron_Data_ROCarray_onefileErrtria_MNM(filename,TI,TW,startt,endt,numer);
    
    results_allHigh_BS(n,1:length(time)) =ROCarea_all_High;
   
    allntrs_High(n,1:2) = [ntrsHigh err_ntrsHigh];
  
end
%%
%%% Population analysis %%%
popresults2_BS = nanmean(results_allHigh_BS); % Calculate average ROC value in each bin
sem2_BS = nanstd(results_allHigh_BS)./sqrt(mean(sum(~isnan(results_allHigh_BS))));

nBest_BS = mean(sum(~isnan(results_allHigh_BS)));

delay_roc_BS  = mean(results_allHigh_BS(:, 32:91), 2);
delay_roc_BS  = delay_roc_BS(find(~isnan(delay_roc_BS)));
sample_roc_BS = mean(results_allHigh_BS(:, 92:101), 2);
sample_roc_BS  = sample_roc_BS(find(~isnan(sample_roc_BS)));

%%
clearvars -except results_allHigh_NS popresults1_NS sem1_NS nBest_NS BrainArea results_allHigh_BS popresults2_BS ...
    sem2_BS nBest_BS time delay_roc_NS sample_roc_NS delay_roc_BS sample_roc_BS TW

%%

colors2 = {[0.4 0.9 1], [0.6 0.9 1]};%dark blue and light blue
colors1 = {[1, 0.6, 0.6], [1, 0.8, 0.8]}; % red
h = figure;
hold on;
% para.legends = {'High','Low'};
% legend(para.legends)
[ha_1] = shadedplot(time+TW,popresults1_NS+sem1_NS,popresults1_NS-sem1_NS,colors1{2},colors1{2});
[ha_2] = shadedplot(time+TW,popresults2_BS+sem2_BS,popresults2_BS-sem2_BS,colors2{2},colors2{2});

plot(time+TW,popresults1_NS,'Color','r','LineWidth',2);
plot(time+TW,popresults2_BS,'Color','b','LineWidth',2);



ha_1(2).DisplayName = ['NS'];
ha_2(2).DisplayName = ['BS'];

line([0 0], [0 1],'color','k','LineStyle','--')
line([.5 .5], [0 1],'color','k','LineStyle','--')
line([3.5 3.5], [0 1],'color','k','LineStyle','--')
% line([4 4], [0 1],'color','k','LineStyle','--')
line([-1 5], [0.5 0.5],'color','k','LineWidth',1)

set(gca, 'XLim',[-1 4], 'YLim',[0.3 0.7])
legend([ha_1(2), ha_2(2)])
title(sprintf([BrainArea,' Neurons (NS=%d,  BS=%d)'], nBest_NS, nBest_BS),'FontSize',15)
xlabel('Time (s)','FontSize',15)
ylabel('Choice probability','FontSize',15)


%% stat testing

disp('NS delay vs 0.5')
[h,p,ci,stats] = ttest(mean(delay_roc_NS,2), 0.5)

disp('NS sample vs 0.5')
[h,p,ci,stats] = ttest(mean(sample_roc_NS,2), 0.5)

disp('BS delay vs 0.5')
[h,p,ci,stats] = ttest(mean(delay_roc_BS,2), 0.5)

disp('BS sample vs 0.5')
[h,p,ci,stats] = ttest(mean(sample_roc_BS,2), 0.5)

disp('NS delay vs BS delay')
[h,p,ci,stats] = ttest2(mean(delay_roc_NS,2), mean(delay_roc_BS,2))

disp('NS sample vs BS sample')
[h,p,ci,stats] = ttest2(mean(sample_roc_NS,2), mean(sample_roc_BS,2))

%%
% disp('NS delay vs 0.5')
% % [h,p,ci,stats] = ttest(results_allHigh_NS, 0.5)
% 
% disp('NS sample vs 0.5')
% [h,p,ci,stats] = ttest(results_allHigh_NS, 0.5)
% 
% disp('BS delay vs 0.5')
% % [h,p,ci,stats] = ttest(results_allHigh_BS, 0.5)
% 
% disp('BS sample vs 0.5')
% [h,p,ci,stats] = ttest(results_allHigh_BS, 0.5)
% 
% disp('NS delay vs BS delay')
% % [h,p,ci,stats] = ttest2(results_allHigh_NS, results_allHigh_BS)
% 
% disp('NS sample vs BS sample')
% [h,p,ci,stats] = ttest2(results_allHigh_NS, results_allHigh_BS)
 

