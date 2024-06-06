clear all; close all; clc; 

BrainArea ='PFC';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NS Neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sig = 1;
NS_BS='NS';
if sig
    [~,~,raw_NS] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\AllSigNeuronsInfo.xlsx',['AllSigNeurons_', BrainArea,'_',NS_BS]);
else
    [~,~,raw_NS] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\NS_BS_Database_20231109a.xlsx',['AllNeurons_', BrainArea,'_',NS_BS]);
end
if strcmp(raw_NS{1,1}, 'Filename')
    raw_NS=raw_NS(2:end, :);
end

%%
%%% Change Parameter %%%
% Area = 'PFC'
TI = 0.05; % sliding time interval(sec) for ROC
TW = 0.5;  % Time Window (sec) for ROC
% TW = 1.5;  % Time Window (sec) for ROC
startt = -1;
% startt = 2.5+TW/2;
endt=2;  %5;
% endt=2.5+TW/2;
numer=1; % use neurons with more than this number of errors >numer
%%%%%%%%%%%%%%%%%%%%%%%%
    
time = startt-TW:TI:endt-TW;
% time = startt:TI:endt;
results_allHigh_NS(1:length(raw_NS),1:length(time))=NaN;

%%
for n = 1:size(raw_NS)
    filename = raw_NS{n,1};%correct trial name
    if strcmpi('MSNG', raw_NS{n, 8})
        [ROCarea_all_High, ntrsHigh, err_ntrsHigh,~, ~, ~] = FS_Neuron_Data_ROCarray_onefileErrtria_MNM(filename,TI,TW,startt,endt,numer);
    else
        [ROCarea_all_High, ntrsHigh, err_ntrsHigh] = ROC_NSBS_together(filename,TI,TW,startt,endt,numer);
    end
    results_allHigh_NS(n,1:length(time)) =ROCarea_all_High;
end
%%
%%% Population analysis %%%

popresults1_NS = nanmean(results_allHigh_NS); % Calculate average ROC value in each bin
sem1_NS = nanstd(results_allHigh_NS)./sqrt(mean(sum(~isnan(results_allHigh_NS))));

nBest_NS = mean(sum(~isnan(results_allHigh_NS)));

delay1_roc_NS = mean(results_allHigh_NS(:, 32:end),2);
delay1_roc_NS = delay1_roc_NS(find(~isnan(delay1_roc_NS)));

%%
% pooled_roc_value_NS = results_allHigh_NS(find(~isnan(results_allHigh_NS)));

%%
clearvars -except results_allHigh_NS popresults1_NS sem1_NS nBest_NS BrainArea delay1_roc_NS delay2_roc_NS TW sig

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BS Neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NS_BS='BS';
if sig
    [~,~,raw_BS] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\AllSigNeuronsInfo.xlsx',['AllSigNeurons_', BrainArea,'_',NS_BS]);
else
    [~,~,raw_BS] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\NS_BS_Database_20231109a.xlsx',['AllNeurons_', BrainArea,'_',NS_BS]);
end

if strcmp(raw_BS{1,1}, 'Filename')
    raw_BS=raw_BS(2:end, :);
end
%%
%%% Change Parameter %%%
% Area = 'PFC'
TI = 0.05; % sliding time interval(sec) for ROC
TW = 0.5;  % Time Window (sec) for ROC
startt = -1;
endt=2; %5;
numer=1; % use neurons with more than this number of errors >numer
%%%%%%%%%%%%%%%%%%%%%%%%
    
time = startt-TW:TI:endt-TW;
% time = startt:TI:endt;
results_allHigh_BS(1:length(raw_BS),1:length(time))=NaN;

%%
for n = 1:size(raw_BS)
    filename = raw_BS{n,1};%correct trial name
    if strcmpi('MSNG', raw_BS{n,8})
        [ROCarea_all_High, ntrsHigh, err_ntrsHigh,~, ~, ~] = FS_Neuron_Data_ROCarray_onefileErrtria_MNM(filename,TI,TW,startt,endt,numer);
    else
        [ROCarea_all_High, ntrsHigh, err_ntrsHigh] = ROC_NSBS_together(filename,TI,TW,startt,endt,numer);
    end
    results_allHigh_BS(n,1:length(time)) =ROCarea_all_High;
end
%%
%%% Population analysis %%%
popresults2_BS = nanmean(results_allHigh_BS); % Calculate average ROC value in each bin
sem2_BS = nanstd(results_allHigh_BS)./sqrt(mean(sum(~isnan(results_allHigh_BS))));

nBest_BS = mean(sum(~isnan(results_allHigh_BS)));

delay1_roc_BS = mean(results_allHigh_BS(:, 32:end),2);
delay1_roc_BS = delay1_roc_BS(find(~isnan(delay1_roc_BS)));

%%
% pooled_roc_value_BS = results_allHigh_BS(find(~isnan(results_allHigh_BS)));

%%
clearvars -except results_allHigh_NS popresults1_NS sem1_NS nBest_NS BrainArea results_allHigh_BS popresults2_BS ...
    sem2_BS nBest_BS time delay1_roc_NS delay2_roc_NS delay1_roc_BS delay2_roc_BS TW

%%

colors2 = {[0.4 0.9 1], [0.6 0.9 1]};%dark blue and light blue
colors1 = {[1, 0.6, 0.6], [1, 0.8, 0.8]}; % red
h = figure;
hold on;
% para.legends = {'High','Low'};
% legend(para.legends)
[ha_1] = shadedplot(time+TW,popresults1_NS+sem1_NS,popresults1_NS-sem1_NS,colors1{2},colors1{2});
% [ha_1] = shadedplot(time,popresults1_NS+sem1_NS,popresults1_NS-sem1_NS,colors1{2},colors1{2});
[ha_2] = shadedplot(time+TW,popresults2_BS+sem2_BS,popresults2_BS-sem2_BS,colors2{2},colors2{2});
% [ha_2] = shadedplot(time,popresults2_BS+sem2_BS,popresults2_BS-sem2_BS,colors2{2},colors2{2});

plot(time+TW,popresults1_NS,'Color','r','LineWidth',2);
% plot(time,popresults1_NS,'Color','r','LineWidth',2);
plot(time+TW,popresults2_BS,'Color','b','LineWidth',2);
% plot(time,popresults2_BS,'Color','b','LineWidth',2);



ha_1(2).DisplayName = ['NS'];
ha_2(2).DisplayName = ['BS'];

line([0 0], [0 1],'color','k','LineStyle','--')
line([.5 .5], [0 1],'color','k','LineStyle','--')
line([2 2], [0 1],'color','k','LineStyle','--')
line([2.5 2.5], [0 1],'color','k','LineStyle','--')
line([4 4], [0 1],'color','k','LineStyle','--')
line([-1 5], [0.5 0.5],'color','k','LineWidth',1)

% set(gca, 'XLim',[-1 2], 'YLim',[0.4 0.6])
set(gca, 'XLim',[-1 2], 'YLim',[0.3 0.7])
legend([ha_1(2), ha_2(2)])
title(sprintf([BrainArea,' Neurons (NS=%d,  BS=%d)'], nBest_NS, nBest_BS),'FontSize',15)
xlabel('Time (s)','FontSize',15)
ylabel('Choice probability','FontSize',15)



%% stat testing
clc 

disp('NS delay1 vs 0.5')
[h,p,ci,stats] = ttest(delay1_roc_NS, 0.5)
meanDifference = mean(delay1_roc_NS) - 0.5;
pooledStdDev = std(delay1_roc_NS);
cohenD = meanDifference / pooledStdDev
% 
% disp('NS delay2 vs 0.5')
% [h,p,ci,stats] = ttest(delay2_roc_NS, 0.5)

disp('BS delay1 vs 0.5')
[h,p,ci,stats] = ttest(delay1_roc_BS, 0.5)
meanDifference = mean(delay1_roc_BS) - 0.5;
pooledStdDev = std(delay1_roc_BS);
cohenD = meanDifference / pooledStdDev
% 
% disp('BS delay2 vs 0.5')
% [h,p,ci,stats] = ttest(delay2_roc_BS, 0.5)

disp('NS delay1 vs BS delay1')
[h,p,ci,stats] = ttest2(delay1_roc_NS, delay1_roc_BS)
% Calculate Cohen's d
meanDifference = mean(delay1_roc_NS) - mean(delay1_roc_BS);
pooledStdDev = sqrt(((std(delay1_roc_NS)^2 + std(delay1_roc_BS)^2) / 2));
cohenD = meanDifference / pooledStdDev
% 
% disp('NS delay2 vs BS delay2')
% [h,p,ci,stats] = ttest2(delay2_roc_NS, mean(delay2_roc_BS,2))


%%
% clc 
% 
% disp('NS pooled vs 0.5')
% [h,p,ci,stats] = ttest(pooled_roc_value_NS, 0.5)
% 
% disp('BS pooled vs 0.5')
% [h,p,ci,stats] = ttest(pooled_roc_value_BS, 0.5)
% 
% disp('NS pooled vs BS pooled')
% [h,p,ci,stats] = ttest2(pooled_roc_value_NS, mean(pooled_roc_value_BS,2))












