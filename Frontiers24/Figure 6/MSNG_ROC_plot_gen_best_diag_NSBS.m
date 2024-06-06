clear all; close all; clc; 

BrainArea='PFC'; % PFC or PPC
NS_BS='NS';

[~,~,raw1] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\NS_BS_Database_20231109a.xlsx',['MSNG_',BrainArea, '_' ,NS_BS]); 
if strcmp(raw1{1,1}, 'Filename')
    raw1=raw1(2:end, :);
end

%% Change Parameter %%%
TI = 0.05; % sliding time interval(sec) for ROC
TW = 0.5; % Time Window (sec) for ROC
startt = -1.0;
endt=5;
time = startt-TW:TI:endt-TW;
results_all_NS(1:length(raw1),1:length(time))=NaN; %for correct
results_all_err(1:length(raw1),1:length(time)) =NaN; %for error

%%
for n = 1:size(raw1)
    filename = [raw1{n,1}, '.mat'];
%     [ROCarea_corr, ntrsBest_corr, ntrsDiag_corr] = ROCclass_best_diag(filename,TI,TW,startt,endt,r1r2,1); %1 for correct 
    [ROCarea_corr, ntrsBest_corr, ntrsDiag_corr] = ROCntr_best_worst(filename,TI,TW,startt,endt); %1 for correct 
    results_all_NS(n,1:length(time)) =ROCarea_corr;
    allntrs_NS(n,1:2) = [ntrsBest_corr ntrsDiag_corr];

    clear ROCarea_corr ntrsBest_corr ntrsDiag_corr
end
%%

popresults1 = nanmean(results_all_NS); % Calculate average ROC value in each bin
sem1 = nanstd(results_all_NS)./sqrt(mean(sum(~isnan(results_all_NS))));
n_NS = mean(sum(~isnan(results_all_NS)));

delay_roc_NS  = mean(results_all_NS(:, 32:91), 2);
delay_roc_NS  = delay_roc_NS(find(~isnan(delay_roc_NS)));
sample_roc_NS = mean(results_all_NS(:, 92:101), 2);
sample_roc_NS  = sample_roc_NS(find(~isnan(sample_roc_NS)));

%%
NS_BS='BS';

[~,~,raw2] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\NS_BS_Database_20231109a.xlsx',['MSNG_',BrainArea, '_' ,NS_BS]); 
if strcmp(raw2{1,1}, 'Filename')
    raw2=raw2(2:end, :);
end

%% Change Parameter %%%
TI = 0.05; % sliding time interval(sec) for ROC
TW = 0.5; % Time Window (sec) for ROC
startt = -1.0;
endt=5;
time = startt-TW:TI:endt-TW;
results_all_BS(1:length(raw2),1:length(time)) =NaN; %for error

%%
for n = 1:size(raw2)
    filename = [raw2{n,1}, '.mat'];
%     [ROCarea_corr, ntrsBest_corr, ntrsDiag_corr] = ROCclass_best_diag(filename,TI,TW,startt,endt,r1r2,1); %1 for correct 
    [ROCarea_corr, ntrsBest_corr, ntrsDiag_corr] = ROCntr_best_worst(filename,TI,TW,startt,endt); %1 for correct 
    results_all_BS(n,1:length(time)) =ROCarea_corr;
    allntrs_corr_BS(n,1:2) = [ntrsBest_corr ntrsDiag_corr];
    
    clear ROCarea_corr ntrsBest_corr ntrsDiag_corr
end
%%

popresults2 = nanmean(results_all_BS); % Calculate average ROC value in each bin
sem2 = nanstd(results_all_BS)./sqrt(mean(sum(~isnan(results_all_BS))));
n_BS = mean(sum(~isnan(results_all_BS)));

delay_roc_BS  = mean(results_all_BS(:, 32:91), 2);
delay_roc_BS  = delay_roc_BS(find(~isnan(delay_roc_BS)));
sample_roc_BS = mean(results_all_BS(:, 92:101), 2);
sample_roc_BS  = sample_roc_BS(find(~isnan(sample_roc_BS)));

%%
colors2 = {[1, 0.6, 0.6]};%red
colors1 = {[0.4 0.9 1]}; %blue
h = figure;
%adding color in between the outlined sections 
line([0 0], [0 1],'Color','k','LineStyle','--');
line([0.5 0.5], [0 1],'Color','k','LineStyle','--');
% line([2 2], [0 1],'Color','k','LineStyle','--');
line([3.5 3.5], [0 1],'Color','k','LineStyle','--');
% line([4 4], [0 1],'Color','k','LineStyle','--'); 
hold on;

shadedplot(time+TW,popresults1+sem1,popresults1-sem1,colors2{1},colors2{1});
shadedplot(time+TW,popresults2+sem2,popresults2-sem2,colors1{1},colors1{1});

p1=plot(time+TW,popresults1,'-r','LineWidth',2);
p2=plot(time+TW,popresults2,'-b','LineWidth',2);

% nNS = num2str(mean(sum(~isnan(results_all_NS))));
% nBS = num2str(mean(sum(~isnan(results_all_BS))));

line([-1 5], [0.5 0.5],'color','k','LineWidth',1)



set(gca, 'XLim',[-1 4], 'YLim',[0.3 0.8]) %1.15*max_max])
% set(gca, 'XLim',[-1 5], 'YLim',[0 1]) %1.15*max_max])
% legend([ha_1(2)])

if strcmp(BrainArea, 'PFC') 
    title(sprintf('MSNG PFC (NS=%d,  BS=%d)', n_NS, n_BS),'FontSize',15)
elseif strcmp(BrainArea, 'PPC')
    title(sprintf('MSNG PPC (NS=%d,  BS=%d)', n_NS, n_BS),'FontSize',15)
end

xlabel('Time(s)','FontSize',15)
ylabel('Area under ROC curve','FontSize',15)
legend([p1 p2],'NS','BS')

% set(h, 'Renderer', 'painters');
% exportgraphics(h, 'ROCntr_PPC_R2_bestVSdia_err.emf', 'ContentType', 'vector');
% print(h, 'ROCntr_PPC_R2_bestVSdia_err.pdf', '-vector', '-bestfit', '-dwinc');


%% stat testing
clc 

disp('NS delay1 vs 0.5')
[h,p,ci,stats] = ttest(mean(delay_roc_NS,2), 0.5)
meanDifference = mean(delay_roc_NS) - 0.5;
pooledStdDev = std(delay_roc_NS);
cohenD = meanDifference / pooledStdDev

disp('NS delay2 vs 0.5')
[h,p,ci,stats] = ttest(mean(sample_roc_NS,2), 0.5)
meanDifference = mean(sample_roc_NS) - 0.5;
pooledStdDev = std(sample_roc_NS);
cohenD = meanDifference / pooledStdDev

disp('BS delay1 vs 0.5')
[h,p,ci,stats] = ttest(mean(delay_roc_BS,2), 0.5)
meanDifference = mean(delay_roc_BS) - 0.5;
pooledStdDev = std(delay_roc_BS);
cohenD = meanDifference / pooledStdDev

disp('BS delay2 vs 0.5')
[h,p,ci,stats] = ttest(mean(sample_roc_BS,2), 0.5)
meanDifference = mean(sample_roc_BS) - 0.5;
pooledStdDev = std(sample_roc_BS);
cohenD = meanDifference / pooledStdDev

disp('NS delay1 vs BS delay1')
[h,p,ci,stats] = ttest2(mean(delay_roc_NS,2), mean(delay_roc_BS,2))
meanDifference = mean(delay_roc_NS) - mean(delay_roc_BS);
pooledStdDev = sqrt(((std(delay_roc_NS)^2 + std(delay_roc_BS)^2) / 2));
cohenD = meanDifference / pooledStdDev

disp('NS delay2 vs BS delay2')
[h,p,ci,stats] = ttest2(mean(sample_roc_NS,2), mean(sample_roc_BS,2))
meanDifference = mean(sample_roc_NS) - mean(sample_roc_BS);
pooledStdDev = sqrt(((std(sample_roc_NS)^2 + std(sample_roc_BS)^2) / 2));
cohenD = meanDifference / pooledStdDev




