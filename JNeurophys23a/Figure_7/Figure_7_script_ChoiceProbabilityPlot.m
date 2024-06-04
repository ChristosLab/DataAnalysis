clear all; close all;
clc; clearvars;

BrainArea='PPC'; % PFC or PPC
r1r2=2; %1 for Remember 1st and 2 for Remember 2nd

[~,~,raw] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\Distractor Project\Rana Scripts\ODRdistVar Rana\ODRdistVar',[BrainArea, 'sigCD']); 
%%
%%% Change Parameter %%%
% Area = 'PFC'
TI = 0.05; % sliding time interval(sec) for ROC
TW = 0.5; % Time Window (sec) for ROC
startt = -1.0;
endt=5;
numer=0; % use neurons with more than this number of errors >numer
%%%%%%%%%%%%%%%%%%%%%%%%

time = startt:TI:endt;
results_allHigh(1:length(raw),1:length(time))=NaN;
results_allLow(1:length(raw),1:length(time))=NaN;

%%
for n = 1:size(raw)
    filename = [raw{n,1},'_',num2str(raw{n,2})];%correct trial name
    [ROCarea_all_High, ntrsHigh, err_ntrsHigh,ROCarea_all_Low, ntrsLow, err_ntrsLow] = Neuron_Data_ROCarray_onefileErrtria_MNM(filename,TI,TW,startt,endt,numer, r1r2);
    
    results_allHigh(n,1:length(time)) =ROCarea_all_High;
    results_allLow(n,1:length(time)) = ROCarea_all_Low;
   
    allntrs_High(n,1:2) = [ntrsHigh err_ntrsHigh];
    allntrs_Low(n,1:2) =  [ntrsLow err_ntrsLow];
  
end
%%
%%% Population analysis %%%

popresults1 = nanmean(results_allHigh); % Calculate average ROC value in each bin
sem1 = nanstd(results_allHigh)./sqrt(sum(~isnan(results_allHigh)));

popresults2 = nanmean(results_allLow); % Calculate average ROC value in each bin
sem2 = nanstd(results_allLow)./sqrt(sum(~isnan(results_allLow)));

% colors = {'b','r','g','m'};
colors2 = {[51/255 205/255 255/255], [0.678 0.922 1]};%dark blue and light blue
colors1 = {[0.7 0.7 0.7], [0.9 0.9 0.9]}; % dark gray and gray
h = figure;
hold on;
% para.legends = {'High','Low'};
% legend(para.legends)
if r1r2==1
    [ha_1] = shadedplot(time,popresults1+sem1,popresults1-sem1,colors1{1},colors1{1});
    [ha_2] = shadedplot(time,popresults2+sem2,popresults2-sem2,colors1{2},colors1{2});
elseif r1r2==2
    [ha_1] = shadedplot(time,popresults1+sem1,popresults1-sem1,colors2{1},colors2{1});
    [ha_2] = shadedplot(time,popresults2+sem2,popresults2-sem2,colors2{2},colors2{2});
end

plot(time,popresults1,'-k','LineWidth',2);
plot(time,popresults2,'--k','LineWidth',2);

nBest = num2str(mean(sum(~isnan(results_allHigh))));
nDia  = num2str(mean(sum(~isnan(results_allLow))));

ha_1(2).DisplayName = ['Best Cue Condition (n=', nBest, ')'];
ha_2(2).DisplayName = ['Diametric Cue Condition (n=', nDia, ')'];

% % plot(cuetime,popresults_go3,':g','LineWidth',2);
% line([0 0], [0 1],'color','k','LineWidth',1, 'LineStyle', '--')
% line([.5 .5], [0 1],'color','k','LineWidth',1, 'LineStyle', '--')
% line([2 2], [0 1],'color','k','LineWidth',1, 'LineStyle', '--')
% line([2.5 2.5], [0 1],'color','k','LineWidth',1, 'LineStyle', '--')
% line([4 4], [0 1],'color','k','LineWidth',1, 'LineStyle', '--')
line([-1 5], [0.5 0.5],'color','k','LineWidth',1)

%adding color in between the outlined sections 
patch1 = patch([0, 0, 0.5, 0.5], [0, 1, 1 ,0], 'r'); 
set(patch1, 'FaceAlpha', 0.15) 
patch1 = patch([2, 2, 2.5, 2.5], [0, 1, 1 ,0], 'r'); 
set(patch1, 'FaceAlpha', 0.15) 
patch1 = patch([4, 4, 5, 5], [0, 1, 1 ,0], 'b'); 
set(patch1, 'FaceAlpha', 0.15)

set(gca, 'XLim',[-1 4.5], 'YLim',[0.3 0.7]) %1.15*max_max])
% set(gca, 'XLim',[-1 5], 'YLim',[0 1]) %1.15*max_max])
legend([ha_1(2), ha_2(2)])
% legend('Preferred','Nonpreferred')

if strcmp(BrainArea, 'PFC') && r1r2==1
    title('PFC Remember 1st (All Classes)','FontSize',15)
elseif strcmp(BrainArea, 'PPC') && r1r2==1
    title('PPC Remember 1st (All Classes)','FontSize',15)
elseif strcmp(BrainArea, 'PPC') && r1r2==2
    title('PPC Remember 2nd (All Classes)','FontSize',15)
elseif strcmp(BrainArea, 'PFC') && r1r2==2
    title('PFC Remember 2nd (All Classes)','FontSize',15)
end

xlabel('Time(s)','FontSize',15)
ylabel('Choice probability','FontSize',15)

%% one-sample ttest 
delay1_mat = nanmean(results_allHigh(:,32:61),2);
x = delay1_mat(~isnan(delay1_mat));

[h, p, ci, stats] = ttest(x,0.5)

delay2_mat = nanmean(results_allHigh(:,72:101),2);
y = delay2_mat(~isnan(delay2_mat));

[h, p, ci, stats] = ttest(y,0.5)

%
delay1_mat = nanmean(results_allLow(:,32:61),2);
x = delay1_mat(~isnan(delay1_mat));

[h, p, ci, stats] = ttest(x,0.5)

delay2_mat = nanmean(results_allLow(:,72:101),2);
y = delay2_mat(~isnan(delay2_mat));

[h, p, ci, stats] = ttest(y,0.5)





