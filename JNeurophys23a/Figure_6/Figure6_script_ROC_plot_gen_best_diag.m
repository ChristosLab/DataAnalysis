clear all; close all; clc; 

BrainArea='PFC'; % PFC or PPC
r1r2=2; % 1 for Remember 1st and 2 for Remember 2nd

[~,~,raw] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\Distractor Project\Rana Scripts\ODRdistVar Rana\ODRdistVar',[BrainArea, 'sigCD']); 

%% Change Parameter %%%
TI = 0.05; % sliding time interval(sec) for ROC
TW = 0.5; % Time Window (sec) for ROC
startt = -1.0;
endt=5;
time = startt:TI:endt;
results_all_corr(1:length(raw),1:length(time))=NaN; %for correct
results_all_err(1:length(raw),1:length(time)) =NaN; %for error

%%
for n = 1:size(raw)
    filename = [raw{n,1},'_',num2str(raw{n,2})];%correct trial name
    [ROCarea_corr, ntrsBest_corr, ntrsDiag_corr] = ROCclass_best_diag(filename,TI,TW,startt,endt,r1r2,1); %1 for correct 
    [ROCarea_err, ntrsBest_err, ntrsDiag_err] = ROCclass_best_diag(filename,TI,TW,startt,endt,r1r2,0); %0 for error
    
    results_all_corr(n,1:length(time)) =ROCarea_corr;
    results_all_err(n,1:length(time))  =ROCarea_err;
   
    allntrs_corr(n,1:2) = [ntrsBest_corr ntrsDiag_corr];
    allntrs_err(n,1:2) = [ntrsBest_err ntrsDiag_err];
  
end
%%

popresults1 = nanmean(results_all_corr); % Calculate average ROC value in each bin
sem1 = nanstd(results_all_corr)./sqrt(sum(~isnan(results_all_corr)));

popresults2 = nanmean(results_all_err); % Calculate average ROC value in each bin
sem2 = nanstd(results_all_err)./sqrt(sum(~isnan(results_all_err)));

colors2 = {[51/255 205/255 255/255], [0.678 0.922 1]};%dark blue and light blue
colors1 = {[0.7 0.7 0.7], [0.9 0.9 0.9]}; % dark gray and gray
h = figure;
hold on;
% para.legends = {'High','Low'};
% legend(para.legends)
if r1r2==1
    shadedplot(time,popresults1+sem1,popresults1-sem1,colors1{1},colors1{1});
    shadedplot(time,popresults2+sem2,popresults2-sem2,colors1{2},colors1{2});
elseif r1r2==2
    shadedplot(time,popresults1+sem1,popresults1-sem1,colors2{1},colors2{1});
    shadedplot(time,popresults2+sem2,popresults2-sem2,colors2{2},colors2{2});
end

p1=plot(time,popresults1,'-k','LineWidth',2);
p2=plot(time,popresults2,'--k','LineWidth',2);

nCorr = num2str(mean(sum(~isnan(results_all_corr))));
nErr  = num2str(mean(sum(~isnan(results_all_err))));

line([-1 5], [0.5 0.5],'color','k','LineWidth',1)

%adding color in between the outlined sections 
patch1 = patch([0, 0, 0.5, 0.5], [0, 0.22, 0.22 ,0], 'r'); 
set(patch1, 'FaceAlpha', 0.15) 
patch1 = patch([2, 2, 2.5, 2.5], [0, 0.22, 0.22 ,0], 'r'); 
set(patch1, 'FaceAlpha', 0.15) 
patch1 = patch([4, 4, 5, 5], [0, 0.22, 0.22 ,0], 'b'); 
set(patch1, 'FaceAlpha', 0.15)

set(gca, 'XLim',[-1 4.5], 'YLim',[0.2 0.9]) %1.15*max_max])
% set(gca, 'XLim',[-1 5], 'YLim',[0 1]) %1.15*max_max])
% legend([ha_1(2)])

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
ylabel('Area under ROC curve','FontSize',15)
legend([p1 p2],['correct(n=',nCorr,')'], ['error(n=',nErr,')'])

% set(h, 'Renderer', 'painters');
% exportgraphics(h, 'ROCntr_PPC_R2_bestVSdia_err.emf', 'ContentType', 'vector');
% print(h, 'ROCntr_PPC_R2_bestVSdia_err.pdf', '-vector', '-bestfit', '-dwinc');
%% one-sample ttest 

% correct
delay1_mat = nanmean(results_all_corr(:,32:61),2);
x = delay1_mat(~isnan(delay1_mat));

[h, p, ci, stats] = ttest(x,0.5)

delay2_mat = nanmean(results_all_corr(:,72:101),2);
y = delay2_mat(~isnan(delay2_mat));

[h, p, ci, stats] = ttest(y,0.5)

%error
delay1_mat = nanmean(results_all_err(:,32:61),2);
x = delay1_mat(~isnan(delay1_mat));

[h, p, ci, stats] = ttest(x,0.5)

delay2_mat = nanmean(results_all_err(:,72:101),2);
y = delay2_mat(~isnan(delay2_mat));

[h, p, ci, stats] = ttest(y,0.5)

%%
% figure
% subplot 211
% plot(time,nanmean(results_best_err),'-k','LineWidth',2); hold on;
% plot(time,nanmean(results_diag_err),'-b','LineWidth',2); hold on;
% % plot(time,nanmean(results_all_err),'--r','LineWidth',2); hold on;
% legend('Error Best ROC', 'Error Diametric ROC')
% 
% 
% subplot 212
% plot(time,nanmean(results_all_err),'--r','LineWidth',2); hold on;
% legend('Error Best vs Diametric ROC')
% 
% sgtitle('PPC remember 2nd')







