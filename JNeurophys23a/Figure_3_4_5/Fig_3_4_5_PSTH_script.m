%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%  PFC sigCD error and correct trials %%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
[~,~,raw] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\Distractor Project\Rana Scripts\ODRdistVar Rana\ODRdistVar','PFCsigCD');
neuron = 0;
%%
for i = 1:length(raw)
% for i=1:130 %for nonsig neurons
    fn_err = [raw{i,1}, '_', num2str(raw{i,2}), '_err.mat'];
    fn_corr = [raw{i,1}, '_', num2str(raw{i,2}), '.mat'];
    % for nonsig neurons
%     x = raw{i,1};
%     fn_corr = [x(1:8), '_', x(10:13), '.mat'];
%     fn_err = [x(1:8), '_', x(10:13), '_err.mat'];

    MatData_err =load(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\Distractor Project\Extraction\Extracted Data\ODRdistVar_Error\', fn_err]);    

    MatData_corr =load(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\Distractor Project\ODRdis_Var_correct_Data\', fn_corr]);
    max_class_corr = Neuron_Data_Max(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\Distractor Project\ODRdis_Var_correct_Data\', fn_corr]);
    
    if max_class_corr(1) == 1
        Classes=[1 2 3 4 5 11 12 13 14 15 6 7 8 9 10 16 17 18 19 20];
    elseif max_class_corr(1) ==6
        Classes = [6 7 8 9 10 16 17 18 19 20 1 2 3 4 5 11 12 13 14 15];
    else
        continue;
    end

    bin_width = 0.05;  % 50 milliseconds bin
    dur = [-1, 5];
    bin_edges = dur(1):bin_width:dur(2); %
    n_bins_cue = numel(histcounts([], bin_edges));

    cls=0;
    if ~isempty(MatData_err.MatData) && ~isempty(MatData_corr.MatData)
        for n= Classes
            %%%%%%%%%%%%%%%% for error %%%%%%%%%%%%%%%%%%%%%%
            allTS_err = []; allTS_corr = [];
            m_counter_err = 0; m_counter_corr = 0;

            if n<=length(MatData_err.MatData.class) && n<=length(MatData_corr.MatData.class)
                if  ~isempty(MatData_err.MatData.class(n).ntr) && ~isempty(MatData_corr.MatData.class(n).ntr)
                    if any(mean([MatData_err.MatData.class(n).ntr.cuerate])) && any(mean([MatData_corr.MatData.class(n).ntr.cuerate]))
                    %if length(MatData_err.MatData.class(n).ntr)>1 && length(MatData_corr.MatData.class(n).ntr)>1
                        %%%%%%%%%%%%%%%% for error %%%%%%%%%%%%%%%%%%%%%%
                        for tr=1:length(MatData_err.MatData.class(n).ntr)
                            if ~isempty(MatData_err.MatData.class(n).ntr(tr).Cue_onT)
                                try
                                    TS        = MatData_err.MatData.class(n).ntr(tr).TS-MatData_err.MatData.class(n).ntr(tr).Cue_onT;
                                    allTS_err     = [allTS_err TS];
                                    m_counter_err = m_counter_err + 1;
                                    clear TS;
                                catch
                                    disp('No Spike time found!!')
                                end
                            end
                        end
                        ntrs_err    = m_counter_err;
                        %%%%%%%%%%%%%%%% for correct %%%%%%%%%%%%%%%%%%%%%%
                        for tr=1:length(MatData_corr.MatData.class(n).ntr)
                            if ~isempty(MatData_corr.MatData.class(n).ntr(tr).Cue_onT)
                                try
                                    TS        = MatData_corr.MatData.class(n).ntr(tr).TS-MatData_corr.MatData.class(n).ntr(tr).Cue_onT;
                                    allTS_corr     = [allTS_corr TS];
                                    m_counter_corr = m_counter_corr + 1;
                                    clear TS;
                                catch
                                    disp('No spike data found')
                                end
                            end
                        end
                        ntrs_corr    = m_counter_corr;
                    %end
                    end
                end
            end
            
            if ~isempty(allTS_err)
                psth_temp_err(cls + 1,:) = histc(allTS_err,bin_edges)/(bin_width*ntrs_err);
            else
                psth_temp_err(cls + 1,:) = ones(1, length(bin_edges))*nan;
            end

            if ~isempty(allTS_corr)
                psth_temp_corr(cls + 1,:) = histc(allTS_corr,bin_edges)/(bin_width*ntrs_corr);
            else
                psth_temp_corr(cls + 1,:) = ones(1, length(bin_edges))*nan;
            end
            cls = cls +1;
        end
        pfc1_best_dia_err(neuron+1, :)   = psth_temp_err(1, :);
        pfc1_best_para_err(neuron+1, :)  = psth_temp_err(2, :);
        pfc1_best_n_err(neuron+1, :)     = psth_temp_err(3, :);
        pfc1_best_best_err(neuron+1, :)  = psth_temp_err(4, :);
        pfc1_best_null_err(neuron+1, :)  = psth_temp_err(5, :);
        pfc2_best_dia_err(neuron+1, :)   = psth_temp_err(6, :);
        pfc2_best_para_err(neuron+1, :)  = psth_temp_err(7, :);
        pfc2_best_n_err(neuron+1, :)     = psth_temp_err(8, :);
        pfc2_best_best_err(neuron+1, :)  = psth_temp_err(9, :);
        pfc2_null_best_err(neuron+1, :)  = psth_temp_err(10, :);
        pfc1_dia_best_err(neuron+1, :)   = psth_temp_err(11, :);
        pfc1_dia_para_err(neuron+1, :)   = psth_temp_err(12, :);
        pfc1_dia_n_err(neuron+1, :)      = psth_temp_err(13, :);
        pfc1_dia_dia_err(neuron+1, :)    = psth_temp_err(14, :);
        pfc1_dia_null_err(neuron+1, :)   = psth_temp_err(15, :);
        pfc2_dia_best_err(neuron+1, :)   = psth_temp_err(16, :);
        pfc2_dia_para_err(neuron+1, :)   = psth_temp_err(17, :);
        pfc2_dia_n_err(neuron+1, :)      = psth_temp_err(18, :);
        pfc2_dia_dia_err(neuron+1, :)    = psth_temp_err(19, :);
        pfc2_null_dia_err(neuron+1, :)   = psth_temp_err(20, :);


        pfc1_best_dia_corr(neuron+1, :)   = psth_temp_corr(1, :);
        pfc1_best_para_corr(neuron+1, :)  = psth_temp_corr(2, :);
        pfc1_best_n_corr(neuron+1, :)     = psth_temp_corr(3, :);
        pfc1_best_best_corr(neuron+1, :)  = psth_temp_corr(4, :);
        pfc1_best_null_corr(neuron+1, :)  = psth_temp_corr(5, :);
        pfc2_best_dia_corr(neuron+1, :)   = psth_temp_corr(6, :);
        pfc2_best_para_corr(neuron+1, :)  = psth_temp_corr(7, :);
        pfc2_best_n_corr(neuron+1, :)     = psth_temp_corr(8, :);
        pfc2_best_best_corr(neuron+1, :)  = psth_temp_corr(9, :);
        pfc2_null_best_corr(neuron+1, :)  = psth_temp_corr(10, :);
        pfc1_dia_best_corr(neuron+1, :)   = psth_temp_corr(11, :);
        pfc1_dia_para_corr(neuron+1, :)   = psth_temp_corr(12, :);
        pfc1_dia_n_corr(neuron+1, :)      = psth_temp_corr(13, :);
        pfc1_dia_dia_corr(neuron+1, :)    = psth_temp_corr(14, :);
        pfc1_dia_null_corr(neuron+1, :)   = psth_temp_corr(15, :);
        pfc2_dia_best_corr(neuron+1, :)   = psth_temp_corr(16, :);
        pfc2_dia_para_corr(neuron+1, :)   = psth_temp_corr(17, :);
        pfc2_dia_n_corr(neuron+1, :)      = psth_temp_corr(18, :);
        pfc2_dia_dia_corr(neuron+1, :)    = psth_temp_corr(19, :);
        pfc2_null_dia_corr(neuron+1, :)   = psth_temp_corr(20, :);
        


        neuron = neuron+1;
        clear psth_temp_err;
        clear psth_temp_corr;

    else
        disp('Empty MatData File!!!');
    end
    clear MatData_err;
    clear MatData_corr;
    
end
%%
pfc1_best_corr = [pfc1_best_dia_corr; pfc1_best_para_corr; pfc1_best_n_corr; pfc1_best_best_corr; pfc1_best_null_corr];
pfc1_best_err = [pfc1_best_dia_err; pfc1_best_para_err; pfc1_best_n_err; pfc1_best_best_err; pfc1_best_null_err];

pfc1_dia_corr = [pfc1_dia_best_corr; pfc1_dia_para_corr; pfc1_dia_n_corr; pfc1_dia_dia_corr; pfc1_dia_null_corr];
pfc1_dia_err = [pfc1_dia_best_err; pfc1_dia_para_err; pfc1_dia_n_err; pfc1_dia_dia_err; pfc1_dia_null_err];

pfc2_best_corr = [pfc2_best_dia_corr; pfc2_best_para_corr; pfc2_best_n_corr; pfc2_best_best_corr; pfc2_null_best_corr];
pfc2_best_err = [pfc2_best_dia_err; pfc2_best_para_err; pfc2_best_n_err; pfc2_best_best_err; pfc2_null_best_err];

pfc2_dia_corr = [pfc2_dia_best_corr; pfc2_dia_para_corr; pfc2_dia_n_corr; pfc2_dia_dia_corr; pfc2_null_dia_corr];
pfc2_dia_err = [pfc2_dia_best_err; pfc2_dia_para_err; pfc2_dia_n_err; pfc2_dia_dia_err; pfc2_null_dia_err];


%%
PFC1_best_corr   = nanmean(pfc1_best_corr);
PFC1_best_corr_sem = nanstd(pfc1_best_corr)./sqrt(sum(~isnan(pfc1_best_corr)));

PFC1_best_err   = nanmean(pfc1_best_err);
PFC1_best_err_sem = nanstd(pfc1_best_err)./sqrt(sum(~isnan(pfc1_best_err)));

PFC1_dia_corr   = nanmean(pfc1_dia_corr);
PFC1_dia_corr_sem = nanstd(pfc1_dia_corr)./sqrt(sum(~isnan(pfc1_dia_corr)));

PFC1_dia_err   = nanmean(pfc1_dia_err);
PFC1_dia_err_sem = nanstd(pfc1_dia_err)./sqrt(sum(~isnan(pfc1_dia_err)));

PFC2_best_corr   = nanmean(pfc2_best_corr);
PFC2_best_corr_sem = nanstd(pfc2_best_corr)./sqrt(sum(~isnan(pfc2_best_corr)));

PFC2_best_err   = nanmean(pfc2_best_err);
PFC2_best_err_sem = nanstd(pfc2_best_err)./sqrt(sum(~isnan(pfc2_best_err)));

PFC2_dia_corr   = nanmean(pfc2_dia_corr);
PFC2_dia_corr_sem = nanstd(pfc2_dia_corr)./sqrt(sum(~isnan(pfc2_dia_corr)));

PFC2_dia_err   = nanmean(pfc2_dia_err);
PFC2_dia_err_sem = nanstd(pfc2_dia_err)./sqrt(sum(~isnan(pfc2_dia_err)));
%%
figure
time = bin_edges;
colors2 = {[51/255 205/255 255/255], [0.678 0.922 1]};%dark blue and light blue
colors1 = {[0.7 0.7 0.7], [0.9 0.9 0.9]}; % dark gray and gray
popresults1 = smooth(PFC1_best_corr, 5)';
popresults2 = smooth(PFC1_best_err, 5)';
subplot(221)

[ha_2] = shadedplot(time,popresults2+PFC1_best_err_sem,popresults2-PFC1_best_err_sem,colors1{2},colors1{2});
[ha_1] = shadedplot(time,popresults1+PFC1_best_corr_sem,popresults1-PFC1_best_corr_sem,colors1{1},colors1{1});

plot(time, popresults1, 'LineWidth', 2,'Color', [0 0 0])
hold on
plot(time, popresults2, 'LineWidth', 2,'Color', [0 0 0], 'LineStyle','--')

ha_1(2).DisplayName = 'Correct';
ha_2(2).DisplayName = 'Error';

%HIGHLIGHTING SPECIFIC SECTIONS OF THE PSTH TIMELINE
patch1 = patch([0, 0, 0.5, 0.5], [0, 45, 45 ,0], 'r'); 
set(patch1, 'FaceAlpha', 0.15) 
patch1 = patch([2, 2, 2.5, 2.5], [0, 45, 45 ,0], 'r'); 
set(patch1, 'FaceAlpha', 0.15) 
patch1 = patch([4, 4, 5, 5], [0, 45, 45 ,0], 'b'); 
set(patch1, 'FaceAlpha', 0.15) 
legend([ha_1(2), ha_2(2)])
title(sprintf('PFC Remember 1st (Best) (Nsolid=%d, Ndashed=%d)', mean(sum(~isnan(pfc1_best_corr))), mean(sum(~isnan(pfc1_best_err)))))
xlim([-1 4.5])
ylim([0 45])
xlabel('Time(s)')
ylabel('Discharge Rate(sp/s)')


%
subplot(222)

popresults1 = smooth(PFC1_dia_corr, 5)';
popresults2 = smooth(PFC1_dia_err, 5)';

[ha_2] = shadedplot(time,popresults2+PFC1_dia_err_sem,popresults2-PFC1_dia_err_sem,colors1{2},colors1{2});
[ha_1] = shadedplot(time,popresults1+PFC1_dia_corr_sem,popresults1-PFC1_dia_corr_sem,colors1{1},colors1{1});

plot(time, popresults1, 'LineWidth', 2,'Color', [0 0 0])
hold on
plot(time, popresults2, 'LineWidth', 2,'Color', [0 0 0], 'LineStyle','--')

%adding color in between the outlined sections 
patch1 = patch([0, 0, 0.5, 0.5], [0, 45, 45 ,0], 'r'); 
set(patch1, 'FaceAlpha', 0.15) 
patch1 = patch([2, 2, 2.5, 2.5], [0, 45, 45 ,0], 'r'); 
set(patch1, 'FaceAlpha', 0.15) 
patch1 = patch([4, 4, 5, 5], [0, 45, 45 ,0], 'b'); 
set(patch1, 'FaceAlpha', 0.15) 

title(sprintf('PFC Remember 1st (Diametric) (Nsolid=%d, Ndashed=%d)', mean(sum(~isnan(pfc1_dia_corr))), mean(sum(~isnan(pfc1_dia_err)))))
xlim([-1 4.5])
ylim([0 45])
xlabel('Time(s)')
ylabel('Discharge Rate(sp/s)')
ha_1(2).DisplayName = 'Correct';
ha_2(2).DisplayName = 'Error';
legend([ha_1(2), ha_2(2)])


%

subplot(223)

popresults1 = smooth(PFC2_best_corr, 5)';
popresults2 = smooth(PFC2_best_err, 5)';

[ha_2] = shadedplot(time,popresults2+PFC2_best_err_sem,popresults2-PFC2_best_err_sem,colors2{2},colors2{2});
[ha_1] = shadedplot(time,popresults1+PFC2_best_corr_sem,popresults1-PFC2_best_corr_sem,colors2{1},colors2{1});

plot(time, popresults1, 'LineWidth', 2,'Color', [0 0 0])
hold on
plot(time, popresults2, 'LineWidth', 2,'Color', [0 0 0], 'LineStyle','--')
hold on
line([0 0], [0 45],'color','k')
line([0.5 0.5], [0 45],'color','k')
%second stimulus presentation
line([2 2], [0 45],'color','k')
line([2.5 2.5], [0 45],'color','k')
%choice period 
line([4 4], [0 45],'color','k')
line([5 5], [0 45],'color','k')
%adding color in between the outlined sections 
patch1 = patch([0, 0, 0.5, 0.5], [0, 45, 45 ,0], 'r'); 
set(patch1, 'FaceAlpha', 0.15) 
patch1 = patch([2, 2, 2.5, 2.5], [0, 45, 45 ,0], 'r'); 
set(patch1, 'FaceAlpha', 0.15) 
patch1 = patch([4, 4, 5, 5], [0, 45, 45 ,0], 'b'); 
set(patch1, 'FaceAlpha', 0.15) 

title(sprintf('PFC Remember 2nd (Best) (Nsolid=%d, Ndashed=%d)', mean(sum(~isnan(pfc2_best_corr))), mean(sum(~isnan(pfc2_best_err)))))
xlim([-1 4.5])
ylim([0 45])
xlabel('Time(s)')
ylabel('Discharge Rate(sp/s)')
ha_1(2).DisplayName = 'Correct';
ha_2(2).DisplayName = 'Error';
legend([ha_1(2), ha_2(2)])

%

subplot(224)

popresults1 = smooth(PFC2_dia_corr, 5)';
popresults2 = smooth(PFC2_dia_err, 5)';

[ha_2] = shadedplot(time,popresults2+PFC2_dia_err_sem,popresults2-PFC2_dia_err_sem,colors2{2},colors2{2});
[ha_1] = shadedplot(time,popresults1+PFC2_dia_corr_sem,popresults1-PFC2_dia_corr_sem,colors2{1},colors2{1});

plot(time, popresults1, 'LineWidth', 2,'Color', [0 0 0])
hold on
plot(time, popresults2, 'LineWidth', 2,'Color', [0 0 0], 'LineStyle','--')
hold on
line([0 0], [0 45],'color','k')
line([0.5 0.5], [0 45],'color','k')
%second stimulus presentation
line([2 2], [0 45],'color','k')
line([2.5 2.5], [0 45],'color','k')
%choice period 
line([4 4], [0 45],'color','k')
line([5 5], [0 45],'color','k')
%adding color in between the outlined sections 
patch1 = patch([0, 0, 0.5, 0.5], [0, 45, 45 ,0], 'r'); 
set(patch1, 'FaceAlpha', 0.15) 
patch1 = patch([2, 2, 2.5, 2.5], [0, 45, 45 ,0], 'r'); 
set(patch1, 'FaceAlpha', 0.15) 
patch1 = patch([4, 4, 5, 5], [0, 45, 45 ,0], 'b'); 
set(patch1, 'FaceAlpha', 0.15) 

title(sprintf('PFC Remember 2nd (Diametric) (Nsolid=%d, Ndashed=%d)', mean(sum(~isnan(pfc2_dia_corr))), mean(sum(~isnan(pfc2_dia_err)))))
xlim([-1 4.5])
ylim([0 45])
xlabel('Time(s)')
ylabel('Discharge Rate(sp/s)')
ha_1(2).DisplayName = 'Correct';
ha_2(2).DisplayName = 'Error';
legend([ha_1(2), ha_2(2)])

%%
clearvars
%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%  PFC nonsigCD error and correct trials %%%%%%%%%%%%%%%%%%%

[~,~,raw] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\Distractor Project\Rana Scripts\ODRdistVar Rana\ODRdistVar','PFCnonsigCD');
neuron = 0;
% for i = 1:length(raw)
for i=1:130
%     fn_err = [raw{i,1}, '_', num2str(raw{i,2}), '_err.mat'];
%     fn_corr = [raw{i,1}, '_', num2str(raw{i,2}), '.mat'];

    x = raw{i,1};
    fn_corr = [x(1:8), '_', x(10:13), '.mat'];
    fn_err = [x(1:8), '_', x(10:13), '_err.mat'];

    MatData_err =load(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\Distractor Project\Extraction\Extracted Data\ODRdistVar_Error\', fn_err]);

    MatData_corr =load(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\Distractor Project\ODRdis_Var_correct_Data\', fn_corr]);
    max_class_corr = Neuron_Data_Max(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\Distractor Project\ODRdis_Var_correct_Data\', fn_corr]);
    
    if max_class_corr(1) == 1
        Classes=[1 2 3 4 5 11 12 13 14 15 6 7 8 9 10 16 17 18 19 20];
    elseif max_class_corr(1) ==6
        Classes = [6 7 8 9 10 16 17 18 19 20 1 2 3 4 5 11 12 13 14 15];
    else
        continue;
    end

    bin_width = 0.05;  % 50 milliseconds bin
    dur = [-1, 5];
    bin_edges = dur(1):bin_width:dur(2); %
    n_bins_cue = numel(histcounts([], bin_edges));

    cls=0;
    if ~isempty(MatData_err.MatData) && ~isempty(MatData_corr.MatData)
        for n= Classes
            %%%%%%%%%%%%%%%% for error %%%%%%%%%%%%%%%%%%%%%%
            allTS_err = []; allTS_corr = [];
            m_counter_err = 0; m_counter_corr = 0;

            if n<=length(MatData_err.MatData.class) && n<=length(MatData_corr.MatData.class)
                if  ~isempty(MatData_err.MatData.class(n).ntr) && ~isempty(MatData_corr.MatData.class(n).ntr)
                    if any(mean([MatData_err.MatData.class(n).ntr.cuerate])) && any(mean([MatData_corr.MatData.class(n).ntr.cuerate]))
                    %if length(MatData_err.MatData.class(n).ntr)>1 && length(MatData_corr.MatData.class(n).ntr)>1
                        %%%%%%%%%%%%%%%% for error %%%%%%%%%%%%%%%%%%%%%%
                        for tr=1:length(MatData_err.MatData.class(n).ntr)
                            if ~isempty(MatData_err.MatData.class(n).ntr(tr).Cue_onT)
                                try
                                    TS        = MatData_err.MatData.class(n).ntr(tr).TS-MatData_err.MatData.class(n).ntr(tr).Cue_onT;
                                    allTS_err     = [allTS_err TS];
                                    m_counter_err = m_counter_err + 1;
                                    clear TS;
                                catch
                                    disp('No Spike time found!!')
                                end
                            end
                        end
                        ntrs_err    = m_counter_err;
                        %%%%%%%%%%%%%%%% for correct %%%%%%%%%%%%%%%%%%%%%%
                        for tr=1:length(MatData_corr.MatData.class(n).ntr)
                            if ~isempty(MatData_corr.MatData.class(n).ntr(tr).Cue_onT)
                                try
                                    TS        = MatData_corr.MatData.class(n).ntr(tr).TS-MatData_corr.MatData.class(n).ntr(tr).Cue_onT;
                                    allTS_corr     = [allTS_corr TS];
                                    m_counter_corr = m_counter_corr + 1;
                                    clear TS;
                                catch
                                    disp('No spike data found')
                                end
                            end
                        end
                        ntrs_corr    = m_counter_corr;
                    %end
                    end
                end
            end
            
            if ~isempty(allTS_err)
                psth_temp_err(cls + 1,:) = histc(allTS_err,bin_edges)/(bin_width*ntrs_err);
            else
                psth_temp_err(cls + 1,:) = ones(1, length(bin_edges))*nan;
            end

            if ~isempty(allTS_corr)
                psth_temp_corr(cls + 1,:) = histc(allTS_corr,bin_edges)/(bin_width*ntrs_corr);
            else
                psth_temp_corr(cls + 1,:) = ones(1, length(bin_edges))*nan;
            end
            cls = cls +1;
        end
        ppc1_best_dia_err(neuron+1, :)   = psth_temp_err(1, :);
        ppc1_best_para_err(neuron+1, :)  = psth_temp_err(2, :);
        ppc1_best_n_err(neuron+1, :)     = psth_temp_err(3, :);
        ppc1_best_best_err(neuron+1, :)  = psth_temp_err(4, :);
        ppc1_best_null_err(neuron+1, :)  = psth_temp_err(5, :);
        ppc2_best_dia_err(neuron+1, :)   = psth_temp_err(6, :);
        ppc2_best_para_err(neuron+1, :)  = psth_temp_err(7, :);
        ppc2_best_n_err(neuron+1, :)     = psth_temp_err(8, :);
        ppc2_best_best_err(neuron+1, :)  = psth_temp_err(9, :);
        ppc2_null_best_err(neuron+1, :)  = psth_temp_err(10, :);
        ppc1_dia_best_err(neuron+1, :)   = psth_temp_err(11, :);
        ppc1_dia_para_err(neuron+1, :)   = psth_temp_err(12, :);
        ppc1_dia_n_err(neuron+1, :)      = psth_temp_err(13, :);
        ppc1_dia_dia_err(neuron+1, :)    = psth_temp_err(14, :);
        ppc1_dia_null_err(neuron+1, :)   = psth_temp_err(15, :);
        ppc2_dia_best_err(neuron+1, :)   = psth_temp_err(16, :);
        ppc2_dia_para_err(neuron+1, :)   = psth_temp_err(17, :);
        ppc2_dia_n_err(neuron+1, :)      = psth_temp_err(18, :);
        ppc2_dia_dia_err(neuron+1, :)    = psth_temp_err(19, :);
        ppc2_null_dia_err(neuron+1, :)   = psth_temp_err(20, :);


        ppc1_best_dia_corr(neuron+1, :)   = psth_temp_corr(1, :);
        ppc1_best_para_corr(neuron+1, :)  = psth_temp_corr(2, :);
        ppc1_best_n_corr(neuron+1, :)     = psth_temp_corr(3, :);
        ppc1_best_best_corr(neuron+1, :)  = psth_temp_corr(4, :);
        ppc1_best_null_corr(neuron+1, :)  = psth_temp_corr(5, :);
        ppc2_best_dia_corr(neuron+1, :)   = psth_temp_corr(6, :);
        ppc2_best_para_corr(neuron+1, :)  = psth_temp_corr(7, :);
        ppc2_best_n_corr(neuron+1, :)     = psth_temp_corr(8, :);
        ppc2_best_best_corr(neuron+1, :)  = psth_temp_corr(9, :);
        ppc2_null_best_corr(neuron+1, :)  = psth_temp_corr(10, :);
        ppc1_dia_best_corr(neuron+1, :)   = psth_temp_corr(11, :);
        ppc1_dia_para_corr(neuron+1, :)   = psth_temp_corr(12, :);
        ppc1_dia_n_corr(neuron+1, :)      = psth_temp_corr(13, :);
        ppc1_dia_dia_corr(neuron+1, :)    = psth_temp_corr(14, :);
        ppc1_dia_null_corr(neuron+1, :)   = psth_temp_corr(15, :);
        ppc2_dia_best_corr(neuron+1, :)   = psth_temp_corr(16, :);
        ppc2_dia_para_corr(neuron+1, :)   = psth_temp_corr(17, :);
        ppc2_dia_n_corr(neuron+1, :)      = psth_temp_corr(18, :);
        ppc2_dia_dia_corr(neuron+1, :)    = psth_temp_corr(19, :);
        ppc2_null_dia_corr(neuron+1, :)   = psth_temp_corr(20, :);
        


        neuron = neuron+1;
        clear psth_temp_err;
        clear psth_temp_corr;

    else
        disp('Empty MatData File!!!');
    end
    clear MatData_err;
    clear MatData_corr;
    
end
%%
ppc1_best_corr = [ppc1_best_dia_corr; ppc1_best_para_corr; ppc1_best_n_corr; ppc1_best_best_corr; ppc1_best_null_corr];
ppc1_best_err = [ppc1_best_dia_err; ppc1_best_para_err; ppc1_best_n_err; ppc1_best_best_err; ppc1_best_null_err];

ppc1_dia_corr = [ppc1_dia_dia_corr; ppc1_dia_para_corr; ppc1_dia_n_corr; ppc1_dia_best_corr; ppc1_dia_null_corr];
ppc1_dia_err = [ppc1_dia_dia_err; ppc1_dia_para_err; ppc1_dia_n_err; ppc1_dia_best_err; ppc1_dia_null_err];

ppc2_best_corr = [ppc2_best_dia_corr; ppc2_best_para_corr; ppc2_best_n_corr; ppc2_best_best_corr; ppc2_null_best_corr];
ppc2_best_err = [ppc2_best_dia_err; ppc2_best_para_err; ppc2_best_n_err; ppc2_best_best_err; ppc2_null_best_err];

ppc2_dia_corr = [ppc2_dia_dia_corr; ppc2_dia_para_corr; ppc2_dia_n_corr; ppc2_dia_best_corr; ppc2_null_dia_corr];
ppc2_dia_err = [ppc2_dia_dia_err; ppc2_dia_para_err; ppc2_dia_n_err; ppc2_dia_best_err; ppc2_null_dia_err];


%%
PPC1_best_corr   = nanmean(ppc1_best_corr);
PPC1_best_corr_sem = nanstd(ppc1_best_corr)./sqrt(sum(~isnan(ppc1_best_corr)));

PPC1_best_err   = nanmean(ppc1_best_err);
PPC1_best_err_sem = nanstd(ppc1_best_err)./sqrt(sum(~isnan(ppc1_best_err)));

PPC1_dia_corr   = nanmean(ppc1_dia_corr);
PPC1_dia_corr_sem = nanstd(ppc1_dia_corr)./sqrt(sum(~isnan(ppc1_dia_corr)));

PPC1_dia_err   = nanmean(ppc1_dia_err);
PPC1_dia_err_sem = nanstd(ppc1_dia_err)./sqrt(sum(~isnan(ppc1_dia_err)));

PPC2_best_corr   = nanmean(ppc2_best_corr);
PPC2_best_corr_sem = nanstd(ppc2_best_corr)./sqrt(sum(~isnan(ppc2_best_corr)));

PPC2_best_err   = nanmean(ppc2_best_err);
PPC2_best_err_sem = nanstd(ppc2_best_err)./sqrt(sum(~isnan(ppc2_best_err)));

PPC2_dia_corr   = nanmean(ppc2_dia_corr);
PPC2_dia_corr_sem = nanstd(ppc2_dia_corr)./sqrt(sum(~isnan(ppc2_dia_corr)));

PPC2_dia_err   = nanmean(ppc2_dia_err);
PPC2_dia_err_sem = nanstd(ppc2_dia_err)./sqrt(sum(~isnan(ppc2_dia_err)));
%%
figure
time = bin_edges;
colors2 = {[51/255 205/255 255/255], [0.678 0.922 1]};%dark blue and light blue
colors1 = {[0.7 0.7 0.7], [0.9 0.9 0.9]}; % dark gray and gray
popresults1 = smooth(PPC1_best_corr, 5)';
popresults2 = smooth(PPC1_best_err, 5)';
subplot(221)

[ha_2] = shadedplot(time,popresults2+PPC1_best_err_sem,popresults2-PPC1_best_err_sem,colors1{2},colors1{2});
[ha_1] = shadedplot(time,popresults1+PPC1_best_corr_sem,popresults1-PPC1_best_corr_sem,colors1{1},colors1{1});

plot(time, popresults1, 'LineWidth', 2,'Color', [0 0 0])
hold on
plot(time, popresults2, 'LineWidth', 2,'Color', [0 0 0], 'LineStyle','--')

ha_1(2).DisplayName = 'Correct';
ha_2(2).DisplayName = 'Error';

%HIGHLIGHTING SPECIFIC SECTIONS OF THE PSTH TIMELINE
patch1 = patch([0, 0, 0.5, 0.5], [0, 45, 45 ,0], 'r'); 
set(patch1, 'FaceAlpha', 0.15) 
patch1 = patch([2, 2, 2.5, 2.5], [0, 45, 45 ,0], 'r'); 
set(patch1, 'FaceAlpha', 0.15) 
patch1 = patch([4, 4, 5, 5], [0, 45, 45 ,0], 'b'); 
set(patch1, 'FaceAlpha', 0.15) 
legend([ha_1(2), ha_2(2)])
title(sprintf('PPC Remember 1st (Best) (Nsolid=%d, Ndashed=%d)', mean(sum(~isnan(ppc1_best_corr))), mean(sum(~isnan(ppc1_best_err)))))
xlim([-1 4.5])
ylim([0 45])
xlabel('Time(s)')
ylabel('Discharge Rate(sp/s)')


%
subplot(222)

popresults1 = smooth(PPC1_dia_corr, 5)';
popresults2 = smooth(PPC1_dia_err, 5)';

[ha_2] = shadedplot(time,popresults2+PPC1_dia_err_sem,popresults2-PPC1_dia_err_sem,colors1{2},colors1{2});
[ha_1] = shadedplot(time,popresults1+PPC1_dia_corr_sem,popresults1-PPC1_dia_corr_sem,colors1{1},colors1{1});

plot(time, popresults1, 'LineWidth', 2,'Color', [0 0 0])
hold on
plot(time, popresults2, 'LineWidth', 2,'Color', [0 0 0], 'LineStyle','--')

%adding color in between the outlined sections 
patch1 = patch([0, 0, 0.5, 0.5], [0, 45, 45 ,0], 'r'); 
set(patch1, 'FaceAlpha', 0.15) 
patch1 = patch([2, 2, 2.5, 2.5], [0, 45, 45 ,0], 'r'); 
set(patch1, 'FaceAlpha', 0.15) 
patch1 = patch([4, 4, 5, 5], [0, 45, 45 ,0], 'b'); 
set(patch1, 'FaceAlpha', 0.15) 

title(sprintf('PPC Remember 1st (Diametric) (Nsolid=%d, Ndashed=%d)', mean(sum(~isnan(ppc1_dia_corr))), mean(sum(~isnan(ppc1_dia_err)))))
xlim([-1 4.5])
ylim([0 45])
xlabel('Time(s)')
ylabel('Discharge Rate(sp/s)')
ha_1(2).DisplayName = 'Correct';
ha_2(2).DisplayName = 'Error';
legend([ha_1(2), ha_2(2)])


%

subplot(223)

popresults1 = smooth(PPC2_best_corr, 5)';
popresults2 = smooth(PPC2_best_err, 5)';

[ha_2] = shadedplot(time,popresults2+PPC2_best_err_sem,popresults2-PPC2_best_err_sem,colors2{2},colors2{2});
[ha_1] = shadedplot(time,popresults1+PPC2_best_corr_sem,popresults1-PPC2_best_corr_sem,colors2{1},colors2{1});

plot(time, popresults1, 'LineWidth', 2,'Color', [0 0 0])
hold on
plot(time, popresults2, 'LineWidth', 2,'Color', [0 0 0], 'LineStyle','--')
hold on
line([0 0], [0 45],'color','k')
line([0.5 0.5], [0 45],'color','k')
%second stimulus presentation
line([2 2], [0 45],'color','k')
line([2.5 2.5], [0 45],'color','k')
%choice period 
line([4 4], [0 45],'color','k')
line([5 5], [0 45],'color','k')
%adding color in between the outlined sections 
patch1 = patch([0, 0, 0.5, 0.5], [0, 45, 45 ,0], 'r'); 
set(patch1, 'FaceAlpha', 0.15) 
patch1 = patch([2, 2, 2.5, 2.5], [0, 45, 45 ,0], 'r'); 
set(patch1, 'FaceAlpha', 0.15) 
patch1 = patch([4, 4, 5, 5], [0, 45, 45 ,0], 'b'); 
set(patch1, 'FaceAlpha', 0.15) 

title(sprintf('PPC Remember 2nd (Best) (Nsolid=%d, Ndashed=%d)', mean(sum(~isnan(ppc2_best_corr))), mean(sum(~isnan(ppc2_best_err)))))
xlim([-1 4.5])
ylim([0 45])
xlabel('Time(s)')
ylabel('Discharge Rate(sp/s)')
ha_1(2).DisplayName = 'Correct';
ha_2(2).DisplayName = 'Error';
legend([ha_1(2), ha_2(2)])

%

subplot(224)

popresults1 = smooth(PPC2_dia_corr, 5)';
popresults2 = smooth(PPC2_dia_err, 5)';

[ha_2] = shadedplot(time,popresults2+PPC2_dia_err_sem,popresults2-PPC2_dia_err_sem,colors2{2},colors2{2});
[ha_1] = shadedplot(time,popresults1+PPC2_dia_corr_sem,popresults1-PPC2_dia_corr_sem,colors2{1},colors2{1});

plot(time, popresults1, 'LineWidth', 2,'Color', [0 0 0])
hold on
plot(time, popresults2, 'LineWidth', 2,'Color', [0 0 0], 'LineStyle','--')
hold on
line([0 0], [0 45],'color','k')
line([0.5 0.5], [0 45],'color','k')
%second stimulus presentation
line([2 2], [0 45],'color','k')
line([2.5 2.5], [0 45],'color','k')
%choice period 
line([4 4], [0 45],'color','k')
line([5 5], [0 45],'color','k')
%adding color in between the outlined sections 
patch1 = patch([0, 0, 0.5, 0.5], [0, 45, 45 ,0], 'r'); 
set(patch1, 'FaceAlpha', 0.15) 
patch1 = patch([2, 2, 2.5, 2.5], [0, 45, 45 ,0], 'r'); 
set(patch1, 'FaceAlpha', 0.15) 
patch1 = patch([4, 4, 5, 5], [0, 45, 45 ,0], 'b'); 
set(patch1, 'FaceAlpha', 0.15) 

title(sprintf('PPC Remember 2nd (Diametric) (Nsolid=%d, Ndashed=%d)', mean(sum(~isnan(ppc2_dia_corr))), mean(sum(~isnan(ppc2_dia_err)))))
xlim([-1 4.5])
ylim([0 45])
xlabel('Time(s)')
ylabel('Discharge Rate(sp/s)')
ha_1(2).DisplayName = 'Correct';
ha_2(2).DisplayName = 'Error';
legend([ha_1(2), ha_2(2)])









%%
clearvars;
%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%  PPC sigCD error and correct trials %%%%%%%%%%%%%%%%%%%

[~,~,raw] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\Distractor Project\Rana Scripts\ODRdistVar Rana\ODRdistVar','PPCsigCD');
neuron = 0;
%%
for i = 1:length(raw)
    fn_err = [raw{i,1}, '_', num2str(raw{i,2}), '_err.mat'];
    MatData_err =load(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\Distractor Project\Extraction\Extracted Data\ODRdistVar_Error\', fn_err]);
    
    fn_corr = [raw{i,1}, '_', num2str(raw{i,2}), '.mat'];
    MatData_corr =load(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\Distractor Project\ODRdis_Var_correct_Data\', fn_corr]);
    max_class_corr = Neuron_Data_Max(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\Distractor Project\ODRdis_Var_correct_Data\', fn_corr]);
    
    if max_class_corr(1) == 1
        Classes=[1 2 3 4 5 11 12 13 14 15 6 7 8 9 10 16 17 18 19 20];
    elseif max_class_corr(1) ==6
        Classes = [6 7 8 9 10 16 17 18 19 20 1 2 3 4 5 11 12 13 14 15];
    else
        continue;
    end

    bin_width = 0.05;  % 50 milliseconds bin
    dur = [-1, 5];
    bin_edges = dur(1):bin_width:dur(2); %
    n_bins_cue = numel(histcounts([], bin_edges));

    cls=0;
    if ~isempty(MatData_err.MatData) && ~isempty(MatData_corr.MatData)
        for n= Classes
            %%%%%%%%%%%%%%%% for error %%%%%%%%%%%%%%%%%%%%%%
            allTS_err = []; allTS_corr = [];
            m_counter_err = 0; m_counter_corr = 0;

            if n<=length(MatData_err.MatData.class) && n<=length(MatData_corr.MatData.class)
                if  ~isempty(MatData_err.MatData.class(n).ntr) && ~isempty(MatData_corr.MatData.class(n).ntr)
                    if any(mean([MatData_err.MatData.class(n).ntr.cuerate])) && any(mean([MatData_corr.MatData.class(n).ntr.cuerate]))
                    %if length(MatData_err.MatData.class(n).ntr)>1 && length(MatData_corr.MatData.class(n).ntr)>1
                        %%%%%%%%%%%%%%%% for error %%%%%%%%%%%%%%%%%%%%%%
                        for tr=1:length(MatData_err.MatData.class(n).ntr)
                            if ~isempty(MatData_err.MatData.class(n).ntr(tr).Cue_onT)
                                try
                                    TS        = MatData_err.MatData.class(n).ntr(tr).TS-MatData_err.MatData.class(n).ntr(tr).Cue_onT;
                                    allTS_err     = [allTS_err TS];
                                    m_counter_err = m_counter_err + 1;
                                    clear TS;
                                catch
                                    disp('No Spike time found!!')
                                end
                            end
                        end
                        ntrs_err    = m_counter_err;
                        %%%%%%%%%%%%%%%% for correct %%%%%%%%%%%%%%%%%%%%%%
                        for tr=1:length(MatData_corr.MatData.class(n).ntr)
                            if ~isempty(MatData_corr.MatData.class(n).ntr(tr).Cue_onT)
                                try
                                    TS        = MatData_corr.MatData.class(n).ntr(tr).TS-MatData_corr.MatData.class(n).ntr(tr).Cue_onT;
                                    allTS_corr     = [allTS_corr TS];
                                    m_counter_corr = m_counter_corr + 1;
                                    clear TS;
                                catch
                                    disp('No spike data found')
                                end
                            end
                        end
                        ntrs_corr    = m_counter_corr;
                    %end
                    end
                end
            end
            
            if ~isempty(allTS_err)
                psth_temp_err(cls + 1,:) = histc(allTS_err,bin_edges)/(bin_width*ntrs_err);
            else
                psth_temp_err(cls + 1,:) = ones(1, length(bin_edges))*nan;
            end

            if ~isempty(allTS_corr)
                psth_temp_corr(cls + 1,:) = histc(allTS_corr,bin_edges)/(bin_width*ntrs_corr);
            else
                psth_temp_corr(cls + 1,:) = ones(1, length(bin_edges))*nan;
            end
            cls = cls +1;
        end
        ppc1_best_dia_err(neuron+1, :)   = psth_temp_err(1, :);
        ppc1_best_para_err(neuron+1, :)  = psth_temp_err(2, :);
        ppc1_best_n_err(neuron+1, :)     = psth_temp_err(3, :);
        ppc1_best_best_err(neuron+1, :)  = psth_temp_err(4, :);
        ppc1_best_null_err(neuron+1, :)  = psth_temp_err(5, :);
        ppc2_best_dia_err(neuron+1, :)   = psth_temp_err(6, :);
        ppc2_best_para_err(neuron+1, :)  = psth_temp_err(7, :);
        ppc2_best_n_err(neuron+1, :)     = psth_temp_err(8, :);
        ppc2_best_best_err(neuron+1, :)  = psth_temp_err(9, :);
        ppc2_null_best_err(neuron+1, :)  = psth_temp_err(10, :);
        ppc1_dia_best_err(neuron+1, :)   = psth_temp_err(11, :);
        ppc1_dia_para_err(neuron+1, :)   = psth_temp_err(12, :);
        ppc1_dia_n_err(neuron+1, :)      = psth_temp_err(13, :);
        ppc1_dia_dia_err(neuron+1, :)    = psth_temp_err(14, :);
        ppc1_dia_null_err(neuron+1, :)   = psth_temp_err(15, :);
        ppc2_dia_best_err(neuron+1, :)   = psth_temp_err(16, :);
        ppc2_dia_para_err(neuron+1, :)   = psth_temp_err(17, :);
        ppc2_dia_n_err(neuron+1, :)      = psth_temp_err(18, :);
        ppc2_dia_dia_err(neuron+1, :)    = psth_temp_err(19, :);
        ppc2_null_dia_err(neuron+1, :)   = psth_temp_err(20, :);


        ppc1_best_dia_corr(neuron+1, :)   = psth_temp_corr(1, :);
        ppc1_best_para_corr(neuron+1, :)  = psth_temp_corr(2, :);
        ppc1_best_n_corr(neuron+1, :)     = psth_temp_corr(3, :);
        ppc1_best_best_corr(neuron+1, :)  = psth_temp_corr(4, :);
        ppc1_best_null_corr(neuron+1, :)  = psth_temp_corr(5, :);
        ppc2_best_dia_corr(neuron+1, :)   = psth_temp_corr(6, :);
        ppc2_best_para_corr(neuron+1, :)  = psth_temp_corr(7, :);
        ppc2_best_n_corr(neuron+1, :)     = psth_temp_corr(8, :);
        ppc2_best_best_corr(neuron+1, :)  = psth_temp_corr(9, :);
        ppc2_null_best_corr(neuron+1, :)  = psth_temp_corr(10, :);
        ppc1_dia_best_corr(neuron+1, :)   = psth_temp_corr(11, :);
        ppc1_dia_para_corr(neuron+1, :)   = psth_temp_corr(12, :);
        ppc1_dia_n_corr(neuron+1, :)      = psth_temp_corr(13, :);
        ppc1_dia_dia_corr(neuron+1, :)    = psth_temp_corr(14, :);
        ppc1_dia_null_corr(neuron+1, :)   = psth_temp_corr(15, :);
        ppc2_dia_best_corr(neuron+1, :)   = psth_temp_corr(16, :);
        ppc2_dia_para_corr(neuron+1, :)   = psth_temp_corr(17, :);
        ppc2_dia_n_corr(neuron+1, :)      = psth_temp_corr(18, :);
        ppc2_dia_dia_corr(neuron+1, :)    = psth_temp_corr(19, :);
        ppc2_null_dia_corr(neuron+1, :)   = psth_temp_corr(20, :);
        


        neuron = neuron+1;
        clear psth_temp_err;
        clear psth_temp_corr;

    else
        disp('Empty MatData File!!!');
    end
    clear MatData_err;
    clear MatData_corr;
    
end
%%
ppc1_best_corr = [ppc1_best_dia_corr; ppc1_best_para_corr; ppc1_best_n_corr; ppc1_best_best_corr; ppc1_best_null_corr];
ppc1_best_err = [ppc1_best_dia_err; ppc1_best_para_err; ppc1_best_n_err; ppc1_best_best_err; ppc1_best_null_err];

ppc1_dia_corr = [ppc1_dia_dia_corr; ppc1_dia_para_corr; ppc1_dia_n_corr; ppc1_dia_best_corr; ppc1_dia_null_corr];
ppc1_dia_err = [ppc1_dia_dia_err; ppc1_dia_para_err; ppc1_dia_n_err; ppc1_dia_best_err; ppc1_dia_null_err];

ppc2_best_corr = [ppc2_best_dia_corr; ppc2_best_para_corr; ppc2_best_n_corr; ppc2_best_best_corr; ppc2_null_best_corr];
ppc2_best_err = [ppc2_best_dia_err; ppc2_best_para_err; ppc2_best_n_err; ppc2_best_best_err; ppc2_null_best_err];

ppc2_dia_corr = [ppc2_dia_dia_corr; ppc2_dia_para_corr; ppc2_dia_n_corr; ppc2_dia_best_corr; ppc2_null_dia_corr];
ppc2_dia_err = [ppc2_dia_dia_err; ppc2_dia_para_err; ppc2_dia_n_err; ppc2_dia_best_err; ppc2_null_dia_err];

% ppc2_best_err = [ppc2_best_dia_err];
% ppc2_dia_err  = [ppc2_dia_best_err];
%%
PPC1_best_corr   = nanmean(ppc1_best_corr);
PPC1_best_corr_sem = nanstd(ppc1_best_corr)./sqrt(sum(~isnan(ppc1_best_corr)));

PPC1_best_err   = nanmean(ppc1_best_err);
PPC1_best_err_sem = nanstd(ppc1_best_err)./sqrt(sum(~isnan(ppc1_best_err)));

PPC1_dia_corr   = nanmean(ppc1_dia_corr);
PPC1_dia_corr_sem = nanstd(ppc1_dia_corr)./sqrt(sum(~isnan(ppc1_dia_corr)));

PPC1_dia_err   = nanmean(ppc1_dia_err);
PPC1_dia_err_sem = nanstd(ppc1_dia_err)./sqrt(sum(~isnan(ppc1_dia_err)));

PPC2_best_corr   = nanmean(ppc2_best_corr);
PPC2_best_corr_sem = nanstd(ppc2_best_corr)./sqrt(sum(~isnan(ppc2_best_corr)));

PPC2_best_err   = nanmean(ppc2_best_err);
PPC2_best_err_sem = nanstd(ppc2_best_err)./sqrt(sum(~isnan(ppc2_best_err)));

PPC2_dia_corr   = nanmean(ppc2_dia_corr);
PPC2_dia_corr_sem = nanstd(ppc2_dia_corr)./sqrt(sum(~isnan(ppc2_dia_corr)));

PPC2_dia_err   = nanmean(ppc2_dia_err);
PPC2_dia_err_sem = nanstd(ppc2_dia_err)./sqrt(sum(~isnan(ppc2_dia_err)));
%%
figure
time = bin_edges;
colors2 = {[51/255 205/255 255/255], [0.678 0.922 1]};%dark blue and light blue
colors1 = {[0.7 0.7 0.7], [0.9 0.9 0.9]}; % dark gray and gray
%%
popresults1 = smooth(PPC1_best_corr, 5)';
popresults2 = smooth(PPC1_best_err, 5)';
subplot(221)

[ha_2] = shadedplot(time,popresults2+PPC1_best_err_sem,popresults2-PPC1_best_err_sem,colors1{2},colors1{2});
[ha_1] = shadedplot(time,popresults1+PPC1_best_corr_sem,popresults1-PPC1_best_corr_sem,colors1{1},colors1{1});

plot(time, popresults1, 'LineWidth', 2,'Color', [0 0 0])
hold on
plot(time, popresults2, 'LineWidth', 2,'Color', [0 0 0], 'LineStyle','--')

ha_1(2).DisplayName = 'Correct';
ha_2(2).DisplayName = 'Error';

%HIGHLIGHTING SPECIFIC SECTIONS OF THE PSTH TIMELINE
patch1 = patch([0, 0, 0.5, 0.5], [0, 45, 45 ,0], 'r'); 
set(patch1, 'FaceAlpha', 0.15) 
patch1 = patch([2, 2, 2.5, 2.5], [0, 45, 45 ,0], 'r'); 
set(patch1, 'FaceAlpha', 0.15) 
patch1 = patch([4, 4, 5, 5], [0, 45, 45 ,0], 'b'); 
set(patch1, 'FaceAlpha', 0.15) 
legend([ha_1(2), ha_2(2)])
title(sprintf('PPC Remember 1st (Best) (Nsolid=%d, Ndashed=%d)', mean(sum(~isnan(ppc1_best_corr))), mean(sum(~isnan(ppc1_best_err)))))
xlim([-1 4.5])
ylim([0 45])
xlabel('Time(s)')
ylabel('Discharge Rate(sp/s)')


%
subplot(222)

popresults1 = smooth(PPC1_dia_corr, 5)';
popresults2 = smooth(PPC1_dia_err, 5)';

[ha_2] = shadedplot(time,popresults2+PPC1_dia_err_sem,popresults2-PPC1_dia_err_sem,colors1{2},colors1{2});
[ha_1] = shadedplot(time,popresults1+PPC1_dia_corr_sem,popresults1-PPC1_dia_corr_sem,colors1{1},colors1{1});

plot(time, popresults1, 'LineWidth', 2,'Color', [0 0 0])
hold on
plot(time, popresults2, 'LineWidth', 2,'Color', [0 0 0], 'LineStyle','--')

%adding color in between the outlined sections 
patch1 = patch([0, 0, 0.5, 0.5], [0, 45, 45 ,0], 'r'); 
set(patch1, 'FaceAlpha', 0.15) 
patch1 = patch([2, 2, 2.5, 2.5], [0, 45, 45 ,0], 'r'); 
set(patch1, 'FaceAlpha', 0.15) 
patch1 = patch([4, 4, 5, 5], [0, 45, 45 ,0], 'b'); 
set(patch1, 'FaceAlpha', 0.15) 

title(sprintf('PPC Remember 1st (Diametric) (Nsolid=%d, Ndashed=%d)', mean(sum(~isnan(ppc1_dia_corr))), mean(sum(~isnan(ppc1_dia_err)))))
xlim([-1 4.5])
ylim([0 45])
xlabel('Time(s)')
ylabel('Discharge Rate(sp/s)')
ha_1(2).DisplayName = 'Correct';
ha_2(2).DisplayName = 'Error';
legend([ha_1(2), ha_2(2)])


%

subplot(223)

popresults1 = smooth(PPC2_best_corr, 5)';
popresults2 = smooth(PPC2_best_err, 5)';

[ha_2] = shadedplot(time,popresults2+PPC2_best_err_sem,popresults2-PPC2_best_err_sem,colors2{2},colors2{2});
[ha_1] = shadedplot(time,popresults1+PPC2_best_corr_sem,popresults1-PPC2_best_corr_sem,colors2{1},colors2{1});

plot(time, popresults1, 'LineWidth', 2,'Color', [0 0 0])
hold on
plot(time, popresults2, 'LineWidth', 2,'Color', [0 0 0], 'LineStyle','--')
hold on
line([0 0], [0 45],'color','k')
line([0.5 0.5], [0 45],'color','k')
%second stimulus presentation
line([2 2], [0 45],'color','k')
line([2.5 2.5], [0 45],'color','k')
%choice period 
line([4 4], [0 45],'color','k')
line([5 5], [0 45],'color','k')
%adding color in between the outlined sections 
patch1 = patch([0, 0, 0.5, 0.5], [0, 45, 45 ,0], 'r'); 
set(patch1, 'FaceAlpha', 0.15) 
patch1 = patch([2, 2, 2.5, 2.5], [0, 45, 45 ,0], 'r'); 
set(patch1, 'FaceAlpha', 0.15) 
patch1 = patch([4, 4, 5, 5], [0, 45, 45 ,0], 'b'); 
set(patch1, 'FaceAlpha', 0.15) 

title(sprintf('PPC Remember 2nd (Best) (Nsolid=%d, Ndashed=%d)', mean(sum(~isnan(ppc2_best_corr))), mean(sum(~isnan(ppc2_best_err)))))
xlim([-1 4.5])
ylim([0 45])
xlabel('Time(s)')
ylabel('Discharge Rate(sp/s)')
ha_1(2).DisplayName = 'Correct';
ha_2(2).DisplayName = 'Error';
legend([ha_1(2), ha_2(2)])

%

subplot(224)

popresults1 = smooth(PPC2_dia_corr, 5)';
popresults2 = smooth(PPC2_dia_err, 5)';

[ha_2] = shadedplot(time,popresults2+PPC2_dia_err_sem,popresults2-PPC2_dia_err_sem,colors2{2},colors2{2});
[ha_1] = shadedplot(time,popresults1+PPC2_dia_corr_sem,popresults1-PPC2_dia_corr_sem,colors2{1},colors2{1});

plot(time, popresults1, 'LineWidth', 2,'Color', [0 0 0])
hold on
plot(time, popresults2, 'LineWidth', 2,'Color', [0 0 0], 'LineStyle','--')
hold on
line([0 0], [0 45],'color','k')
line([0.5 0.5], [0 45],'color','k')
%second stimulus presentation
line([2 2], [0 45],'color','k')
line([2.5 2.5], [0 45],'color','k')
%choice period 
line([4 4], [0 45],'color','k')
line([5 5], [0 45],'color','k')
%adding color in between the outlined sections 
patch1 = patch([0, 0, 0.5, 0.5], [0, 45, 45 ,0], 'r'); 
set(patch1, 'FaceAlpha', 0.15) 
patch1 = patch([2, 2, 2.5, 2.5], [0, 45, 45 ,0], 'r'); 
set(patch1, 'FaceAlpha', 0.15) 
patch1 = patch([4, 4, 5, 5], [0, 45, 45 ,0], 'b'); 
set(patch1, 'FaceAlpha', 0.15) 

title(sprintf('PPC Remember 2nd (Diametric) (Nsolid=%d, Ndashed=%d)', mean(sum(~isnan(ppc2_dia_corr))), mean(sum(~isnan(ppc2_dia_err)))))
xlim([-1 4.5])
ylim([0 45])
xlabel('Time(s)')
ylabel('Discharge Rate(sp/s)')
ha_1(2).DisplayName = 'Correct';
ha_2(2).DisplayName = 'Error';
legend([ha_1(2), ha_2(2)])

%% to check PPC R2 best vs dia error
% ppc2_best_err = [ppc2_null_best_err];
% ppc2_dia_err  = [ppc2_null_dia_err];
% PPC2_best_err   = nanmean(ppc2_best_err);
% PPC2_best_err_sem = nanstd(ppc2_best_err)./sqrt(sum(~isnan(ppc2_best_err)));
% PPC2_dia_err   = nanmean(ppc2_dia_err);
% PPC2_dia_err_sem = nanstd(ppc2_dia_err)./sqrt(sum(~isnan(ppc2_dia_err)));
% 
% popresults1 = smooth(PPC2_best_err, 5)';
% popresults2 = smooth(PPC2_dia_err, 5)';
% 
% [ha_2] = shadedplot(time,popresults2+PPC2_dia_err_sem,popresults2-PPC2_dia_err_sem,colors2{2},colors2{2});
% [ha_1] = shadedplot(time,popresults1+PPC2_best_err_sem,popresults1-PPC2_best_err_sem,colors2{1},colors2{1});
% 
% plot(time, popresults1, 'LineWidth', 2,'Color', [0 0 0])
% hold on
% plot(time, popresults2, 'LineWidth', 2,'Color', [0 0 0], 'LineStyle','--')
% hold on
% line([0 0], [0 45],'color','k')
% line([0.5 0.5], [0 45],'color','k')
% %second stimulus presentation
% line([2 2], [0 45],'color','k')
% line([2.5 2.5], [0 45],'color','k')
% %choice period 
% line([4 4], [0 45],'color','k')
% line([5 5], [0 45],'color','k')
% %adding color in between the outlined sections 
% patch1 = patch([0, 0, 0.5, 0.5], [0, 45, 45 ,0], 'r'); 
% set(patch1, 'FaceAlpha', 0.15) 
% patch1 = patch([2, 2, 2.5, 2.5], [0, 45, 45 ,0], 'r'); 
% set(patch1, 'FaceAlpha', 0.15) 
% patch1 = patch([4, 4, 5, 5], [0, 45, 45 ,0], 'b'); 
% set(patch1, 'FaceAlpha', 0.15) 
% 
% title(sprintf('PPC Remember 2nd (Best VS Diametric) (Nsolid=%d, Ndashed=%d)', mean(sum(~isnan(ppc2_best_err))), mean(sum(~isnan(ppc2_dia_err)))))
% xlim([-1 4.5])
% ylim([0 45])
% xlabel('Time(s)')
% ylabel('Discharge Rate(sp/s)')
% ha_1(2).DisplayName = 'Best Error';
% ha_2(2).DisplayName = 'Diametric Error';
% legend([ha_1(2), ha_2(2)])
%% 
clearvars
%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%  PPC nonsigCD error and correct trials %%%%%%%%%%%%%%%%%%%

[~,~,raw] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\Distractor Project\Rana Scripts\ODRdistVar Rana\ODRdistVar','PPCnonsigCD');
neuron = 0;
for i = 1:156
    x = raw{i,1};
    fn_err = [x(1:8), '_', x(10:13), '_err.mat'];
    MatData_err =load(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\Distractor Project\Extraction\Extracted Data\ODRdistVar_Error\', fn_err]);
%     max_class_err = Neuron_Data_Max(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\Distractor Project\Extraction\Extracted Data\ODRdistVar_Error\', fn_err]);
    
    fn_corr = [x(1:8), '_', x(10:13), '.mat'];
    MatData_corr =load(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\Distractor Project\ODRdis_Var_correct_Data\', fn_corr]);
    max_class_corr = Neuron_Data_Max(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\Distractor Project\ODRdis_Var_correct_Data\', fn_corr]);
    
    if max_class_corr(1) == 1
        Classes=[1 2 3 4 5 11 12 13 14 15 6 7 8 9 10 16 17 18 19 20];
    elseif max_class_corr(1) ==6
        Classes = [6 7 8 9 10 16 17 18 19 20 1 2 3 4 5 11 12 13 14 15];
    else
        continue;
    end

    bin_width = 0.05;  % 50 milliseconds bin
    dur = [-1, 5];
    bin_edges = dur(1):bin_width:dur(2); %
    n_bins_cue = numel(histcounts([], bin_edges));

    cls=0;
    if ~isempty(MatData_err.MatData) && ~isempty(MatData_corr.MatData)
        for n= Classes
            %%%%%%%%%%%%%%%% for error %%%%%%%%%%%%%%%%%%%%%%
            allTS_err = []; allTS_corr = [];
            m_counter_err = 0; m_counter_corr = 0;

            if n<=length(MatData_err.MatData.class) && n<=length(MatData_corr.MatData.class)
                if  ~isempty(MatData_err.MatData.class(n).ntr) && ~isempty(MatData_corr.MatData.class(n).ntr)
                    if any(mean([MatData_err.MatData.class(n).ntr.cuerate])) && any(mean([MatData_corr.MatData.class(n).ntr.cuerate]))
                    %if length(MatData_err.MatData.class(n).ntr)>1 && length(MatData_corr.MatData.class(n).ntr)>1
                        %%%%%%%%%%%%%%%% for error %%%%%%%%%%%%%%%%%%%%%%
                        for tr=1:length(MatData_err.MatData.class(n).ntr)
                            if ~isempty(MatData_err.MatData.class(n).ntr(tr).Cue_onT)
                                try
                                    TS        = MatData_err.MatData.class(n).ntr(tr).TS-MatData_err.MatData.class(n).ntr(tr).Cue_onT;
                                    allTS_err     = [allTS_err TS];
                                    m_counter_err = m_counter_err + 1;
                                    clear TS;
                                catch
                                    disp('No Spike time found!!')
                                end
                            end
                        end
                        ntrs_err    = m_counter_err;
                        %%%%%%%%%%%%%%%% for correct %%%%%%%%%%%%%%%%%%%%%%
                        for tr=1:length(MatData_corr.MatData.class(n).ntr)
                            if ~isempty(MatData_corr.MatData.class(n).ntr(tr).Cue_onT)
                                try
                                    TS        = MatData_corr.MatData.class(n).ntr(tr).TS-MatData_corr.MatData.class(n).ntr(tr).Cue_onT;
                                    allTS_corr     = [allTS_corr TS];
                                    m_counter_corr = m_counter_corr + 1;
                                    clear TS;
                                catch
                                    disp('No spike data found')
                                end
                            end
                        end
                        ntrs_corr    = m_counter_corr;
                    %end
                    end
                end
            end
            
            if ~isempty(allTS_err)
                psth_temp_err(cls + 1,:) = histc(allTS_err,bin_edges)/(bin_width*ntrs_err);
            else
                psth_temp_err(cls + 1,:) = ones(1, length(bin_edges))*nan;
            end

            if ~isempty(allTS_corr)
                psth_temp_corr(cls + 1,:) = histc(allTS_corr,bin_edges)/(bin_width*ntrs_corr);
            else
                psth_temp_corr(cls + 1,:) = ones(1, length(bin_edges))*nan;
            end
            cls = cls +1;
        end
        ppc1_best_dia_err(neuron+1, :)   = psth_temp_err(1, :);
        ppc1_best_para_err(neuron+1, :)  = psth_temp_err(2, :);
        ppc1_best_n_err(neuron+1, :)     = psth_temp_err(3, :);
        ppc1_best_best_err(neuron+1, :)  = psth_temp_err(4, :);
        ppc1_best_null_err(neuron+1, :)  = psth_temp_err(5, :);
        ppc2_best_dia_err(neuron+1, :)   = psth_temp_err(6, :);
        ppc2_best_para_err(neuron+1, :)  = psth_temp_err(7, :);
        ppc2_best_n_err(neuron+1, :)     = psth_temp_err(8, :);
        ppc2_best_best_err(neuron+1, :)  = psth_temp_err(9, :);
        ppc2_null_best_err(neuron+1, :)  = psth_temp_err(10, :);
        ppc1_dia_best_err(neuron+1, :)   = psth_temp_err(11, :);
        ppc1_dia_para_err(neuron+1, :)   = psth_temp_err(12, :);
        ppc1_dia_n_err(neuron+1, :)      = psth_temp_err(13, :);
        ppc1_dia_dia_err(neuron+1, :)    = psth_temp_err(14, :);
        ppc1_dia_null_err(neuron+1, :)   = psth_temp_err(15, :);
        ppc2_dia_best_err(neuron+1, :)   = psth_temp_err(16, :);
        ppc2_dia_para_err(neuron+1, :)   = psth_temp_err(17, :);
        ppc2_dia_n_err(neuron+1, :)      = psth_temp_err(18, :);
        ppc2_dia_dia_err(neuron+1, :)    = psth_temp_err(19, :);
        ppc2_null_dia_err(neuron+1, :)   = psth_temp_err(20, :);


        ppc1_best_dia_corr(neuron+1, :)   = psth_temp_corr(1, :);
        ppc1_best_para_corr(neuron+1, :)  = psth_temp_corr(2, :);
        ppc1_best_n_corr(neuron+1, :)     = psth_temp_corr(3, :);
        ppc1_best_best_corr(neuron+1, :)  = psth_temp_corr(4, :);
        ppc1_best_null_corr(neuron+1, :)  = psth_temp_corr(5, :);
        ppc2_best_dia_corr(neuron+1, :)   = psth_temp_corr(6, :);
        ppc2_best_para_corr(neuron+1, :)  = psth_temp_corr(7, :);
        ppc2_best_n_corr(neuron+1, :)     = psth_temp_corr(8, :);
        ppc2_best_best_corr(neuron+1, :)  = psth_temp_corr(9, :);
        ppc2_null_best_corr(neuron+1, :)  = psth_temp_corr(10, :);
        ppc1_dia_best_corr(neuron+1, :)   = psth_temp_corr(11, :);
        ppc1_dia_para_corr(neuron+1, :)   = psth_temp_corr(12, :);
        ppc1_dia_n_corr(neuron+1, :)      = psth_temp_corr(13, :);
        ppc1_dia_dia_corr(neuron+1, :)    = psth_temp_corr(14, :);
        ppc1_dia_null_corr(neuron+1, :)   = psth_temp_corr(15, :);
        ppc2_dia_best_corr(neuron+1, :)   = psth_temp_corr(16, :);
        ppc2_dia_para_corr(neuron+1, :)   = psth_temp_corr(17, :);
        ppc2_dia_n_corr(neuron+1, :)      = psth_temp_corr(18, :);
        ppc2_dia_dia_corr(neuron+1, :)    = psth_temp_corr(19, :);
        ppc2_null_dia_corr(neuron+1, :)   = psth_temp_corr(20, :);
        


        neuron = neuron+1;
        clear psth_temp_err;
        clear psth_temp_corr;

    else
        disp('Empty MatData File!!!');
    end
    clear MatData_err;
    clear MatData_corr;
    
end
%%
ppc1_best_corr = [ppc1_best_dia_corr; ppc1_best_para_corr; ppc1_best_n_corr; ppc1_best_best_corr; ppc1_best_null_corr];
ppc1_best_err = [ppc1_best_dia_err; ppc1_best_para_err; ppc1_best_n_err; ppc1_best_best_err; ppc1_best_null_err];

ppc1_dia_corr = [ppc1_dia_dia_corr; ppc1_dia_para_corr; ppc1_dia_n_corr; ppc1_dia_best_corr; ppc1_dia_null_corr];
ppc1_dia_err = [ppc1_dia_dia_err; ppc1_dia_para_err; ppc1_dia_n_err; ppc1_dia_best_err; ppc1_dia_null_err];

ppc2_best_corr = [ppc2_best_dia_corr; ppc2_best_para_corr; ppc2_best_n_corr; ppc2_best_best_corr; ppc2_null_best_corr];
ppc2_best_err = [ppc2_best_dia_err; ppc2_best_para_err; ppc2_best_n_err; ppc2_best_best_err; ppc2_null_best_err];

ppc2_dia_corr = [ppc2_dia_dia_corr; ppc2_dia_para_corr; ppc2_dia_n_corr; ppc2_dia_best_corr; ppc2_null_dia_corr];
ppc2_dia_err = [ppc2_dia_dia_err; ppc2_dia_para_err; ppc2_dia_n_err; ppc2_dia_best_err; ppc2_null_dia_err];


%%
PPC1_best_corr   = nanmean(ppc1_best_corr);
PPC1_best_corr_sem = nanstd(ppc1_best_corr)./sqrt(sum(~isnan(ppc1_best_corr)));

PPC1_best_err   = nanmean(ppc1_best_err);
PPC1_best_err_sem = nanstd(ppc1_best_err)./sqrt(sum(~isnan(ppc1_best_err)));

PPC1_dia_corr   = nanmean(ppc1_dia_corr);
PPC1_dia_corr_sem = nanstd(ppc1_dia_corr)./sqrt(sum(~isnan(ppc1_dia_corr)));

PPC1_dia_err   = nanmean(ppc1_dia_err);
PPC1_dia_err_sem = nanstd(ppc1_dia_err)./sqrt(sum(~isnan(ppc1_dia_err)));

PPC2_best_corr   = nanmean(ppc2_best_corr);
PPC2_best_corr_sem = nanstd(ppc2_best_corr)./sqrt(sum(~isnan(ppc2_best_corr)));

PPC2_best_err   = nanmean(ppc2_best_err);
PPC2_best_err_sem = nanstd(ppc2_best_err)./sqrt(sum(~isnan(ppc2_best_err)));

PPC2_dia_corr   = nanmean(ppc2_dia_corr);
PPC2_dia_corr_sem = nanstd(ppc2_dia_corr)./sqrt(sum(~isnan(ppc2_dia_corr)));

PPC2_dia_err   = nanmean(ppc2_dia_err);
PPC2_dia_err_sem = nanstd(ppc2_dia_err)./sqrt(sum(~isnan(ppc2_dia_err)));
%%
figure
time = bin_edges;
colors2 = {[51/255 205/255 255/255], [0.678 0.922 1]};%dark blue and light blue
colors1 = {[0.7 0.7 0.7], [0.9 0.9 0.9]}; % dark gray and gray
popresults1 = smooth(PPC1_best_corr, 5)';
popresults2 = smooth(PPC1_best_err, 5)';
subplot(221)

[ha_2] = shadedplot(time,popresults2+PPC1_best_err_sem,popresults2-PPC1_best_err_sem,colors1{2},colors1{2});
[ha_1] = shadedplot(time,popresults1+PPC1_best_corr_sem,popresults1-PPC1_best_corr_sem,colors1{1},colors1{1});

plot(time, popresults1, 'LineWidth', 2,'Color', [0 0 0])
hold on
plot(time, popresults2, 'LineWidth', 2,'Color', [0 0 0], 'LineStyle','--')

ha_1(2).DisplayName = 'Correct';
ha_2(2).DisplayName = 'Error';

%HIGHLIGHTING SPECIFIC SECTIONS OF THE PSTH TIMELINE
patch1 = patch([0, 0, 0.5, 0.5], [0, 45, 45 ,0], 'r'); 
set(patch1, 'FaceAlpha', 0.15) 
patch1 = patch([2, 2, 2.5, 2.5], [0, 45, 45 ,0], 'r'); 
set(patch1, 'FaceAlpha', 0.15) 
patch1 = patch([4, 4, 5, 5], [0, 45, 45 ,0], 'b'); 
set(patch1, 'FaceAlpha', 0.15) 
legend([ha_1(2), ha_2(2)])
title(sprintf('PPC Remember 1st (Best) (Nsolid=%d, Ndashed=%d)', mean(sum(~isnan(ppc1_best_corr))), mean(sum(~isnan(ppc1_best_err)))))
xlim([-1 4.5])
ylim([0 45])
xlabel('Time(s)')
ylabel('Discharge Rate(sp/s)')


%
subplot(222)

popresults1 = smooth(PPC1_dia_corr, 5)';
popresults2 = smooth(PPC1_dia_err, 5)';

[ha_2] = shadedplot(time,popresults2+PPC1_dia_err_sem,popresults2-PPC1_dia_err_sem,colors1{2},colors1{2});
[ha_1] = shadedplot(time,popresults1+PPC1_dia_corr_sem,popresults1-PPC1_dia_corr_sem,colors1{1},colors1{1});

plot(time, popresults1, 'LineWidth', 2,'Color', [0 0 0])
hold on
plot(time, popresults2, 'LineWidth', 2,'Color', [0 0 0], 'LineStyle','--')

%adding color in between the outlined sections 
patch1 = patch([0, 0, 0.5, 0.5], [0, 45, 45 ,0], 'r'); 
set(patch1, 'FaceAlpha', 0.15) 
patch1 = patch([2, 2, 2.5, 2.5], [0, 45, 45 ,0], 'r'); 
set(patch1, 'FaceAlpha', 0.15) 
patch1 = patch([4, 4, 5, 5], [0, 45, 45 ,0], 'b'); 
set(patch1, 'FaceAlpha', 0.15) 

title(sprintf('PPC Remember 1st (Diametric) (Nsolid=%d, Ndashed=%d)', mean(sum(~isnan(ppc1_dia_corr))), mean(sum(~isnan(ppc1_dia_err)))))
xlim([-1 4.5])
ylim([0 45])
xlabel('Time(s)')
ylabel('Discharge Rate(sp/s)')
ha_1(2).DisplayName = 'Correct';
ha_2(2).DisplayName = 'Error';
legend([ha_1(2), ha_2(2)])


%

subplot(223)

popresults1 = smooth(PPC2_best_corr, 5)';
popresults2 = smooth(PPC2_best_err, 5)';

[ha_2] = shadedplot(time,popresults2+PPC2_best_err_sem,popresults2-PPC2_best_err_sem,colors2{2},colors2{2});
[ha_1] = shadedplot(time,popresults1+PPC2_best_corr_sem,popresults1-PPC2_best_corr_sem,colors2{1},colors2{1});

plot(time, popresults1, 'LineWidth', 2,'Color', [0 0 0])
hold on
plot(time, popresults2, 'LineWidth', 2,'Color', [0 0 0], 'LineStyle','--')
hold on
line([0 0], [0 45],'color','k')
line([0.5 0.5], [0 45],'color','k')
%second stimulus presentation
line([2 2], [0 45],'color','k')
line([2.5 2.5], [0 45],'color','k')
%choice period 
line([4 4], [0 45],'color','k')
line([5 5], [0 45],'color','k')
%adding color in between the outlined sections 
patch1 = patch([0, 0, 0.5, 0.5], [0, 45, 45 ,0], 'r'); 
set(patch1, 'FaceAlpha', 0.15) 
patch1 = patch([2, 2, 2.5, 2.5], [0, 45, 45 ,0], 'r'); 
set(patch1, 'FaceAlpha', 0.15) 
patch1 = patch([4, 4, 5, 5], [0, 45, 45 ,0], 'b'); 
set(patch1, 'FaceAlpha', 0.15) 

title(sprintf('PPC Remember 2nd (Best) (Nsolid=%d, Ndashed=%d)', mean(sum(~isnan(ppc2_best_corr))), mean(sum(~isnan(ppc2_best_err)))))
xlim([-1 4.5])
ylim([0 45])
xlabel('Time(s)')
ylabel('Discharge Rate(sp/s)')
ha_1(2).DisplayName = 'Correct';
ha_2(2).DisplayName = 'Error';
legend([ha_1(2), ha_2(2)])

%

subplot(224)

popresults1 = smooth(PPC2_dia_corr, 5)';
popresults2 = smooth(PPC2_dia_err, 5)';

[ha_2] = shadedplot(time,popresults2+PPC2_dia_err_sem,popresults2-PPC2_dia_err_sem,colors2{2},colors2{2});
[ha_1] = shadedplot(time,popresults1+PPC2_dia_corr_sem,popresults1-PPC2_dia_corr_sem,colors2{1},colors2{1});

plot(time, popresults1, 'LineWidth', 2,'Color', [0 0 0])
hold on
plot(time, popresults2, 'LineWidth', 2,'Color', [0 0 0], 'LineStyle','--')
hold on
line([0 0], [0 45],'color','k')
line([0.5 0.5], [0 45],'color','k')
%second stimulus presentation
line([2 2], [0 45],'color','k')
line([2.5 2.5], [0 45],'color','k')
%choice period 
line([4 4], [0 45],'color','k')
line([5 5], [0 45],'color','k')
%adding color in between the outlined sections 
patch1 = patch([0, 0, 0.5, 0.5], [0, 45, 45 ,0], 'r'); 
set(patch1, 'FaceAlpha', 0.15) 
patch1 = patch([2, 2, 2.5, 2.5], [0, 45, 45 ,0], 'r'); 
set(patch1, 'FaceAlpha', 0.15) 
patch1 = patch([4, 4, 5, 5], [0, 45, 45 ,0], 'b'); 
set(patch1, 'FaceAlpha', 0.15) 

title(sprintf('PPC Remember 2nd (Diametric) (Nsolid=%d, Ndashed=%d)', mean(sum(~isnan(ppc2_dia_corr))), mean(sum(~isnan(ppc2_dia_err)))))
xlim([-1 4.5])
ylim([0 45])
xlabel('Time(s)')
ylabel('Discharge Rate(sp/s)')
ha_1(2).DisplayName = 'Correct';
ha_2(2).DisplayName = 'Error';
legend([ha_1(2), ha_2(2)])










