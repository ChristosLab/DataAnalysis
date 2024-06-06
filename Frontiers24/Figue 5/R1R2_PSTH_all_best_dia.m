clc; clear all; close all;

NS_BS='NS';
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%  PFC all neurons %%%%%%%%%%%%%%%%%%%

[~,~,raw1] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\NS_BS_Database_20231109a.xlsx',['R1R2_PFC_' ,NS_BS]);
% [~,~,raw1] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\AllSigNeuronsInfo.xlsx',['R1R2_PFC_' ,NS_BS]);
if strcmp(raw1{1,1}, 'Filename')
    raw1=raw1(2:end, :);
end
%%
neuron = 0;
for i = 1:length(raw1)
    i
    fn_corr = [raw1{i,1}, '.mat'];  %for new format 
    MatData_corr =load(['ALLDataCorrErr\', fn_corr]);
    max_class_corr = Neuron_Data_Max(fn_corr);
    max_class_corr_cuerate = max_class_corr(2);
    max_class_corr_fixrate = max_class_corr(3);
    max_class_corr_cdrate = max_class_corr(4);
    if ~isempty(max_class_corr)
        if fn_corr(8)=='2'
            if max_class_corr(1) == 1
                Classes=[1 2 3 4 5 6 7 8 9 10];
            elseif max_class_corr(1) ==6
                Classes = [6 7 8 9 10 1 2 3 4 5];
            else
                continue;
            end
        elseif fn_corr(8)=='1'
            if max_class_corr(1) <=4
                Classes=[max_class_corr(1) max_class_corr(1)+4];
            elseif max_class_corr(1) >4 && max_class_corr(1)<=8
                Classes=[max_class_corr(1) max_class_corr(1)-4 ];
            else
                disp(fn_corr)
                continue;
            end 
        end
%         disp(max_class_corr(1))
    else
        disp(fn_corr)
    end

    bin_width = 0.05;  % 50 milliseconds bin
    dur = [-1-bin_width, 5-bin_width];
    bin_edges = dur(1):bin_width:dur(2); %
    n_bins_cue = numel(histcounts([], bin_edges));

    cls=0;
    if ~isempty(MatData_corr.MatData)
        for n= Classes
            %%%%%%%%%%%%%%%% for error %%%%%%%%%%%%%%%%%%%%%%
            allTS_corr = [];
            m_counter_corr = 0;

            if n<=length(MatData_corr.MatData.class) && ~isempty(MatData_corr.MatData.class(n).ntr)
                if any(mean([MatData_corr.MatData.class(n).ntr.cuerate]))
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
                end
            else
                disp("Class does not exist")
            end

            if ~isempty(allTS_corr)
                psth_temp_corr(cls + 1,:) = histc(allTS_corr,bin_edges)/(bin_width*ntrs_corr);
            else
                psth_temp_corr(cls + 1,:) = ones(1, length(bin_edges))*nan;
            end
            cls = cls +1;
        end
        if fn_corr(8)=='2'
            pfc1_best_dia_corr(neuron+1, :)   = psth_temp_corr(1, :);
            pfc1_best_para_corr(neuron+1, :)  = psth_temp_corr(2, :);
            pfc1_best_n_corr(neuron+1, :)     = psth_temp_corr(3, :);
            pfc1_best_best_corr(neuron+1, :)  = psth_temp_corr(4, :);
            pfc1_best_null_corr(neuron+1, :)  = psth_temp_corr(5, :);
            pfc1_dia_best_corr(neuron+1, :)   = psth_temp_corr(6, :);
            pfc1_dia_para_corr(neuron+1, :)   = psth_temp_corr(7, :);
            pfc1_dia_n_corr(neuron+1, :)      = psth_temp_corr(8, :);
            pfc1_dia_dia_corr(neuron+1, :)    = psth_temp_corr(9, :);
            pfc1_dia_null_corr(neuron+1, :)   = psth_temp_corr(10, :);
            pfc1_best_cuerate(neuron+1,:)     = max_class_corr_cuerate;
            pfc1_best_fixrate(neuron+1,:)     = max_class_corr_fixrate;
            pfc1_best_cdrate(neuron+1,:)      = max_class_corr_cdrate;
        elseif fn_corr(8)=='1'
            pfc1_best_dia_corr(neuron+1, :)   = psth_temp_corr(1, :);
            pfc1_dia_best_corr(neuron+1, :)   = psth_temp_corr(2, :);
            pfc1_best_cuerate(neuron+1,:)     = max_class_corr_cuerate;
            pfc1_best_fixrate(neuron+1,:)     = max_class_corr_fixrate;
            pfc1_best_cdrate(neuron+1,:)     = max_class_corr_cdrate;
        end

        neuron = neuron+1;
%         disp(neuron)
        clear psth_temp_corr;

    else
        disp('Empty MatData File!!!');
    end
    clear MatData_corr;
    
end
%%
pfc1_best_corr = [pfc1_best_dia_corr; pfc1_best_para_corr; pfc1_best_n_corr; pfc1_best_best_corr; pfc1_best_null_corr];
pfc1_dia_corr = [pfc1_dia_best_corr; pfc1_dia_para_corr; pfc1_dia_n_corr; pfc1_dia_dia_corr; pfc1_dia_null_corr];


%%
pfc1_best_corr     = pfc1_best_corr(~all(pfc1_best_corr == 0, 2), :);
PFC1_best_corr     = nanmean(pfc1_best_corr);
PFC1_best_corr_sem = nanstd(pfc1_best_corr)./sqrt(mean(sum(~isnan(pfc1_best_corr))));
PFC1_best_n        = mean(sum(~isnan(pfc1_best_corr)));

R1R2_PFC_NS_allbestcue = nanmean(pfc1_best_corr(:, 23:32), 2);
R1R2_PFC_NS_allbestcue = R1R2_PFC_NS_allbestcue(find(~isnan(R1R2_PFC_NS_allbestcue)));
R1R2_PFC_NS_allbestcd = nanmean(pfc1_best_corr(:, 33:62), 2);
R1R2_PFC_NS_allbestcd = R1R2_PFC_NS_allbestcd(find(~isnan(R1R2_PFC_NS_allbestcd)));

pfc1_dia_corr     = pfc1_dia_corr(~all(pfc1_dia_corr == 0, 2), :);
PFC1_dia_corr     = nanmean(pfc1_dia_corr);
PFC1_dia_corr_sem = nanstd(pfc1_dia_corr)./sqrt(mean(sum(~isnan(pfc1_dia_corr))));
PFC1_dia_n        = mean(sum(~isnan(pfc1_dia_corr)));

R1R2_PFC_NS_cuerates = pfc1_best_cuerate(find(~isnan(pfc1_best_cuerate)));
R1R2_PFC_NS_fixrates = pfc1_best_fixrate(find(~isnan(pfc1_best_fixrate)));
R1R2_PFC_NS_cdrates  = pfc1_best_cdrate(find(~isnan(pfc1_best_cdrate)));

clearvars -except PFC1_best_corr PFC1_best_corr_sem PFC1_dia_corr PFC1_dia_corr_sem NS_BS raw1 PFC1_best_n ...
    PFC1_dia_n R1R2_PFC_NS_cuerates R1R2_PFC_NS_fixrates R1R2_PFC_NS_cdrates R1R2_PFC_NS_allbestcue R1R2_PFC_NS_allbestcd

%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%  PPC all neurons %%%%%%%%%%%%%%%%%%%

% [~,~,raw2] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\NS_BS_Database_20231109a.xlsx',['R1R2_PPC_' ,NS_BS]);
[~,~,raw2] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\AllSigNeuronsInfo.xlsx',['R1R2_PPC_' ,NS_BS]);
if strcmp(raw2{1,1}, 'Filename')
    raw2=raw2(2:end, :);
end
%%
neuron = 0;
for i = 1:length(raw2)
    i
    fn_corr = [raw2{i,1}, '.mat'];  %for new format 
    MatData_corr =load(['ALLDataCorrErr\', fn_corr]);
    max_class_corr = Neuron_Data_Max(fn_corr);
    max_class_corr_cuerate = max_class_corr(2);
    max_class_corr_fixrate = max_class_corr(3);
    max_class_corr_cdrate = max_class_corr(4);
    if ~isempty(max_class_corr)
        if fn_corr(8)=='2'
            if max_class_corr(1) == 1
                Classes=[1 2 3 4 5 6 7 8 9 10];
            elseif max_class_corr(1) ==6
                Classes = [6 7 8 9 10 1 2 3 4 5];
            else
                continue;
            end
        elseif fn_corr(8)=='1'
            if max_class_corr(1) <=4
                Classes=[max_class_corr(1) max_class_corr(1)+4];
            elseif max_class_corr(1) >4 && max_class_corr(1)<=8
                Classes=[max_class_corr(1) max_class_corr(1)-4 ];
            else
                disp(fn_corr)
                continue;
            end 
        end
%         disp(max_class_corr(1))
    else
        disp(fn_corr)
    end

    bin_width = 0.05;  % 50 milliseconds bin
    dur = [-1-bin_width, 5-bin_width];
    bin_edges = dur(1):bin_width:dur(2); %
    n_bins_cue = numel(histcounts([], bin_edges));

    cls=0;
    if ~isempty(MatData_corr.MatData)
        for n= Classes
            %%%%%%%%%%%%%%%% for error %%%%%%%%%%%%%%%%%%%%%%
            allTS_corr = [];
            m_counter_corr = 0;

            if n<=length(MatData_corr.MatData.class)
                if  ~isempty(MatData_corr.MatData.class(n).ntr)
                    if any(mean([MatData_corr.MatData.class(n).ntr.cuerate]))
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
                    end
                end
            end

            if ~isempty(allTS_corr)
                psth_temp_corr(cls + 1,:) = histc(allTS_corr,bin_edges)/(bin_width*ntrs_corr);
            else
                psth_temp_corr(cls + 1,:) = ones(1, length(bin_edges))*nan;
            end
            cls = cls +1;
        end
        if fn_corr(8)=='2'
            ppc1_best_dia_corr(neuron+1, :)   = psth_temp_corr(1, :);
            ppc1_best_para_corr(neuron+1, :)  = psth_temp_corr(2, :);
            ppc1_best_n_corr(neuron+1, :)     = psth_temp_corr(3, :);
            ppc1_best_best_corr(neuron+1, :)  = psth_temp_corr(4, :);
            ppc1_best_null_corr(neuron+1, :)  = psth_temp_corr(5, :);
            ppc1_dia_best_corr(neuron+1, :)   = psth_temp_corr(6, :);
            ppc1_dia_para_corr(neuron+1, :)   = psth_temp_corr(7, :);
            ppc1_dia_n_corr(neuron+1, :)      = psth_temp_corr(8, :);
            ppc1_dia_dia_corr(neuron+1, :)    = psth_temp_corr(9, :);
            ppc1_dia_null_corr(neuron+1, :)   = psth_temp_corr(10, :);
            ppc1_best_cuerate(neuron+1,:)     = max_class_corr_cuerate;
            ppc1_best_fixrate(neuron+1,:)     = max_class_corr_fixrate;
            ppc1_best_cdrate(neuron+1,:)     = max_class_corr_cdrate;
        elseif fn_corr(8)=='1'
            ppc1_best_dia_corr(neuron+1, :)   = psth_temp_corr(1, :);
            ppc1_dia_best_corr(neuron+1, :)   = psth_temp_corr(2, :);
            ppc1_best_cuerate(neuron+1,:)     = max_class_corr_cuerate;
            ppc1_best_fixrate(neuron+1,:)     = max_class_corr_fixrate;
            ppc1_best_cdrate(neuron+1,:)     = max_class_corr_cdrate;
        end

        neuron = neuron+1;
%         disp(neuron)
        clear psth_temp_corr;

    else
        disp('Empty MatData File!!!');
    end
    clear MatData_corr;
    
end
%%
ppc1_best_corr = [ppc1_best_dia_corr; ppc1_best_para_corr; ppc1_best_n_corr; ppc1_best_best_corr; ppc1_best_null_corr];
ppc1_dia_corr = [ppc1_dia_best_corr; ppc1_dia_para_corr; ppc1_dia_n_corr; ppc1_dia_dia_corr; ppc1_dia_null_corr];

%%
ppc1_best_corr     = ppc1_best_corr(~all(ppc1_best_corr == 0, 2), :);
PPC1_best_corr     = nanmean(ppc1_best_corr);
PPC1_best_corr_sem = nanstd(ppc1_best_corr)./sqrt(mean(sum(~isnan(ppc1_best_corr))));
PPC1_best_n        = mean(sum(~isnan(ppc1_best_corr)));

R1R2_PPC_NS_allbestcue = nanmean(ppc1_best_corr(:, 23:32), 2);
R1R2_PPC_NS_allbestcue = R1R2_PPC_NS_allbestcue(find(~isnan(R1R2_PPC_NS_allbestcue)));
R1R2_PPC_NS_allbestcd = nanmean(ppc1_best_corr(:, 33:62), 2);
R1R2_PPC_NS_allbestcd = R1R2_PPC_NS_allbestcd(find(~isnan(R1R2_PPC_NS_allbestcd)));

ppc1_dia_corr     = ppc1_dia_corr(~all(ppc1_dia_corr == 0, 2), :);
PPC1_dia_corr     = nanmean(ppc1_dia_corr);
PPC1_dia_corr_sem = nanstd(ppc1_dia_corr)./sqrt(mean(sum(~isnan(ppc1_dia_corr))));
PPC1_dia_n        = mean(sum(~isnan(ppc1_dia_corr)));

R1R2_PPC_NS_cuerates = ppc1_best_cuerate(find(~isnan(ppc1_best_cuerate)));
R1R2_PPC_NS_fixrates = ppc1_best_fixrate(find(~isnan(ppc1_best_fixrate)));
R1R2_PPC_NS_cdrates  = ppc1_best_cdrate(find(~isnan(ppc1_best_cdrate)));
%%
figure
time = bin_edges+bin_width;
if strcmpi(NS_BS, 'BS')
    color='b';
    colors = {[0.3010 0.7450 0.9330], [0 1 1]};
elseif strcmpi(NS_BS, 'NS')
    color='r';
    colors = {[1, 0.6, 0.6], [1, 0.8, 0.8]};
end

popresults1 = smooth(PFC1_best_corr, 1)';
popresults2 = smooth(PFC1_dia_corr, 1)';
subplot(221)

[ha_1] = shadedplot(time,popresults1+PFC1_best_corr_sem,popresults1-PFC1_best_corr_sem,colors{1},colors{1});
[ha_2] = shadedplot(time,popresults2+PFC1_dia_corr_sem,popresults2-PFC1_dia_corr_sem,colors{2},colors{2});

plot(time, popresults2, 'LineWidth', 2,'Color', color, 'LineStyle','--')
hold on
plot(time, popresults1, 'LineWidth', 2,'Color', color)
hold on

ha_1(2).DisplayName = 'Best';
ha_2(2).DisplayName = 'Diametric';

%HIGHLIGHTING SPECIFIC SECTIONS OF THE PSTH TIMELINE
line([0 0], [0  100],'Color','k','LineStyle','--');
line([0.5 0.5], [0  100],'Color','k','LineStyle','--');
line([2 2], [0  100],'Color','k','LineStyle','--');
line([2.5 2.5], [0  100],'Color','k','LineStyle','--');
line([4 4], [0  100],'Color','k','LineStyle','--');

legend([ha_1(2), ha_2(2)])

% title(sprintf('PFC Remember 1st (NS=%d, nSolid=%d, nDashed=%d)',length(raw1), PFC1_best_n, PFC1_dia_n))
title(sprintf('PFC Remember 1st (NS=%d)',length(raw1)))
xlim([-1 4.5])
ylim([0 40])
xlabel('Time (s)', 'FontSize', 15)
ylabel('Discharge Rate (sp/s)', 'FontSize', 15)
% grid on
% box off

%
subplot(223)
popresults1 = smooth(PPC1_best_corr, 1)';
popresults2 = smooth(PPC1_dia_corr, 1)';

[ha_1] = shadedplot(time,popresults1+PPC1_best_corr_sem,popresults1-PPC1_best_corr_sem,colors{1},colors{1});
[ha_2] = shadedplot(time,popresults2+PPC1_dia_corr_sem,popresults2-PPC1_dia_corr_sem,colors{2},colors{2});

plot(time, popresults2, 'LineWidth', 2,'Color', color, 'LineStyle','--')
hold on
plot(time, popresults1, 'LineWidth', 2,'Color', color)
hold on

ha_1(2).DisplayName = 'Best';
ha_2(2).DisplayName = 'Diametric';

%HIGHLIGHTING SPECIFIC SECTIONS OF THE PSTH TIMELINE
line([0 0], [0  100],'Color','k','LineStyle','--');
line([0.5 0.5], [0  100],'Color','k','LineStyle','--');
line([2 2], [0  100],'Color','k','LineStyle','--');
line([2.5 2.5], [0  100],'Color','k','LineStyle','--');
line([4 4], [0  100],'Color','k','LineStyle','--');

legend([ha_1(2), ha_2(2)])

% title(sprintf('PPC Remember 1st (NS=%d, nSolid=%d, nDashed=%d)',length(raw2), PPC1_best_n, PPC1_dia_n))
title(sprintf('PPC Remember 1st (NS=%d)',length(raw2)))
xlim([-1 4.5])
ylim([0 40])
xlabel('Time (s)', 'FontSize', 15)
ylabel('Discharge Rate (sp/s)', 'FontSize', 15)
% grid on
% box off

%%


clc; clearvars -except R1R2_PFC_NS_cuerates R1R2_PPC_NS_cuerates R1R2_PPC_NS_fixrates R1R2_PFC_NS_fixrates ...
    R1R2_PFC_NS_cdrates R1R2_PPC_NS_cdrates R1R2_PFC_NS_allbestcue R1R2_PFC_NS_allbestcd R1R2_PPC_NS_allbestcue R1R2_PPC_NS_allbestcd

NS_BS='BS';
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%  PFC all neurons %%%%%%%%%%%%%%%%%%%

% [~,~,raw1] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\NS_BS_Database_20231109a.xlsx',['R1R2_PFC_' ,NS_BS]);
[~,~,raw1] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\AllSigNeuronsInfo.xlsx',['R1R2_PFC_' ,NS_BS]);
if strcmp(raw1{1,1}, 'Filename')
    raw1=raw1(2:end, :);
end
%%
neuron = 0;
for i = 1:length(raw1)
    i
    fn_corr = [raw1{i,1}, '.mat'];  %for new format 
    MatData_corr =load(['ALLDataCorrErr\', fn_corr]);
    max_class_corr = Neuron_Data_Max(fn_corr);
    max_class_corr_cuerate = max_class_corr(2);
    max_class_corr_fixrate = max_class_corr(3);
    max_class_corr_cdrate = max_class_corr(4);
    if ~isempty(max_class_corr)
        if fn_corr(8)=='2'
            if max_class_corr(1) == 1
                Classes=[1 2 3 4 5 6 7 8 9 10];
            elseif max_class_corr(1) ==6
                Classes = [6 7 8 9 10 1 2 3 4 5];
            else
                continue;
            end
        elseif fn_corr(8)=='1'
            if max_class_corr(1) <=4
                Classes=[max_class_corr(1) max_class_corr(1)+4];
            elseif max_class_corr(1) >4 && max_class_corr(1)<=8
                Classes=[max_class_corr(1) max_class_corr(1)-4 ];
            else
                disp(fn_corr)
                continue;
            end 
        end
%         disp(max_class_corr(1))
    else
        disp(fn_corr)
    end

    bin_width = 0.05;  % 50 milliseconds bin
    dur = [-1-bin_width, 5-bin_width];
    bin_edges = dur(1):bin_width:dur(2); %
    n_bins_cue = numel(histcounts([], bin_edges));

    cls=0;
    if ~isempty(MatData_corr.MatData)
        for n= Classes
            %%%%%%%%%%%%%%%% for error %%%%%%%%%%%%%%%%%%%%%%
            allTS_corr = [];
            m_counter_corr = 0;

            if n<=length(MatData_corr.MatData.class)
                if  ~isempty(MatData_corr.MatData.class(n).ntr)
                    if any(mean([MatData_corr.MatData.class(n).ntr.cuerate]))
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
                    end
                end
            end

            if ~isempty(allTS_corr)
                psth_temp_corr(cls + 1,:) = histc(allTS_corr,bin_edges)/(bin_width*ntrs_corr);
            else
                psth_temp_corr(cls + 1,:) = ones(1, length(bin_edges))*nan;
            end
            cls = cls +1;
        end
        if fn_corr(8)=='2'
            pfc1_best_dia_corr(neuron+1, :)   = psth_temp_corr(1, :);
            pfc1_best_para_corr(neuron+1, :)  = psth_temp_corr(2, :);
            pfc1_best_n_corr(neuron+1, :)     = psth_temp_corr(3, :);
            pfc1_best_best_corr(neuron+1, :)  = psth_temp_corr(4, :);
            pfc1_best_null_corr(neuron+1, :)  = psth_temp_corr(5, :);
            pfc1_dia_best_corr(neuron+1, :)   = psth_temp_corr(6, :);
            pfc1_dia_para_corr(neuron+1, :)   = psth_temp_corr(7, :);
            pfc1_dia_n_corr(neuron+1, :)      = psth_temp_corr(8, :);
            pfc1_dia_dia_corr(neuron+1, :)    = psth_temp_corr(9, :);
            pfc1_dia_null_corr(neuron+1, :)   = psth_temp_corr(10, :);
            pfc1_best_cuerate(neuron+1,:)     = max_class_corr_cuerate;
            pfc1_best_fixrate(neuron+1,:)     = max_class_corr_fixrate;
            pfc1_best_cdrate(neuron+1,:)      = max_class_corr_cdrate;
        elseif fn_corr(8)=='1'
            pfc1_best_dia_corr(neuron+1, :)   = psth_temp_corr(1, :);
            pfc1_dia_best_corr(neuron+1, :)   = psth_temp_corr(2, :);
            pfc1_best_cuerate(neuron+1,:)     = max_class_corr_cuerate;
            pfc1_best_fixrate(neuron+1,:)     = max_class_corr_fixrate;
            pfc1_best_cdrate(neuron+1,:)      = max_class_corr_cdrate;
        end

        neuron = neuron+1;
%         disp(neuron)
        clear psth_temp_corr;

    else
        disp('Empty MatData File!!!');
    end
    clear MatData_corr;
    
end
%%
pfc1_best_corr = [pfc1_best_dia_corr; pfc1_best_para_corr; pfc1_best_n_corr; pfc1_best_best_corr; pfc1_best_null_corr];
pfc1_dia_corr = [pfc1_dia_best_corr; pfc1_dia_para_corr; pfc1_dia_n_corr; pfc1_dia_dia_corr; pfc1_dia_null_corr];


%%
pfc1_best_corr     = pfc1_best_corr(~all(pfc1_best_corr == 0, 2), :);
PFC1_best_corr     = nanmean(pfc1_best_corr);
PFC1_best_corr_sem = nanstd(pfc1_best_corr)./sqrt(mean(sum(~isnan(pfc1_best_corr))));
PFC1_best_n        = mean(sum(~isnan(pfc1_best_corr)));

R1R2_PFC_BS_allbestcue = nanmean(pfc1_best_corr(:, 23:32), 2);
R1R2_PFC_BS_allbestcue = R1R2_PFC_BS_allbestcue(find(~isnan(R1R2_PFC_BS_allbestcue)));
R1R2_PFC_BS_allbestcd = nanmean(pfc1_best_corr(:, 33:62), 2);
R1R2_PFC_BS_allbestcd = R1R2_PFC_BS_allbestcd(find(~isnan(R1R2_PFC_BS_allbestcd)));

pfc1_dia_corr     = pfc1_dia_corr(~all(pfc1_dia_corr == 0, 2), :);
PFC1_dia_corr     = nanmean(pfc1_dia_corr);
PFC1_dia_corr_sem = nanstd(pfc1_dia_corr)./sqrt(mean(sum(~isnan(pfc1_dia_corr))));
PFC1_dia_n        = mean(sum(~isnan(pfc1_dia_corr)));

R1R2_PFC_BS_cuerates = pfc1_best_cuerate(find(~isnan(pfc1_best_cuerate)));
R1R2_PFC_BS_fixrates = pfc1_best_fixrate(find(~isnan(pfc1_best_fixrate)));
R1R2_PFC_BS_cdrates  = pfc1_best_cdrate(find(~isnan(pfc1_best_cdrate)));

clearvars -except PFC1_best_corr PFC1_best_corr_sem PFC1_dia_corr PFC1_dia_corr_sem NS_BS raw1 PFC1_best_n PFC1_dia_n ...
    R1R2_PFC_BS_cuerates R1R2_PFC_NS_cuerates R1R2_PPC_NS_cuerates R1R2_PFC_BS_fixrates R1R2_PFC_NS_fixrates ...
    R1R2_PPC_NS_fixrates R1R2_PFC_NS_cdrates R1R2_PPC_NS_cdrates R1R2_PFC_BS_cdrates R1R2_PFC_NS_allbestcue R1R2_PFC_NS_allbestcd  ...
    R1R2_PPC_NS_allbestcue R1R2_PPC_NS_allbestcd R1R2_PFC_BS_allbestcue R1R2_PFC_BS_allbestcd

%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%  PPC all neurons %%%%%%%%%%%%%%%%%%%

% [~,~,raw2] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\NS_BS_Database_20231109a.xlsx',['R1R2_PPC_' ,NS_BS]);
[~,~,raw2] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\AllSigNeuronsInfo.xlsx',['R1R2_PPC_' ,NS_BS]);
if strcmp(raw2{1,1}, 'Filename')
    raw2=raw2(2:end, :);
end
%%
neuron = 0;
for i = 1:length(raw2)
    i
    fn_corr = [raw2{i,1}, '.mat'];  %for new format 
    MatData_corr =load(['ALLDataCorrErr\', fn_corr]);
    max_class_corr = Neuron_Data_Max(fn_corr);
    max_class_corr_cuerate = max_class_corr(2);
    max_class_corr_fixrate = max_class_corr(3);
    max_class_corr_cdrate = max_class_corr(4);
    if ~isempty(max_class_corr)
        if fn_corr(8)=='2'
            if max_class_corr(1) == 1
                Classes=[1 2 3 4 5 6 7 8 9 10];
            elseif max_class_corr(1) ==6
                Classes = [6 7 8 9 10 1 2 3 4 5];
            else
                continue;
            end
        elseif fn_corr(8)=='1'
            if max_class_corr(1) <=4
                Classes=[max_class_corr(1) max_class_corr(1)+4];
            elseif max_class_corr(1) >4 && max_class_corr(1)<=8
                Classes=[max_class_corr(1) max_class_corr(1)-4 ];
            else
                disp(fn_corr)
                continue;
            end 
        end
%         disp(max_class_corr(1))
    else
        disp(fn_corr)
    end

    bin_width = 0.05;  % 50 milliseconds bin
    dur = [-1-bin_width, 5-bin_width];
    bin_edges = dur(1):bin_width:dur(2); %
    n_bins_cue = numel(histcounts([], bin_edges));

    cls=0;
    if ~isempty(MatData_corr.MatData)
        for n= Classes
            %%%%%%%%%%%%%%%% for error %%%%%%%%%%%%%%%%%%%%%%
            allTS_corr = [];
            m_counter_corr = 0;

            if n<=length(MatData_corr.MatData.class)
                if  ~isempty(MatData_corr.MatData.class(n).ntr)
                    if any(mean([MatData_corr.MatData.class(n).ntr.cuerate]))
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
                    end
                end
            end

            if ~isempty(allTS_corr)
                psth_temp_corr(cls + 1,:) = histc(allTS_corr,bin_edges)/(bin_width*ntrs_corr);
            else
                psth_temp_corr(cls + 1,:) = ones(1, length(bin_edges))*nan;
            end
            cls = cls +1;
        end
        if fn_corr(8)=='2'
            ppc1_best_dia_corr(neuron+1, :)   = psth_temp_corr(1, :);
            ppc1_best_para_corr(neuron+1, :)  = psth_temp_corr(2, :);
            ppc1_best_n_corr(neuron+1, :)     = psth_temp_corr(3, :);
            ppc1_best_best_corr(neuron+1, :)  = psth_temp_corr(4, :);
            ppc1_best_null_corr(neuron+1, :)  = psth_temp_corr(5, :);
            ppc1_dia_best_corr(neuron+1, :)   = psth_temp_corr(6, :);
            ppc1_dia_para_corr(neuron+1, :)   = psth_temp_corr(7, :);
            ppc1_dia_n_corr(neuron+1, :)      = psth_temp_corr(8, :);
            ppc1_dia_dia_corr(neuron+1, :)    = psth_temp_corr(9, :);
            ppc1_dia_null_corr(neuron+1, :)   = psth_temp_corr(10, :);
            ppc1_best_cuerate(neuron+1,:)     = max_class_corr_cuerate;
            ppc1_best_fixrate(neuron+1,:)     = max_class_corr_fixrate;
            ppc1_best_cdrate(neuron+1,:)     = max_class_corr_cdrate;
        elseif fn_corr(8)=='1'
            ppc1_best_dia_corr(neuron+1, :)   = psth_temp_corr(1, :);
            ppc1_dia_best_corr(neuron+1, :)   = psth_temp_corr(2, :);
            ppc1_best_cuerate(neuron+1,:)     = max_class_corr_cuerate;
            ppc1_best_fixrate(neuron+1,:)     = max_class_corr_fixrate;
            ppc1_best_cdrate(neuron+1,:)     = max_class_corr_cdrate;
        end

        neuron = neuron+1;
%         disp(neuron)
        clear psth_temp_corr;

    else
        disp('Empty MatData File!!!');
    end
    clear MatData_corr;
    
end
%%
ppc1_best_corr = [ppc1_best_dia_corr; ppc1_best_para_corr; ppc1_best_n_corr; ppc1_best_best_corr; ppc1_best_null_corr];
ppc1_dia_corr = [ppc1_dia_best_corr; ppc1_dia_para_corr; ppc1_dia_n_corr; ppc1_dia_dia_corr; ppc1_dia_null_corr];

%%
ppc1_best_corr     = ppc1_best_corr(~all(ppc1_best_corr == 0, 2), :);
PPC1_best_corr     = nanmean(ppc1_best_corr);
PPC1_best_corr_sem = nanstd(ppc1_best_corr)./sqrt(mean(sum(~isnan(ppc1_best_corr))));
PPC1_best_n        = mean(sum(~isnan(ppc1_best_corr)));

R1R2_PPC_BS_allbestcue = nanmean(ppc1_best_corr(:, 23:32), 2);
R1R2_PPC_BS_allbestcue = R1R2_PPC_BS_allbestcue(find(~isnan(R1R2_PPC_BS_allbestcue)));
R1R2_PPC_BS_allbestcd = nanmean(ppc1_best_corr(:, 33:62), 2);
R1R2_PPC_BS_allbestcd = R1R2_PPC_BS_allbestcd(find(~isnan(R1R2_PPC_BS_allbestcd)));

ppc1_dia_corr     = ppc1_dia_corr(~all(ppc1_dia_corr == 0, 2), :);
PPC1_dia_corr     = nanmean(ppc1_dia_corr);
PPC1_dia_corr_sem = nanstd(ppc1_dia_corr)./sqrt(mean(sum(~isnan(ppc1_dia_corr))));
PPC1_dia_n        = mean(sum(~isnan(ppc1_dia_corr)));

R1R2_PPC_BS_cuerates = ppc1_best_cuerate(find(~isnan(ppc1_best_cuerate)));
R1R2_PPC_BS_fixrates = ppc1_best_fixrate(find(~isnan(ppc1_best_fixrate)));
R1R2_PPC_BS_cdrates  = ppc1_best_cdrate(find(~isnan(ppc1_best_cdrate)));

%%
time = bin_edges+bin_width;
if strcmpi(NS_BS, 'BS')
    color='b';
    colors = {[0.3010 0.7450 0.9330], [0 1 1]};
elseif strcmpi(NS_BS, 'NS')
    color='r';
    colors = {[1, 0.6, 0.6], [1, 0.8, 0.8]};
end

popresults1 = smooth(PFC1_best_corr, 1)';
popresults2 = smooth(PFC1_dia_corr, 1)';
subplot(222)

[ha_1] = shadedplot(time,popresults1+PFC1_best_corr_sem,popresults1-PFC1_best_corr_sem,colors{1},colors{1});
[ha_2] = shadedplot(time,popresults2+PFC1_dia_corr_sem,popresults2-PFC1_dia_corr_sem,colors{2},colors{2});

plot(time, popresults2, 'LineWidth', 2,'Color', color, 'LineStyle','--')
hold on
plot(time, popresults1, 'LineWidth', 2,'Color', color)
hold on

ha_1(2).DisplayName = 'Best';
ha_2(2).DisplayName = 'Diametric';

%HIGHLIGHTING SPECIFIC SECTIONS OF THE PSTH TIMELINE
line([0 0], [0  100],'Color','k','LineStyle','--');
line([0.5 0.5], [0  100],'Color','k','LineStyle','--');
line([2 2], [0  100],'Color','k','LineStyle','--');
line([2.5 2.5], [0  100],'Color','k','LineStyle','--');
line([4 4], [0  100],'Color','k','LineStyle','--');

legend([ha_1(2), ha_2(2)])

% title(sprintf('PFC Remember 1st (BS=%d, nSolid=%d, nDashed=%d)',length(raw1), PFC1_best_n, PFC1_dia_n))
title(sprintf('PFC Remember 1st (BS=%d)',length(raw1)))
xlim([-1 4.5])
ylim([0 40])
xlabel('Time (s)', 'FontSize', 15)
ylabel('Discharge Rate (sp/s)', 'FontSize', 15)
% grid on
% box off

%
subplot(224)
popresults1 = smooth(PPC1_best_corr, 1)';
popresults2 = smooth(PPC1_dia_corr, 1)';

[ha_1] = shadedplot(time,popresults1+PPC1_best_corr_sem,popresults1-PPC1_best_corr_sem,colors{1},colors{1});
[ha_2] = shadedplot(time,popresults2+PPC1_dia_corr_sem,popresults2-PPC1_dia_corr_sem,colors{2},colors{2});

plot(time, popresults2, 'LineWidth', 2,'Color', color, 'LineStyle','--')
hold on
plot(time, popresults1, 'LineWidth', 2,'Color', color)
hold on

ha_1(2).DisplayName = 'Best';
ha_2(2).DisplayName = 'Diametric';

%HIGHLIGHTING SPECIFIC SECTIONS OF THE PSTH TIMELINE
line([0 0], [0  100],'Color','k','LineStyle','--');
line([0.5 0.5], [0  100],'Color','k','LineStyle','--');
line([2 2], [0  100],'Color','k','LineStyle','--');
line([2.5 2.5], [0  100],'Color','k','LineStyle','--');
line([4 4], [0  100],'Color','k','LineStyle','--');

legend([ha_1(2), ha_2(2)])

% title(sprintf('PPC Remember 1st (BS=%d, nSolid=%d, nDashed=%d)',length(raw2), PPC1_best_n, PPC1_dia_n))
title(sprintf('PPC Remember 1st (BS=%d)',length(raw2)))
xlim([-1 4.5])
ylim([0  40])
xlabel('Time (s)', 'FontSize', 15)
ylabel('Discharge Rate (sp/s)', 'FontSize', 15)
% grid on
% box off

%
sgtitle('R1R2')

%%
clearvars -except R1R2_PFC_BS_cuerates R1R2_PFC_NS_cuerates R1R2_PPC_NS_cuerates R1R2_PFC_BS_fixrates ...
    R1R2_PFC_NS_fixrates R1R2_PPC_NS_fixrates R1R2_PPC_BS_cuerates R1R2_PPC_BS_fixrates ...
    R1R2_PFC_NS_cdrates R1R2_PPC_NS_cdrates R1R2_PFC_BS_cdrates R1R2_PPC_BS_cdrates R1R2_PFC_NS_allbestcue R1R2_PFC_NS_allbestcd ...
    R1R2_PPC_NS_allbestcue R1R2_PPC_NS_allbestcd R1R2_PFC_BS_allbestcue R1R2_PFC_BS_allbestcd R1R2_PPC_BS_allbestcue R1R2_PPC_BS_allbestcd 

%%

% Perform Wilcoxon rank-sum test
[p, h, stats] = ranksum(R1R2_PFC_NS_fixrates', R1R2_PFC_BS_fixrates')
[p, h, stats] = ranksum(R1R2_PPC_NS_fixrates', R1R2_PPC_BS_fixrates')

[p, h, stats] = ranksum(R1R2_PFC_NS_cuerates', R1R2_PFC_BS_cuerates')
[h, p, ci, stats] = ttest2(R1R2_PFC_NS_cuerates', R1R2_PFC_BS_cuerates')
meanDifference = mean(R1R2_PFC_NS_cuerates) - mean(R1R2_PFC_BS_cuerates);
pooledStdDev = sqrt(((std(R1R2_PFC_NS_cuerates)^2 + std(R1R2_PFC_BS_cuerates)^2) / 2));
cohenD = meanDifference / pooledStdDev


[p, h, stats] = ranksum(R1R2_PPC_NS_cuerates', R1R2_PPC_BS_cuerates')
[h, p, ci, stats] = ttest2(R1R2_PPC_NS_cuerates', R1R2_PPC_BS_cuerates')
meanDifference = mean(R1R2_PPC_NS_cuerates) - mean(R1R2_PPC_BS_cuerates);
pooledStdDev = sqrt(((std(R1R2_PPC_NS_cuerates)^2 + std(R1R2_PPC_BS_cuerates)^2) / 2));
cohenD = meanDifference / pooledStdDev

[p, h, stats] = ranksum(R1R2_PFC_NS_cdrates', R1R2_PFC_BS_cdrates')
[h, p, ci, stats] = ttest2(R1R2_PFC_NS_cdrates', R1R2_PFC_BS_cdrates')
meanDifference = mean(R1R2_PFC_NS_cdrates) - mean(R1R2_PFC_BS_cdrates);
pooledStdDev = sqrt(((std(R1R2_PFC_NS_cdrates)^2 + std(R1R2_PFC_BS_cdrates)^2) / 2));
cohenD = meanDifference / pooledStdDev

[p, h, stats] = ranksum(R1R2_PPC_NS_cdrates', R1R2_PPC_BS_cdrates')
[h, p, ci, stats] = ttest2(R1R2_PPC_NS_cdrates', R1R2_PPC_BS_cdrates')
meanDifference = mean(R1R2_PPC_NS_cdrates) - mean(R1R2_PPC_BS_cdrates);
pooledStdDev = sqrt(((std(R1R2_PPC_NS_cdrates)^2 + std(R1R2_PPC_BS_cdrates)^2) / 2));
cohenD = meanDifference / pooledStdDev

%%
[p, h, stats] = ranksum(R1R2_PFC_NS_allbestcue', R1R2_PFC_BS_allbestcue')
[h, p, ci, stats] = ttest2(R1R2_PFC_NS_allbestcue', R1R2_PFC_BS_allbestcue')
meanDifference = mean(R1R2_PFC_NS_allbestcue) - mean(R1R2_PFC_BS_allbestcue);
pooledStdDev = sqrt(((std(R1R2_PFC_NS_allbestcue)^2 + std(R1R2_PFC_BS_allbestcue)^2) / 2));
cohenD = meanDifference / pooledStdDev

[p, h, stats] = ranksum(R1R2_PPC_NS_allbestcue', R1R2_PPC_BS_allbestcue')
[h, p, ci, stats] = ttest2(R1R2_PPC_NS_allbestcue', R1R2_PPC_BS_allbestcue')
meanDifference = mean(R1R2_PPC_NS_allbestcue) - mean(R1R2_PPC_BS_allbestcue);
pooledStdDev = sqrt(((std(R1R2_PPC_NS_allbestcue)^2 + std(R1R2_PPC_BS_allbestcue)^2) / 2));
cohenD = meanDifference / pooledStdDev

[p, h, stats] = ranksum(R1R2_PFC_NS_allbestcd', R1R2_PFC_BS_allbestcd')
[h, p, ci, stats] = ttest2(R1R2_PFC_NS_allbestcd', R1R2_PFC_BS_allbestcd')
meanDifference = mean(R1R2_PFC_NS_allbestcd) - mean(R1R2_PFC_BS_allbestcd);
pooledStdDev = sqrt(((std(R1R2_PFC_NS_allbestcd)^2 + std(R1R2_PFC_BS_allbestcd)^2) / 2));
cohenD = meanDifference / pooledStdDev


[p, h, stats] = ranksum(R1R2_PPC_NS_allbestcd', R1R2_PPC_BS_allbestcd')
[h, p, ci, stats] = ttest2(R1R2_PPC_NS_allbestcd', R1R2_PPC_BS_allbestcd')
meanDifference = mean(R1R2_PPC_NS_allbestcd) - mean(R1R2_PPC_BS_allbestcd);
pooledStdDev = sqrt(((std(R1R2_PPC_NS_allbestcd)^2 + std(R1R2_PPC_BS_allbestcd)^2) / 2));
cohenD = meanDifference / pooledStdDev



%%
% fig1 = open('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\Figure5\R1R2_PSTH_sig.fig');
% set(fig1, 'Renderer', 'painters');
% exportgraphics(fig1, 'C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\Figure5\R1R2_PSTH_sig.emf', 'ContentType', 'vector');
% print(fig1, 'C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\Figure5\R1R2_PSTH_sig.pdf', '-vector', '-bestfit', '-dwinc');
% 
% 





