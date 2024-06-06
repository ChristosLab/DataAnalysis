clc; clear all; close all;

NS_BS='NS';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  PFC all neurons  %%%%%%%%%%%%%%%%%%%%

% [~,~,raw1] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\NS_BS_Database_20231109a.xlsx',['MSNG_PFC_', NS_BS]);
[~,~,raw1] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\AllSigNeuronsInfo.xlsx',['MSNG_PFC_', NS_BS]);
if strcmp(raw1{1,1}, 'Filename')
    raw1=raw1(2:end, :);
end
%%
neuron = 0; no_max_classes=0; no_min_classes=0;
for i = 1:length(raw1)
    fn_corr = [raw1{i,1}, '.mat'];
    try
        MatData_corr =load(['ALLDataCorrErr\', fn_corr]);
        disp(fn_corr)
    catch
        try
            fn_corr=[fn_corr(1:9), '1', fn_corr(10:end)];
            MatData_corr =load(['ALLDataCorrErr\', fn_corr]);
            disp(fn_corr)
        catch
            disp('Neuron file not found')
        end
    end
    [max_class_corr, max_class_corr_cuerate, fixrate, max_class_corr_cdrate, distant_class_corr] = MSNG_best_distant(fn_corr);

    if rem(max_class_corr,2)==0
        max_classes=[max_class_corr-1, max_class_corr];
    elseif rem(max_class_corr,2)~=0 && max_class_corr~=17
        max_classes=[max_class_corr, max_class_corr+1];
    elseif max_class_corr==17
        max_classes=max_class_corr;
    end

    if length( distant_class_corr)==1
         distant_classes=[ distant_class_corr,  distant_class_corr+1];
    elseif length( distant_class_corr)==2
         distant_classes=[ distant_class_corr(1),  distant_class_corr(1)+1, distant_class_corr(2),  distant_class_corr(2)+1];
    end
    
    no_max_classes = no_max_classes+length(max_classes);
    no_min_classes = no_min_classes+length(distant_classes);

    fprintf('BestCue=%d,   DistantCue=%d\n',max_class_corr(1), distant_class_corr(1));

    if isempty(max_class_corr) || isempty(distant_class_corr)
        disp('No Best and/or worst cue found')
        continue
    end

    bin_width = 0.05;  % 50 milliseconds bin
    dur = [-1-bin_width, 5-bin_width];
    bin_edges = dur(1):bin_width:dur(2); %
    n_bins_cue = numel(histcounts([], bin_edges));

    cls=0;
    if ~isempty(MatData_corr.MatData)
        for n= [max_classes, distant_classes]
            allTS_corr = [];
            m_counter_corr = 0;
            if n<=length(MatData_corr.MatData.class) && ~isempty(MatData_corr.MatData.class(n).ntr)
                for tr=1:length(MatData_corr.MatData.class(n).ntr)
                    if ~isempty(MatData_corr.MatData.class(n).ntr(tr).Cue_onT)
                        try
                            TS             = MatData_corr.MatData.class(n).ntr(tr).TS-MatData_corr.MatData.class(n).ntr(tr).Cue_onT;
                            allTS_corr     = [allTS_corr TS];
                            m_counter_corr = m_counter_corr + 1;
                            clear TS;
                        catch
                            disp('No spike data found')
                        end
                    end
                end
                ntrs_corr    = m_counter_corr;
            else
                disp('Class does not exist')
            end

            if ~isempty(allTS_corr)
                psth_temp_corr(cls + 1,:) = histc(allTS_corr,bin_edges)/(bin_width*ntrs_corr);
            else
                psth_temp_corr(cls + 1,:) = ones(1, length(bin_edges))*nan;
            end
            cls = cls +1;
        end
        if length(max_classes)==2 && length(distant_classes)==2
            pfc1_bestcue_m(neuron+1, :)     = psth_temp_corr(1, :);
            pfc1_bestcue_nm(neuron+1, :)    = psth_temp_corr(2, :);
            pfc1_distantcue_m(neuron+1, :)  = psth_temp_corr(3, :);
            pfc1_distantcue_nm(neuron+1, :) = psth_temp_corr(4, :);
            pfc1_bestcue_rate(neuron+1, :)  = max_class_corr_cuerate;
            pfc1_bestfix_rate(neuron+1, :)  = fixrate;
            pfc1_bestcd_rate(neuron+1, :)  = max_class_corr_cdrate;

        elseif length(max_classes)==1 && length(distant_classes)==4 && max_classes == 17
            pfc1_bestcue_m(neuron+1, :)     = psth_temp_corr(1, :);
            pfc1_distantcue_m(neuron+1, :)  = nanmean([psth_temp_corr(2, :); psth_temp_corr(4, :)]);
            pfc1_distantcue_nm(neuron+1, :) = nanmean([psth_temp_corr(3, :); psth_temp_corr(5, :)]);
            pfc1_bestcue_rate(neuron+1, :)  = max_class_corr_cuerate;
            pfc1_bestfix_rate(neuron+1, :)  = fixrate;
            pfc1_bestcd_rate(neuron+1, :)  = max_class_corr_cdrate;
        end


        neuron = neuron+1;
        disp(neuron)
        clear psth_temp_corr;

    else
        disp('Empty MatData File!!!');
    end
    clear MatData_corr;
    
end

%%
pfc1_bestcue     =[pfc1_bestcue_m; pfc1_bestcue_nm];
pfc1_bestcue     = pfc1_bestcue(~all(pfc1_bestcue == 0, 2), :);
PFC1_bestcue     = nanmean(pfc1_bestcue);
PFC1_bestcue_sem = nanstd(pfc1_bestcue)./sqrt(mean(sum(~isnan(pfc1_bestcue))));
PFC1_best_n      = mean(sum(~isnan(pfc1_bestcue)));

MSNG_PFC_NS_allbestcue = nanmean(pfc1_bestcue(:, 23:32), 2);
MSNG_PFC_NS_allbestcue = MSNG_PFC_NS_allbestcue(find(~isnan(MSNG_PFC_NS_allbestcue)));
MSNG_PFC_NS_allbestcd = nanmean(pfc1_bestcue(:, 33:92), 2);
MSNG_PFC_NS_allbestcd = MSNG_PFC_NS_allbestcd(find(~isnan(MSNG_PFC_NS_allbestcd)));

pfc1_distantcue     =[pfc1_distantcue_m; pfc1_distantcue_nm];
pfc1_distantcue     = pfc1_distantcue(~all(pfc1_distantcue == 0, 2), :);
PFC1_distantcue     = nanmean(pfc1_distantcue);
PFC1_distantcue_sem = nanstd(pfc1_distantcue)./sqrt(mean(sum(~isnan(pfc1_distantcue))));
PFC1_distant_n      = mean(sum(~isnan(pfc1_distantcue)));

MSNG_PFC_NS_cuerates = pfc1_bestcue_rate(find(~isnan(pfc1_bestcue_rate)));
MSNG_PFC_NS_fixrates = pfc1_bestfix_rate(find(~isnan(pfc1_bestfix_rate)));
MSNG_PFC_NS_cdrates  = pfc1_bestcd_rate(find(~isnan(pfc1_bestcd_rate)));

%%
clearvars -except PFC1_distantcue_sem PFC1_bestcue PFC1_bestcue_sem PFC1_distantcue raw1 NS_BS ...
    PFC1_best_n PFC1_distant_n MSNG_PFC_NS_cuerates MSNG_PFC_NS_fixrates MSNG_PFC_NS_cdrates MSNG_PFC_NS_allbestcue MSNG_PFC_NS_allbestcd
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  PPC all neurons  %%%%%%%%%%%%%%%%%%%%

% [~,~,raw2] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\NS_BS_Database_20231109a.xlsx',['MSNG_PPC_', NS_BS]);
[~,~,raw2] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\AllSigNeuronsInfo.xlsx',['MSNG_PPC_', NS_BS]);
if strcmp(raw2{1,1}, 'Filename')
    raw2=raw2(2:end, :);
end
%%
neuron = 0; no_max_classes=0; no_min_classes=0;
for i = 1:length(raw2)
    fn_corr = [raw2{i,1}, '.mat'];
    try
        MatData_corr =load(['ALLDataCorrErr\', fn_corr]);
        disp(fn_corr)
    catch
        try
            fn_corr=[fn_corr(1:9), '1', fn_corr(10:end)];
            MatData_corr =load(['ALLDataCorrErr\', fn_corr]);
            disp(fn_corr)
        catch
            disp('Neuron file not found')
        end
    end
    [max_class_corr, max_class_corr_cuerate, fixrate, max_class_corr_cdrate, distant_class_corr] = MSNG_best_distant(fn_corr);

    if rem(max_class_corr,2)==0
        max_classes=[max_class_corr-1, max_class_corr];
    elseif rem(max_class_corr,2)~=0 && max_class_corr~=17
        max_classes=[max_class_corr, max_class_corr+1];
    elseif max_class_corr==17
        max_classes=max_class_corr;
    end

    if length( distant_class_corr)==1
         distant_classes=[ distant_class_corr,  distant_class_corr+1];
    elseif length( distant_class_corr)==2
         distant_classes=[ distant_class_corr(1),  distant_class_corr(1)+1, distant_class_corr(2),  distant_class_corr(2)+1];
    end
    
    no_max_classes = no_max_classes+length(max_classes);
    no_min_classes = no_min_classes+length(distant_classes);

    fprintf('BestCue=%d,   DistantCue=%d\n',max_class_corr(1), distant_class_corr(1));

    if isempty(max_class_corr) || isempty(distant_class_corr)
        disp('No Best and/or worst cue found')
        continue
    end

    bin_width = 0.05;  % 50 milliseconds bin
    dur = [-1-bin_width, 5-bin_width];
    bin_edges = dur(1):bin_width:dur(2); %
    n_bins_cue = numel(histcounts([], bin_edges));

    cls=0;
    if ~isempty(MatData_corr.MatData)
        for n= [max_classes, distant_classes]
            allTS_corr = [];
            m_counter_corr = 0;
            if  n<=length(MatData_corr.MatData.class) && ~isempty(MatData_corr.MatData.class(n).ntr)
                for tr=1:length(MatData_corr.MatData.class(n).ntr)
                    if ~isempty(MatData_corr.MatData.class(n).ntr(tr).Cue_onT)
                        try
                            TS             = MatData_corr.MatData.class(n).ntr(tr).TS-MatData_corr.MatData.class(n).ntr(tr).Cue_onT;
                            allTS_corr     = [allTS_corr TS];
                            m_counter_corr = m_counter_corr + 1;
                            clear TS;
                        catch
                            disp('No spike data found')
                        end
                    end
                end
                ntrs_corr    = m_counter_corr;
            else
                disp('Class does not exist')
            end

            if ~isempty(allTS_corr)
                psth_temp_corr(cls + 1,:) = histc(allTS_corr,bin_edges)/(bin_width*ntrs_corr);
            else
                psth_temp_corr(cls + 1,:) = ones(1, length(bin_edges))*nan;
            end
            cls = cls +1;
        end
        if length(max_classes)==2 && length(distant_classes)==2
            ppc1_bestcue_m(neuron+1, :)     = psth_temp_corr(1, :);
            ppc1_bestcue_nm(neuron+1, :)    = psth_temp_corr(2, :);
            ppc1_distantcue_m(neuron+1, :)  = psth_temp_corr(3, :);
            ppc1_distantcue_nm(neuron+1, :) = psth_temp_corr(4, :);
            ppc1_bestcue_rate(neuron+1, :)  = max_class_corr_cuerate;
            ppc1_bestfix_rate(neuron+1, :)  = fixrate;
            ppc1_bestcd_rate(neuron+1, :)  = max_class_corr_cdrate;

        elseif length(max_classes)==1 && length(distant_classes)==4 && max_classes == 17
            ppc1_bestcue_m(neuron+1, :)     = psth_temp_corr(1, :);
            ppc1_distantcue_m(neuron+1, :)  = nanmean([psth_temp_corr(2, :); psth_temp_corr(4, :)]);
            ppc1_distantcue_nm(neuron+1, :) = nanmean([psth_temp_corr(3, :); psth_temp_corr(5, :)]);
            ppc1_bestcue_rate(neuron+1, :)  = max_class_corr_cuerate;
            ppc1_bestfix_rate(neuron+1, :)  = fixrate;
            ppc1_bestcd_rate(neuron+1, :)  = max_class_corr_cdrate;
        end


        neuron = neuron+1;
        disp(neuron)
        clear psth_temp_corr;

    else
        disp('Empty MatData File!!!');
    end
    clear MatData_corr;
    
end

%%
ppc1_bestcue     =[ppc1_bestcue_m; ppc1_bestcue_nm];
ppc1_bestcue     = ppc1_bestcue(~all(ppc1_bestcue == 0, 2), :);
PPC1_bestcue     = nanmean(ppc1_bestcue);
PPC1_bestcue_sem = nanstd(ppc1_bestcue)./sqrt(mean(sum(~isnan(ppc1_bestcue))));
PPC1_best_n      = mean(sum(~isnan(ppc1_bestcue)));

MSNG_PPC_NS_allbestcue = nanmean(ppc1_bestcue(:, 23:32), 2);
MSNG_PPC_NS_allbestcue = MSNG_PPC_NS_allbestcue(find(~isnan(MSNG_PPC_NS_allbestcue)));
MSNG_PPC_NS_allbestcd = nanmean(ppc1_bestcue(:, 33:92), 2);
MSNG_PPC_NS_allbestcd = MSNG_PPC_NS_allbestcd(find(~isnan(MSNG_PPC_NS_allbestcd)));

ppc1_distantcue     =[ppc1_distantcue_m; ppc1_distantcue_nm];
ppc1_distantcue     = ppc1_distantcue(~all(ppc1_distantcue == 0, 2), :);
PPC1_distantcue     = nanmean(ppc1_distantcue);
PPC1_distantcue_sem = nanstd(ppc1_distantcue)./sqrt(mean(sum(~isnan(ppc1_distantcue))));
PPC1_distant_n      = mean(sum(~isnan(ppc1_distantcue)));

MSNG_PPC_NS_cuerates = ppc1_bestcue_rate(find(~isnan(ppc1_bestcue_rate)));
MSNG_PPC_NS_fixrates = ppc1_bestfix_rate(find(~isnan(ppc1_bestfix_rate)));
MSNG_PPC_NS_cdrates  = ppc1_bestcd_rate(find(~isnan(ppc1_bestcd_rate)));

%%
figure
time = bin_edges+bin_width;
if strcmpi(NS_BS, 'BS')
    color='b';
%     colors = {[0.4 0.9 1], [0.8 0.8 1]};
    colors = {[0.3010 0.7450 0.9330], [0 1 1]};
elseif strcmpi(NS_BS, 'NS')
    color='r';
    colors = {[1, 0.6, 0.6], [1, 0.8, 0.8]};
end

popresults1 = smooth(PFC1_bestcue, 1)';
popresults2 = smooth(PFC1_distantcue, 1)';
subplot(221)

[ha_1] = shadedplot(time,popresults1+PFC1_bestcue_sem,popresults1-PFC1_bestcue_sem,colors{1},colors{1});
[ha_2] = shadedplot(time,popresults2+PFC1_distantcue_sem,popresults2-PFC1_distantcue_sem,colors{2},colors{2});

plot(time, popresults2, 'LineWidth', 2,'Color', color, 'LineStyle','--')
hold on
plot(time, popresults1, 'LineWidth', 2,'Color', color)
hold on

ha_1(2).DisplayName = 'Best Cue';
ha_2(2).DisplayName = 'Most distant Cue';

%HIGHLIGHTING SPECIFIC SECTIONS OF THE PSTH TIMELINE
line([0 0], [0 100],'Color','k','LineStyle','--');
line([0.5 0.5], [0 100],'Color','k','LineStyle','--');
line([3.5 3.5], [0 100],'Color','k','LineStyle','--');
line([4 4], [0 100],'Color','k','LineStyle','--');

legend([ha_1(2), ha_2(2)])

% title(sprintf('PFC (NS=%d, nSolid=%d, nDashed=%d)',length(raw1), PFC1_best_n, PFC1_distant_n))
title(sprintf('PFC (NS=%d)',length(raw1)))
xlim([-1 4.5])
ylim([0 20])
xlabel('Time (s)', 'FontSize', 15)
ylabel('Discharge Rate (sp/s)', 'FontSize', 15)
% grid on
% box off

%
subplot(223)

popresults1 = smooth(PPC1_bestcue, 1)';
popresults2 = smooth(PPC1_distantcue, 1)';

[ha_1] = shadedplot(time,popresults1+PPC1_bestcue_sem,popresults1-PPC1_bestcue_sem,colors{1},colors{1});
[ha_2] = shadedplot(time,popresults2+PPC1_distantcue_sem,popresults2-PPC1_distantcue_sem,colors{2},colors{2});

plot(time, popresults2, 'LineWidth', 2,'Color', color, 'LineStyle','--')
hold on
plot(time, popresults1, 'LineWidth', 2,'Color', color)
hold on

ha_1(2).DisplayName = 'Best Cue';
ha_2(2).DisplayName = 'Most distant Cue';

%HIGHLIGHTING SPECIFIC SECTIONS OF THE PSTH TIMELINE
line([0 0], [0 100],'Color','k','LineStyle','--');
line([0.5 0.5], [0  100],'Color','k','LineStyle','--');
line([3.5 3.5], [0  100],'Color','k','LineStyle','--');
line([4 4], [0  100],'Color','k','LineStyle','--');

legend([ha_1(2), ha_2(2)])

% title(sprintf('PPC (NS=%d, nSolid=%d, nDashed=%d)',length(raw2), PPC1_best_n, PPC1_distant_n))
title(sprintf('PPC (NS=%d)',length(raw2)))
xlim([-1 4.5])
ylim([0 20])
xlabel('Time (s)', 'FontSize', 15)
ylabel('Discharge Rate (sp/s)', 'FontSize', 15)
% grid on
% box off

% pause
%%
clearvars -except MSNG_PFC_NS_cuerates MSNG_PFC_NS_fixrates MSNG_PPC_NS_cuerates MSNG_PPC_NS_fixrates MSNG_PFC_NS_cdrates MSNG_PPC_NS_cdrates ...
    MSNG_PFC_NS_allbestcue MSNG_PFC_NS_allbestcd MSNG_PPC_NS_allbestcue MSNG_PPC_NS_allbestcd

% clc; clear all; close all;

NS_BS='BS';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  PFC all neurons  %%%%%%%%%%%%%%%%%%%%

% [~,~,raw1] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\NS_BS_Database_20231109a.xlsx',['MSNG_PFC_', NS_BS]);
[~,~,raw1] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\AllSigNeuronsInfo.xlsx',['MSNG_PFC_', NS_BS]);
if strcmp(raw1{1,1}, 'Filename')
    raw1=raw1(2:end, :);
end
%%
neuron = 0; no_max_classes=0; no_min_classes=0;
for i = 1:length(raw1)
    fn_corr = [raw1{i,1}, '.mat'];
    try
        MatData_corr =load(['ALLDataCorrErr\', fn_corr]);
        disp(fn_corr)
    catch
        try
            fn_corr=[fn_corr(1:9), '1', fn_corr(10:end)];
            MatData_corr =load(['ALLDataCorrErr\', fn_corr]);
            disp(fn_corr)
        catch
            disp('Neuron file not found')
        end
    end
    [max_class_corr, max_class_corr_cuerate, fixrate, max_class_corr_cdrate, distant_class_corr] = MSNG_best_distant(fn_corr);

    if rem(max_class_corr,2)==0
        max_classes=[max_class_corr-1, max_class_corr];
    elseif rem(max_class_corr,2)~=0 && max_class_corr~=17
        max_classes=[max_class_corr, max_class_corr+1];
    elseif max_class_corr==17
        max_classes=max_class_corr;
    end

    if length( distant_class_corr)==1
         distant_classes=[ distant_class_corr,  distant_class_corr+1];
    elseif length( distant_class_corr)==2
         distant_classes=[ distant_class_corr(1),  distant_class_corr(1)+1, distant_class_corr(2),  distant_class_corr(2)+1];
    end
    
    no_max_classes = no_max_classes+length(max_classes);
    no_min_classes = no_min_classes+length(distant_classes);

    fprintf('BestCue=%d,   DistantCue=%d\n',max_class_corr(1), distant_class_corr(1));

    if isempty(max_class_corr) || isempty(distant_class_corr)
        disp('No Best and/or worst cue found')
        continue
    end

    bin_width = 0.05;  % 50 milliseconds bin
    dur = [-1-bin_width, 5-bin_width];
    bin_edges = dur(1):bin_width:dur(2); %
    n_bins_cue = numel(histcounts([], bin_edges));

    cls=0;
    if ~isempty(MatData_corr.MatData)
        for n= [max_classes, distant_classes]
            allTS_corr = [];
            m_counter_corr = 0;
            if  n<=length(MatData_corr.MatData.class) && ~isempty(MatData_corr.MatData.class(n).ntr)
                for tr=1:length(MatData_corr.MatData.class(n).ntr)
                    if ~isempty(MatData_corr.MatData.class(n).ntr(tr).Cue_onT)
                        try
                            TS             = MatData_corr.MatData.class(n).ntr(tr).TS-MatData_corr.MatData.class(n).ntr(tr).Cue_onT;
                            allTS_corr     = [allTS_corr TS];
                            m_counter_corr = m_counter_corr + 1;
                            clear TS;
                        catch
                            disp('No spike data found')
                        end
                    end
                end
                ntrs_corr    = m_counter_corr;
            else
                disp('Class does not exist')
            end

            if ~isempty(allTS_corr)
                psth_temp_corr(cls + 1,:) = histc(allTS_corr,bin_edges)/(bin_width*ntrs_corr);
            else
                psth_temp_corr(cls + 1,:) = ones(1, length(bin_edges))*nan;
            end
            cls = cls +1;
        end
        if length(max_classes)==2 && length(distant_classes)==2
            pfc1_bestcue_m(neuron+1, :)     = psth_temp_corr(1, :);
            pfc1_bestcue_nm(neuron+1, :)    = psth_temp_corr(2, :);
            pfc1_distantcue_m(neuron+1, :)  = psth_temp_corr(3, :);
            pfc1_distantcue_nm(neuron+1, :) = psth_temp_corr(4, :);
            pfc1_bestcue_rate(neuron+1, :)  = max_class_corr_cuerate;
            pfc1_bestfix_rate(neuron+1, :)  = fixrate;
            pfc1_bestcd_rate(neuron+1, :)  = max_class_corr_cdrate;

        elseif length(max_classes)==1 && length(distant_classes)==4 && max_classes == 17
            pfc1_bestcue_m(neuron+1, :)     = psth_temp_corr(1, :);
            pfc1_distantcue_m(neuron+1, :)  = nanmean([psth_temp_corr(2, :); psth_temp_corr(4, :)]);
            pfc1_distantcue_nm(neuron+1, :) = nanmean([psth_temp_corr(3, :); psth_temp_corr(5, :)]);
            pfc1_bestcue_rate(neuron+1, :)  = max_class_corr_cuerate;
            pfc1_bestfix_rate(neuron+1, :)  = fixrate;
            pfc1_bestcd_rate(neuron+1, :)  = max_class_corr_cdrate;
        end


        neuron = neuron+1;
        disp(neuron)
        clear psth_temp_corr;

    else
        disp('Empty MatData File!!!');
    end
    clear MatData_corr;
    
end

%%
pfc1_bestcue     =[pfc1_bestcue_m; pfc1_bestcue_nm];
pfc1_bestcue     = pfc1_bestcue(~all(pfc1_bestcue == 0, 2), :);
PFC1_bestcue     = nanmean(pfc1_bestcue);
PFC1_bestcue_sem = nanstd(pfc1_bestcue)./sqrt(mean(sum(~isnan(pfc1_bestcue))));
PFC1_best_n      = mean(sum(~isnan(pfc1_bestcue)));

MSNG_PFC_BS_allbestcue = nanmean(pfc1_bestcue(:, 23:32), 2);
MSNG_PFC_BS_allbestcue = MSNG_PFC_BS_allbestcue(find(~isnan(MSNG_PFC_BS_allbestcue)));
MSNG_PFC_BS_allbestcd = nanmean(pfc1_bestcue(:, 33:92), 2);
MSNG_PFC_BS_allbestcd = MSNG_PFC_BS_allbestcd(find(~isnan(MSNG_PFC_BS_allbestcd)));

pfc1_distantcue     =[pfc1_distantcue_m; pfc1_distantcue_nm];
pfc1_distantcue     = pfc1_distantcue(~all(pfc1_distantcue == 0, 2), :);
PFC1_distantcue     = nanmean(pfc1_distantcue);
PFC1_distantcue_sem = nanstd(pfc1_distantcue)./sqrt(mean(sum(~isnan(pfc1_distantcue))));
PFC1_distant_n      = mean(sum(~isnan(pfc1_distantcue)));

MSNG_PFC_BS_cuerates = pfc1_bestcue_rate(find(~isnan(pfc1_bestcue_rate)));
MSNG_PFC_BS_fixrates = pfc1_bestfix_rate(find(~isnan(pfc1_bestfix_rate)));
MSNG_PFC_BS_cdrates  = pfc1_bestcd_rate(find(~isnan(pfc1_bestcd_rate)));

%%
clearvars -except PFC1_distantcue_sem PFC1_bestcue PFC1_bestcue_sem PFC1_distantcue raw1 NS_BS PFC1_best_n ...
    PFC1_distant_n MSNG_PFC_NS_cuerates MSNG_PPC_NS_cuerates MSNG_PFC_BS_cuerates MSNG_PFC_NS_fixrates ...
    MSNG_PPC_NS_fixrates MSNG_PFC_BS_fixrates MSNG_PFC_NS_cdrates MSNG_PPC_NS_cdrates MSNG_PFC_BS_cdrates ...
    MSNG_PFC_NS_allbestcue MSNG_PFC_NS_allbestcd MSNG_PPC_NS_allbestcue MSNG_PPC_NS_allbestcd MSNG_PFC_BS_allbestcue MSNG_PFC_BS_allbestcd
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  PPC all neurons  %%%%%%%%%%%%%%%%%%%%

% [~,~,raw2] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\NS_BS_Database_20231109a.xlsx',['MSNG_PPC_', NS_BS]);
[~,~,raw2] = xlsread('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\AllSigNeuronsInfo.xlsx',['MSNG_PPC_', NS_BS]);
if strcmp(raw2{1,1}, 'Filename')
    raw2=raw2(2:end, :);
end
%%
neuron = 0; no_max_classes=0; no_min_classes=0;
for i = 1:length(raw2)
    fn_corr = [raw2{i,1}, '.mat'];
    try
        MatData_corr =load(['ALLDataCorrErr\', fn_corr]);
        disp(fn_corr)
    catch
        try
            fn_corr=[fn_corr(1:9), '1', fn_corr(10:end)];
            MatData_corr =load(['ALLDataCorrErr\', fn_corr]);
            disp(fn_corr)
        catch
            disp('Neuron file not found')
        end
    end
    [max_class_corr, max_class_corr_cuerate, fixrate, max_class_corr_cdrate, distant_class_corr] = MSNG_best_distant(fn_corr);

    if rem(max_class_corr,2)==0
        max_classes=[max_class_corr-1, max_class_corr];
    elseif rem(max_class_corr,2)~=0 && max_class_corr~=17
        max_classes=[max_class_corr, max_class_corr+1];
    elseif max_class_corr==17
        max_classes=max_class_corr;
    end

    if length( distant_class_corr)==1
         distant_classes=[ distant_class_corr,  distant_class_corr+1];
    elseif length( distant_class_corr)==2
         distant_classes=[ distant_class_corr(1),  distant_class_corr(1)+1, distant_class_corr(2),  distant_class_corr(2)+1];
    end
    
    no_max_classes = no_max_classes+length(max_classes);
    no_min_classes = no_min_classes+length(distant_classes);

    fprintf('BestCue=%d,   DistantCue=%d\n',max_class_corr(1), distant_class_corr(1));

    if isempty(max_class_corr) || isempty(distant_class_corr)
        disp('No Best and/or worst cue found')
        continue
    end

    bin_width = 0.05;  % 50 milliseconds bin
    dur = [-1-bin_width, 5-bin_width];
    bin_edges = dur(1):bin_width:dur(2); %
    n_bins_cue = numel(histcounts([], bin_edges));

    cls=0;
    if ~isempty(MatData_corr.MatData)
        for n= [max_classes, distant_classes]
            allTS_corr = [];
            m_counter_corr = 0;
            if  n<=length(MatData_corr.MatData.class) && ~isempty(MatData_corr.MatData.class(n).ntr)
                for tr=1:length(MatData_corr.MatData.class(n).ntr)
                    if ~isempty(MatData_corr.MatData.class(n).ntr(tr).Cue_onT)
                        try
                            TS             = MatData_corr.MatData.class(n).ntr(tr).TS-MatData_corr.MatData.class(n).ntr(tr).Cue_onT;
                            allTS_corr     = [allTS_corr TS];
                            m_counter_corr = m_counter_corr + 1;
                            clear TS;
                        catch
                            disp('No spike data found')
                        end
                    end
                end
                ntrs_corr    = m_counter_corr;
            else
                disp('Class does not exist')
            end

            if ~isempty(allTS_corr)
                psth_temp_corr(cls + 1,:) = histc(allTS_corr,bin_edges)/(bin_width*ntrs_corr);
            else
                psth_temp_corr(cls + 1,:) = ones(1, length(bin_edges))*nan;
            end
            cls = cls +1;
        end
        if length(max_classes)==2 && length(distant_classes)==2
            ppc1_bestcue_m(neuron+1, :)     = psth_temp_corr(1, :);
            ppc1_bestcue_nm(neuron+1, :)    = psth_temp_corr(2, :);
            ppc1_distantcue_m(neuron+1, :)  = psth_temp_corr(3, :);
            ppc1_distantcue_nm(neuron+1, :) = psth_temp_corr(4, :);
            ppc1_bestcue_rate(neuron+1, :)  = max_class_corr_cuerate;
            ppc1_bestfix_rate(neuron+1, :)  = fixrate;
            ppc1_bestcd_rate(neuron+1, :)   = max_class_corr_cdrate;

        elseif length(max_classes)==1 && length(distant_classes)==4 && max_classes == 17
            ppc1_bestcue_m(neuron+1, :)     = psth_temp_corr(1, :);
            ppc1_distantcue_m(neuron+1, :)  = nanmean([psth_temp_corr(2, :); psth_temp_corr(4, :)]);
            ppc1_distantcue_nm(neuron+1, :) = nanmean([psth_temp_corr(3, :); psth_temp_corr(5, :)]);
            ppc1_bestcue_rate(neuron+1, :)  = max_class_corr_cuerate;
            ppc1_bestfix_rate(neuron+1, :)  = fixrate;
            ppc1_bestcd_rate(neuron+1, :)   = max_class_corr_cdrate;
        end


        neuron = neuron+1;
        disp(neuron)
        clear psth_temp_corr;

    else
        disp('Empty MatData File!!!');
    end
    clear MatData_corr;
    
end

%%
ppc1_bestcue     =[ppc1_bestcue_m; ppc1_bestcue_nm];
ppc1_bestcue     = ppc1_bestcue(~all(ppc1_bestcue == 0, 2), :);
PPC1_bestcue     = nanmean(ppc1_bestcue);
PPC1_bestcue_sem = nanstd(ppc1_bestcue)./sqrt(mean(sum(~isnan(ppc1_bestcue))));
PPC1_best_n      = mean(sum(~isnan(ppc1_bestcue)));

MSNG_PPC_BS_allbestcue = nanmean(ppc1_bestcue(:, 23:32), 2);
MSNG_PPC_BS_allbestcue = MSNG_PPC_BS_allbestcue(find(~isnan(MSNG_PPC_BS_allbestcue)));
MSNG_PPC_BS_allbestcd = nanmean(ppc1_bestcue(:, 33:92), 2);
MSNG_PPC_BS_allbestcd = MSNG_PPC_BS_allbestcd(find(~isnan(MSNG_PPC_BS_allbestcd)));

ppc1_distantcue     =[ppc1_distantcue_m; ppc1_distantcue_nm];
ppc1_distantcue     = ppc1_distantcue(~all(ppc1_distantcue == 0, 2), :);
PPC1_distantcue     = nanmean(ppc1_distantcue);
PPC1_distantcue_sem = nanstd(ppc1_distantcue)./sqrt(mean(sum(~isnan(ppc1_distantcue))));
PPC1_distant_n      = mean(sum(~isnan(ppc1_distantcue)));

MSNG_PPC_BS_cuerates = ppc1_bestcue_rate(find(~isnan(ppc1_bestcue_rate)));
MSNG_PPC_BS_fixrates = ppc1_bestfix_rate(find(~isnan(ppc1_bestfix_rate)));
MSNG_PPC_BS_cdrates  = ppc1_bestcd_rate(find(~isnan(ppc1_bestcd_rate)));


%%
time = bin_edges+bin_width;
if strcmpi(NS_BS, 'BS')
    color='b';
%     colors = {[0.4 0.9 1], [0.8 0.8 1]};
    colors = {[0.3010 0.7450 0.9330], [0 1 1]};
elseif strcmpi(NS_BS, 'NS')
    color='r';
    colors = {[1, 0.6, 0.6], [1, 0.8, 0.8]};
end

popresults1 = smooth(PFC1_bestcue, 1)';
popresults2 = smooth(PFC1_distantcue, 1)';
subplot(222)

[ha_1] = shadedplot(time,popresults1+PFC1_bestcue_sem,popresults1-PFC1_bestcue_sem,colors{1},colors{1});
[ha_2] = shadedplot(time,popresults2+PFC1_distantcue_sem,popresults2-PFC1_distantcue_sem,colors{2},colors{2});

plot(time, popresults2, 'LineWidth', 2,'Color', color, 'LineStyle','--')
hold on
plot(time, popresults1, 'LineWidth', 2,'Color', color)
hold on

ha_1(2).DisplayName = 'Best Cue';
ha_2(2).DisplayName = 'Most distant Cue';

%HIGHLIGHTING SPECIFIC SECTIONS OF THE PSTH TIMELINE
line([0 0], [0  100],'Color','k','LineStyle','--');
line([0.5 0.5], [0  100],'Color','k','LineStyle','--');
line([3.5 3.5], [0  100],'Color','k','LineStyle','--');
line([4 4], [0  100],'Color','k','LineStyle','--');

legend([ha_1(2), ha_2(2)])

% title(sprintf('PFC (BS=%d, nSolid=%d, nDashed=%d)',length(raw1), PFC1_best_n, PFC1_distant_n))
title(sprintf('PFC (BS=%d)',length(raw1)))
xlim([-1 4.5])
ylim([0 20])
xlabel('Time (s)', 'FontSize', 15)
ylabel('Discharge Rate (sp/s)', 'FontSize', 15)
% grid on
% box off

%
subplot(224)

popresults1 = smooth(PPC1_bestcue, 1)';
popresults2 = smooth(PPC1_distantcue, 1)';

[ha_1] = shadedplot(time,popresults1+PPC1_bestcue_sem,popresults1-PPC1_bestcue_sem,colors{1},colors{1});
[ha_2] = shadedplot(time,popresults2+PPC1_distantcue_sem,popresults2-PPC1_distantcue_sem,colors{2},colors{2});

plot(time, popresults2, 'LineWidth', 2,'Color', color, 'LineStyle','--')
hold on
plot(time, popresults1, 'LineWidth', 2,'Color', color)
hold on

ha_1(2).DisplayName = 'Best Cue';
ha_2(2).DisplayName = 'Most distant Cue';

%HIGHLIGHTING SPECIFIC SECTIONS OF THE PSTH TIMELINE
line([0 0], [0  100],'Color','k','LineStyle','--');
line([0.5 0.5], [0  100],'Color','k','LineStyle','--');
line([3.5 3.5], [0  100],'Color','k','LineStyle','--');
line([4 4], [0  100],'Color','k','LineStyle','--');

legend([ha_1(2), ha_2(2)])

% title(sprintf('PPC (BS=%d, nSolid=%d, nDashed=%d)',length(raw2), PPC1_best_n, PPC1_distant_n))
title(sprintf('PPC (BS=%d)',length(raw2)))
xlim([-1 4.5])
ylim([0 20])
xlabel('Time (s)', 'FontSize', 15)
ylabel('Discharge Rate (sp/s)', 'FontSize', 15)
% grid on
% box off
sgtitle('MSNG')


%%
clearvars -except MSNG_PFC_NS_cuerates MSNG_PPC_NS_cuerates MSNG_PFC_BS_cuerates MSNG_PPC_BS_cuerates ...
    MSNG_PFC_NS_fixrates MSNG_PPC_NS_fixrates MSNG_PFC_BS_fixrates MSNG_PPC_BS_fixrates MSNG_PFC_NS_cdrates MSNG_PPC_NS_cdrates ...
    MSNG_PFC_BS_cdrates MSNG_PPC_BS_cdrates MSNG_PFC_NS_allbestcue MSNG_PFC_NS_allbestcd MSNG_PPC_NS_allbestcue MSNG_PPC_NS_allbestcd ...
    MSNG_PFC_BS_allbestcue MSNG_PFC_BS_allbestcd MSNG_PPC_BS_allbestcue MSNG_PPC_BS_allbestcd
%%

% Perform Wilcoxon rank-sum test
[p, h, stats] = ranksum(MSNG_PFC_NS_fixrates', MSNG_PFC_BS_fixrates')
[p, h, stats] = ranksum(MSNG_PPC_NS_fixrates', MSNG_PPC_BS_fixrates')

[p, h, stats] = ranksum(MSNG_PFC_NS_cuerates', MSNG_PFC_BS_cuerates')
[h, p, ci, stats] = ttest2(MSNG_PFC_NS_cuerates', MSNG_PFC_BS_cuerates')
meanDifference = mean(MSNG_PFC_NS_cuerates) - mean(MSNG_PFC_BS_cuerates);
pooledStdDev = sqrt(((std(MSNG_PFC_NS_cuerates)^2 + std(MSNG_PFC_BS_cuerates)^2) / 2));
cohenD = meanDifference / pooledStdDev


[p, h, stats] = ranksum(MSNG_PPC_NS_cuerates', MSNG_PPC_BS_cuerates')
[h, p, ci, stats] = ttest2(MSNG_PPC_NS_cuerates', MSNG_PPC_BS_cuerates')
meanDifference = mean(MSNG_PPC_NS_cuerates) - mean(MSNG_PPC_BS_cuerates);
pooledStdDev = sqrt(((std(MSNG_PPC_NS_cuerates)^2 + std(MSNG_PPC_BS_cuerates)^2) / 2));
cohenD = meanDifference / pooledStdDev

[p, h, stats] = ranksum(MSNG_PFC_NS_cdrates', MSNG_PFC_BS_cdrates')
[h, p, ci, stats] = ttest2(MSNG_PFC_NS_cdrates', MSNG_PFC_BS_cdrates')
meanDifference = mean(MSNG_PFC_NS_cdrates) - mean(MSNG_PFC_BS_cdrates);
pooledStdDev = sqrt(((std(MSNG_PFC_NS_cdrates)^2 + std(MSNG_PFC_BS_cdrates)^2) / 2));
cohenD = meanDifference / pooledStdDev


[p, h, stats] = ranksum(MSNG_PPC_NS_cdrates', MSNG_PPC_BS_cdrates')
[h, p, ci, stats] = ttest2(MSNG_PPC_NS_cdrates', MSNG_PPC_BS_cdrates')
meanDifference = mean(MSNG_PPC_NS_cdrates) - mean(MSNG_PPC_BS_cdrates);
pooledStdDev = sqrt(((std(MSNG_PPC_NS_cdrates)^2 + std(MSNG_PPC_BS_cdrates)^2) / 2));
cohenD = meanDifference / pooledStdDev

%%
[p, h, stats] = ranksum(MSNG_PFC_NS_allbestcue', MSNG_PFC_BS_allbestcue')
[p, h, ci, stats] = ttest2(MSNG_PFC_NS_allbestcue', MSNG_PFC_BS_allbestcue')
meanDifference = mean(MSNG_PFC_NS_allbestcue) - mean(MSNG_PFC_BS_allbestcue);
pooledStdDev = sqrt(((std(MSNG_PFC_NS_allbestcue)^2 + std(MSNG_PFC_BS_allbestcue)^2) / 2));
cohenD = meanDifference / pooledStdDev

[p, h, stats] = ranksum(MSNG_PPC_NS_allbestcue', MSNG_PPC_BS_allbestcue')
[p, h, ci, stats] = ttest2(MSNG_PPC_NS_allbestcue', MSNG_PPC_BS_allbestcue')
meanDifference = mean(MSNG_PPC_NS_allbestcue) - mean(MSNG_PPC_BS_allbestcue);
pooledStdDev = sqrt(((std(MSNG_PPC_NS_allbestcue)^2 + std(MSNG_PPC_BS_allbestcue)^2) / 2));
cohenD = meanDifference / pooledStdDev

[p, h, stats] = ranksum(MSNG_PFC_NS_allbestcd', MSNG_PFC_BS_allbestcd')
[p, h, ci, stats] = ttest2(MSNG_PFC_NS_allbestcd', MSNG_PFC_BS_allbestcd')
meanDifference = mean(MSNG_PFC_NS_allbestcd) - mean(MSNG_PFC_BS_allbestcd);
pooledStdDev = sqrt(((std(MSNG_PFC_NS_allbestcd)^2 + std(MSNG_PFC_BS_allbestcd)^2) / 2));
cohenD = meanDifference / pooledStdDev


[p, h, stats] = ranksum(MSNG_PPC_NS_allbestcd', MSNG_PPC_BS_allbestcd')
[p, h, ci, stats] = ttest2(MSNG_PPC_NS_allbestcd', MSNG_PPC_BS_allbestcd')
meanDifference = mean(MSNG_PPC_NS_allbestcd) - mean(MSNG_PPC_BS_allbestcd);
pooledStdDev = sqrt(((std(MSNG_PPC_NS_allbestcd)^2 + std(MSNG_PPC_BS_allbestcd)^2) / 2));
cohenD = meanDifference / pooledStdDev






%%
% fig1 = open('C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\Figure4\MSNG_PSTH_sig.fig');
% set(fig1, 'Renderer', 'painters');
% exportgraphics(fig1, 'C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\Figure4\MSNG_PSTH_sig.emf', 'ContentType', 'vector');
% print(fig1, 'C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\Figure4\MSNG_PSTH_sig.pdf', '-vector', '-bestfit', '-dwinc');
% 




