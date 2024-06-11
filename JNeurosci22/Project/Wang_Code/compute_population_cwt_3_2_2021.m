clear
close all
zw_setpath
%%
fname_ = 'complete_cwt_repo_3_2_2-21.mat';
load(fullfile(project_dir, output_database, fname_), 'cwt_repo');
fname_ = 'lfp_tbl.mat';
load(fullfile(project_dir, output_database, fname_), 'lfp_tbl');
%%  CWT parameters
fs = 500;
cue_dur = [-1 + 1/fs, 4]; % Signal length 
sac_dur = [-3 + 1/fs, 2];
normalizer = 4;
kernel_flag = 3;
down_sample = 10;
f_range = 2:2:128;
avg_method = 1; %   0 - arithmetic; 1 - geometric
%%
cue_n_sample = (diff(cue_dur)*fs + 1)/down_sample;
sac_n_sample = (diff(sac_dur)*fs + 1)/down_sample;
%%
target_frs = [4,8; 8,16; 16, 32; 32, numel(f_range)];
baseline_bin = [26, 50];
%%
temp_cwt_cue = zeros([numel(cwt_repo), numel(f_range), cue_n_sample]);
temp_cwt_sac = zeros([numel(cwt_repo), numel(f_range), sac_n_sample]);
for i = 1:numel(cwt_repo)
    if lfp_tbl.task_id(i) == 2
        temp_cwt_cue(i, :, :) = nan(numel(f_range), cue_n_sample);
        temp_cwt_sac(i, :, :) = nan(numel(f_range), sac_n_sample);
        continue %  Skipping anti-saccade sessions
    end
    %   CueOn aligned
    if isempty(cwt_repo(i).class(1).cue_cwt)
        temp_cwt_cue(i, :, :) = nan(numel(f_range), cue_n_sample);
    else
        single_cwt_cue = zeros(0, numel(f_range), cue_n_sample);
        for j = 1:numel(cwt_repo(i).class)
            single_cwt_cue = [single_cwt_cue; cwt_repo(i).class(j).cue_cwt];
        end
        temp_cwt_cue(i, :, :) = nanmean(single_cwt_cue, 1);
    end
    %   SaccadeOn aligned
    if isempty(cwt_repo(i).class(1).sac_cwt)
        temp_cwt_sac(i, :, :) = nan(numel(f_range), sac_n_sample);
    else
        single_cwt_sac = zeros(0, numel(f_range), sac_n_sample);
        for j = 1:numel(cwt_repo(i).class)
            single_cwt_sac = [single_cwt_sac; cwt_repo(i).class(j).sac_cwt];
        end
        temp_cwt_sac(i, :, :) = nanmean(single_cwt_sac, 1);
    end
    i
end
%%
% save(fullfile(project_dir, output_database, 'temp_cwt.mat'), 'temp_cwt');
%%
figure(); plot(mean(squeeze(nanmean(temp_cwt_cue(t1, 12:17, :), 1)), 1));ylim([6, 8]*10^-3);
figure(); plot(mean(squeeze(nanmean(temp_cwt_cue(t2, 12:17, :), 1)), 1));ylim([6, 8]*10^-3);
%%
figure(); plot(mean(squeeze(nanmean(temp_cwt_cue(t1, 32:64, :), 1)), 1));ylim([.7, 1.4]*10^-3);
figure(); plot(mean(squeeze(nanmean(temp_cwt_cue(t2, 32:64, :), 1)), 1));ylim([.7, 1.4]*10^-3);
%%
temp_cwt_cue_b = zeros([numel(cwt_repo), numel(f_range), cue_n_sample]);
temp_cwt_sac_b = zeros([numel(cwt_repo), numel(f_range), sac_n_sample]);
for i = 1:numel(cwt_repo)
    if lfp_tbl.task_id(i) == 2
        temp_cwt_cue_b(i, :, :) = nan(numel(f_range), cue_n_sample);
        temp_cwt_sac_b(i, :, :) = nan(numel(f_range), sac_n_sample);
        continue
    end
    %   Cue aligned
    if isempty(cwt_repo(i).class(1).cue_cwt)
        temp_cwt_cue_b(i, :, :) = nan(numel(f_range), cue_n_sample);
    else
        single_cwt_cue = zeros(0, numel(f_range), cue_n_sample);
        for j = 1:numel(cwt_repo(i).class)
            single_cwt_cue = [single_cwt_cue; cwt_repo(i).class(j).cue_cwt./mean(cwt_repo(i).class(j).cue_cwt(:, :, baseline_bin(1):baseline_bin(2)), 3)];
        end
        temp_cwt_cue_b(i, :, :) = nanmean(single_cwt_cue, 1);
    end
    %   Sac aligned
    if isempty(cwt_repo(i).class(1).sac_cwt)
        temp_cwt_sac_b(i, :, :) = nan(numel(f_range), sac_n_sample);
    else
        single_cwt_sac = zeros(0, numel(f_range), sac_n_sample);
        for j = 1:numel(cwt_repo(i).class)
            %   NOTICE: CUE aligned baseline used for SAC aligned signal
            single_cwt_sac = [single_cwt_sac; cwt_repo(i).class(j).sac_cwt./mean(cwt_repo(i).class(j).cue_cwt(:, :, baseline_bin(1):baseline_bin(2)), 3)];
        end
        temp_cwt_sac_b(i, :, :) = nanmean(single_cwt_sac, 1);
    end

    i
end
%%
figure(); plot(mean(squeeze(nanmean(temp_cwt_cue_b(t1, 12:17, :), 1)), 1))%;ylim([6, 8]*10^-3);
figure(); plot(mean(squeeze(nanmean(temp_cwt_cue_b(t2, 12:17, :), 1)), 1))%;ylim([6, 8]*10^-3);
%%
figure(); plot(mean(squeeze(nanmean(temp_cwt_cue_b([site_tuning_cat{1, [1,3,4], :, 1}], 32:64, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
hold on
plot(mean(squeeze(nanmean(temp_cwt_cue_b([site_tuning_cat{1, [1,3,4], :, 2}], 32:64, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
plot(mean(squeeze(nanmean(temp_cwt_cue_b([site_tuning_cat{2, [1,3,4], :, 1}], 32:64, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
plot(mean(squeeze(nanmean(temp_cwt_cue_b([site_tuning_cat{2, [1,3,4], :, 2}], 32:64, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
%%
figure(); plot(mean(squeeze(nanmean(temp_cwt_sac_b([site_tuning_cat{1, [1,3,4], :, 1}], 32:64, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
hold on
plot(mean(squeeze(nanmean(temp_cwt_sac_b([site_tuning_cat{1, [1,3,4], :, 2}], 32:64, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
plot(mean(squeeze(nanmean(temp_cwt_sac_b([site_tuning_cat{2, [1,3,4], :, 1}], 32:64, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
plot(mean(squeeze(nanmean(temp_cwt_sac_b([site_tuning_cat{2, [1,3,4], :, 2}], 32:64, :), 1)), 1))%;ylim([.7, 1.4]*10^-3);
%%
temp_baseline = zeros([numel(cwt_repo), numel(f_range)]);
for i = 1:numel(cwt_repo)
    if lfp_tbl.task_id(i) == 2
        temp_baseline(i, :) = nan(numel(f_range), 1);
    end
    if isempty(cwt_repo(i).class(1).cue_cwt)
        temp_baseline(i, :) = nan(numel(f_range), 1);
    else
        single_cwt_cue = zeros(0, numel(f_range));
        for j = 1:numel(cwt_repo(i).class)
        single_cwt_cue = [single_cwt_cue; mean(cwt_repo(i).class(j).cue_cwt(:, :, baseline_bin(1):baseline_bin(2)), 3)];
        end
        temp_baseline(i, :) = nanmean(single_cwt_cue, 1);
    end
    i
end
%%
fname_ = 'temp_cwt_3_2_2021.mat';
save(fullfile(project_dir, output_database, fname_), 'temp_baseline', 'temp_cwt_cue_b', 'temp_cwt_sac_b', 'temp_cwt_cue', 'temp_cwt_sac');
%%
figure(); plot(mean(squeeze(nanmean(temp_cwt_cue(t1, 12:17, :)./temp_baseline(t1, 12:17), 1)), 1));%ylim([6, 8]*10^-3);
hold on
plot(mean(squeeze(nanmean(temp_cwt_cue(t2, 12:17, :)./temp_baseline(t2, 12:17), 1)), 1));%ylim([6, 8]*10^-3);
figure(); plot(mean(squeeze(nanmean(temp_cwt_cue(t1, 32:64, :)./temp_baseline(t1, 32:64), 1)), 1));%ylim([6, 8]*10^-3);
hold on
plot(mean(squeeze(nanmean(temp_cwt_cue(t2, 32:64, :)./temp_baseline(t2, 32:64), 1)), 1));%ylim([6, 8]*10^-3)
figure(); plot(mean(squeeze(nanmean(temp_cwt_cue(t1, 8:16, :)./temp_baseline(t1, 8:16), 1)), 1));%ylim([6, 8]*10^-3);
hold on
plot(mean(squeeze(nanmean(temp_cwt_cue(t2, 8:16, :)./temp_baseline(t2, 8:16), 1)), 1));%ylim([6, 8]*10^-3)
%%
