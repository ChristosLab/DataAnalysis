%%
bin_width = 0.05;  % 50-millisecond bin
bin_width_for_pev = 0.25; % 250-millisecond bin for PEV
step = 0.02; %  20 millisecond step matching that of SPD
t_range_cue = [-1, 4];
t_range_sac = [-2, 3];
bin_edges_cue = t_range_cue(1):step:t_range_cue(2);
bin_edges_sac = t_range_sac(1):step:t_range_sac(2);
n_bins_cue = numel(histcounts([], bin_edges_cue));
n_bins_sac = numel(histcounts([], bin_edges_sac));
%%
neuron_list_ = 1:size(dodr_neuron_tbl, 1);
for i = neuron_list_
    if numel(dodr_neuron_repo(i).class) < 1
        continue
    end
    i
    for j = 1:numel(dodr_neuron_repo(i).class)
        if isempty(dodr_neuron_repo(i).class(j).ntr)
            continue
        end
        psth_cue_         = zeros(numel(dodr_neuron_repo(i).class(j).ntr), n_bins_cue);
        psth_sac_         = zeros(numel(dodr_neuron_repo(i).class(j).ntr), n_bins_sac);
        psth_cue_for_pev_ = zeros(numel(dodr_neuron_repo(i).class(j).ntr), n_bins_cue);
        psth_cue_for_sac_ = zeros(numel(dodr_neuron_repo(i).class(j).ntr), n_bins_sac);
        for m = 1:numel(dodr_neuron_repo(i).class(j).ntr)
            TS_cue_ = dodr_neuron_repo(i).class(j).ntr(m).TS - dodr_neuron_repo(i).class(j).ntr(m).Cue_onT;
            psth_cue_(m, :)         = zw_spike_time_to_psth(TS_cue_, bin_width, step, t_range_cue);
            psth_cue_for_pev_(m, :) = zw_spike_time_to_psth(TS_cue_, bin_width_for_pev, step, t_range_cue);
            if isfield(dodr_neuron_repo(i).class(j).ntr(m), 'Saccade_onT')
                if ~isempty(dodr_neuron_repo(i).class(j).ntr(m).Saccade_onT)
                    TS_sac_ = dodr_neuron_repo(i).class(j).ntr(m).TS - dodr_neuron_repo(i).class(j).ntr(m).Saccade_onT;
                    psth_sac_(m, :)         = zw_spike_time_to_psth(TS_sac_, bin_width, step, t_range_sac);
                    psth_sac_for_pev_(m, :) = zw_spike_time_to_psth(TS_sac_, bin_width_for_pev, step, t_range_sac);
                else
                    psth_sac_(m, :) = NaN;
                    psth_sac_for_pev_(m, :) = NaN;
                end
            else
                psth_sac_(m, :) = NaN;
                psth_sac_for_pev_(m, :) = NaN;
            end
        end
        dodr_neuron_repo(i).class(j).psth_cue_upsampled          = psth_cue_;
        dodr_neuron_repo(i).class(j).psth_cue_upsampled_for_pev  = psth_cue_for_pev_;
        dodr_neuron_repo(i).class(j).psth_sac_upsampled          = psth_sac_;
        dodr_neuron_repo(i).class(j).psth_sac_upsampled_for_pev  = psth_sac_for_pev_;
    end
end
%%
fname_ = 'dodr_neuron_repo_w_psth.mat';
save(fullfile(project_dir, output_database, fname_), 'dodr_neuron_repo');
%%
for i = neuron_list_
    for j = 1:8
%         temp_psth(i, j, :)                 = nanmean(dodr_neuron_repo(i).class(j).psth_cue_upsampled - mean(dodr_neuron_repo(i).class(j).psth_cue_upsampled(:, 26:50), 2), 1);
        temp_psth_1s_baseline(i, j, :)     = nanmean(dodr_neuron_repo(i).class(j).psth_cue_upsampled - mean(dodr_neuron_repo(i).class(j).psth_cue_upsampled(:, 1:50), 2), 1);
        temp_psth_raw(i, j, :)             = nanmean(dodr_neuron_repo(i).class(j).psth_cue_upsampled, 1);
%         temp_psth_sac(i, j, :)             = nanmean(dodr_neuron_repo(i).class(j).psth_sac_upsampled - mean(dodr_neuron_repo(i).class(j).psth_cue_upsampled(:, 26:50), 2), 1);
        temp_psth_1s_baseline_sac(i, j, :) = nanmean(dodr_neuron_repo(i).class(j).psth_sac_upsampled - mean(dodr_neuron_repo(i).class(j).psth_cue_upsampled(:, 1:50), 2), 1);
        temp_psth_raw_sac(i, j, :)         = nanmean(dodr_neuron_repo(i).class(j).psth_sac_upsampled, 1);
    end
end
%%
for i = 1:size(temp_psth_raw, 1)
%   First Stimulus
    [y_, i_] = max(sum(temp_psth_raw(i, :,51:100), 3));
    dodr_neuron_best_class_raw_target(i) = i_;
%   Second Stimulus
    [y_, i_] = max(sum(temp_psth_raw(i, :,101:150), 3));
    dodr_neuron_best_class_raw_distractor(i) = i_;
end
%%
best_psth_raw_target         = nan(numel(dodr_neuron_repo), n_bins_cue);
best_psth_raw_sac_target     = nan(numel(dodr_neuron_repo), n_bins_sac);
best_psth_raw_distractor     = nan(numel(dodr_neuron_repo), n_bins_cue);
best_psth_raw_sac_distractor = nan(numel(dodr_neuron_repo), n_bins_sac);
for i = neuron_list_
    best_psth_raw_target(i, :)         = temp_psth_raw(i, dodr_neuron_best_class_raw_target(i), :);
    best_psth_raw_sac_target(i, :)     = temp_psth_raw_sac(i, dodr_neuron_best_class_raw_target(i), :);
    best_psth_raw_distractor(i, :)     = temp_psth_raw(i, dodr_neuron_best_class_raw_distractor(i), :);
    best_psth_raw_sac_distractor(i, :) = temp_psth_raw_sac(i, dodr_neuron_best_class_raw_distractor(i), :);
end
%%
for i = 1:size(temp_psth_1s_baseline, 1)
%   First Stimulus
    [y_, i_] = max(sum(temp_psth_1s_baseline(i, :,51:100), 3));
    dodr_neuron_best_class_1s_baseline_target(i) = i_;
%   Second Stimulus
    [y_, i_] = max(sum(temp_psth_1s_baseline(i, :,101:150), 3));
    dodr_neuron_best_class_1s_baseline_distractor(i) = i_;
end
%%
best_psth_1s_baseline_target         = nan(numel(dodr_neuron_repo), n_bins_cue);
best_psth_1s_baseline_sac_target     = nan(numel(dodr_neuron_repo), n_bins_sac);
best_psth_1s_baseline_distractor     = nan(numel(dodr_neuron_repo), n_bins_cue);
best_psth_1s_baseline_sac_distractor = nan(numel(dodr_neuron_repo), n_bins_sac);
for i = neuron_list_
    best_psth_1s_baseline_target(i, :)         = temp_psth_1s_baseline(i, dodr_neuron_best_class_raw_target(i), :);
    best_psth_1s_baseline_sac_target(i, :)     = temp_psth_1s_baseline(i, dodr_neuron_best_class_raw_target(i), :);
    best_psth_1s_baseline_distractor(i, :)     = temp_psth_1s_baseline(i, dodr_neuron_best_class_raw_distractor(i), :);
    best_psth_1s_baseline_sac_distractor(i, :) = temp_psth_1s_baseline(i, dodr_neuron_best_class_raw_distractor(i), :);
end
%%
best_psth_1s_baseline_opp         = nan(numel(dodr_neuron_repo), n_bins_cue);
best_psth_1s_baseline_sac_opp     = nan(numel(dodr_neuron_repo), n_bins_sac);
for i = neuron_list_
    current_class    = dodr_neuron_best_class_raw_target(i);
    current_class_op = mod(current_class + 4, 8);
    if current_class_op == 0
        current_class_op = 8;
    end
    best_psth_1s_baseline_opp(i, :)         = temp_psth_1s_baseline(i, current_class_op, :);
    best_psth_1s_baseline_sac_opp(i, :)     = temp_psth_1s_baseline_sac(i, current_class_op, :);
end
%%
best_psth_raw_opp         = nan(numel(dodr_neuron_repo), n_bins_cue);
best_psth_raw_sac_opp     = nan(numel(dodr_neuron_repo), n_bins_sac);
for i = neuron_list_
    current_class    = dodr_neuron_best_class_raw_target(i);
    current_class_op = mod(current_class + 4, 8);
    if current_class_op == 0
        current_class_op = 8;
    end
    best_psth_raw_opp(i, :)         = temp_psth_raw(i, current_class_op, :);
    best_psth_raw_sac_opp(i, :)     = temp_psth_raw_sac(i, current_class_op, :);
end
%%
save(fullfile(project_dir, output_database, 'dodr_best_psth.mat'), ...
    'dodr_neuron_best_class_raw_target', 'dodr_neuron_best_class_raw_distractor', ...
    'dodr_neuron_best_class_1s_baseline_target', 'dodr_neuron_best_class_1s_baseline_distractor', ...
    'temp_psth_raw', 'temp_psth_raw_sac', ...
    'best_psth_raw_target', 'best_psth_raw_sac_target', ...
    'best_psth_raw_distractor', 'best_psth_raw_sac_distractor', ...
    'best_psth_raw_opp', 'best_psth_raw_sac_opp', ...
    'temp_psth_1s_baseline', 'temp_psth_1s_baseline_sac', ...
    'best_psth_1s_baseline_target','best_psth_1s_baseline_sac_target', ...
    'best_psth_1s_baseline_opp','best_psth_1s_baseline_sac_opp', ...
    'best_psth_1s_baseline_distractor','best_psth_1s_baseline_sac_distractor' ...
    );
%%
bin_width = 0.02;  % 20 milliseconds bin
bin_edges_cue = -1:bin_width:3; %   [-1, 3] around cue onset, should work for both tasks
n_bins_cue = numel(histcounts([], bin_edges_cue));
%%  neuron class ANOVA
% neuron_analysis_output_sac = struct([]);
for i = neuron_list_
        tic
        n_class_ = numel(dodr_neuron_repo(i).class);
        if isempty(dodr_neuron_repo(i).class)
            continue
        end
        psth_cue_ = {dodr_neuron_repo(i).class.psth_cue_upsampled_for_pev};
        psth_sac_ = {dodr_neuron_repo(i).class.psth_sac_upsampled_for_pev};
        for j = 1:n_bins_cue
            anova_y_cue_ = [];
            anova_label_cue_ = [];
            anova_y_sac_ = [];
            anova_label_sac_ = [];
            for m = 1:n_class_
                if m > 8
                    class_id_ = mod(m, 8);
                else
                    class_id_ = m;
                end
                anova_y_cue_ = [anova_y_cue_, psth_cue_{m}(:, j)'];
                anova_label_cue_ = [anova_label_cue_, class_id_ + zeros(size(psth_cue_{m}(:, j)'))];
                anova_y_sac_ = [anova_y_sac_, psth_sac_{m}(:, j)'];
                anova_label_sac_ = [anova_label_sac_, class_id_ + zeros(size(psth_sac_{m}(:, j)'))];
            end
            [~, out_tb_, out_stats_] = anova1(anova_y_cue_, anova_label_cue_, 'off');
            neuron_analysis_output_cue(i, j).tb = out_tb_;
            neuron_analysis_output_cue(i, j).stats = out_stats_;
            [~, out_tb_, out_stats_] = anova1(anova_y_sac_, anova_label_sac_, 'off');
            neuron_analysis_output_sac(i, j).tb = out_tb_;
            neuron_analysis_output_sac(i, j).stats = out_stats_;
        end
        toc
end
%%  Save neuron_analysis_output_cue
fname_ = 'dodr_neuron_analysis_output';
save(fullfile(project_dir, output_database, fname_), 'neuron_analysis_output_cue');
%%  Compute PEV
neuron_pev_cue = nan(size(neuron_analysis_output_cue));
neuron_pev_sac = nan(size(neuron_analysis_output_sac));
for i = neuron_list_
%     i
        for j = 1:size(neuron_analysis_output_cue, 2)
            if isempty(neuron_analysis_output_cue(i, j).tb)
                        neuron_pev_cue(i, :) = NaN;        
        neuron_pev_cue(i, :) = NaN;
            else
                neuron_pev_cue(i, j) = zw_pev_v2(neuron_analysis_output_cue(i, j).tb);
                if isnan(neuron_pev_cue(i, j))
                    i
                    j
                end
            end
        end
        for j = 1:size(neuron_analysis_output_sac, 2)
            neuron_pev_sac(i, j) = zw_pev(neuron_analysis_output_sac(i, j).tb);
        end
end
%%  Save PEV
fname_ = 'dodr_neuron_pev';
save(fullfile(project_dir, output_database, fname_), 'neuron_pev_cue', 'neuron_pev_sac');
%%
responsive = any(table2array(dodr_neuron_tbl(:, 9:12)), 2);
%%  Quick plot checks
figure; 
hold on
plot(mean(best_psth_raw_target(dodr_neuron_tbl.stage == 1, :), 1))
plot(mean(best_psth_raw_target(dodr_neuron_tbl.stage == 2, :), 1))
%% 
figure; 
hold on
plot(mean(best_psth_1s_baseline_target(dodr_neuron_tbl.stage == 1, :), 1))
plot(mean(best_psth_1s_baseline_target(dodr_neuron_tbl.stage == 2, :), 1))
%% 
figure; 
hold on
plot(mean(best_psth_raw_distractor(dodr_neuron_tbl.stage == 1, :), 1))
plot(mean(best_psth_raw_distractor(dodr_neuron_tbl.stage == 2, :), 1))
%% 
figure; 
hold on
plot(mean(best_psth_1s_baseline_distractor(dodr_neuron_tbl.stage == 1, :), 1))
plot(mean(best_psth_1s_baseline_distractor(dodr_neuron_tbl.stage == 2, :), 1))
%%
figure; 
hold on
plot(mean(best_psth_1s_baseline_target(find(neuron_epoch_delay_tuning_flag(dodr_neuron_tbl.stage == 1)), :), 1))
plot(mean(best_psth_1s_baseline_target(find(neuron_epoch_delay_tuning_flag(dodr_neuron_tbl.stage == 2)), :), 1))
%%
figure; 
hold on
plot(mean(best_psth_1s_baseline_target(find(~neuron_epoch_delay_tuning_flag(dodr_neuron_tbl.stage == 1)), :), 1))
plot(mean(best_psth_1s_baseline_target(find(~neuron_epoch_delay_tuning_flag(dodr_neuron_tbl.stage == 2)), :), 1))
%%
figure
imagesc(squeeze(nanmean(temp_cwt_cue_b(dodr_n_y,:,:), 1)) - 1, [-1,1]);
colormap(clm_)
figure
imagesc(squeeze(nanmean(temp_cwt_cue_b(dodr_n_a,:,:), 1)) - 1, [-1,1]);
colormap(clm_)