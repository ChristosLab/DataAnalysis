fname_ = 'upsample_neuron_repo.mat';
load(fullfile(project_dir, output_database, fname_), 'neuron_repo');
load(fullfile(project_dir, output_database, 'neuron_site_categories_nonsig.mat'));
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
neuron_list_ = [neuron_tuning_cat_nonsig{:, [1,3,4], :}, neuron_tuning_cat{:, [1,3,4], :}];
for i = neuron_list_ % All LFP matched sig and non-sig neurons
    if numel(neuron_repo(i).class) < 1
        continue
    end
    i
    for j = 1:numel(neuron_repo(i).class)
        if isempty(neuron_repo(i).class(j).ntr)
            continue
        end
        psth_cue_         = zeros(numel(neuron_repo(i).class(j).ntr), n_bins_cue);
        psth_sac_         = zeros(numel(neuron_repo(i).class(j).ntr), n_bins_sac);
        psth_cue_for_pev_ = zeros(numel(neuron_repo(i).class(j).ntr), n_bins_cue);
        psth_cue_for_sac_ = zeros(numel(neuron_repo(i).class(j).ntr), n_bins_sac);
        for m = 1:numel(neuron_repo(i).class(j).ntr)
            TS_cue_ = neuron_repo(i).class(j).ntr(m).TS - neuron_repo(i).class(j).ntr(m).Cue_onT;
            psth_cue_(m, :)         = zw_spike_time_to_psth(TS_cue_, bin_width, step, t_range_cue);
            psth_cue_for_pev_(m, :) = zw_spike_time_to_psth(TS_cue_, bin_width_for_pev, step, t_range_cue);
            if isfield(neuron_repo(i).class(j).ntr(m), 'Saccade_onT')
                if ~isempty(neuron_repo(i).class(j).ntr(m).Saccade_onT)
                    TS_sac_ = neuron_repo(i).class(j).ntr(m).TS - neuron_repo(i).class(j).ntr(m).Saccade_onT;
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
        neuron_repo(i).class(j).psth_cue_upsampled          = psth_cue_;
        neuron_repo(i).class(j).psth_cue_upsampled_for_pev  = psth_cue_for_pev_;
        neuron_repo(i).class(j).psth_sac_upsampled          = psth_sac_;
        neuron_repo(i).class(j).psth_sac_upsampled_for_pev  = psth_sac_for_pev_;
    end
end
%%
fname_ = 'upsample_neuron_repo_3_11_2021.mat';
save(fullfile(project_dir, output_database, fname_), 'neuron_repo');
%%
load(fullfile(project_dir, output_database, 'best_psth_upsampled.mat'), ...
    'temp_psth', 'best_psth', 'neuron_best_class', ...
    'temp_psth_1s_baseline', 'best_psth_1s_baseline', 'neuron_best_class_1s_baseline', ...
    'temp_psth_raw', 'best_psth_raw', 'neuron_best_class_raw');
%%
%%
for i = neuron_list_
    if ~(neuron_tbl.task_id(i) == 1) || ~(neuron_tbl.PFC(i) == 1)
        temp_psth_raw(i, :) = nan;
        temp_psth(i, :) = nan;
        continue
    end
    for j = 1:8
        temp_psth(i, j, :)                 = nanmean(neuron_repo(i).class(j).psth_cue_upsampled - mean(neuron_repo(i).class(j).psth_cue_upsampled(:, 26:50), 2), 1);
        temp_psth_1s_baseline(i, j, :)     = nanmean(neuron_repo(i).class(j).psth_cue_upsampled - mean(neuron_repo(i).class(j).psth_cue_upsampled(:, 1:50), 2), 1);
        temp_psth_raw(i, j, :)             = nanmean(neuron_repo(i).class(j).psth_cue_upsampled, 1);
        temp_psth_sac(i, j, :)             = nanmean(neuron_repo(i).class(j).psth_sac_upsampled - mean(neuron_repo(i).class(j).psth_cue_upsampled(:, 26:50), 2), 1);
        temp_psth_1s_baseline_sac(i, j, :) = nanmean(neuron_repo(i).class(j).psth_sac_upsampled - mean(neuron_repo(i).class(j).psth_cue_upsampled(:, 1:50), 2), 1);
        temp_psth_raw_sac(i, j, :)         = nanmean(neuron_repo(i).class(j).psth_sac_upsampled, 1);
    end
end
%%
for i = 1:size(temp_psth_raw, 1)
    [y_, i_] = max(sum(temp_psth_raw(i, :, 51:150), 3));
    neuron_best_class_raw(i) = i_;
    [y_, i_] = min(sum(temp_psth_raw(i, :, 51:150), 3));
    neuron_worst_class_raw(i) = i_;
end
%%
best_psth_raw     = nan(numel(neuron_repo), n_bins_cue);
worst_psth_raw     = nan(numel(neuron_repo), n_bins_cue);
best_psth_raw_sac = nan(numel(neuron_repo), n_bins_sac);
for i = neuron_list_
    best_psth_raw(i, :)     = temp_psth_raw(i, neuron_best_class_raw(i), :);
    best_psth_raw_sac(i, :) = temp_psth_raw_sac(i, neuron_best_class_raw(i), :);
    worst_psth_raw(i, :)     = temp_psth_raw(i, neuron_worst_class_raw(i), :);
end
%%
for i = 1:size(temp_psth, 1)
    [y_, i_] = max(sum(temp_psth(i, :,51:150), 3));
    neuron_best_class(i) = i_;
end
%%
best_psth     = nan(numel(neuron_repo), n_bins_cue);
best_psth_sac = nan(numel(neuron_repo), n_bins_sac);
for i = neuron_list_
    best_psth(i, :)     = temp_psth(i, neuron_best_class(i), :);
    best_psth_sac(i, :) = temp_psth_sac(i, neuron_best_class(i), :);
end
%%
for i = 1:size(temp_psth, 1)
    [y_, i_] = max(sum(temp_psth_1s_baseline(i, :, 51:150), 3));
    neuron_best_class_1s_baseline(i) = i_;
    [y_, i_] = min(sum(temp_psth_1s_baseline(i, :, 51:150), 3));
    neuron_worst_class_1s_baseline(i) = i_;

end
%%
best_psth_1s_baseline     = nan(numel(neuron_repo), n_bins_cue);
best_psth_1s_baseline_sac = nan(numel(neuron_repo), n_bins_sac);
worst_psth_1s_baseline     = nan(numel(neuron_repo), n_bins_cue);
for i = neuron_list_
    best_psth_1s_baseline(i, :)     = temp_psth(i, neuron_best_class_1s_baseline(i), :);
    best_psth_1s_baseline_sac(i, :) = temp_psth_sac(i, neuron_best_class_1s_baseline(i), :);
    worst_psth_1s_baseline(i, :)     = temp_psth(i, neuron_worst_class_1s_baseline(i), :);
end
%%
save(fullfile(project_dir, output_database, 'best_psth_upsampled_3_11_2021.mat'), ...
    'temp_psth', 'temp_psth_sac', 'best_psth', 'best_psth_sac', 'neuron_best_class', ...
    'temp_psth_1s_baseline', 'temp_psth_1s_baseline_sac', 'best_psth_1s_baseline','best_psth_1s_baseline_sac', 'neuron_best_class_1s_baseline', ...
    'temp_psth_raw', 'temp_psth_raw_sac', 'best_psth_raw','best_psth_raw_sac', 'neuron_best_class_raw');
%%
save(fullfile(project_dir, output_database, 'worst_psth_upsampled_12_20_2021.mat'), ...
    'worst_psth_1s_baseline', 'neuron_worst_class_1s_baseline', ...
    'worst_psth_raw', 'neuron_worst_class_raw');
%%
bin_width = 0.02;  % 20 milliseconds bin
bin_edges_cue = -1:bin_width:3; %   [-1, 3] around cue onset, should work for both tasks
n_bins_cue = numel(histcounts([], bin_edges_cue));
%%  Load neuron_analysis_output_cue (sig)
fname_ = 'neuron_analysis_output_cue_ps_upsampled';
load(fullfile(project_dir, output_database, fname_), 'neuron_analysis_output_cue');
%%  neuron class ANOVA
% neuron_analysis_output_sac = struct([]);
for i = neuron_list_
        tic
        n_class_ = numel(neuron_repo(i).class);
        if isempty(neuron_repo(i).class)
            continue
        end
        psth_cue_ = {neuron_repo(i).class.psth_cue_upsampled_for_pev};
        psth_sac_ = {neuron_repo(i).class.psth_sac_upsampled_for_pev};
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
fname_ = 'neuron_analysis_output_cue_ps_upsampled_3_11_2021';
save(fullfile(project_dir, output_database, fname_), 'neuron_analysis_output_cue');
%%  Load PEV (sig)
fname_ = 'neuron_pev_cue_ps_upsampled';
load(fullfile(project_dir, output_database, fname_), 'neuron_pev_cue');
%%  Compute PEV
neuron_pev_cue = nan(size(neuron_analysis_output_cue));
neuron_pev_sac = nan(size(neuron_analysis_output_sac));
for i = neuron_list_
%     i
    if neuron_tbl.task_id(i) == 1
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
    else
        neuron_pev_cue(i, :) = NaN;        
        neuron_pev_sac(i, :) = NaN;
    end
end
%%  Save PEV
fname_ = 'neuron_pev_ps_upsampled_3_11_2021';
save(fullfile(project_dir, output_database, fname_), 'neuron_pev_cue', 'neuron_pev_sac');
%%
%   Output data files updated 3/4/2021 -ZW