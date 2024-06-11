fname_ = 'upsample_neuron_repo.mat';
load(fullfile(project_dir, output_database, fname_), 'neuron_repo');
%%
bin_width = 0.05;  % 50 milliseconds bin
step = 0.02; %  20 millisecond step matching that of SPD
t_range_cue = [-1, 3];
t_range_sac = [-3, 1];
bin_edges_cue = -1:step:3;
bin_edges_sac = -3:step:1;
n_bins_cue = numel(histcounts([], bin_edges_cue));
n_bins_sac = numel(histcounts([], bin_edges_sac));
for i = [neuron_tuning_cat_nonsig{:, [1,3,4], :}] % All LFP matched non-sig neurons
    if numel(neuron_repo(i).class) < 1
        continue
    end
    i
    for j = 1:numel(neuron_repo(i).class)
        if isempty(neuron_repo(i).class(j).ntr)
            continue
        end
        psth_cue_ = zeros(numel(neuron_repo(i).class(j).ntr), n_bins_cue);
        psth_sac_ = zeros(numel(neuron_repo(i).class(j).ntr), n_bins_sac);
        for m = 1:numel(neuron_repo(i).class(j).ntr)
            TS_cue_ = neuron_repo(i).class(j).ntr(m).TS - neuron_repo(i).class(j).ntr(m).Cue_onT;
            psth_cue_(m, :) = zw_spike_time_to_psth(TS_cue_, bin_width, step, t_range_cue);
            if isfield(neuron_repo(i).class(j).ntr(m), 'Saccade_onT')
                if ~isempty(neuron_repo(i).class(j).ntr(m).Saccade_onT)
                    TS_sac_ = neuron_repo(i).class(j).ntr(m).TS - neuron_repo(i).class(j).ntr(m).Saccade_onT;
                    psth_sac_(m, :) = zw_spike_time_to_psth(TS_sac_, bin_width, step, t_range_sac);
                else
                    psth_sac_(m, :) = NaN;
                end
            else
                psth_sac_(m, :) = NaN;
            end
        end
        neuron_repo(i).class(j).psth_cue_upsampled = psth_cue_;
        neuron_repo(i).class(j).psth_sac_upsampled = psth_sac_;
    end
end
%%
fname_ = 'upsample_neuron_repo.mat';
save(fullfile(project_dir, output_database, fname_), 'neuron_repo');
%%
load(fullfile(project_dir, output_database, 'best_psth_upsampled.mat'), ...
    'temp_psth', 'best_psth', 'neuron_best_class', ...
    'temp_psth_1s_baseline', 'best_psth_1s_baseline', 'neuron_best_class_1s_baseline', ...
    'temp_psth_raw', 'best_psth_raw', 'neuron_best_class_raw');
%%
%%
for i = [neuron_tuning_cat_nonsig{:, [1,3,4], :}]
    if ~(neuron_tbl.task_id(i) == 1) || ~(neuron_tbl.PFC(i) == 1)
        temp_psth_raw(i, :) = nan;
        temp_psth(i, :) = nan;
        continue
    end
    for j = 1:8
        temp_psth(i, j, :) = nanmean(neuron_repo(i).class(j).psth_cue_upsampled - mean(neuron_repo(i).class(j).psth_cue_upsampled(:, 26:50), 2), 1);
        temp_psth_1s_baseline(i, j, :) = nanmean(neuron_repo(i).class(j).psth_cue_upsampled - mean(neuron_repo(i).class(j).psth_cue_upsampled(:, 1:50), 2), 1);
        temp_psth_raw(i, j, :) = nanmean(neuron_repo(i).class(j).psth_cue_upsampled, 1);
    end
end
%%
for i = 1:size(temp_psth_raw, 1)
    [y_, i_] = max(sum(temp_psth_raw(i, :, 51:150), 3));
    neuron_best_class_raw(i) = i_;
end
%%
for i = 1:size(best_psth_raw, 1)
    best_psth_raw(i, :) = temp_psth_raw(i, neuron_best_class_raw(i), :);
end
%%
for i = 1:size(temp_psth, 1)
    [y_, i_] = max(sum(temp_psth(i, :,51:150), 3));
    neuron_best_class(i) = i_;
end
%%
for i = 1:size(best_psth, 1)
    best_psth(i, :) = temp_psth(i, neuron_best_class(i), :);
end
%%
for i = 1:size(temp_psth, 1)
    [y_, i_] = max(sum(temp_psth_1s_baseline(i, :, 51:150), 3));
    neuron_best_class_1s_baseline(i) = i_;
end
%%
for i = 1:size(best_psth, 1)
    best_psth_1s_baseline(i, :) = temp_psth(i, neuron_best_class_1s_baseline(i), :);
end
%%
save(fullfile(project_dir, output_database, 'best_psth_upsampled.mat'), ...
    'temp_psth', 'best_psth', 'neuron_best_class', ...
    'temp_psth_1s_baseline', 'best_psth_1s_baseline', 'neuron_best_class_1s_baseline', ...
    'temp_psth_raw', 'best_psth_raw', 'neuron_best_class_raw');
%%
bin_width = 0.02;  % 20 milliseconds bin
bin_edges_cue = -1:bin_width:3; %   [-1, 3] around cue onset, should work for both tasks
n_bins_cue = numel(histcounts([], bin_edges_cue));
%%  Load neuron_analysis_output_cue (sig)
fname_ = 'neuron_analysis_output_cue_ps_upsampled';
load(fullfile(project_dir, output_database, fname_), 'neuron_analysis_output_cue');
%%  neuron class ANOVA
% neuron_analysis_output_sac = struct([]);
for i = [neuron_tuning_cat_nonsig{:, [1,3,4], :}]
    if ~(neuron_tbl.task_id(i) == 1) %  Omit neuron if anti-saccade
        continue
    else
        tic
        n_class_ = numel(neuron_repo(i).class);
        if isempty(neuron_repo(i).class)
            continue
        end
        psth_cue_ = {neuron_repo(i).class.psth_cue_upsampled};
%         psth_sac_ = {neuron_repo(i).class.psth_sac};
        for j = 1:n_bins_cue
            anova_y_cue_ = [];
            anova_label_cue_ = [];
%             anova_y_sac_ = [];
%             anova_label_sac_ = [];
            for m = 1:n_class_
                if m > 8
                    class_id_ = mod(m, 8);
                else
                    class_id_ = m;
                end
                anova_y_cue_ = [anova_y_cue_, psth_cue_{m}(:, j)'];
                anova_label_cue_ = [anova_label_cue_, class_id_ + zeros(size(psth_cue_{m}(:, j)'))];
%                 anova_y_sac_ = [anova_y_, psth_sac_{m}(:, j)'];
%                 anova_label_sac_ = [anova_label_, class_id_ + zeros(size(psth_sac_{m}(:, j)'))];
            end
            [~, out_tb_, out_stats_] = anova1(anova_y_cue_, anova_label_cue_, 'off');
            neuron_analysis_output_cue(i, j).tb = out_tb_;
            neuron_analysis_output_cue(i, j).stats = out_stats_;
%             [~, out_tb_, out_stats_] = anova1(anova_y_sac_, anova_label_sac_, 'off');
%             neuron_analysis_output_sac(i, j).tb = out_tb_;
%             neuron_analysis_output_sac(i, j).stats = out_stats_;
        end
        toc
    end
end
%%  Save neuron_analysis_output_cue
fname_ = 'neuron_analysis_output_cue_ps_upsampled';
save(fullfile(project_dir, output_database, fname_), 'neuron_analysis_output_cue');
%%  Load PEV (sig)
fname_ = 'neuron_pev_cue_ps_upsampled';
load(fullfile(project_dir, output_database, fname_), 'neuron_pev_cue');
%%  Compute PEV
% neuron_pev_sac = zeros(size(neuron_analysis_output_sac));
for i = [neuron_tuning_cat_nonsig{:, [1,3,4], :}]
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
%         for j = 1:size(neuron_analysis_output_sac, 2)
%             neuron_pev_sac(i, j) = zw_pev(neuron_analysis_output_sac(i, j).tb);
%         end
    else
        neuron_pev_cue(i, :) = NaN;        
%         neuron_pev_sac(i, :) = NaN;
    end
end
%%  Save PEV
fname_ = 'neuron_pev_cue_ps_upsampled';
save(fullfile(project_dir, output_database, fname_), 'neuron_pev_cue');
%%
%   Output data files updated 3/4/2021 -ZW