%%  Whole delay period neuron tuning
fname_ = 'complete_neuron_repo_w_psth.mat';
load(fullfile(project_dir, output_database, fname_));
%%
% n_bins_lfp = (diff(sac_dur)*fs + 1)/down_sample;
neuron_epoch_analysis_output_cue = struct([]);
% lfp_analysis_output_sac = struct([]);
tic
for i = 1:numel(neuron_repo)
    if ~(neuron_tbl.task_id(i) == 1) %  Omit neuron if anti-saccade
        neuron_epoch_analysis_output_cue(i).tb = [];
        neuron_epoch_analysis_output_cue(i).stats = [];
        continue
    end
    if isempty(neuron_repo(i).class) %  Skip neuron if missing data
        neuron_epoch_analysis_output_cue(i).tb = [];
        neuron_epoch_analysis_output_cue(i).stats = [];
        continue
    end
    n_class_ = numel(neuron_repo(i).class);
    cue_neuron_epoch_ = {neuron_repo(i).class.psth_cue};
    %     sac_cwt_ = {cwt_repo(i).class.sac_cwt};
    %         for j = 1:n_bins_lfp
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
        anova_y_cue_ = [anova_y_cue_, sum(cue_neuron_epoch_{m}(:, 31:60), 2)'];
        anova_label_cue_ = [anova_label_cue_, class_id_ + zeros(size(sum(cue_neuron_epoch_{m}(:, 31:60), 2)'))];
        %                 anova_y_sac_ = [anova_y_, sac_cwt_{m}(:, j)'];
        %                 anova_label_sac_ = [anova_label_, class_id_ + zeros(size(sac_cwt_{m}(:, j)'))];
    end
    [~, out_tb_, out_stats_] = anova1(anova_y_cue_, anova_label_cue_, 'off');
    neuron_epoch_analysis_output_cue(i).tb = out_tb_;
    neuron_epoch_analysis_output_cue(i).stats = out_stats_;
    %         end
    toc
end
%%
n_counter_ = 0;
neuron_epoch_delay_tuning_flag = zeros(size(neuron_repo));
for i = 1:size(neuron_epoch_delay_tuning_flag, 2)
    if ~isempty(neuron_epoch_analysis_output_cue(i).tb)
        n_counter_ = n_counter_ + 1;
        neuron_epoch_delay_tuning_flag(i) = neuron_epoch_analysis_output_cue(i).tb{2, 6} < 0.05;
        temp_p_(i) = neuron_epoch_analysis_output_cue(i).tb{2, 6};
    end
end
%%
fname_ = 'neuron_epoch_delay_tuning_flag.mat';
save(fullfile(project_dir, output_database, fname_), 'neuron_epoch_delay_tuning_flag', 'neuron_epoch_analysis_output_cue');
