%%  Whole delay period LFP tuning
% n_bins_lfp = (diff(sac_dur)*fs + 1)/down_sample;
tfr_analysis_output_cue = struct([]);
% lfp_analysis_output_sac = struct([]);
for i = 1:numel(cwt_repo)
    if lfp_tbl.task_id(i) == 2 %  Omit neuron if not included in LFP sites
        continue
    end
    tic
    n_class_ = numel(cwt_repo(i).class);
    cue_tfr_ = {baseline_normalized_tfr_repo(i).class.cue_tfr};
%     sac_cwt_ = {cwt_repo(i).class.sac_cwt};
    for k = 1:n_bands
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
                anova_y_cue_ = [anova_y_cue_, sum(cue_tfr_{m}{k}(:, 76:150), 2)'];
                anova_label_cue_ = [anova_label_cue_, class_id_ + zeros(size(sum(cue_tfr_{m}{k}(:, 76:150), 2)'))];
%                 anova_y_sac_ = [anova_y_, sac_cwt_{m}(:, j)'];
%                 anova_label_sac_ = [anova_label_, class_id_ + zeros(size(sac_cwt_{m}(:, j)'))];
            end
            [~, out_tb_, out_stats_] = anova1(anova_y_cue_, anova_label_cue_, 'off');
            tfr_analysis_output_cue(i, k).tb = out_tb_;
            tfr_analysis_output_cue(i, k).stats = out_stats_;
%         end
    end
    toc
end
%%
tfr_delay_tuning_flag = zeros(size(tfr_analysis_output_cue));
for i = 1:size(tfr_delay_tuning_flag, 1)
    for j = 1:size(tfr_delay_tuning_flag, 2)
        if ~isempty(tfr_analysis_output_cue(i,j).tb)
        tfr_delay_tuning_flag(i, j) = tfr_analysis_output_cue(i,j).tb{2, 6} < 0.05;
        end
    end
end
%%
fname_ = 'tfr_delay_tuning_flag.mat';
save(fullfile(project_dir, output_database, fname_), 'tfr_delay_tuning_flag');
