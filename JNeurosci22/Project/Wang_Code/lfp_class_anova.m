%%  lfp class ANOVA
%%
target_frs = [4,8; 8,16; 16, 32; 32, 64];
n_bands = size(target_frs, 1);
%%
tfr_repo = struct();
tfr_repo.class.cue_tfr = cell(n_bands, 1);
% tfr_repo.class.sac_tfr = cell(8,1);
%%  Compute tfr
for i = 1:numel(cwt_repo)
    if lfp_tbl.task_id(i) == 2 %  Omit neuron if not included in LFP sites
        continue
    end
    tic
    n_class_ = numel(cwt_repo(i).class);
    cue_cwt_ = {cwt_repo(i).class.cue_cwt};
%     sac_cwt_ = {cwt_repo(i).class.sac_cwt};
    for k = 1:n_bands
            for m = 1:n_class_
               temp_tfr_ = sum(cwt_repo(i).class(m).cue_cwt(:, target_frs(k, 1):target_frs(k, 2), :), 2);
               temp_tfr_ = reshape(temp_tfr_, [numel(temp_tfr_)/200, 200]);
                tfr_repo(i).class(m).cue_tfr{k} = temp_tfr_;
%                 tfr_repo(i).class(m).sac_tfr{k} = squeeze(sum(cwt_repo(i).class(m).sac_cwt(:, target_frs(k, 1):target_frs(k, 2), :), 2));
            end
            %             [~, out_tb_, out_stats_] = anova1(anova_y_, anova_label_, 'off');
            %             neuron_analysis_output_cue(i, j).tb = out_tb_;
            %             neuron_analysis_output_cue(i, j).stats = out_stats_;
        toc
    end
end
%%  ANOVA
n_bins_lfp = (diff(sac_dur)*fs + 1)/down_sample;
lfp_analysis_output_cue = struct([]);
% lfp_analysis_output_sac = struct([]);
for i = 1:numel(cwt_repo)
    if lfp_tbl.task_id(i) == 2 %  Omit neuron if not included in LFP sites
        continue
    end
    tic
    n_class_ = numel(cwt_repo(i).class);
    cue_tfr_ = {tfr_repo(i).class.cue_tfr};
%     sac_cwt_ = {cwt_repo(i).class.sac_cwt};
    for k = 1:n_bands
        for j = 1:n_bins_lfp
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
                anova_y_cue_ = [anova_y_cue_, cue_tfr_{m}{k}(:, j)'];
                anova_label_cue_ = [anova_label_cue_, class_id_ + zeros(size(cue_tfr_{m}{k}(:, j)'))];
%                 anova_y_sac_ = [anova_y_, sac_cwt_{m}(:, j)'];
%                 anova_label_sac_ = [anova_label_, class_id_ + zeros(size(sac_cwt_{m}(:, j)'))];
            end
            [~, out_tb_, out_stats_] = anova1(anova_y_cue_, anova_label_cue_, 'off');
            lfp_analysis_output_cue(i, j, k).tb = out_tb_;
            lfp_analysis_output_cue(i, j, k).stats = out_stats_;
        end
    end
    toc
end
%%  Compute PEV
lfp_pev_cue = zeros(size(lfp_analysis_output_cue));
% lfp_pev_sac = zeros(size(lfp_analysis_output_sac));
tic
for i = 1:size(lfp_analysis_output_cue, 1)
    if lfp_tbl.task_id(i) == 2
        continue
    end
    for j = 1:size(lfp_analysis_output_cue, 2)
        for k = 1:size(lfp_analysis_output_cue, 3)
            try
                lfp_pev_cue(i, j, k) = zw_pev(lfp_analysis_output_cue(i, j, k).tb);
            catch
                i
                lfp_pev_cue(i, j, k) = nan;
            end
        end
    end
    %         for j = 1:size(neuron_analysis_output_sac, 2)
    %             lfp_pev_sac(i, j) = zw_pev(neuron_analysis_output_sac(i, j).tb);
    %         end
%     lfp_pev_cue(i, :) = NaN;
    %         lfp_pev_sac(i, :) = NaN;
    toc
end
%%
fname_ = 'lfp_pve_raw.mat';
save(fullfile(project_dir, output_database, fname_), 'lfp_pev_cue', 'lfp_analysis_output_cue');