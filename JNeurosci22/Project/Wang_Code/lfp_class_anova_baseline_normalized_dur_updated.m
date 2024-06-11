%%  lfp class ANOVA
%%
zw_setpath;
%%
fname_ = 'complete_cwt_repo_3_2_2-21';
load(fullfile(project_dir, output_database, fname_));
%%
target_frs = [4,8; 8,16; 16, 32; 32, 64];
n_bands = size(target_frs, 1);
baseline_bin = [26, 50];
%%  CWT parameters
fs = 500;
cue_dur = [-1 + 1/fs, 4]; % Signal length 
normalizer = 4;
kernel_flag = 3;
down_sample = 10;
f_range = 2:2:128;
baseline_bin = [26, 50];
avg_method = 1; %   0 - arithmetic; 1 - geometric
%%
baseline_normalized_tfr_repo = struct();
baseline_normalized_tfr_repo.class.cue_tfr = cell(n_bands, 1);
% baseline_normalized_tfr_repo.class.sac_tfr = cell(8,1);
%%  Compute tfr
for i = [t{:}]
    tic
    n_class_ = numel(cwt_repo(i).class);
    cue_cwt_ = {cwt_repo(i).class.cue_cwt};
%     sac_cwt_ = {cwt_repo(i).class.sac_cwt};
    for k = 1:n_bands
            for m = 1:n_class_
               temp_tfr_ = sum(cwt_repo(i).class(m).cue_cwt(:, target_frs(k, 1):target_frs(k, 2), :), 2)...
                   ./nanmean(sum(cwt_repo(i).class(m).cue_cwt(:, target_frs(k, 1):target_frs(k, 2), baseline_bin(1):baseline_bin(2)), 2), 3);
               temp_tfr_ = reshape(temp_tfr_, [numel(temp_tfr_)/250, 250]);
                baseline_normalized_tfr_repo(i).class(m).cue_tfr{k} = temp_tfr_;
%                 baseline_normalized_tfr_repo(i).class(m).sac_tfr{k} = squeeze(sum(cwt_repo(i).class(m).sac_cwt(:, target_frs(k, 1):target_frs(k, 2), :), 2));
            end
            %             [~, out_tb_, out_stats_] = anova1(anova_y_, anova_label_, 'off');
            %             neuron_analysis_output_cue(i, j).tb = out_tb_;
            %             neuron_analysis_output_cue(i, j).stats = out_stats_;
    end
    toc
end
%%
fname_ = 'baseline_normalized_tfr_repo.mat';
save(fullfile(project_dir, output_database, fname_), 'baseline_normalized_tfr_repo');
%%  ANOVA
n_bins_lfp = (diff(cue_dur)*fs + 1)/down_sample;
lfp_analysis_output_cue = struct([]);
% lfp_analysis_output_sac = struct([]);
for i = [t{:}]
    tic
    n_class_ = numel(baseline_normalized_tfr_repo(i).class);
    cue_tfr_ = {baseline_normalized_tfr_repo(i).class.cue_tfr};
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
for i = [t{:}]
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
fname_ = 'lfp_pve_baseline_normalized_dur_updated.mat';
save(fullfile(project_dir, output_database, fname_), 'lfp_pev_cue');
