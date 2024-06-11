%%  Whole delay period LFP modulation
% n_bins_lfp = (diff(sac_dur)*fs + 1)/down_sample;
tfr_modulation_output_cue = struct([]);
% lfp_analysis_output_sac = struct([]);
for i = 1:numel(cwt_repo)
    if ~ismember(i,[site_tuning_cat{:}]) %  Omit neuron if not included in LFP sites
        continue
    end
    tic
    n_class_ = numel(cwt_repo(i).class);
    cue_tfr_ = {baseline_normalized_tfr_repo(i).class.cue_tfr};
%     sac_cwt_ = {cwt_repo(i).class.sac_cwt};
    for k = 1:n_bands
%         for j = 1:n_bins_lfp
            pow_baseline_ = [];
            pow_cue_ = [];
            pow_delay_ = [];
            p_ = [];
            mean_ = [];
            for m = 1:n_class_
                if m > 8
                    class_id_ = mod(m, 8);
                else
                    class_id_ = m;
                end
                
                pow_baseline_ = nanmean(cue_tfr_{m}{k}(:, 26:50), 2)';
                pow_cue_ = nanmean(cue_tfr_{m}{k}(:, 51:75), 2)';
                pow_delay_ = nanmean(cue_tfr_{m}{k}(:, 76:150), 2)';
                pow_total_ = nanmean(cue_tfr_{m}{k}(:, 51:150), 2)';
                p__ = zeros(3, 1);
                [~, p__(1)] = ttest(pow_baseline_, pow_cue_);
                [~, p__(2)] = ttest(pow_baseline_, pow_delay_); 
                [~, p__(3)] = ttest(pow_baseline_, pow_total_);
                p_ = [p_, p__];
%                 p_ = [p_, [signrank(pow_baseline_, pow_cue_); signrank(pow_baseline_, pow_delay_); signrank(pow_baseline_, pow_total_)]];
                mean_ = [mean_, [nanmean(pow_cue_); nanmean(pow_delay_); nanmean(pow_total_)]];
            end
            tfr_modulation_output_cue(i, k).p = p_;
            tfr_modulation_output_cue(i, k).mean = mean_;
%         end
    end
    toc
end
%%
tfr_cue_modulation_flag = zeros([size(tfr_modulation_output_cue), 2]);
tfr_delay_modulation_flag = zeros([size(tfr_modulation_output_cue), 2]);
tfr_total_modulation_flag = zeros([size(tfr_modulation_output_cue), 2]);

for i = 1:size(tfr_delay_modulation_flag, 1)
    for j = 1:size(tfr_delay_modulation_flag, 2)
        if ~isempty(tfr_modulation_output_cue(i,j).p)
        tfr_cue_modulation_flag(i, j, 1) = any((tfr_modulation_output_cue(i,j).p(1, :) < 0.05/8).*(tfr_modulation_output_cue(i,j).mean(1, :) > 1));
        tfr_delay_modulation_flag(i, j, 1) = any((tfr_modulation_output_cue(i,j).p(2, :) < 0.05/8).*(tfr_modulation_output_cue(i,j).mean(2, :) > 1));
        tfr_total_modulation_flag(i, j, 1) = any((tfr_modulation_output_cue(i,j).p(3, :) < 0.05/8).*(tfr_modulation_output_cue(i,j).mean(3, :) > 1));
        tfr_cue_modulation_flag(i, j, 2) = any((tfr_modulation_output_cue(i,j).p(1, :) < 0.05/8).*(tfr_modulation_output_cue(i,j).mean(1, :) < 1));
        tfr_delay_modulation_flag(i, j, 2) = any((tfr_modulation_output_cue(i,j).p(2, :) < 0.05/8).*(tfr_modulation_output_cue(i,j).mean(2, :) < 1));
        tfr_total_modulation_flag(i, j, 2) = any((tfr_modulation_output_cue(i,j).p(3, :) < 0.05/8).*(tfr_modulation_output_cue(i,j).mean(3, :) < 1));
% 
%         tfr_cue_modulation_flag(i, j, 1) = any((tfr_modulation_output_cue(i,j).p(1, :) < 0.05).*(tfr_modulation_output_cue(i,j).mean(1, :) > 1));
%         tfr_delay_modulation_flag(i, j, 1) = any((tfr_modulation_output_cue(i,j).p(2, :) < 0.05).*(tfr_modulation_output_cue(i,j).mean(2, :) > 1));
%         tfr_total_modulation_flag(i, j, 1) = any((tfr_modulation_output_cue(i,j).p(3, :) < 0.05).*(tfr_modulation_output_cue(i,j).mean(3, :) > 1));
%         tfr_cue_modulation_flag(i, j, 2) = any((tfr_modulation_output_cue(i,j).p(1, :) < 0.05).*(tfr_modulation_output_cue(i,j).mean(1, :) < 1));
%         tfr_delay_modulation_flag(i, j, 2) = any((tfr_modulation_output_cue(i,j).p(2, :) < 0.05).*(tfr_modulation_output_cue(i,j).mean(2, :) < 1));
%         tfr_total_modulation_flag(i, j, 2) = any((tfr_modulation_output_cue(i,j).p(3, :) < 0.05).*(tfr_modulation_output_cue(i,j).mean(3, :) < 1));
        end
    end
end
%%
n_young_gamma_tune_delay = sum(tfr_delay_modulation_flag([t{1, [1,3,4]}], 3, 1));
n_adult_gamma_tune_delay = sum(tfr_delay_modulation_flag([t{2, [1,3,4]}], 3, 1));

[p_, h_] = fishertest([n_young_gamma_tune_delay, n_adult_gamma_tune_delay; numel([t{1, [1,3,4]}]) - n_young_gamma_tune_delay, numel([t{2, [1,3,4]}]) - n_adult_gamma_tune_delay])

%%  Tuning numbers
n_young_gamma_tune_delay = sum(tfr_delay_modulation_flag([t{1, [1,3,4]}], 3, 1));
n_adult_gamma_tune_delay = sum(tfr_delay_modulation_flag([t{2, [1,3,4]}], 3, 1));

[p_, h_] = fishertest([n_young_gamma_tune_delay, n_adult_gamma_tune_delay; numel([t{1, [1,3,4]}]) - n_young_gamma_tune_delay, numel([t{2, [1,3,4]}]) - n_adult_gamma_tune_delay])

n_young_gamma_tune_cue = sum(tfr_cue_modulation_flag([t{1, [1,3,4]}], 3, 1));
n_adult_gamma_tune_cue = sum(tfr_cue_modulation_flag([t{2, [1,3,4]}], 3, 1));
[p_, h_] = fishertest([n_young_gamma_tune_cue, n_adult_gamma_tune_cue; numel([t{1, [1,3,4]}]) - n_young_gamma_tune_cue, numel([t{2, [1,3,4]}]) - n_adult_gamma_tune_cue])

n_young_gamma_tune_total = sum(tfr_total_modulation_flag([t{1, [1,3,4]}], 3, 1));
n_adult_gamma_tune_total = sum(tfr_total_modulation_flag([t{2, [1,3,4]}], 3, 1));
[p_, h_] = fishertest([n_young_gamma_tune_total, n_adult_gamma_tune_total; numel([t{1, [1,3,4]}]) - n_young_gamma_tune_total, numel([t{2, [1,3,4]}]) - n_adult_gamma_tune_total])
%%
% fname_ = 'tfr_modulation.mat';
% save(fullfile(project_dir, output_database, fname_), 'tfr_total_modulation_flag', 'tfr_cue_modulation_flag', 'tfr_delay_modulation_flag', 'tfr_modulation_output_cue');