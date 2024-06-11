%%  Whole delay period LFP modulation
% n_bins_lfp = (diff(sac_dur)*fs + 1)/down_sample;
tfr_modulation_all_class = struct([]);
% lfp_analysis_output_sac = struct([]);
for i = 1:numel(cwt_repo)
    if ~ismember(i,[site_tuning_cat{:}]) %  Omit neuron if not included in LFP sites
        for k = 1:n_bands
            tfr_modulation_all_class(i, k).p = [];
            tfr_modulation_all_class(i, k).mean = [];
        end
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
            pow_total_ = [];
            p_ = [];
            mean_ = [];
            for m = 1:n_class_
                if m > 8
                    class_id_ = mod(m, 8);
                else
                    class_id_ = m;
                end
                pow_baseline_ = [pow_baseline_, nanmean(cue_tfr_{m}{k}(:, 26:50), 2)'];
                pow_cue_ = [pow_cue_, nanmean(cue_tfr_{m}{k}(:, 51:75), 2)'];
                pow_delay_ = [pow_delay_, nanmean(cue_tfr_{m}{k}(:, 76:150), 2)'];
                pow_total_ = [pow_total_, nanmean(cue_tfr_{m}{k}(:, 51:150), 2)'];
                mean_ = [mean_, [nanmean(pow_cue_); nanmean(pow_delay_); nanmean(pow_total_)]];
            end
            p__ = zeros(3, 1);
            [~, p__(1)] = ttest(pow_baseline_, pow_cue_);
            [~, p__(2)] = ttest(pow_baseline_, pow_delay_);
            [~, p__(3)] = ttest(pow_baseline_, pow_total_);
            tfr_modulation_all_class(i, k).p = p__;
            tfr_modulation_all_class(i, k).mean = mean_;
            %         end
    end
    toc
end
%%
tfr_cue_modulation_flag_all_class = zeros([size(tfr_modulation_all_class), 2]);
tfr_delay_modulation_flag_all_class = zeros([size(tfr_modulation_all_class), 2]);
tfr_total_modulation_flag_all_class = zeros([size(tfr_modulation_all_class), 2]);

for i = 1:size(tfr_total_modulation_flag_all_class, 1)
    for j = 1:size(tfr_total_modulation_flag_all_class, 2)
        if ~isempty(tfr_modulation_all_class(i,j).p)
        tfr_cue_modulation_flag_all_class(i, j, 1) = any((tfr_modulation_all_class(i,j).p(1, :) < 0.05).*(tfr_modulation_all_class(i,j).mean(1, :) > 1));
        tfr_delay_modulation_flag_all_class(i, j, 1) = any((tfr_modulation_all_class(i,j).p(2, :) < 0.05).*(tfr_modulation_all_class(i,j).mean(2, :) > 1));
        tfr_total_modulation_flag_all_class(i, j, 1) = any((tfr_modulation_all_class(i,j).p(3, :) < 0.05).*(tfr_modulation_all_class(i,j).mean(3, :) > 1));
        tfr_cue_modulation_flag_all_class(i, j, 2) = any((tfr_modulation_all_class(i,j).p(1, :) < 0.05).*(tfr_modulation_all_class(i,j).mean(1, :) < 1));
        tfr_delay_modulation_flag_all_class(i, j, 2) = any((tfr_modulation_all_class(i,j).p(2, :) < 0.05).*(tfr_modulation_all_class(i,j).mean(2, :) < 1));
        tfr_total_modulation_flag_all_class(i, j, 2) = any((tfr_modulation_all_class(i,j).p(3, :) < 0.05).*(tfr_modulation_all_class(i,j).mean(3, :) < 1));
% 
        end
    end
end
%%
%%  Tuning numbers
n_young_gamma_tune_delay = sum(tfr_delay_modulation_flag_all_class([t{1, [1,3,4]}], 3, 1));
n_adult_gamma_tune_delay = sum(tfr_delay_modulation_flag_all_class([t{2, [1,3,4]}], 3, 1));
[p_, h_] = fishertest([n_young_gamma_tune_delay, n_adult_gamma_tune_delay; numel([t{1, [1,3,4]}]) - n_young_gamma_tune_delay, numel([t{2, [1,3,4]}]) - n_adult_gamma_tune_delay])
n_young_high_gamma_tune_delay = sum(tfr_delay_modulation_flag_all_class([t{1, [1,3,4]}], 4, 1));
n_adult_high_gamma_tune_delay = sum(tfr_delay_modulation_flag_all_class([t{2, [1,3,4]}], 4, 1));
[p_, h_] = fishertest([n_young_high_gamma_tune_delay, n_adult_high_gamma_tune_delay; numel([t{1, [1,3,4]}]) - n_young_high_gamma_tune_delay, numel([t{2, [1,3,4]}]) - n_adult_high_gamma_tune_delay])
%%
%%  Modulated delay power compare
site_mod_cat = cell(2, 4, 2);
for i_stage = 1:2
    for i_monk = 1:4
        %               stage,   monkey, delay gamma mod
        site_mod_cat{i_stage, i_monk, 1} = intersect(t{i_stage, i_monk}, find(any(mapping_mat.*~tfr_delay_modulation_flag_all_class(:, 3, 1), 2))');
        site_mod_cat{i_stage, i_monk, 2} = intersect(t{i_stage, i_monk}, find(any(mapping_mat.*tfr_delay_modulation_flag_all_class(:, 3, 1), 2))');
        site_mod_cat{i_stage, i_monk, 3} = intersect(t{i_stage, i_monk}, find(any(mapping_mat.*~tfr_delay_modulation_flag_all_class(:, 4, 1), 2))');
        site_mod_cat{i_stage, i_monk, 4} = intersect(t{i_stage, i_monk}, find(any(mapping_mat.*tfr_delay_modulation_flag_all_class(:, 4, 1), 2))');        
    end
end
%%
clc
band_i = 3;
set_gamma_power_mod_adolescent = mean(lfp_mod_cue([site_mod_cat{1, [1, 3, 4], 2}], 76:150, band_i), 2);
set_gamma_power_mod_adult = mean(lfp_mod_cue([site_mod_cat{2, [1, 3, 4], 2}], 76:150, band_i), 2);
[h_, p_, ~, stats] = ttest2(set_gamma_power_mod_adolescent, set_gamma_power_mod_adult)
% 
band_i = 4;
set_gamma_power_mod_adolescent = mean(lfp_mod_cue([site_mod_cat{1, [1, 3, 4], 4}], 76:150, band_i), 2);
set_gamma_power_mod_adult = mean(lfp_mod_cue([site_mod_cat{2, [1, 3, 4], 4}], 76:150, band_i), 2);
[h_, p_, ~, stats] = ttest2(set_gamma_power_mod_adolescent, set_gamma_power_mod_adult)