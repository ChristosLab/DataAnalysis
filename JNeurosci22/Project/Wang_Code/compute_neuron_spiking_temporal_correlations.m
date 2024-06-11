cor_window_ = 76:150;
% cor_window_ = 51:75;

lfp_mod_psth_corr = nan(size(mapping_mat, 2), 8, size(target_frs, 1));
figure
for i = 1:size(mapping_mat, 2) %    looping through
    if any(mapping_mat(:, i))
        if isempty(neuron_repo(i).class)
            continue
        end
        t_cwt_r_ = cwt_repo(find(mapping_mat(:, i)));
        for k = 1:numel(t_cwt_r_.class)
            upsampled_psth = spline(1/20:1/20:4, neuron_repo(i).class(k).psth_cue, 1/50:1/50:4);
            upsampled_psth = upsampled_psth(:, cor_window_);
            for j = 1:size(lfp_mod_cue, 3) %    Bands
                current_cwt_ = squeeze(sum(t_cwt_r_.class(k).cue_cwt(:, target_frs(j, 1):target_frs(j, 2), cor_window_), 2));
                %   Exclude potential single trials with NaN SPD
                valid_trials = find(~any(isnan(current_cwt_), 2));
                current_cwt_ = current_cwt_(valid_trials, :);
                lfp_mod_psth_corr(i, k, j) = corr(reshape(upsampled_psth(valid_trials, :), numel(upsampled_psth(valid_trials, :)), 1), reshape(current_cwt_, numel(current_cwt_), 1), 'type', 'Spearman');
                lfp_mod_psth_corr(i, k, j)
                subplot(2,1,1)
                yyaxis left
%                 plot(reshape(upsampled_psth(valid_trials, :), numel(upsampled_psth(valid_trials, :)), 1))
                plot(upsampled_psth(valid_trials(1), :));
                yyaxis right
%                 plot(reshape(current_cwt_(valid_trials, :), numel(current_cwt_(valid_trials, :)), 1))
                plot(current_cwt_(valid_trials(1), :));
                xlabel(num2str(j))
                subplot(2,1,2)
                plot(xcorr(reshape(current_cwt_, numel(current_cwt_), 1), reshape(upsampled_psth(valid_trials, :), numel(upsampled_psth(valid_trials, :)), 1), 'unbiased'))
                pause;
%                 current_corr_r_ = zeros(size(upsampled_psth, 1), 1);
%                 for l = 1:size(upsampled_psth, 1)
%                     r_mat = corrcoef(upsampled_psth(l, cor_window_), current_cwt_(l, cor_window_));
%                     current_corr_r_(l) = r_mat(1, 2);
%                 end
%                 lfp_mod_psth_corr(i, k, j) = nanmean(current_corr_r_); %    Mean across trials
            end
        end
    end
%     i
end
%%
best_lfp_mod_psth_corr = zeros(size(lfp_mod_psth_corr,1), size(lfp_mod_psth_corr,3));
for i = 1:size(best_lfp_mod_psth_corr, 1)
    best_lfp_mod_psth_corr(i, :) = lfp_mod_psth_corr(i, neuron_best_class_1s_baseline(i), :);
end
%%
figure; histogram(best_lfp_mod_psth_corr([neuron_tuning_cat{1, [1,3,4], :, 2}], 1), -1:0.05:1, 'Normalization', 'probability')
hold on; histogram(best_lfp_mod_psth_corr([neuron_tuning_cat{2, [1,3,4], :, 2}], 1), -1:0.05:1, 'Normalization', 'probability')

figure; histogram(best_lfp_mod_psth_corr([neuron_tuning_cat{1, [1,3,4], :, 2}], 2), -1:0.05:1, 'Normalization', 'probability')
hold on; histogram(best_lfp_mod_psth_corr([neuron_tuning_cat{2, [1,3,4], :, 2}], 2), -1:0.05:1, 'Normalization', 'probability')

figure; histogram(best_lfp_mod_psth_corr([neuron_tuning_cat{1, [1,3,4], :, 2}], 3), -1:0.05:1, 'Normalization', 'probability')
hold on; histogram(best_lfp_mod_psth_corr([neuron_tuning_cat{2, [1,3,4], :, 2}], 3), -1:0.05:1, 'Normalization', 'probability')
%%
figure; histogram(best_lfp_mod_psth_corr([neuron_tuning_cat{1, [1,3,4], :, 2}], 4), -1:0.05:1, 'Normalization', 'probability', 'FaceColor', 'b', 'FaceAlpha', 1)
hold on; histogram(best_lfp_mod_psth_corr([neuron_tuning_cat{2, [1,3,4], :, 2}], 4), -1:0.05:1, 'Normalization', 'probability', 'FaceColor', 'r', 'FaceAlpha', 1)
%%
hist_plot({best_lfp_mod_psth_corr([neuron_tuning_cat{1, [1,3,4], :, 2}], 4), best_lfp_mod_psth_corr([neuron_tuning_cat{2, [1,3,4], :, 2}], 4)}, ...
    -.7:0.05:.7, {'Adolescent', 'Adult'}, {'b', 'r'}, 'probability');
%%
hist_plot({best_lfp_mod_psth_corr([neuron_tuning_cat{1, [1,3,4], :, 2}], 3), best_lfp_mod_psth_corr([neuron_tuning_cat{2, [1,3,4], :, 2}], 3)}, ...
    -.7:0.05:.7, {'Adolescent', 'Adult'}, {'b', 'r'}, 'probability');
%%
hist_plot({best_lfp_mod_psth_corr([neuron_tuning_cat{1, [1,3,4], :, 2}], 2), best_lfp_mod_psth_corr([neuron_tuning_cat{2, [1,3,4], :, 2}], 2)}, ...
    -.7:0.05:.7, {'Adolescent', 'Adult'}, {'b', 'r'}, 'probability');
%%
hist_plot({best_lfp_mod_psth_corr([neuron_tuning_cat{1, [1,3,4], :, 2}], 1), best_lfp_mod_psth_corr([neuron_tuning_cat{2, [1,3,4], :, 2}], 1)}, ...
    -.7:0.05:.7, {'Adolescent', 'Adult'}, {'b', 'r'}, 'probability');