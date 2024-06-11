fname_ = 'upsample_neuron_repo_3_11_2021.mat';
load(fullfile(project_dir, output_database, fname_), 'neuron_repo');
fname_ = 'complete_cwt_repo_3_2_2-21.mat';
load(fullfile(project_dir, output_database, fname_), 'cwt_repo');
%%
cor_window_ = 76:150;
lfp_mod_psth_corr = nan(size(mapping_mat, 2), 8, 64);
%%

for i = 1:size(mapping_mat, 2) %    looping through
    if any(mapping_mat(:, i))
        if isempty(neuron_repo(i).class)
            continue
        end
        t_cwt_r_ = cwt_repo(find(mapping_mat(:, i)));
        for k = 1:numel(t_cwt_r_.class)
            upsampled_psth = neuron_repo(i).class(k).psth_cue_upsampled_for_pev; %  50-ms boxcar kernel
            upsampled_psth = upsampled_psth(:, cor_window_);
            for j = 1:64 %    Bands
                current_cwt_ = squeeze(t_cwt_r_.class(k).cue_cwt(:, j, cor_window_));
                %   Exclude potential single trials with NaN SPD
                valid_trials = find(~any(isnan(current_cwt_), 2));
                current_cwt_ = current_cwt_(valid_trials, :);
%                 lfp_mod_psth_corr(i, k, j) = corr(reshape(upsampled_psth(valid_trials, :), numel(upsampled_psth(valid_trials, :)), 1), reshape(current_cwt_, numel(current_cwt_), 1), 'type', 'Spearman');
                lfp_mod_psth_corr(i, k, j) = corr(reshape(upsampled_psth(valid_trials, :), numel(upsampled_psth(valid_trials, :)), 1), reshape(current_cwt_, numel(current_cwt_), 1), 'type', 'Pearson');
                if sum(reshape(upsampled_psth(valid_trials, :), numel(upsampled_psth(valid_trials, :)), 1)) == 0
                    lfp_mod_psth_corr(i, k, j) = 0;
                end 
                if isnan(lfp_mod_psth_corr(i, k, j))
                    i
                    j
                    k                    
                end
%                 subplot(2,1,1)
%                 yyaxis left
%                 plot(reshape(upsampled_psth(valid_trials, :), numel(upsampled_psth(valid_trials, :)), 1))
%                 plot(upsampled_psth(valid_trials(1), :));
%                 yyaxis right
%                 plot(reshape(current_cwt_(valid_trials, :), numel(current_cwt_(valid_trials, :)), 1))
%                 plot(current_cwt_(valid_trials(1), :));
%                 xlabel(num2str(j))
%                 subplot(2,1,2)
%                 plot(xcorr(reshape(current_cwt_, numel(current_cwt_), 1), reshape(upsampled_psth(valid_trials, :), numel(upsampled_psth(valid_trials, :)), 1), 'unbiased'))
%                 pause;
%                                 current_corr_r_ = zeros(size(upsampled_psth, 1), 1);
%                                 for l = 1:size(upsampled_psth, 1)
%                                     r_mat = corrcoef(upsampled_psth(l, cor_window_), current_cwt_(l, cor_window_));
%                                     current_corr_r_(l) = r_mat(1, 2);
%                                 end
%                                 lfp_mod_psth_corr(i, k, j) = nanmean(current_corr_r_); %    Mean across trials
            end
        end
    end
    i
end
%%
best_lfp_mod_psth_corr = zeros(size(lfp_mod_psth_corr,1), size(lfp_mod_psth_corr,3));
for i = 1:size(best_lfp_mod_psth_corr, 1)
    best_lfp_mod_psth_corr(i, :) = lfp_mod_psth_corr(i, neuron_best_class_1s_baseline(i), :);
end
%%
best_lfp_mod_psth_corr_spearman = best_lfp_mod_psth_corr;
%%
for i = 1:64
[p_temp_comp(i), ~] = ranksum(best_lfp_mod_psth_corr([neuron_tuning_cat{1, :, :, :}], i), best_lfp_mod_psth_corr([neuron_tuning_cat{2, :, :, :}], i));
p_temp_adolescent(i) = signrank(best_lfp_mod_psth_corr([neuron_tuning_cat{1, :, :, :}], i));
p_temp_adult(i) = signrank(best_lfp_mod_psth_corr([neuron_tuning_cat{2, :, :, :}], i));
end
%%
figure
hold on
plot(2:2:128, mean(best_lfp_mod_psth_corr([neuron_tuning_cat{1, :, :, :}], :), 1), 'LineWidth', 1.2, 'Color', 'b')
fill(...
    [f_range, fliplr(f_range)], ...
    [nanmean(inputs{i}, 1) + nanstd(inputs{i}, 1)/sqrt(size(inputs{i}, 1)), fliplr(nanmean(inputs{i}, 1) - nanstd(inputs{i}, 1)/sqrt(size(inputs{i}, 1)))], color_map{i}, 'FaceAlpha', 0.15, 'EdgeAlpha', 0.15);
plot(2:2:128, mean(best_lfp_mod_psth_corr([neuron_tuning_cat{2, :, :, :}], :), 1), 'LineWidth', 1.2, 'Color', 'r')
[~, h_] = bonf_holm(p_temp_adolescent, 0.05);
plot(f_range(h_), h_(h_) - 0.99, '.', 'MarkerSize', 10, 'Color', [0.2, 0.2, 0.6]);
[~, h_] = bonf_holm(p_temp_adult, 0.05);
plot(f_range(h_), h_(h_) - 1.01, '.', 'MarkerSize', 10, 'Color', [0.6, 0.2, 0.2]);
[~, h_] = bonf_holm(p_temp_comp, 0.05);
plot(f_range(h_), h_(h_) - 0.7, '*k');%%
xlabel('Frequency (Hz)')
ylabel('Correlation coefficient')
xlim([2, 128])
ylim([-0.05, 0.3])
set_plot_poster
%%
%%
for i = 1:64
[p_temp_comp(i), ~] = ranksum(best_lfp_mod_psth_corr_spearman([neuron_tuning_cat{1, :, :, :}], i), best_lfp_mod_psth_corr_spearman([neuron_tuning_cat{2, :, :, :}], i));
p_temp_adolescent(i) = signrank(best_lfp_mod_psth_corr_spearman([neuron_tuning_cat{1, :, :, :}], i));
p_temp_adult(i) = signrank(best_lfp_mod_psth_corr_spearman([neuron_tuning_cat{2, :, :, :}], i));
end
%%
h = zeros(3, numel(f_range));
[~, h_] = bonf_holm(p_temp_adolescent, 0.05);
h(1, :) = h_;
[~, h_] = bonf_holm(p_temp_adult, 0.05);
h(2, :) = h_;
[~, h_] = bonf_holm(p_temp_comp, 0.05);
h(3, :) = h_;
corr_plot({best_lfp_mod_psth_corr_spearman([neuron_tuning_cat{1, :, :, :}], :), best_lfp_mod_psth_corr_spearman([neuron_tuning_cat{2, :, :, :}], :)}, f_range, ...
    h, {'Adolescent', 'Adult'}, {'b', 'r'});
fig_name = 'temp_cor_fr_lfp';
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
savefig(gcf, fullfile(project_dir, fig_lib, fig_name));
%%