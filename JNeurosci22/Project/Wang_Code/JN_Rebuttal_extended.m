%%  Matching behavior numbers between stages
figure('Units', 'inches', 'Position', [2,2,4.5, 4.5]);
hold on;
sites_to_sort_y = [t{1, 3}]';
sites_to_sort_a = [t{2, 3}]';
% sites_to_sort_y = [t{1, 4}]';
% sites_to_sort_a = [t{2, 4}]';
[~, i_] = sort([behavior_repo(sites_to_sort_y).s_rate], 'descend');
sites_to_sort_y = sites_to_sort_y(i_);
[~, i_] = sort([behavior_repo(sites_to_sort_a).s_rate]);
sites_to_sort_a = sites_to_sort_a(i_);
[pc_y, miu_y] = zw_cum_mean([behavior_repo(sites_to_sort_y).s_rate]);
plot(100*pc_y, 100*miu_y, '.b', 'DisplayName', 'Young average', 'MarkerSize', 8);
plot(100*pc_y, 100*[behavior_repo(sites_to_sort_y).s_rate], '--b', 'DisplayName', 'Young single session', 'LineWidth', 1)
[pc_a, miu_a] = zw_cum_mean(sort([behavior_repo(sites_to_sort_a).s_rate]));
plot(100*pc_a, 100*miu_a, '.r', 'DisplayName', 'Adult average', 'MarkerSize', 8);
plot(100*pc_a, 100*[behavior_repo(sites_to_sort_a).s_rate], '--r', 'DisplayName', 'Adult single session', 'LineWidth', 1)
%   Knowing that the cum mean is monotinic, loop by pc_
pc_edge_ = zeros(2, 2);
for j = 1:(numel(pc_y) - 1)
    for i = 1:(numel(pc_a) - 1)
        if (miu_a(i) <= miu_y(j)) && (miu_a(i + 1) >= miu_y(j + 1) && pc_a(i) <= pc_y(j + 1) && pc_a(i + 1) >= pc_y(j))
            pc_edge_ = [j, j + 1; i, i + 1];
        end
    end
end
yl_ = ylim();
plot(100*mean(pc_a(pc_edge_(2, :))) + [0, 0], yl_, 'k', 'DisplayName', 'Criterion');
xtickformat('percentage');
ytickformat('percentage');
ylabel('Percentage correct choice');
xlabel('Percentile of sites recorded')
legend('Location', 'south');
set_standard_figure(12, 'bold', 1.6);
fig_name = 'performance_matching';
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dtiff', '-r400');
%%  Performance matched LFP mod
title_st = {'Alpha (8-16 Hz)', 'Beta (16-32 Hz)', 'Gamma (32-64 Hz)', 'High-Gamma (64-128 Hz)'};
for band_i_ = 1:4
    % mod_plot({lfp_mod_cue(sites_to_sort_y(1:pc_edge_(1,1)), : , band_i_), ...
    %     lfp_mod_cue(sites_to_sort_a(1:pc_edge_(2,1)), : , band_i_)}, ...
    %     linspace(-1,4,250), ...
    %     {'Ado matched', 'Adu matched'}, {'b', 'r'}, [], 12, [3.2, 3.2])
    %
    % mod_plot({lfp_mod_cue(sites_to_sort_y(1:pc_edge_(1,1)), : , band_i_), ...
    %     lfp_mod_cue(sites_to_sort_a(1:pc_edge_(2,1)), : , band_i_), ...
    %     lfp_mod_cue(sites_to_sort_y(pc_edge_(1,1):end), : , band_i_), ...
    %     lfp_mod_cue(sites_to_sort_a(pc_edge_(2,1):end), : , band_i_)}, ...
    %     linspace(-1,4,250), ...
    %     {'Ado matched', 'Adu matched', 'Ado rest', 'Adu rest'}, {'b', 'r', [0, 0, 0.6], [0.6, 0, 0]}, [], 12, [3.2, 3.2]);
    %
    mod_plot({lfp_mod_cue(sites_to_sort_y(1:pc_edge_(1,1)), : , band_i_), ...
        lfp_mod_cue(sites_to_sort_a(1:pc_edge_(2,1)), : , band_i_), ...
        lfp_mod_cue(sites_to_sort_y, : , band_i_), ...
        lfp_mod_cue(sites_to_sort_a, : , band_i_)}, ...
        linspace(-1,4,250), ...
        {'Ado matched', 'Adu matched', 'Ado all', 'Adu all'}, {'b', 'r', [0.2, 0.2, 0.6], [0.6, 0.2, 0.2]});
    xlim([-1, 2.5]);
    title(title_st{band_i_})
    xlabel('Time from cue onset (s)')
    if band_i_ ~= 3
        legend off
    end
    fig_name = sprintf('y_vs_a_performance_matched_band_%d_subj_%d', band_i_, 3);
    print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
    % print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dtiff', '-r400');
end
%%
%%  Matching behavior numbers between stages
figure('Units', 'inches', 'Position', [2,2,4.5, 4.5]);
hold on;
sites_to_sort_y = [t{1, 4}]';
sites_to_sort_a = [t{2, 4}]';
[~, i_] = sort([behavior_repo(sites_to_sort_y).s_rate], 'descend');
sites_to_sort_y = sites_to_sort_y(i_);
[~, i_] = sort([behavior_repo(sites_to_sort_a).s_rate]);
sites_to_sort_a = sites_to_sort_a(i_);
[pc_y, miu_y] = zw_cum_mean([behavior_repo(sites_to_sort_y).s_rate]);
plot(100*pc_y, 100*miu_y, '.b', 'DisplayName', 'Young average', 'MarkerSize', 8);
plot(100*pc_y, 100*[behavior_repo(sites_to_sort_y).s_rate], '--b', 'DisplayName', 'Young single session', 'LineWidth', 1)
[pc_a, miu_a] = zw_cum_mean(sort([behavior_repo(sites_to_sort_a).s_rate]));
plot(100*pc_a, 100*miu_a, '.r', 'DisplayName', 'Adult average', 'MarkerSize', 8);
plot(100*pc_a, 100*[behavior_repo(sites_to_sort_a).s_rate], '--r', 'DisplayName', 'Adult single session', 'LineWidth', 1)
%   Knowing that the cum mean is monotinic, loop by pc_
pc_edge_ = zeros(2, 2);
for j = 1:(numel(pc_y) - 1)
    for i = 1:(numel(pc_a) - 1)
        if (miu_a(i) <= miu_y(j)) && (miu_a(i + 1) >= miu_y(j + 1) && pc_a(i) <= pc_y(j + 1) && pc_a(i + 1) >= pc_y(j))
            pc_edge_ = [j, j + 1; i, i + 1];
        end
    end
end
yl_ = ylim();
plot(100*mean(pc_a(pc_edge_(2, :))) + [0, 0], yl_, 'k', 'DisplayName', 'Criterion');
xtickformat('percentage');
ytickformat('percentage');
ylabel('Percentage correct choice');
xlabel('Percentile of sites recorded')
legend('Location', 'south');
set_standard_figure(12, 'bold', 1.6);
fig_name = 'performance_matching';
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dtiff', '-r400');
%%  Performance matched LFP mod
title_st = {'Alpha (8-16 Hz)', 'Beta (16-32 Hz)', 'Gamma (32-64 Hz)', 'High-Gamma (64-128 Hz)'};
for band_i_ = 1:4
    % mod_plot({lfp_mod_cue(sites_to_sort_y(1:pc_edge_(1,1)), : , band_i_), ...
    %     lfp_mod_cue(sites_to_sort_a(1:pc_edge_(2,1)), : , band_i_)}, ...
    %     linspace(-1,4,250), ...
    %     {'Ado matched', 'Adu matched'}, {'b', 'r'}, [], 12, [3.2, 3.2])
    %
    % mod_plot({lfp_mod_cue(sites_to_sort_y(1:pc_edge_(1,1)), : , band_i_), ...
    %     lfp_mod_cue(sites_to_sort_a(1:pc_edge_(2,1)), : , band_i_), ...
    %     lfp_mod_cue(sites_to_sort_y(pc_edge_(1,1):end), : , band_i_), ...
    %     lfp_mod_cue(sites_to_sort_a(pc_edge_(2,1):end), : , band_i_)}, ...
    %     linspace(-1,4,250), ...
    %     {'Ado matched', 'Adu matched', 'Ado rest', 'Adu rest'}, {'b', 'r', [0, 0, 0.6], [0.6, 0, 0]}, [], 12, [3.2, 3.2]);
    %
    mod_plot({lfp_mod_cue(sites_to_sort_y(1:pc_edge_(1,1)), : , band_i_), ...
        lfp_mod_cue(sites_to_sort_a(1:pc_edge_(2,1)), : , band_i_), ...
        lfp_mod_cue(sites_to_sort_y, : , band_i_), ...
        lfp_mod_cue(sites_to_sort_a, : , band_i_)}, ...
        linspace(-1,4,250), ...
        {'Ado matched', 'Adu matched', 'Ado all', 'Adu all'}, {'b', 'r', [0.2, 0.2, 0.6], [0.6, 0.2, 0.2]});
    xlim([-1, 2.5]);
    title(title_st{band_i_})
    xlabel('Time from cue onset (s)')
    if band_i_ ~= 3
        legend off
    end
    fig_name = sprintf('y_vs_a_performance_matched_band_%d_subj_%d', band_i_, 4);
    print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
    % print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dtiff', '-r400');
end
%%
%%  Spiking firing rate (1000ms baseline) by spiking tuning
psth_plot({worst_psth_1s_baseline([neuron_tuning_cat{1, [1,3,4], :, 2}], :), ...
    worst_psth_1s_baseline([neuron_tuning_cat{1, [1,3,4], :, 1}], :), ...
    worst_psth_1s_baseline([neuron_tuning_cat{2, [1,3,4], :, 2}], :), ...
    worst_psth_1s_baseline([neuron_tuning_cat{2, [1,3,4], :, 1}], :)}, ...
    linspace(-1,4,250), ...
    {'Adolescent in', 'Adolescent non-in', 'Adult in', 'Adult non-in'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
%     {'Adolescent informative', 'Adolescent non-informative', 'Adult informative', 'Adult non-informative'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
xlim([-1. 2.5])
ylim([-5, 10])
title('Evoked Firing Rate')
xlabel('Time from cue onset (s)')
fig_name = ['worst_1s_baseline_spiking_firing_draft_-1000_2500'];
legend off
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
savefig(gcf, fullfile(project_dir, fig_lib, fig_name))
%%  Mean of neuronal worst psth differences: adolescent-vs-adult by inf-vs-noninf
set_neuron_worstpsth_inf_adolescent = mean(worst_psth_1s_baseline([neuron_mod_cat{1, [1,3,4], 2}], 51:75), 2);
set_neuron_worstpsth_inf_adult = mean(worst_psth_1s_baseline([neuron_mod_cat{2, [1,3,4], 2}], 51:75), 2);
set_neuron_worstpsth_noninf_adolescent = mean(worst_psth_1s_baseline([neuron_mod_cat{1, [1,3,4], 1}], 51:75), 2);
set_neuron_worstpsth_noninf_adult = mean(worst_psth_1s_baseline([neuron_mod_cat{2, [1,3,4], 1}], 51:75), 2);
% [h_, p_] = ttest2(set_neuron_pev_mod_adolescent, set_neuron_pev_mod_adult)
%%  ANOVA neuronal worst psth: Stage-by-neuron inf
anova_input_cell = {set_neuron_worstpsth_noninf_adolescent', set_neuron_worstpsth_noninf_adult'; ...
    set_neuron_worstpsth_inf_adolescent', set_neuron_worstpsth_inf_adult'};
zw_anova_from_cell(anova_input_cell)
%%  Updating mod_psth_temp_cor with median and percentiles
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
fname_ = 'best_lfp_mod_psth_corr_spearman';
save(fullfile(project_dir, output_database, fname_), 'best_lfp_mod_psth_corr', 'p_temp_adolescent', 'p_temp_adult', 'p_temp_comp');
%%
h = zeros(3, numel(f_range));
[~, h_] = bonf_holm(p_temp_adolescent, 0.05);
h(1, :) = h_;
[~, h_] = bonf_holm(p_temp_adult, 0.05);
h(2, :) = h_;
[~, h_] = bonf_holm(p_temp_comp, 0.05);
h(3, :) = h_;
median_corr_plot({best_lfp_mod_psth_corr_spearman([neuron_tuning_cat{1, :, :, :}], :), best_lfp_mod_psth_corr_spearman([neuron_tuning_cat{2, :, :, :}], :)}, f_range, ...
    h, {'Adolescent', 'Adult'}, {'b', 'r'});
fig_name = 'temp_cor_fr_lfp_median';
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
savefig(gcf, fullfile(project_dir, fig_lib, fig_name));
%%
