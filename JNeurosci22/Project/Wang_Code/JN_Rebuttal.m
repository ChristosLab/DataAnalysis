%%  JN comments rebuttal
clear
zw_setpath;
%%
fname_ = 'odr_behavior_repo';
load(fullfile(project_dir, output_database, fname_));
%
fname_ = 'tfr_modulation.mat';
load(fullfile(project_dir, output_database, fname_), 'tfr_total_modulation_flag', 'tfr_cue_modulation_flag', 'tfr_delay_modulation_flag', 'tfr_modulation_output_cue');
fname_ = 'neuron_site_categories.mat';
load(fullfile(project_dir, output_database, fname_), 'site_mod_cat', 'site_tuning_cat', 'neuron_mod_cat', 'neuron_tuning_cat', 'mapping_mat', 'responsive');
fname_ = 'lfp_neuron_matching.mat';
load(fullfile(project_dir, output_database, fname_), 't', 'n', 'mapping_mat');
load(fullfile(project_dir, output_database, 'lfp_tbl.mat'))
load(fullfile(project_dir, output_database, 'neuron_tbl_w_pfc.mat'))
%%  Reviwer 1 - Comment 1
%% gamma delay - behavior
figure
hold on
plot(mean([lfp_mod_cue([t{1, :}]', 76:150, 3)], 2), [behavior_repo([t{1, :}]').s_rate], 'r*')
plot(mean([lfp_mod_cue([t{2, :}]'', 76:150, 3)], 2), [behavior_repo([t{2, :}]').s_rate], 'b*')
% yl_ = ylim();
% plot(mean(mean([lfp_mod_cue([t{1, :}]', 76:150, 3)], 2)) + [0, 0], yl_, 'r');
% plot(mean(mean([lfp_mod_cue([t{2, :}]', 76:150, 3)], 2)) + [0, 0], yl_, 'b');
%% Cue power TO behavior -- by subject
target_subject_ = [1, 3, 4];
monkey_id = {'ind', 'jac', 'lem', 'ken'};
ylabel_name_ = {'Alpha', 'Beta', 'Gamma', 'High gamma'};
target_epoch_ = 51:75;
for k = 1:4
    figure('Units', 'inches', 'Position', [2,2,8,3])
    for i = 1:numel(target_subject_)
        subplot(1, 3, i);
        hold on
        x1 = [behavior_repo([t{1, target_subject_(i)}]').s_rate]';
        x2 = [behavior_repo([t{2, target_subject_(i)}]').s_rate]';
        y1 = mean([lfp_mod_cue([t{1, target_subject_(i)}]', target_epoch_, k)], 2);
        y2 = mean([lfp_mod_cue([t{2, target_subject_(i)}]', target_epoch_, k)], 2);
        plot(100*x1, y1, 'b*');
        plot(100*x2, y2, 'r*');
        [rho1, p1] = corr(x1, y1, 'Type', 'spearman');
        %     [rho2, p2] = corr(x2, y2, 'Type', 'spearman');
        text(70, 0.9, sprintf('rho = %.2f\np = %.3f', rho1, p1));
        xlim(100*[0.65, 1]);
        ylim([0.8, 1.5]);
        xtickformat('percentage');
        xlabel('Percentage correct');
        ylabel(sprintf('%s power', ylabel_name_{k}));
        title(monkey_id{target_subject_(i)});
        set_standard_figure(12, 'bold', 1.6);
    end
    fig_name = sprintf('cue_pow_vs_beh_corr_%s', ylabel_name_{k});
    print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
    print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dtiff', '-r400');
end

%% Delay power TO behavior -- by subject
target_subject_ = [1, 3, 4];
monkey_id = {'ind', 'jac', 'lem', 'ken'};
ylabel_name_ = {'Alpha', 'Beta', 'Gamma', 'High gamma'};
target_epoch_ = 76:150;
for k = 1:4
    figure('Units', 'inches', 'Position', [2,2,8,3])
    for i = 1:numel(target_subject_)
        subplot(1, 3, i);
        hold on
        x1 = [behavior_repo([t{1, target_subject_(i)}]').s_rate]';
        x2 = [behavior_repo([t{2, target_subject_(i)}]').s_rate]';
        y1 = mean([lfp_mod_cue([t{1, target_subject_(i)}]', target_epoch_, k)], 2);
        y2 = mean([lfp_mod_cue([t{2, target_subject_(i)}]', target_epoch_, k)], 2);
        plot(100*x1, y1, 'b*');
        plot(100*x2, y2, 'r*');
        [rho1, p1] = corr(x1, y1, 'Type', 'spearman');
        %     [rho2, p2] = corr(x2, y2, 'Type', 'spearman');
        text(70, 0.9, sprintf('rho = %.2f\np = %.3f', rho1, p1));
        xlim(100*[0.65, 1]);
        ylim([0.8, 1.5]);
        xtickformat('percentage');
        xlabel('Percentage correct');
        ylabel(sprintf('%s power', ylabel_name_{k}));
        title(monkey_id{target_subject_(i)});
        set_standard_figure(12, 'bold', 1.6);
    end
    fig_name = sprintf('delay_pow_vs_beh_corr_%s', ylabel_name_{k});
    print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
    print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dtiff', '-r400');
end
%%  Session order TO delay gamma -- by subject
figure('Units', 'inches', 'Position', [2,2,14, 3])
for i = [1, 3, 4]
    subplot(1, 4, i);
    hold on
    plot(lfp_tbl.session_id(([t{1, i}])), mean([lfp_mod_cue([t{1, i}]', 76:150, 3)], 2), 'b*')
    plot(lfp_tbl.session_id(([t{2, i}])), mean([lfp_mod_cue([t{2, i}]', 76:150, 3)], 2), 'r*')
    ylim([0.8, 1.5]);
    %     xlim([0.8, 1.5]);
    set_standard_figure
end
%%  Session order TO behavior -- by subject
figure('Units', 'inches', 'Position', [2,2,14, 3])
for i = [1, 3, 4]
    subplot(1, 4, i);
    hold on
    plot(lfp_tbl.session_id(([t{1, i}])), [behavior_repo([t{1, i}]').s_rate], 'b*')
    plot(lfp_tbl.session_id(([t{2, i}])), [behavior_repo([t{2, i}]').s_rate], 'r*')
    ylim([0.65, 1]);
    %     xlim([0.8, 1.5]);
    set_standard_figure
end

%%  Moment-by-moment corr beteeen gamma power and behavior
for_plot_ = zeros(100, 2, 4);
for i = 1:size(for_plot_, 1)
    for j = 1:2
        for k = [1, 3, 4]
            c_ = corrcoef(mean([lfp_mod_cue([t{j, k}]', 50 + i, 3)], 2), [behavior_repo([t{j, k}]').s_rate]);
            %             c_ = corrcoef(mean([lfp_mod_cue([t{j, k}]', 76:(76+i - 1), 3)], 2), [behavior_repo([t{j, k}]').s_rate]);
            for_plot_(i, j, k) = c_(2, 1);
        end
    end
end
%%
figure('Units', 'inches', 'Position', [2,2,10, 3])
for i = [1, 3, 4]
    subplot(1, 4, i)
    hold on
    plot(for_plot_(:, 1, i), 'b');
    plot(for_plot_(:, 2, i), 'r');
    ylim([-0.8, .25])
    set_standard_figure
end
%%  Session order TO delay firing rate -- by subject
figure('Units', 'inches', 'Position', [2,2,14, 3])
for i = [1, 3, 4]
    subplot(1, 4, i);
    xdata = {neuron_tbl.session_id([n{1, i}]), ...
        neuron_tbl.session_id([n{2, i}])};
    ydata = {mean([best_psth_raw([n{1, i}]', 1:50)], 2), ...
        mean([best_psth_raw([n{2, i}]', 1:50)], 2)};
    hold on
    plot(xdata{1}, ydata{1}, 'b*');
    plot(xdata{2}, ydata{2}, 'r*');
    xl_ = xlim();
    plot(xl_, [0, 0] + mean(ydata{1}), 'b');
    plot(xl_, [0, 0] + mean(ydata{2}), 'r');
    ylim([-10, 80]);
    set_standard_figure
end
%%  Session order TO delay firing rate -- by subject
figure('Units', 'inches', 'Position', [2,2,6, 3])
for i = [1, 3, 4]
    subplot(1, 4, i);
    xdata = {neuron_tbl.session_id([n{1, i}]), ...
        neuron_tbl.session_id([n{2, i}])};
    ydata = {mean([best_psth_raw([n{1, i}]', 76:150)], 2)./mean([best_psth_raw([n{1, i}]', 1:50)], 2), ...
        mean([best_psth_raw([n{2, i}]', 76:150)], 2)./mean([best_psth_raw([n{2, i}]', 1:50)], 2)};
    hold on
    plot(xdata{1}, ydata{1}, 'b*');
    plot(xdata{2}, ydata{2}, 'r*');
    xl_ = xlim();
    plot(xl_, [0, 0] + mean(ydata{1}), 'b');
    plot(xl_, [0, 0] + mean(ydata{2}), 'r');
    %     ylim([0, ])
    set_standard_figure(12, 'bold', 1.6);
end
%%  Matching behavior numbers between stages
figure('Units', 'inches', 'Position', [2,2,4.5, 4.5]);
hold on;
% sites_to_sort_y = [t{1, 3}]';
% sites_to_sort_a = [t{2, 3}]';
% sites_to_sort_y = [t{1, 4}]';
% sites_to_sort_a = [t{2, 4}]';
sites_to_sort_y = [t{1, :}]';
sites_to_sort_a = [t{2, :}]';
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
    fig_name = sprintf('y_vs_a_performance_matched_band_%d', band_i_);
    print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
    % print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dtiff', '-r400');
end
%%  LFP mod by subject
monkey_id = {'IN', 'JA', 'LE', 'KE'};
target_subject_ = [1, 3, 4];
n_bands = 4;
title_st = {'Alpha (8-16 Hz)', 'Beta (16-32 Hz)', 'Gamma (32-64 Hz)', 'High-Gamma (64-128 Hz)'};
for k = 1:numel(target_subject_)
    for i = 1:4
        
        mod_plot({lfp_mod_cue([t{1, target_subject_(k)}] , :, i), lfp_mod_cue([t{2, target_subject_(k)}] , :, i)}, ...
            linspace(-1,4,250), {'Adolescent', 'Adult'}, {'b' 'r'}, [], 8, [2, 2])
        title([title_st{i}, ' ', monkey_id{target_subject_(k)}]);
        xlabel('Time from cue onset (s)')
        xlim([-1, 2.5]);
        if i ~= 3
            legend off
        end
        fig_name = ['lfp_pow_band_', num2str(i), '_cue_subject_', num2str(target_subject_(k))];
        print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
%         print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r300');
        % savefig(gcf, fullfile(project_dir, fig_lib, fig_name));
    end
end
%%  Reviewer 1 - Comment 2
%%  Error vs. correct ANOVA
%  Mean of neuronal PEV differences: adolescent_gamma vs adult_gamma
set_neuron_pev_mod_adolescent = mean(neuron_pev_cue([neuron_mod_cat{1, [1,3,4], 2}], 31:60), 2);
set_neuron_pev_mod_adult = mean(neuron_pev_cue([neuron_mod_cat{2, [1,3,4], 2}], 31:60), 2);
set_neuron_pev_non_mod_adolescent = mean(neuron_pev_cue([neuron_mod_cat{1, [1,3,4], 1}], 31:60), 2);
set_neuron_pev_non_mod_adult = mean(neuron_pev_cue([neuron_mod_cat{2, [1,3,4], 1}], 31:60), 2);
[h_, p_] = ttest2(set_neuron_pev_mod_adolescent, set_neuron_pev_mod_adult)
%  ANOVA Neuron PEV: Stage-by-Gamma mod
anova_input_cell = {set_neuron_pev_non_mod_adolescent', set_neuron_pev_non_mod_adult'; set_neuron_pev_mod_adolescent', set_neuron_pev_mod_adult'};
zw_anova_from_cell(anova_input_cell)

%%  Reviewer 1 - Comment 3
%%  Neuron - gamma best class comparison
relative_class = nan(size(neuron_tbl, 1), 1);
for i = 1:size(neuron_tbl, 1)
    current_lfp_class_ = gamma_power_best_class(find(mapping_mat(:, i)));
    if ~isempty(current_lfp_class_)
        relative_class(i) = get_relative_rad(neuron_best_class_1s_baseline(i), current_lfp_class_);
    end
end
%   Overlaying histogram for neuron-relative-to-lfp class
figure; hold on;
% histogram(relative_class([neuron_tuning_cat{1, :, 2, 2}]), 'Normalization', 'probability', 'FaceColor', 'b');
% histogram(relative_class([neuron_tuning_cat{2, :, 2, 2}]), 'Normalization', 'probability', 'FaceColor', 'r');
%%  Average tuning curves
fname_ = 'baseline_normalized_tfr_repo.mat';
load(fullfile(project_dir, output_database, fname_), 'baseline_normalized_tfr_repo');
%%  LFP tuning curve
set_lfp_mod_by_class_ = cell(0);
target_sites = {...
    [site_tuning_cat{1,:,1,:}], [site_tuning_cat{1,:,2,:}]; ...
    [site_tuning_cat{2,:,1,:}], [site_tuning_cat{2,:,2,:}]};
target_sites_hg = {...
    [site_tuning_cat_hg{1,:,1,:}], [site_tuning_cat_hg{1,:,2,:}]; ...
    [site_tuning_cat_hg{2,:,1,:}], [site_tuning_cat_hg{2,:,2,:}]};
% target_sites = {...
%     [site_tuning_cat{1,:,:,1}], [site_tuning_cat{1,:,:,2}]; ...
%     [site_tuning_cat{2,:,:,1}], [site_tuning_cat{2,:,:,2}]};

target_bins  = 76:150;
for band_i_ = 1:2
for target_stage = 1:2
for i = 1:size(target_sites, 2)
    current_sites_    = target_sites{target_stage, i};
    current_sites_hg_ = target_sites_hg{target_stage, i};
    set_lfp_mod_by_class_{target_stage, i, 1} = zeros(numel(current_sites_), 8);
    set_lfp_best_class_{target_stage, i, 1} = gamma_power_best_class(current_sites_);
    set_lfp_mod_by_class_{target_stage, i, 2} = zeros(numel(current_sites_hg_), 8);
    set_lfp_best_class_{target_stage, i, 2} = high_gamma_power_best_class(current_sites_hg_);
for j = 1:numel(current_sites_)
    for k = 1:8
    set_lfp_mod_by_class_{target_stage, i, 1}(j, k) = nanmean(nanmean(baseline_normalized_tfr_repo(current_sites_(j)).class(k).cue_tfr{3}(:, target_bins)));
    end
end
for j = 1:numel(current_sites_hg_)
    for k = 1:8
    set_lfp_mod_by_class_{target_stage, i, 2}(j, k) = nanmean(nanmean(baseline_normalized_tfr_repo(current_sites_hg_(j)).class(k).cue_tfr{4}(:, target_bins)));
    end
end

% baseline_normalized_tfr_repo;
end
end
end
%%  By stage
tuning_curve_plot(set_lfp_mod_by_class_(1,:), set_lfp_best_class_(1,:), {'Ado non', 'ado sel'}, {[0.2, 0.2, 0.6], [0, 0, 1]})
tuning_curve_plot(set_lfp_mod_by_class_(2,:), set_lfp_best_class_(2,:), {'Adu non', 'adu sel'}, {[0.6, 0.2, 0.2], [1, 0, 0]})
%%  All 4
tuning_curve_plot(set_lfp_mod_by_class_(:, :, 1), set_lfp_best_class_(:, :, 1), {'Ado non', 'Adu non', 'ado sel', 'adu sel'}, {[0.2, 0.2, 0.6], [0.6, 0.2, 0.2], [0, 0, 1], [1, 0, 0]})
title(title_st{3})
ylim([1, 1.85]);
fig_name = 'gamma_tuning_curve_pub';
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
savefig(gcf, fullfile(project_dir, fig_lib, fig_name));
tuning_curve_plot(set_lfp_mod_by_class_(:, :, 2), set_lfp_best_class_(:, :, 2), {'Ado non', 'Adu non', 'ado sel', 'adu sel'}, {[0.2, 0.2, 0.6], [0.6, 0.2, 0.2], [0, 0, 1], [1, 0, 0]})
title(title_st{4})
ylim([1, 1.85]);
fig_name = 'high_gamma_tuning_curve_pub';
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
savefig(gcf, fullfile(project_dir, fig_lib, fig_name));
%%  PSTH tuning_curve
%%
set_psth_by_class_ = cell(0);
target_neurons = {...
    [neuron_tuning_cat{1, :, :, 1}], [neuron_tuning_cat{1, :, :, 2}]; ...
    [neuron_tuning_cat{2, :, :, 1}], [neuron_tuning_cat{2, :, :, 2}]};
target_bins  = 76:150;
for target_stage = 1:2
for i = 1:size(target_neurons, 2)
    current_neurons_ = target_neurons{target_stage, i};
set_psth_by_class_{target_stage, i} = zeros(numel(current_neurons_), 8);
    set_neuron_best_class_{target_stage, i} = neuron_best_class_1s_baseline(current_neurons_);
for j = 1:numel(current_neurons_)
    for k = 1:8
    set_psth_by_class_{target_stage, i}(j, k) = nanmean(nanmean(temp_psth_1s_baseline(current_neurons_(j), k, target_bins)));
    end
end
end
end
%%  By stage
neuron_tuning_curve_plot(set_psth_by_class_(1,:), set_neuron_best_class_(1,:), {'Ado non', 'ado sel'}, {[0.2, 0.2, 0.6], [0, 0, 1]})
neuron_tuning_curve_plot(set_psth_by_class_(2,:), set_neuron_best_class_(2,:), {'Adu non', 'adu sel'}, {[0.6, 0.2, 0.2], [1, 0, 0]})
%%  All 4
neuron_tuning_curve_plot(set_psth_by_class_(:), set_neuron_best_class_(:), {'Ado non', 'Adu non', 'ado sel', 'adu sel'}, {[0.2, 0.2, 0.6], [0.6, 0.2, 0.2], [0, 0, 1], [1, 0, 0]})
title('Spiking')
fig_name = 'neuron_tuning_curve_pub';
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
savefig(gcf, fullfile(project_dir, fig_lib, fig_name));
%%   Relative angle - Polar plot
%   All neurons
% relative_class_counts_y = histcounts(relative_class([neuron_tuning_cat{1, :, :, :}]), -3.5:4.5);
% relative_class_counts_a = histcounts(relative_class([neuron_tuning_cat{2, :, :, :}]), -3.5:4.5);
% class_polar_plot(relative_class_counts_y, relative_class_counts_a);
%   Selective neurons and Selective LFP only 
relative_class_counts_y = histcounts(relative_class([neuron_tuning_cat{1, :, 2, 2}]), -3.5:4.5);
relative_class_counts_a = histcounts(relative_class([neuron_tuning_cat{2, :, 2, 2}]), -3.5:4.5);
class_polar_plot(relative_class_counts_y, relative_class_counts_a, '-');
fig_name = 'selective_tuning_polarplot';
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
savefig(gcf, fullfile(project_dir, fig_lib, fig_name))
%%
relative_class_counts_y(4)
sum(relative_class_counts_y)
relative_class_counts_a(4)
sum(relative_class_counts_a)
%%
[h_, p_, stat_] = fishertest([relative_class_counts_y(4), sum(relative_class_counts_y([1:3, 5:8])); relative_class_counts_a(4), sum(relative_class_counts_a([1:3, 5:8]))])
%%
%%  Reviewer 2 - comment 1
[h_, p_, ~, stat_] = vartest2(mean(lfp_mod_cue([t{2, :}], 51:75, 3), 2), mean(lfp_mod_cue([t{1, :}], 51:75, 3), 2))
[h_, p_, ~, stat_] = vartest2(mean(lfp_mod_cue([t{2, :}], 76:150, 3), 2), mean(lfp_mod_cue([t{1, :}], 76:150, 3), 2))
%%  Reviewer 2 - comment 1 -- reversed
[h_, p_, ~, stat_] = vartest2(mean(lfp_mod_cue([t{1, :}], 51:75, 3), 2), mean(lfp_mod_cue([t{2, :}], 51:75, 3), 2))
[h_, p_, ~, stat_] = vartest2(mean(lfp_mod_cue([t{1, :}], 76:150, 3), 2), mean(lfp_mod_cue([t{2, :}], 76:150, 3), 2))

%%  Adds depth information to lfp_tbl
dep_ = mapping_mat'.*neuron_tbl.Depth;
for i = 1:size(lfp_tbl, 1)
    if(~any(mapping_mat(i, :)))
        lfp_tbl.Depth(i) = nan;
        continue
    end
    current_dep_ = dep_(find(mapping_mat(i, :)), i);
    current_dep_(isnan(current_dep_)) = [];
    if isempty(current_dep_)
        continue
    end
    if all(current_dep_ == current_dep_', 'all')
        lfp_tbl.Depth(i) = current_dep_(1);
    else
        current_dep_
    end
end
%%
figure
hold on
plot(lfp_tbl.Depth([t{1, :}]), mean(lfp_mod_cue([t{1, :}], 76:150, 3), 2), '*')
plot(lfp_tbl.Depth([t{2, :}]), mean(lfp_mod_cue([t{2, :}], 76:150, 3), 2), '*')
%%  Add brain region information to lfp/neuron table
fname_ = fullfile(project_dir, output_database, 'NeuronArea.xlsx');
neuron_area = readtable(fname_, 'sheet', 'Sheet2');
neuron_area.Filename = lower(neuron_area.Filename);
neuron_area.area_code = 1*strcmp(neuron_area.Area, '8') + 2*strcmp(neuron_area.Area, '46');
%%
for i = 1:size(lfp_tbl, 1)
    found_entry_ = find(strcmp(neuron_area.Filename, lfp_tbl.Filename(i)));
    if ~isempty(found_entry_)
        lfp_tbl.area_code(i) = neuron_area.area_code(found_entry_(1));
    else
        lfp_tbl.area_code(i) = 0;
    end
end
%%
set_lfp_ = cell([2, 2]);
out_ = zeros(size(set_lfp_));
in_ = set_lfp_;
target_bins = 76:150;
% target_bins = 51:75;
for i = 1:2
    for j = 1:2
        set_lfp_{i, j} = intersect([t{i, :}], find(lfp_tbl.area_code == j));
        in_{i, j}      = mean(lfp_mod_cue(set_lfp_{i, j}, target_bins, 3), 2)';
        var(mean(lfp_mod_cue(set_lfp_{i, j}, target_bins, 3), 2));
        out_(i,j) = mean(mean(lfp_mod_cue(set_lfp_{i, j}, target_bins, 3)));
    end
end
%% F-test of each region
[h_, p_, ~, stat_] = vartest2([in_{2, 1}], [in_{1, 1}])
%%
[h_, p_, ~, stat_] = vartest2([in_{2, 2}], [in_{1, 2}])
%%
title_st = {'Alpha (8-16 Hz)', 'Beta (16-32 Hz)', 'Gamma (32-64 Hz)', 'High-Gamma (64-128 Hz)'};
for band_i_ = 1:4
    mod_plot({lfp_mod_cue(set_lfp_{1, 1}, : , band_i_), ...
        lfp_mod_cue(set_lfp_{2, 1}, : , band_i_), ...
        lfp_mod_cue(set_lfp_{1, 2}, : , band_i_), ...
        lfp_mod_cue(set_lfp_{2, 2}, : , band_i_)}, ...
        linspace(-1,4,250), ...
        {'Adolescent 8a', 'Adult 8a', 'Adolescent 46', 'Adult 46'}, {'b', 'r', [0.2, 0.2, 0.6], [0.6, 0.2, 0.2]}, [], 8, [2, 2]);
    xlim([-1, 2.5]);
    title(title_st{band_i_})
    xlabel('Time from cue onset (s)')
    if band_i_ ~= 3
        legend off
    end
    fig_name = sprintf('lfp_mod_stage_by_region_band_%d', band_i_);
    print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
%     print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
%     savefig(gcf, fullfile(project_dir, fig_lib, fig_name));
end

%%  Reviewer 2 - Comment 2
fname_ = 'single_trial_gamma_psth_corr_shuffle';
load(fullfile(project_dir, output_database, fname_), ...
    'rho_mat_gamma_psth_adult', 'rho_mat_high_gamma_psth_adult', ...
    'rho_mat_gamma_psth_young', 'rho_mat_high_gamma_psth_young');
%%
nanmedian(diag(rho_mat_gamma_psth_young))
prctile(diag(rho_mat_gamma_psth_young), [25, 75])
nanmedian(diag(rho_mat_gamma_psth_adult))
prctile(diag(rho_mat_gamma_psth_adult), [25, 75])
%%
nanmedian(diag(rho_mat_high_gamma_psth_young))
prctile(diag(rho_mat_high_gamma_psth_young), [25, 75])
nanmedian(diag(rho_mat_high_gamma_psth_adult))
prctile(diag(rho_mat_high_gamma_psth_adult), [25, 75])
%%
%%  Shuffling test for median
fname_ = 'Figure5_3_15_21_workspace.mat';
load(fullfile(project_dir, output_database, fname_), 'young_target_', 'adult_target_');
%%
sort_window = 76:150;
n_it = 100000;
shuffled_rho_hg_y = zeros(1, n_it);
shuffled_rho_hg_a = zeros(1, n_it);
shuffled_rho_g_y = zeros(1, n_it);
shuffled_rho_g_a = zeros(1, n_it);
n_y_ = numel(young_target_);
n_a_ = numel(adult_target_);
for i = 1:n_it
    y_ = randperm(n_y_);
    l_y_ = sub2ind([n_y_, n_y_], 1:n_y_, y_);
    a_ = randperm(n_a_);
    l_a_ = sub2ind([n_a_, n_a_], 1:n_a_, a_);
    shuffled_rho_hg_y(i) = nanmedian(rho_mat_high_gamma_psth_young(l_y_));
    shuffled_rho_hg_a(i) = nanmedian(rho_mat_high_gamma_psth_adult(l_a_));
    shuffled_rho_g_y(i) = nanmedian(rho_mat_gamma_psth_young(l_y_));
    shuffled_rho_g_a(i) = nanmedian(rho_mat_gamma_psth_adult(l_a_));
    i
end
%%
[f_, x_] = ecdf(shuffled_rho_hg_y);
i_ = find(x_ < nanmedian(diag(rho_mat_high_gamma_psth_young)), 1, 'last')
1 - f_(i_)
x_(floor(0.975*numel(x_)))
x_(ceil(0.025*numel(x_)))
%%
[f_, x_] = ecdf(shuffled_rho_hg_a);
i_ = find(x_ < nanmedian(diag(rho_mat_high_gamma_psth_adult)), 1, 'last')
1 - f_(i_)
x_(floor(0.975*numel(x_)))
x_(ceil(0.025*numel(x_)))
%%
[f_, x_] = ecdf(shuffled_rho_g_y);
i_ = find(x_ < nanmedian(diag(rho_mat_gamma_psth_young)), 1, 'last')
1 - f_(i_)
x_(floor(0.975*numel(x_)))
x_(ceil(0.025*numel(x_)))

%%
[f_, x_] = ecdf(shuffled_rho_g_a);
i_ = find(x_ < nanmedian(diag(rho_mat_gamma_psth_adult)), 1, 'last');
1 - f_(i_)
x_(floor(0.975*numel(x_)))
x_(ceil(0.025*numel(x_)))