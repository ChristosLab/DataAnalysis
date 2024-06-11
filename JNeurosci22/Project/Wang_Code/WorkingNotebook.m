figure
imagesc(neuron_pev_cue([neuron_tuning_cat{2, [1,3,4],:, 2}], :), [0, 0.6])
%%  Sort by neuron PEV
target_subset = [neuron_tuning_cat{1, [1,3,4],:, 2}];
[~, peak_timepoint] = max(neuron_pev_cue(target_subset, :), [], 2);
[~, trace_order] = sort(peak_timepoint);
min_maxed_pev_ = (neuron_pev_cue(target_subset(trace_order), :) - min(neuron_pev_cue(target_subset(trace_order), :), [], 2))./(max(neuron_pev_cue(target_subset(trace_order), :), [], 2) - min(neuron_pev_cue(target_subset(trace_order), :), [], 2));
figure
imagesc(min_maxed_pev_, [0, 1]);
colormap('jet')
h = colorbar;
ylabel(h, 'Norm. PEV');
% 
band_i = 3;
mod_to_plot = lfp_mod_cue(map_neuron_to_site(mapping_mat, target_subset(trace_order)), :, band_i) - 1; %    Converting to change from baseline
mod_to_plot = mod_to_plot./max(mod_to_plot, [], 2);
figure
imagesc(mod_to_plot, [-1, 1])
colormap(clm_)
% 
band_i = 4;
mod_to_plot = lfp_mod_cue(map_neuron_to_site(mapping_mat, target_subset(trace_order)), :, band_i) - 1; %    Converting to change from baseline
mod_to_plot = mod_to_plot./max(mod_to_plot, [], 2);
figure
imagesc(mod_to_plot, [-1, 1])
colormap(clm_)
%%
%%  Sort by best_psth_1s_baseline
band_i = 4;
target_subset = [neuron_tuning_cat{1, [1,3,4],:, 2}];
[~, peak_timepoint] = max(best_psth_1s_baseline(target_subset, 76:200), [], 2);
[~, trace_order] = sort(peak_timepoint);
figure
imagesc(best_psth_1s_baseline(target_subset(trace_order), :), [-5, 60])
figure
imagesc(neuron_pev_cue(target_subset(trace_order), :))
figure
imagesc(lfp_mod_cue(map_neuron_to_site(mapping_mat, target_subset(trace_order)), :, band_i), [1, 1.5])
figure
imagesc(lfp_pev_cue(map_neuron_to_site(mapping_mat, target_subset(trace_order)), :, band_i), [0, 0.6])
%%
%%  Sort by best_psth_1s_baseline plot preferred stimulus LFP
band_i = 3;
peak_search_dur = 51:200;
target_subset = [neuron_tuning_cat{1, [1,3,4],:, :}];
[~, peak_timepoint] = max(best_psth_1s_baseline(target_subset, peak_search_dur), [], 2);
[~, trace_order] = sort(peak_timepoint);
figure
imagesc(best_psth_1s_baseline(target_subset(trace_order), :), [-5, 60])
% figure
% imagesc(neuron_pev_cue(target_subset(trace_order), :))
figure
imagesc(lfp_mod_neuron_best(target_subset(trace_order), :, band_i), [1, 2])

%%
%%  Sort by LFP power modulation to preferred stimulus
band_i = 4;
peak_search_dur = 76:150;
target_subset = [neuron_tuning_cat{2, [1,3,4],:, :}];
[~, peak_timepoint] = max(lfp_mod_neuron_best(target_subset, peak_search_dur, band_i), [], 2);
[sorted_peak_timepoint, trace_order] = sort(peak_timepoint);
figure
imagesc(lfp_mod_neuron_best(target_subset(trace_order), :, band_i), [1, 3])
hold on
plot(sorted_peak_timepoint + peak_search_dur(1) - 1, 1:numel(target_subset), 'r');
figure
imagesc(best_psth_1s_baseline(target_subset(trace_order), :), [-5, 40])
hold on
plot(sorted_peak_timepoint + peak_search_dur(1) - 1, 1:numel(target_subset), 'r');
% figure
% imagesc(lfp_pev_cue(target_subset(trace_order), :, band_i))

%%  Sort by LFP power modulation
band_i = 4;
target_subset = map_neuron_to_site(mapping_mat, [neuron_tuning_cat{2, [1,3,4],:, 2}]);
[~, peak_timepoint] = max(lfp_mod_cue(target_subset, :, band_i), [], 2);
[sorted_peak_timepoint, trace_order] = sort(peak_timepoint);
figure
imagesc(lfp_mod_cue(target_subset(trace_order), :, band_i), [1, 2])
hold on
plot(sorted_peak_timepoint, 1:numel(target_subset), 'r');
figure
imagesc(lfp_pev_cue(target_subset(trace_order), :, band_i))

%%  Match firing rate between stages
figure; 
histogram(mean(best_psth_raw([neuron_tuning_cat{1, [1,3,4],:, 2}], 76:150), 2))
hold on
histogram(mean(best_psth_raw([neuron_tuning_cat{2, [1,3,4],:, 2}], 76:150), 2))
figure; 
%%
figure
histogram(mean(best_psth_1s_baseline([neuron_tuning_cat{1, [1,3,4],:, :}], 76:150), 2), -5:1:30)
hold on
histogram(mean(best_psth_1s_baseline([neuron_tuning_cat{2, [1,3,4],:, :}], 76:150), 2), -5:1:30)
%%  Create raw fr matched dataset
adult_psth_raw_ = best_psth_raw([neuron_tuning_cat{2, [1,3,4],:, :}], :);
% Sort mean delay firing rate from low to high
[sorted_adult_psth_raw_delay_mean_ , idx_adult_psth_raw_] = sort(mean(adult_psth_raw_(:, 76:150), 2));
adolescent_psth_raw_ = best_psth_raw([neuron_tuning_cat{1, [1,3,4],:, :}], :);
[sorted_adolescent_psth_raw_delay_mean_ , idx_adolescent_psth_raw_] = sort(mean(adolescent_psth_raw_(:, 76:150), 2));
% Cumulate adolescent mean delay firing high to low (dropping low firing units)
[percentile_adolescent_, means_adolescent_] = zw_cum_mean(flipud(sorted_adolescent_psth_raw_delay_mean_));
% Cumulate adult mean delay firing low to high (dropping high firing units)
[percentile_adult_, means_adult_] = zw_cum_mean(sorted_adult_psth_raw_delay_mean_);
% 
figure
plot(percentile_adolescent_, fliplr(means_adolescent_))
hold on
plot(percentile_adult_, fliplr(means_adult_))
xlabel('Percent neurons dropped')
% Visual inspection selects the 10.6 percentile as crossing point
corssing = 0.106;
[~, i_adolescent_] = min(abs(percentile_adolescent_ - corssing))
[~, i_adult_] = min(abs(percentile_adult_ - (1 - corssing)))
mean(sorted_adolescent_psth_raw_delay_mean_(i_adolescent_:end))
mean(sorted_adult_psth_raw_delay_mean_(1:i_adult_))
[h_, p_] = ttest2(sorted_adolescent_psth_raw_delay_mean_(i_adolescent_:end), sorted_adult_psth_raw_delay_mean_(1:i_adult_))
%%  Create evoked fr matched dataset
adult_psth_raw_ = best_psth_1s_baseline([neuron_tuning_cat{2, [1,3,4],:, :}], :);
% Sort mean delay firing rate from low to high
[sorted_adult_psth_raw_delay_mean_ , idx_adult_psth_raw_] = sort(mean(adult_psth_raw_(:, 76:150), 2));
adolescent_psth_raw_ = best_psth_1s_baseline([neuron_tuning_cat{1, [1,3,4],:, :}], :);
[sorted_adolescent_psth_raw_delay_mean_ , idx_adolescent_psth_raw_] = sort(mean(adolescent_psth_raw_(:, 76:150), 2));
% Cumulate adolescent mean delay firing high to low (dropping low firing units)
[percentile_adolescent_, means_adolescent_] = zw_cum_mean(flipud(sorted_adolescent_psth_raw_delay_mean_));
% Cumulate adult mean delay firing low to high (dropping high firing units)
[percentile_adult_, means_adult_] = zw_cum_mean(sorted_adult_psth_raw_delay_mean_);
% 
figure
plot(percentile_adolescent_, fliplr(means_adolescent_))
hold on
plot(percentile_adult_, fliplr(means_adult_))
xlabel('Percent neurons dropped')
%%
% Visual inspection selects the 10.6 percentile as crossing point
clc
corssing = 0.055;
[~, i_adolescent_] = min(abs(percentile_adolescent_ - corssing))
[~, i_adult_] = min(abs(percentile_adult_ - (1 - corssing)))
mean(sorted_adolescent_psth_raw_delay_mean_(i_adolescent_:end))
mean(sorted_adult_psth_raw_delay_mean_(1:i_adult_))
[h_, p_, ~, stats_] = ttest2(sorted_adolescent_psth_raw_delay_mean_(i_adolescent_:end), sorted_adult_psth_raw_delay_mean_(1:i_adult_))
std(sorted_adolescent_psth_raw_delay_mean_(i_adolescent_:end))
std(sorted_adult_psth_raw_delay_mean_(1:i_adult_))
%%  Extracting indices
adolescent_set = [neuron_tuning_cat{1, [1,3,4],:, :}];
adult_set = [neuron_tuning_cat{2, [1,3,4],:, :}];
redueced_adolescent_set = adolescent_set(idx_adolescent_psth_raw_(i_adolescent_:numel(idx_adolescent_psth_raw_)));
redueced_adult_set = adult_set(idx_adult_psth_raw_(1:i_adult_)); 
%%
fname_ = 'reduced_neuron_sets';
save(fullfile(project_dir, output_database, fname_), 'redueced_adolescent_set', 'redueced_adult_set');
%%
stage_i = 1;
band_i = 1;
figure;
hold on
plot(mean(lfp_mod_cue(map_neuron_to_site(mapping_mat, [neuron_tuning_cat{stage_i, [1,3,4], :, :}]), 76:150, band_i), 2), mean(best_psth_1s_baseline([neuron_tuning_cat{stage_i, [1,3,4], :, :}], 76:150), 2), 'db', 'MarkerSize', 2)
stage_i = 2;
plot(mean(lfp_mod_cue(map_neuron_to_site(mapping_mat, [neuron_tuning_cat{stage_i, [1,3,4], :, :}]), 76:150, band_i), 2), mean(best_psth_1s_baseline([neuron_tuning_cat{stage_i, [1,3,4], :, :}], 76:150), 2), 'dr', 'MarkerSize', 2)
stage_i = 1;
figure;
hold on
plot(mean(lfp_mod_cue(map_neuron_to_site(mapping_mat, redueced_adolescent_set), 76:150, band_i), 2), mean(best_psth_1s_baseline(redueced_adolescent_set, 76:150), 2), 'db', 'MarkerSize', 2)
stage_i = 2;
plot(mean(lfp_mod_cue(map_neuron_to_site(mapping_mat, redueced_adult_set), 76:150, band_i), 2), mean(best_psth_1s_baseline(redueced_adult_set, 76:150), 2), 'dr', 'MarkerSize', 2)
%%
stage_i = 1;
band_i = 2;
figure;
hold on
plot(mean(lfp_mod_cue(map_neuron_to_site(mapping_mat, [neuron_tuning_cat{stage_i, [1,3,4], :, :}]), 76:150, band_i), 2), mean(best_psth_1s_baseline([neuron_tuning_cat{stage_i, [1,3,4], :, :}], 76:150), 2), 'db', 'MarkerSize', 2)
stage_i = 2;
plot(mean(lfp_mod_cue(map_neuron_to_site(mapping_mat, [neuron_tuning_cat{stage_i, [1,3,4], :, :}]), 76:150, band_i), 2), mean(best_psth_1s_baseline([neuron_tuning_cat{stage_i, [1,3,4], :, :}], 76:150), 2), 'dr', 'MarkerSize', 2)
stage_i = 1;
figure;
hold on
plot(mean(lfp_mod_cue(map_neuron_to_site(mapping_mat, redueced_adolescent_set), 76:150, band_i), 2), mean(best_psth_1s_baseline(redueced_adolescent_set, 76:150), 2), 'db', 'MarkerSize', 2)
stage_i = 2;
plot(mean(lfp_mod_cue(map_neuron_to_site(mapping_mat, redueced_adult_set), 76:150, band_i), 2), mean(best_psth_1s_baseline(redueced_adult_set, 76:150), 2), 'dr', 'MarkerSize', 2)
%%
stage_i = 1;
band_i = 3;
figure;
fig_name = ['delay_evoked_firing_vs_gamma_mod'];
hold on
plot(mean(lfp_mod_cue(map_neuron_to_site(mapping_mat, [neuron_tuning_cat{stage_i, [1,3,4], :, :}]), 76:150, band_i), 2) - 1, mean(best_psth_1s_baseline([neuron_tuning_cat{stage_i, [1,3,4], :, :}], 76:150), 2), 'db', 'MarkerSize', 2, ...
    'DisplayName', 'Adolescent');
stage_i = 2;
plot(mean(lfp_mod_cue(map_neuron_to_site(mapping_mat, [neuron_tuning_cat{stage_i, [1,3,4], :, :}]), 76:150, band_i), 2) - 1, mean(best_psth_1s_baseline([neuron_tuning_cat{stage_i, [1,3,4], :, :}], 76:150), 2), 'dr', 'MarkerSize', 2, ...
    'DisplayName', 'Adult');
stage_i = 1;
ylabel('Evoked firing rate (spikes/s)')
xlabel('Gamma power change from baseline')
ylim([-10, 75])
title('Delay period')
legend
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1);
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
figure;
fig_name = ['delay_evoked_firing_vs_gamma_mod_reduced_set'];
hold on
plot(mean(lfp_mod_cue(map_neuron_to_site(mapping_mat, redueced_adolescent_set), 76:150, band_i), 2) - 1, mean(best_psth_1s_baseline(redueced_adolescent_set, 76:150), 2), 'db', 'MarkerSize', 2, ...
    'DisplayName', 'Adolescent')
stage_i = 2;
plot(mean(lfp_mod_cue(map_neuron_to_site(mapping_mat, redueced_adult_set), 76:150, band_i), 2) - 1, mean(best_psth_1s_baseline(redueced_adult_set, 76:150), 2), 'dr', 'MarkerSize', 2, ...
    'DisplayName', 'Adult')
ylabel('Evoked firing rate (spikes/s)')
xlabel('Gamma power change from baseline')
ylim([-10, 75])
title('Delay period (firing rate matched)')
legend
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1);
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
%%
stage_i = 1;
band_i = 4;
figure;
fig_name = ['delay_evoked_firing_vs_high_gamma_mod'];
hold on
plot(mean(lfp_mod_cue(map_neuron_to_site(mapping_mat, [neuron_tuning_cat{stage_i, [1,3,4], :, :}]), 76:150, band_i), 2) - 1, mean(best_psth_1s_baseline([neuron_tuning_cat{stage_i, [1,3,4], :, :}], 76:150), 2), 'db', 'MarkerSize', 2, ...
    'DisplayName', 'Adolescent');
stage_i = 2;
plot(mean(lfp_mod_cue(map_neuron_to_site(mapping_mat, [neuron_tuning_cat{stage_i, [1,3,4], :, :}]), 76:150, band_i), 2) - 1, mean(best_psth_1s_baseline([neuron_tuning_cat{stage_i, [1,3,4], :, :}], 76:150), 2), 'dr', 'MarkerSize', 2, ...
    'DisplayName', 'Adult');
stage_i = 1;
ylabel('Evoked firing rate (spikes/s)')
xlabel('High gamma power change from baseline')
ylim([-10, 75])
title('Delay period')
legend
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1);
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
figure;
fig_name = ['delay_evoked_firing_vs_high_gamma_mod_reduced_set'];
hold on
plot(mean(lfp_mod_cue(map_neuron_to_site(mapping_mat, redueced_adolescent_set), 76:150, band_i), 2) - 1, mean(best_psth_1s_baseline(redueced_adolescent_set, 76:150), 2), 'db', 'MarkerSize', 2, ...
    'DisplayName', 'Adolescent')
stage_i = 2;
plot(mean(lfp_mod_cue(map_neuron_to_site(mapping_mat, redueced_adult_set), 76:150, band_i), 2) - 1, mean(best_psth_1s_baseline(redueced_adult_set, 76:150), 2), 'dr', 'MarkerSize', 2, ...
    'DisplayName', 'Adult')
ylabel('Evoked firing rate (spikes/s)')
xlabel('High gamma power change from baseline')
ylim([-10, 75])
title('Delay period (firing rate matched)')
legend
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1);
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
%%
figure; 
hold on
histogram(mean(lfp_pev_cue([site_tuning_cat{1,[1,3,4], :, :}], 76:150, 3), 2), 'Normalization', 'probability');
histogram(mean(neuron_pev_cue([neuron_tuning_cat{1,[1,3,4],:, :}], 76:150), 2), 'Normalization', 'probability');
%%
figure; 
hold on
histogram(mean(lfp_pev_cue([site_tuning_cat{2,[1,3,4],:, :}], 76:150, 3), 2), 'Normalization', 'probability');
histogram(mean(neuron_pev_cue([neuron_tuning_cat{2,[1,3,4],:, :}], 76:150), 2), 'Normalization', 'probability');
%%
figure; 
hold on
histogram(mean(lfp_pev_cue([site_tuning_cat{1,[1,3,4],:, :}], 76:150, 3), 2), 'Normalization', 'probability');
histogram(mean(lfp_pev_cue([site_tuning_cat{2,[1,3,4],:, :}], 76:150, 3), 2), 'Normalization', 'probability');
%%
figure; 
hold on
histogram(mean(neuron_pev_cue([neuron_tuning_cat{1,[1,3,4],:, :}], 76:150), 2), 'Normalization', 'probability');
histogram(mean(neuron_pev_cue([neuron_tuning_cat{2,[1,3,4],:, :}], 76:150), 2), 'Normalization', 'probability');
%%
temp_set_ = [site_tuning_cat{2,[1,3,4],2, 2}];
figure
for i = temp_set_
    plot(lfp_pev_cue(i, :, 3))
    pause
end
%%
