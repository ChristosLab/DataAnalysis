young_target_ = find(neuron_tbl_by_trial.responsive.*(neuron_tbl_by_trial.stage == 1).*keep_trial);
adult_target_ = find(neuron_tbl_by_trial.responsive.*(neuron_tbl_by_trial.stage == 2).*keep_trial);
sort_window = 76:150;
min_maxed_psth_       = (psth_cue_ - min(psth_cue_(:, sort_window), [], 2))./(max(psth_cue_(:, sort_window), [], 2) - min(psth_cue_(:, sort_window), [], 2));
evoked_psth_ = psth_cue_ - mean(psth_cue_(:, 1:50), 2) ;
min_maxed_evoked_psth_ = evoked_psth_./max(abs(evoked_psth_(:, sort_window)), [], 2);
min_maxed_gamma_      = (cue_gamma_ - min(cue_gamma_(:, sort_window), [], 2))./(max(cue_gamma_(:, sort_window), [], 2) - min(cue_gamma_(:, sort_window), [], 2));
min_maxed_cue_gamma_mod_ = cue_gamma_mod_./max(abs(cue_gamma_mod_(:, sort_window)), [], 2);
min_maxed_high_gamma_ = (cue_high_gamma_ - min(cue_high_gamma_(:, sort_window), [], 2))./(max(cue_high_gamma_(:, sort_window), [], 2) - min(cue_high_gamma_(:, sort_window), [], 2));
min_maxed_all_gamma_  = (cue_all_gamma_ - min(cue_all_gamma_(:, sort_window), [], 2))./(max(cue_all_gamma_(:, sort_window), [], 2) - min(cue_all_gamma_(:, sort_window), [], 2));
%%
save(fullfile(project_dir, output_database, 'Figure5_3_15_21_workspace'));
%%
[rho_mat_high_gamma_psth_young, p_mat_high_gamma_psth_young] = corr(cue_high_gamma_mod_(young_target_, sort_window)', psth_cue_(young_target_, sort_window)', 'type', 'Spearman');
[rho_mat_high_gamma_psth_adult, p_mat_high_gamma_psth_adult] = corr(cue_high_gamma_mod_(adult_target_, sort_window)', psth_cue_(adult_target_, sort_window)', 'type', 'Spearman');
[rho_mat_gamma_psth_young, p_mat_gamma_psth_young] = corr(cue_gamma_mod_(young_target_, sort_window)', psth_cue_(young_target_, sort_window)', 'type', 'Spearman');
[rho_mat_gamma_psth_adult, p_mat_gamma_psth_adult] = corr(cue_gamma_mod_(adult_target_, sort_window)', psth_cue_(adult_target_, sort_window)', 'type', 'Spearman');
%%  Shuffling test
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
    shuffled_rho_hg_y(i) = nanmean(rho_mat_high_gamma_psth_young(l_y_));
    shuffled_rho_hg_a(i) = nanmean(rho_mat_high_gamma_psth_adult(l_a_));
    shuffled_rho_g_y(i) = nanmean(rho_mat_gamma_psth_young(l_y_));
    shuffled_rho_g_a(i) = nanmean(rho_mat_gamma_psth_adult(l_a_));
    i
end
%%
[f_, x_] = ecdf(shuffled_rho_hg_y);
i_ = find(x_ < nanmean(diag(rho_mat_high_gamma_psth_young)), 1, 'last')
1 - f_(i_)
x_(floor(0.975*numel(x_)))
x_(ceil(0.025*numel(x_)))
%%
[f_, x_] = ecdf(shuffled_rho_hg_a);
i_ = find(x_ < nanmean(diag(rho_mat_high_gamma_psth_adult)), 1, 'last')
1 - f_(i_)
x_(floor(0.975*numel(x_)))
x_(ceil(0.025*numel(x_)))
%%
[f_, x_] = ecdf(shuffled_rho_g_y);
i_ = find(x_ < nanmean(diag(rho_mat_gamma_psth_young)), 1, 'last')
1 - f_(i_)
x_(floor(0.975*numel(x_)))
x_(ceil(0.025*numel(x_)))

%%
[f_, x_] = ecdf(shuffled_rho_g_a);
i_ = find(x_ < nanmean(diag(rho_mat_gamma_psth_adult)), 1, 'last');
1 - f_(i_)
x_(floor(0.975*numel(x_)))
x_(ceil(0.025*numel(x_)))
%%
fname_ = 'single_trial_gamma_psth_corr_shuffle';
save(fullfile(project_dir, output_database, fname_), ...
    'rho_mat_high_gamma_psth_young', 'rho_mat_high_gamma_psth_adult', ...
    'rho_mat_gamma_psth_young', 'rho_mat_gamma_psth_adult', ...
    'shuffled_rho_hg_y', 'shuffled_rho_hg_a', 'shuffled_rho_g_y', 'shuffled_rho_g_a');
%%
%%
%     'rho_mat_high_gamma_psth_young', 'rho_mat_high_gamma_psth_adult', ...
%     'rho_mat_gamma_psth_young', 'rho_mat_gamma_psth_adult', ...
%     'shuffled_rho_hg_y', 'shuffled_rho_hg_a', 'shuffled_rho_g_y', 'shuffled_rho_g_a');
%%  Do the sorting YOUNG
peak_timepoint = zeros(numel(young_target_), 1);
for i = 1:numel(young_target_)
%     [~, i_] = max(min_maxed_gamma_(young_target_(i), sort_window));
%     [~, i_] = max(min_maxed_cue_gamma_mod_(young_target_(i), sort_window));
        [~, i_] = max(cue_gamma_mod_(young_target_(i), sort_window));
%     [~, i_] = max(min_maxed_high_gamma_(young_target_(i), sort_window));
%     [~, i_] = max(min_maxed_all_gamma_(young_target_(i), sort_window));
%     [~, i_] = max(min_maxed_psth_(adult_target_(i) sort_window));        
    peak_timepoint(i) = i_;
%         peak_timepoint(i) = zw_center_of_mass(min_maxed_pev_(i, 76:150));

end
[sorted_peak_timepoint, trace_order] = sort(peak_timepoint);
peak_timepoint_gamma_young = peak_timepoint;
%%
%   Gamma mod YOUNG
figure('Unit', 'inches', 'Position', [2, 2, 3.5, 2.5])
imagesc([-1, 4], 1:numel(trace_order), min_maxed_cue_gamma_mod_(young_target_(trace_order), :), [-1, 1]);
% imagesc(cue_gamma_mod_(young_target_(trace_order), :), [-.6, .6]);
% colormap('jet')
h = colorbar;
colormap(clm_);
ylabel(h, 'Norm. gamma (32-64 Hz) modulation');
h.FontSize = 8;
% plot(sorted_peak_timepoint+75, 1:numel(target_subset), 'r');
% 
ylabel('Trial #')
xlim([-1, 3])
xlabel('Time from cue onset (s)')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',10,'FontWeight','Bold', 'LineWidth', 1);
fig_name = ['young_gamma_mod_trace_gamma_sort_3_15_21'];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
% 
%   PSTH YOUNG
% 
figure('Unit', 'inches', 'Position', [2, 2, 3.5, 2.5])
imagesc([-1, 4], 1:numel(trace_order), min_maxed_psth_(young_target_(trace_order), :), [0, 1]);
% imagesc(psth_cue_(young_target_(trace_order), :), [0, 100]);
% imagesc(evoked_psth_(young_target_(trace_order), :), [0, 15]);
% imagesc(min_maxed_evoked_psth_(young_target_(trace_order), :), [0, 1]);
% ylabel(h, 'Norm. evoked firing rate');
% colormap('jet')
h = colorbar;
ylabel(h, 'Norm. firing rate');
h.FontSize = 8;
xlim([-1, 3])
xlabel('Time from cue onset (s)')
% plot(sorted_peak_timepoint+75, 1:numel(target_subset), 'r');
% 
ylabel('Trial #')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',10,'FontWeight','Bold', 'LineWidth', 1);
fig_name = ['young_psth_trace_gamma_sort_3_15_21'];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');

%%  Do the sorting ADULT
peak_timepoint = zeros(numel(adult_target_), 1);
for i = 1:numel(adult_target_)
%     [~, i_] = max(min_maxed_gamma_(young_target_(i), sort_window));
%     [~, i_] = max(min_maxed_cue_gamma_mod_(young_target_(i), sort_window));
        [~, i_] = max(cue_gamma_mod_(adult_target_(i), sort_window));
%     [~, i_] = max(min_maxed_high_gamma_(young_target_(i), sort_window));
%     [~, i_] = max(min_maxed_all_gamma_(young_target_(i), sort_window));
%     [~, i_] = max(min_maxed_psth_(adult_target_(i) sort_window));        
    peak_timepoint(i) = i_;
%         peak_timepoint(i) = zw_center_of_mass(min_maxed_pev_(i, 76:150));

end
[sorted_peak_timepoint, trace_order] = sort(peak_timepoint);
peak_timepoint_gamma_adult = peak_timepoint;
%%
%   Gamma mod ADULT
figure('Unit', 'inches', 'Position', [2, 2, 3.5, 2.5])
imagesc([-1, 4], 1:numel(trace_order), min_maxed_cue_gamma_mod_(adult_target_(trace_order), :), [-1, 1]);
% imagesc(cue_gamma_mod_(adult_target_(trace_order), :), [-.6, .6]);
% colormap('jet')
h = colorbar;
colormap(clm_);
ylabel(h, 'Norm. gamma (32-64 Hz) modulation');
h.FontSize = 8;
% plot(sorted_peak_timepoint+75, 1:numel(target_subset), 'r');
% 
ylabel('Trial #')
xlim([-1, 3])
xlabel('Time from cue onset (s)')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',10,'FontWeight','Bold', 'LineWidth', 1);
fig_name = ['adult_gamma_mod_trace_gamma_sort_3_15_21'];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
% 
%   PSTH ADULT
% 
figure('Unit', 'inches', 'Position', [2, 2, 3.5, 2.5])
imagesc([-1, 4], 1:numel(trace_order), min_maxed_psth_(adult_target_(trace_order), :), [0, 1]);
% imagesc(psth_cue_(adult_target_(trace_order), :), [0, 100]);
% imagesc(evoked_psth_(adult_target_(trace_order), :), [0, 15]);
% imagesc(min_maxed_evoked_psth_(adult_target_(trace_order), :), [0, 1]);
% ylabel(h, 'Norm. evoked firing rate');
% colormap('jet')
h = colorbar;
ylabel(h, 'Norm. firing rate');
h.FontSize = 8;
xlim([-1, 3])
xlabel('Time from cue onset (s)')
% plot(sorted_peak_timepoint+75, 1:numel(target_subset), 'r');
% 
ylabel('Trial #')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',10,'FontWeight','Bold', 'LineWidth', 1);
fig_name = ['adult_psth_trace_gamma_sort_3_15_21'];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
%%