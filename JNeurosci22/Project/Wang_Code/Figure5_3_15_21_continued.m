%%  Do the sorting YOUNG
peak_timepoint = zeros(numel(young_target_), 1);
for i = 1:numel(young_target_)
%     [~, i_] = max(min_maxed_gamma_(young_target_(i), sort_window));
%     [~, i_] = max(min_maxed_cue_gamma_mod_(young_target_(i), sort_window));
%         [~, i_] = max(cue_gamma_mod_(young_target_(i), sort_window));
    [~, i_] = max(min_maxed_cue_high_gamma_mod_(young_target_(i), sort_window));
%     [~, i_] = max(min_maxed_all_gamma_(young_target_(i), sort_window));
%     [~, i_] = max(min_maxed_psth_(adult_target_(i) sort_window));        
    peak_timepoint(i) = i_;
%         peak_timepoint(i) = zw_center_of_mass(min_maxed_pev_(i, 76:150));

end
[sorted_peak_timepoint, trace_order] = sort(peak_timepoint);
peak_timepoint_hi_gamma_young = peak_timepoint;
%%
%   High Gamma mod YOUNG
figure('Unit', 'inches', 'Position', [2, 2, 3.5, 2.5])
imagesc([-1, 4], 1:numel(trace_order), min_maxed_cue_high_gamma_mod_(young_target_(trace_order), :), [-1, 1]);
% imagesc(cue_gamma_mod_(young_target_(trace_order), :), [-.6, .6]);
% colormap('jet')
h = colorbar;
colormap(clm_);
ylabel(h, 'Norm. high gamma (64-128 Hz) modulation');
h.FontSize = 8;
% plot(sorted_peak_timepoint+75, 1:numel(target_subset), 'r');
% 
ylabel('Trial #')
xlim([-1, 3])
xlabel('Time from cue onset (s)')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',10,'FontWeight','Bold', 'LineWidth', 1);
fig_name = ['young_high_gamma_mod_trace_high_gamma_sort_3_15_21'];
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
fig_name = ['young_psth_trace_high_gamma_sort_3_15_21'];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');

%%  Do the sorting ADULT
peak_timepoint = zeros(numel(adult_target_), 1);
for i = 1:numel(adult_target_)
%     [~, i_] = max(min_maxed_gamma_(young_target_(i), sort_window));
%     [~, i_] = max(min_maxed_cue_gamma_mod_(young_target_(i), sort_window));
%         [~, i_] = max(cue_gamma_mod_(adult_target_(i), sort_window));
    [~, i_] = max(min_maxed_cue_high_gamma_mod_(adult_target_(i), sort_window));
%     [~, i_] = max(min_maxed_all_gamma_(young_target_(i), sort_window));
%     [~, i_] = max(min_maxed_psth_(adult_target_(i) sort_window));        
    peak_timepoint(i) = i_;
%         peak_timepoint(i) = zw_center_of_mass(min_maxed_pev_(i, 76:150));

end
[sorted_peak_timepoint, trace_order] = sort(peak_timepoint);
peak_timepoint_hi_gamma_adult = peak_timepoint;
%%
%   Gamma mod ADULT
figure('Unit', 'inches', 'Position', [2, 2, 3.5, 2.5])
imagesc([-1, 4], 1:numel(trace_order), min_maxed_cue_high_gamma_mod_(adult_target_(trace_order), :), [-1, 1]);
% imagesc(cue_gamma_mod_(adult_target_(trace_order), :), [-.6, .6]);
% colormap('jet')
h = colorbar;
colormap(clm_);
ylabel(h, 'Norm. high gamma (64-128 Hz) modulation');
h.FontSize = 8;
% plot(sorted_peak_timepoint+75, 1:numel(target_subset), 'r');
% 
ylabel('Trial #')
xlim([-1, 3])
xlabel('Time from cue onset (s)')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',10,'FontWeight','Bold', 'LineWidth', 1);
fig_name = ['adult_high_gamma_mod_trace_high_gamma_sort_3_15_21'];
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
fig_name = ['adult_psth_trace_high_gamma_sort_3_15_21'];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
%%
peak_timepoint_gamma_adult(peak_timepoint_gamma_adult < 2) = [];
peak_timepoint_gamma_adult(peak_timepoint_gamma_adult >74) = [];
peak_timepoint_gamma_young(peak_timepoint_gamma_young < 2) = [];
peak_timepoint_gamma_young(peak_timepoint_gamma_young >74) = [];
peak_timepoint_hi_gamma_adult(peak_timepoint_hi_gamma_adult < 2) = [];
peak_timepoint_hi_gamma_adult(peak_timepoint_hi_gamma_adult > 74) = [];
peak_timepoint_hi_gamma_young(peak_timepoint_hi_gamma_young < 2) = [];
peak_timepoint_hi_gamma_young(peak_timepoint_hi_gamma_young > 74) = [];
%%
[h_, p_, stats] = kstest2(peak_timepoint_gamma_adult, peak_timepoint_gamma_young)
[h_, p_, stats] = kstest2(peak_timepoint_hi_gamma_adult, peak_timepoint_hi_gamma_young)
