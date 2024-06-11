close all
clear
zw_setpath
%%  Loading existing data
%   Loads: cwt trial-averaged cwt for each session X electrode (temp_cwt)
fname_ = 'temp_cwt.mat';
load(fullfile(project_dir, output_database, fname_), 'temp_cwt');
% 
%   Loads: 
fname_ = 'tfr_baseline.mat';
load(fullfile(project_dir, output_database, fname_), 'temp_baseline', 'temp_cwt_b');
%
%   Loads: SPIKING PEV calculated for each relevant neuron (neuron_pev_cue)
fname_ = 'neuron_pev_cue_ps.mat';
load(fullfile(project_dir, output_database, fname_), 'neuron_pev_cue');
% 
%   Loads: CWT PEV calculated for each site at 4 bands (lfp_pev_cue)
fname_ = 'lfp_pve_baseline_normalized.mat';
load(fullfile(project_dir, output_database, fname_), 'lfp_pev_cue');
% 
%   Loads: Trial-average PSTH (best_psth_raw, best_psth) for the best class
%   (neuron_best_class, neuron_best_class_raw) out of 8 that are either
%   before or after subtracting baseline; Also the trial-averge PSTH for
%   all classes (temp_psth, temp_psth_raw)
%   CREATED IN: COMPUTE_POPULATION_PSTH_12_4
fname_ = 'best_psth.mat';
load(fullfile(project_dir, output_database, fname_), ...
    'temp_psth', 'best_psth', 'neuron_best_class', ...
    'temp_psth_1s_baseline', 'best_psth_1s_baseline', 'neuron_best_class_1s_baseline', ...
    'temp_psth_raw', 'best_psth_raw', 'neuron_best_class_raw');% 
%   Loads: 
%%  Loading pre-computed flags
% 
%   Loads: 1) the matrix (mapping_mat) that matches LFP to SPIKING indices
%   and   2) the basic grouping infomation (Monkey X stage) for LFP (t) and
%   for   3) SPIKING (n)
fname_ = 'lfp_neuron_matching.mat';
load(fullfile(project_dir, output_database, fname_), 't', 'n', 'mapping_mat');
%
fname_ = 'tfr_delay_tuning_flag.mat';
load(fullfile(project_dir, output_database, fname_), 'tfr_delay_tuning_flag');
%
%   Loads: Band power significance increase from baseline in different
%   epochs ('tfr_cue_modulation_flag', 'tfr_delay_modulation_flag'
%   'tfr_total_modulation_flag' not computed
%   CREATED IN: WHOLE_DELAY_PERIOD_LFP_MODULATION
fname_ = 'tfr_modulation.mat';
load(fullfile(project_dir, output_database, fname_), 'tfr_cue_modulation_flag', 'tfr_delay_modulation_flag', 'tfr_modulation_output_cue');
% 
%   Loads: Neuron tuning flag in the delay period by ANOVA
%   (neuron_epoch_delay_tuning_flag)
%   Not loaded: Full stats ('neuron_epoch_analysis_output_cue')
fname_ = 'neuron_epoch_delay_tuning_flag.mat';
load(fullfile(project_dir, output_database, fname_), 'neuron_epoch_delay_tuning_flag');
% 
% 
%%  Variable cleaning
% Take care of accidental nan PEV
neuron_pev_cue(isnan(neuron_pev_cue)) = 0;
%%  Compute more needed variables
% Computes: Trial and class average LFP modulation for quicker plotting
lfp_mod_cue = zeros(size(lfp_pev_cue));
for i = 1:size(lfp_mod_cue, 3)
    lfp_mod_cue(:, :, i) = squeeze(nanmean(temp_cwt(:, target_frs(i, 1):target_frs(i, 2), :)./temp_baseline(:, target_frs(i, 1):target_frs(i, 2)), 2));
end
% 
% Computes the best class for delay period gamma power
gamma_power_best_class = zeros([size(tfr_modulation_output_cue, 1), 1]);
best_gamma_power = zeros([size(tfr_modulation_output_cue, 1), 1]);
for i = 1:numel(gamma_power_best_class)
    if ~isempty(tfr_modulation_output_cue(i, 3).mean)
        [y_, i_] = max(tfr_modulation_output_cue(i, 3).mean(2, :));
        gamma_power_best_class(i) = i_;
        best_gamma_power(i) = y_;
    end
end
zw_anova_from_cell({best_gamma_power([site_tuning_cat{1, [1,3,4], :, 1}])', ...
    best_gamma_power([site_tuning_cat{1, [1,3,4], :, 2}])'; ...
    best_gamma_power([site_tuning_cat{2, [1,3,4], :, 1}])', ...
    best_gamma_power([site_tuning_cat{2, [1,3,4], :, 2}])'})
% Conclusion: best class gamma power (instead of the cross-calss average)
% still shows the same trend (adolescent > adult; info > non-info sites)
%%  Extact lfp modulation corresponding to each neuron's preferred stimulus
lfp_mod_neuron_best = nan(size(neuron_best_class_1s_baseline, 1), size(lfp_mod_cue, 2), size(lfp_mod_cue, 3));
for i = [neuron_tuning_cat{:, [1,3,4], :}]
    target_class = neuron_best_class_1s_baseline(i);
    for j = 1:size(target_frs, 1)
        current_cwt_ = cwt_repo(map_neuron_to_site(mapping_mat, i)).class(target_class).cue_cwt;
        temp_mod_ = squeeze(nanmean(current_cwt_./nanmean(current_cwt_(:, :, 26:50), 3), 1));
        lfp_mod_neuron_best(i, :, j) = squeeze(nanmean(temp_mod_(target_frs(j, 1):target_frs(j, 2), :), 1));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Spectrogram for adolescent and adult respectively

clm_n = 100;
clm_ = [[zeros(1, clm_n);zeros(1, clm_n); linspace(1, 0, clm_n)]'; [linspace(0, 1, clm_n); zeros(1, clm_n);zeros(1, clm_n)]'];
c_range = [-.51, .51];
figure();
imagesc([-1, 3], [2,128], (squeeze(nanmean(temp_cwt_b([t{1, [1,3,4]}], :, :), 1)) - 1)*100, c_range*100);
set(gca, 'YDir','normal'); 
colormap(clm_)
xlabel('Time from cue onset (s)')
ylabel('Frequency (Hz)')
h = colorbar();
h.Ruler.TickLabelFormat='%g%%';
ylabel(h, '% Power of baseline')
set_plot_12_4();
title('Adolescent')
fig_name = 'spectrogram_change_young_draft'
add_epochline([1 1 1])
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
% 
figure();
imagesc([-1, 3], [2,128], (squeeze(nanmean(temp_cwt_b([t{2, [1,3,4]}], :, :), 1)) - 1)*100, c_range*100);
set(gca, 'YDir','normal'); 
colormap(clm_)
xlabel('Time from cue onset (s)')
ylabel('Frequency (Hz)')
h = colorbar();
h.Ruler.TickLabelFormat='%g%%';
ylabel(h, '% Power of baseline')
set_plot_12_4();
title('Adult')
fig_name = 'spectrogram_change_adult_draft'
add_epochline([1 1 1])
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
%%
%% LFP modulation
for i = 1:n_bands
% mod_plot({squeeze(nanmean(temp_cwt([t13;t14], target_frs(i, 1):target_frs(i, 2), :)./temp_baseline([t13;t14], target_frs(i, 1):target_frs(i, 2)), 2)), ...
%     squeeze(nanmean(temp_cwt([t23;t24], target_frs(i, 1):target_frs(i, 2), :)./temp_baseline([t23;t24], target_frs(i, 1):target_frs(i, 2)), 2))}, ...
%     linspace(-1,3,200), {'Adolescent', 'Adult'}, {'b' 'r'})
mod_plot({squeeze(nanmean(temp_cwt([t{1, [1,3,4]}], target_frs(i, 1):target_frs(i, 2), :)./temp_baseline([t{1, [1,3,4]}], target_frs(i, 1):target_frs(i, 2)), 2)), ...
    squeeze(nanmean(temp_cwt([t{2, [1,3,4]}], target_frs(i, 1):target_frs(i, 2), :)./temp_baseline([t{2, [1,3,4]}], target_frs(i, 1):target_frs(i, 2)), 2))}, ...
    linspace(-1,3,200), {'Adolescent', 'Adult'}, {'b' 'r'})
title(title_st{i})
xlabel('Time from cue onset (s)')
if i ~= 3
    legend off
end
fig_name = ['lfp_pow_band_', num2str(i), '_draft'];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
end
%%  High Gamma power by high gamma power modulation 
band_i_ = 4;
mod_plot({lfp_mod_cue([site_mod_cat{1, [1,3,4], 4}], :, band_i_), ...
    lfp_mod_cue([site_mod_cat{1, [1,3,4], 3}], :, band_i_), ...
    lfp_mod_cue([site_mod_cat{2, [1,3,4], 4}], :, band_i_), ...
    lfp_mod_cue([site_mod_cat{2, [1,3,4], 3}], :, band_i_)}, ...
    linspace(-1,3,200), ...
    {'Adolescent modulated', 'Adolescent non-modulated', 'Adult modulated', 'Adult non-modulated'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title(title_st{band_i_})
legend off
xlabel('Time from cue onset (s)')
fig_name = ['high_gamma_power_modulation_by_gamma_power_modulation_draft'];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
%%  Gamma power by gamma power modulation 
band_i_ = 3;
mod_plot({lfp_mod_cue([site_mod_cat{1, [1,3,4], 2}], :, band_i_), ...
    lfp_mod_cue([site_mod_cat{1, [1,3,4], 1}], :, band_i_), ...
    lfp_mod_cue([site_mod_cat{2, [1,3,4], 2}], :, band_i_), ...
    lfp_mod_cue([site_mod_cat{2, [1,3,4], 1}], :, band_i_)}, ...
    linspace(-1,3,200), ...
    {'Adolescent gamma-modulated', 'Adolescent non-gamma-modulated', 'Adult gamma-modulated', 'Adult non-gamma-modulated'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title(title_st{band_i_})
xlabel('Time from cue onset (s)')
fig_name = ['gamma_power_modulation_by_gamma_power_modulation_draft'];
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
%%  Gamma power modulation by spiking tuning
for band_i_ = 1:4
mod_plot({lfp_mod_cue([site_tuning_cat{1, [1,3,4], :, 2}], :, band_i_), ...
    lfp_mod_cue([site_tuning_cat{1, [1,3,4], :, 1}], :, band_i_), ...
    lfp_mod_cue([site_tuning_cat{2, [1,3,4], :, 2}], :, band_i_), ...
    lfp_mod_cue([site_tuning_cat{2, [1,3,4], :, 1}], :, band_i_)}, ...
    linspace(-1,3,200), ...
    {'Adolescent informative', 'Adolescent non-informative', 'Adult informative', 'Adult non-informative'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title(title_st{band_i_})
xlabel('Time from cue onset (s)')
if band_i_ ~= 2
    legend off
end
fig_name = ['lfp_pow_band_', num2str(band_i_), 'power_modulation_by_spiking_tuning', '_draft'];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
end
%%  High Gamma power modulation by spiking tuning
band_i_ = 4;
mod_plot({lfp_mod_cue([site_tuning_cat{1, [1,3,4], :, 2}], :, band_i_), ...
    lfp_mod_cue([site_tuning_cat{1, [1,3,4], :, 1}], :, band_i_), ...
    lfp_mod_cue([site_tuning_cat{2, [1,3,4], :, 2}], :, band_i_), ...
    lfp_mod_cue([site_tuning_cat{2, [1,3,4], :, 1}], :, band_i_)}, ...
    linspace(-1,3,200), ...
    {'Adolescent informative', 'Adolescent non-informative', 'Adult informative', 'Adult non-informative'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title('High Gamma Power')
xlabel('Time from cue onset (s)')
fig_name = ['high_gamma_power_modulation_by_spiking_tuning_draft'];
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
%%
%%  Spiking firing rate raw
psth_plot({best_psth_raw([neuron_tuning_cat{1, [1,3,4], :, :}], :), ...
    best_psth_raw([neuron_tuning_cat{2, [1,3,4], :, :}], :)}, ...
    linspace(-1,3,200), ...
    {'Adolescent', 'Adult'}, {'b','r'})
title('Absolute Firing Rate')
xlabel('Time from cue onset (s)')
fig_name = ['raw_spiking_firing_draft'];
legend off
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');

%%  Spiking firing rate (1000ms baseline)
psth_plot({best_psth_1s_baseline([neuron_tuning_cat{1, [1,3,4], :, :}], :), ...
    best_psth_1s_baseline([neuron_tuning_cat{2, [1,3,4], :, :}], :)}, ...
    linspace(-1,3,200), ...
    {'Adolescent', 'Adult'}, {'b','r'})
title('Evoked Firing Rate')
xlabel('Time from cue onset (s)')
fig_name = ['evoked_spiking_firing_draft'];
legend off
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');

%%
%%  Spiking firing rate (1000ms baseline) by spiking tuning
psth_plot({best_psth_1s_baseline([neuron_tuning_cat{1, [1,3,4], :, 2}], :), ...
    best_psth_1s_baseline([neuron_tuning_cat{1, [1,3,4], :, 1}], :), ...
    best_psth_1s_baseline([neuron_tuning_cat{2, [1,3,4], :, 2}], :), ...
    best_psth_1s_baseline([neuron_tuning_cat{2, [1,3,4], :, 1}], :)}, ...
    linspace(-1,3,200), ...
    {'Adolescent informative', 'Adolescent non-informative', 'Adult informative', 'Adult non-informative'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title('Evoked Firing Rate')
xlabel('Time from cue onset (s)')
fig_name = ['spiking_firing_rate_by_spiking_tuning_draft'];
legend off
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
%%  Neuron PEV by spiking tuning
pev_plot({neuron_pev_cue([neuron_tuning_cat{1, [1,3,4], :, 2}], :), ...
    neuron_pev_cue([neuron_tuning_cat{1, [1,3,4], :, 1}], :), ...
    neuron_pev_cue([neuron_tuning_cat{2, [1,3,4], :, 2}], :), ...
    neuron_pev_cue([neuron_tuning_cat{2, [1,3,4], :, 1}], :)}, ...
    linspace(-1,3,200), ...
    {'Adolescent informative', 'Adolescent non-informative', 'Adult informative', 'Adult non-informative'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title('Spiking PEV')
xlabel('Time from cue onset (s)')
fig_name = ['neuron_pev_by_spiking_tuning_draft'];
legend off
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
%%
%%  Sort by neuron PEV peaks, then plot 1) Norm. PEV; 2) Norm. gamma power; 3) Norm. high-gamma power
target_subset = [neuron_tuning_cat{2, [1,3,4],:, 2}];
min_maxed_pev_ = (neuron_pev_cue(target_subset, :) - min(neuron_pev_cue(target_subset, :), [], 2))./(max(neuron_pev_cue(target_subset, :), [], 2) - min(neuron_pev_cue(target_subset, :), [], 2));

peak_timepoint = zeros(numel(target_subset), 1);
for i = 1:numel(target_subset)
    [~, i_] = max(min_maxed_pev_(i, 76:200));
    peak_timepoint(i) = i_;
%         peak_timepoint(i) = zw_center_of_mass(min_maxed_pev_(i, 76:150));

end
[sorted_peak_timepoint, trace_order] = sort(peak_timepoint);


figure
imagesc([-1, 3], 1:numel(trace_order), min_maxed_pev_(trace_order, :), [0, 1]);
% colormap('jet')
h = colorbar;
ylabel(h, 'Norm. Spiking PEV');
hold on
% plot(sorted_peak_timepoint+75, 1:numel(target_subset), 'r');

% 
ylabel('Neuron #')
xlabel('Time from cue onset (s)')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1);
fig_name = ['adult_pev_trace_draft'];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
band_i = 3;
mod_to_plot = lfp_mod_cue(map_neuron_to_site(mapping_mat, target_subset(trace_order)), :, band_i) - 1; %    Converting to change from baseline
mod_to_plot = mod_to_plot./max(mod_to_plot, [], 2);
figure
imagesc([-1, 3], 1:numel(trace_order), mod_to_plot, [-1, 1])
colormap(clm_)
h = colorbar;
ylabel(h, 'Norm. Gamma Power');
ylabel('Neuron #')
xlabel('Time from cue onset (s)')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1);
fig_name = ['adult_gamma_mod_trace_draft'];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');

%
band_i = 4;
mod_to_plot = lfp_mod_cue(map_neuron_to_site(mapping_mat, target_subset(trace_order)), :, band_i) - 1; %    Converting to change from baseline
mod_to_plot = mod_to_plot./max(mod_to_plot, [], 2);
figure
imagesc([-1, 3], 1:numel(trace_order), mod_to_plot, [-1, 1])
colormap(clm_)
h = colorbar;
ylabel(h, 'Norm. High Gamma Power');
ylabel('Neuron #')
xlabel('Time from cue onset (s)')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1);
fig_name = ['adult_high_gamma_mod_trace_draft'];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');