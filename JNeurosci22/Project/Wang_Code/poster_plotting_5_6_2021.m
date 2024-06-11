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
% fname_ = 'neuron_pev_cue_ps.mat';
fname_ = 'neuron_pev_ps_upsampled_3_11_2021.mat';
load(fullfile(project_dir, output_database, fname_), 'neuron_pev_cue');
% 
%   Loads: CWT PEV calculated for each site at 4 bands (lfp_pev_cue)
fname_ = 'lfp_pve_baseline_normalized_dur_updated.mat';
load(fullfile(project_dir, output_database, fname_), 'lfp_pev_cue');
% 
%   Loads: Trial-average PSTH (best_psth_raw, best_psth) for the best class
%   (neuron_best_class, neuron_best_class_raw) out of 8 that are either
%   before or after subtracting baseline; Also the trial-averge PSTH for
%   all classes (temp_psth, temp_psth_raw)
%   CREATED IN: COMPUTE_POPULATION_PSTH_12_4
fname_ = 'best_psth_upsampled_3_11_2021.mat';
load(fullfile(project_dir, output_database, fname_), ...
    'temp_psth', 'best_psth', 'neuron_best_class', ...
    'temp_psth_1s_baseline', 'best_psth_1s_baseline', 'neuron_best_class_1s_baseline', ...
    'temp_psth_raw', 'best_psth_raw', 'neuron_best_class_raw');% 
%   Loads within-session averages of cwt 
fname_ = 'temp_cwt_3_2_2021.mat';
load(fullfile(project_dir, output_database, fname_), 'temp_baseline', 'temp_cwt_cue_b', 'temp_cwt_sac_b', 'temp_cwt_cue', 'temp_cwt_sac');
%   Loads behavior summary
fname_ = 'odr_behavior_repo';
load(fullfile(project_dir, output_database, fname_), 'behavior_repo');
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
%   'tfr_modulation_output_cue' not loaded
%   'tfr_total_modulation_flag' not computed
%   CREATED IN: WHOLE_DELAY_PERIOD_LFP_MODULATION
fname_ = 'tfr_modulation.mat';
load(fullfile(project_dir, output_database, fname_), 'tfr_cue_modulation_flag', 'tfr_delay_modulation_flag');
% 
%   Loads: Neuron tuning flag in the delay period by ANOVA
%   (neuron_epoch_delay_tuning_flag)
%   Not loaded: Full stats ('neuron_epoch_analysis_output_cue')
fname_ = 'neuron_epoch_delay_tuning_flag.mat';
load(fullfile(project_dir, output_database, fname_), 'neuron_epoch_delay_tuning_flag');
%
fname_ = 'tfr_modulation.mat';
load(fullfile(project_dir, output_database, fname_), 'tfr_cue_modulation_flag', 'tfr_delay_modulation_flag', 'tfr_modulation_output_cue');
% 
% 
%%  Variable cleaning
% Take care of accidental nan PEV
neuron_pev_cue(isnan(neuron_pev_cue)) = 0;
%%  Compute more needed variables
% Computes: Trial and class average LFP modulation for quicker plotting
% lfp_mod_cue = zeros(size(lfp_pev_cue));
n_bands = 4;
target_frs = [4,8; 8,16; 16, 32; 32, 64];
lfp_mod_cue = zeros(size(temp_cwt_cue, 1), size(temp_cwt_cue, 3),  n_bands);
lfp_mod_sac = zeros(size(temp_cwt_sac, 1), size(temp_cwt_sac, 3),  n_bands);
for i = 1:size(lfp_mod_cue, 3)
    lfp_mod_cue(:, :, i) = squeeze(nanmean(temp_cwt_cue(:, target_frs(i, 1):target_frs(i, 2), :)./temp_baseline(:, target_frs(i, 1):target_frs(i, 2)), 2));
    lfp_mod_sac(:, :, i) = squeeze(nanmean(temp_cwt_sac(:, target_frs(i, 1):target_frs(i, 2), :)./temp_baseline(:, target_frs(i, 1):target_frs(i, 2)), 2));
end
%%
% Computes the best class for delay period beta power decrease
power_best_class = zeros([size(tfr_modulation_output_cue, 1), 1]);
best_power = zeros([size(tfr_modulation_output_cue, 1), 1]);
for i = 1:numel(gamma_power_best_class)
    if ~isempty(tfr_modulation_output_cue(i, 2).mean)
        [y_, i_] = min(tfr_modulation_output_cue(i, 3).mean(2, :)./tfr_modulation_output_cue(i, 2).mean(2, :));
        power_best_class(i) = i_;
        best_power(i) = y_;
    end
end
[~, i_] = max(best_power([t{:}]));
temp_ttt_(i_)
best_power(temp_ttt_(i_))
%%
% Computes the best class for delay period gamma/beta power decrease
beta_power_best_class = zeros([size(tfr_modulation_output_cue, 1), 1]);
best_beta_power = zeros([size(tfr_modulation_output_cue, 1), 1]);
for i = 1:numel(gamma_power_best_class)
    if ~isempty(tfr_modulation_output_cue(i, 2).mean)
        [y_, i_] = min(tfr_modulation_output_cue(i, 2).mean(2, :));
        beta_power_best_class(i) = i_;
        best_beta_power(i) = y_;
    end
end
%%
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
%%  Spectrogram for adolescent and adult respectively (cue aligned)

clm_n = 100;
clm_ = [[zeros(1, clm_n);zeros(1, clm_n); linspace(1, 0, clm_n)]'; [linspace(0, 1, clm_n); zeros(1, clm_n);zeros(1, clm_n)]'];
c_range = [-.51, .51];
figure('Unit', 'inches', 'Position', [2, 2, 4, 3.2]);
imagesc([-1, 4], [2,128], (squeeze(nanmean(temp_cwt_cue_b([t{1, [1,3,4]}], :, :), 1)) - 1)*100, c_range*100);
set(gca, 'YDir','normal'); 
colormap(clm_)
xlabel('Time from cue onset (s)')
ylabel('Frequency (Hz)')
h = colorbar();
h.Ruler.TickLabelFormat='%g%%';
ylabel(h, '% Power of baseline')
set_plot_poster();
title('Adolescent')
fig_name = 'spectrogram_change_young_cue_draft'
add_epochline([1 1 1])
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
% savefig(gcf, fullfile(project_dir, fig_lib, fig_name));

% 
figure('Unit', 'inches', 'Position', [2, 2, 4, 3.2]);
imagesc([-1, 4], [2,128], (squeeze(nanmean(temp_cwt_cue_b([t{2, [1,3,4]}], :, :), 1)) - 1)*100, c_range*100);
set(gca, 'YDir','normal'); 
colormap(clm_)
xlabel('Time from cue onset (s)')
ylabel('Frequency (Hz)')
h = colorbar();
h.Ruler.TickLabelFormat='%g%%';
ylabel(h, '% Power of baseline')
set_plot_poster();
title('Adult')
fig_name = 'spectrogram_change_adult_cue_draft'
add_epochline([1 1 1])
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
% savefig(gcf, fullfile(project_dir, fig_lib, fig_name));
%%
%% LFP modulation (cue aligned)
n_bands = 4;
title_st = {'Alpha (8-16 Hz)', 'Beta (16-32 Hz)', 'Gamma (32-64 Hz)', 'High-Gamma (64-128 Hz)'};
for i = 1:n_bands
% mod_plot({squeeze(nanmean(temp_cwt([t13;t14], target_frs(i, 1):target_frs(i, 2), :)./temp_baseline([t13;t14], target_frs(i, 1):target_frs(i, 2)), 2)), ...
%     squeeze(nanmean(temp_cwt([t23;t24], target_frs(i, 1):target_frs(i, 2), :)./temp_baseline([t23;t24], target_frs(i, 1):target_frs(i, 2)), 2))}, ...
%     linspace(-1,4,250), {'Adolescent', 'Adult'}, {'b' 'r'})
mod_plot({squeeze(nanmean(temp_cwt_cue([t{1, [1,3,4]}], target_frs(i, 1):target_frs(i, 2), :)./temp_baseline([t{1, [1,3,4]}], target_frs(i, 1):target_frs(i, 2)), 2)), ...
    squeeze(nanmean(temp_cwt_cue([t{2, [1,3,4]}], target_frs(i, 1):target_frs(i, 2), :)./temp_baseline([t{2, [1,3,4]}], target_frs(i, 1):target_frs(i, 2)), 2))}, ...
    linspace(-1,4,250), {'Adolescent', 'Adult'}, {'b' 'r'}, [], 12, [3.2, 3.2])
title(title_st{i})
xlabel('Time from cue onset (s)')
if i ~= 3
    legend off
end
fig_name = ['lfp_pow_band_', num2str(i), '_cue_draft'];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
savefig(gcf, fullfile(project_dir, fig_lib, fig_name));
end
%%  High Gamma power by high gamma power modulation (cue aligned)
band_i_ = 4;
mod_plot({lfp_mod_cue([site_mod_cat{1, [1,3,4], 4}], :, band_i_), ...
    lfp_mod_cue([site_mod_cat{1, [1,3,4], 3}], :, band_i_), ...
    lfp_mod_cue([site_mod_cat{2, [1,3,4], 4}], :, band_i_), ...
    lfp_mod_cue([site_mod_cat{2, [1,3,4], 3}], :, band_i_)}, ...
    linspace(-1,4,250), ...
    {'Ado mod', 'Ado non-mod', 'Adu mod', 'Adu non-mod'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]}, ...
    [], 12, [3.2, 3.2])
%     {'Adolescent modulated', 'Adolescent non-modulated', 'Adult modulated', 'Adult non-modulated'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title(title_st{band_i_})
% legend off
xlabel('Time from cue onset (s)')
fig_name = ['high_gamma_power_modulation_by_gamma_power_modulation_cue_draft'];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
savefig(gcf, fullfile(project_dir, fig_lib, fig_name));
%%  Gamma power by gamma power modulation (cue aligned)
band_i_ = 3;
mod_plot({lfp_mod_cue([site_mod_cat{1, [1,3,4], 2}], :, band_i_), ...
    lfp_mod_cue([site_mod_cat{1, [1,3,4], 1}], :, band_i_), ...
    lfp_mod_cue([site_mod_cat{2, [1,3,4], 2}], :, band_i_), ...
    lfp_mod_cue([site_mod_cat{2, [1,3,4], 1}], :, band_i_)}, ...
    linspace(-1,4,250), ...
    {'Ado mod', 'Ado non-mod', 'Adu mod', 'Adu non-mod'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]}, ...
    [], 12, [3.2, 3.2])
%     {'Adolescent modulated', 'Adolescent non-modulated', 'Adult modulated', 'Adult non-modulated'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title(title_st{band_i_})
xlabel('Time from cue onset (s)')
fig_name = ['gamma_power_modulation_by_gamma_power_modulation_cue_draft'];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
savefig(gcf, fullfile(project_dir, fig_lib, fig_name));
% 
%%  High Gamma power by high gamma power modulation (sac aligned)
band_i_ = 4;
mod_plot_sac({lfp_mod_sac([site_mod_cat{1, [1,3,4], 4}], :, band_i_), ...
    lfp_mod_sac([site_mod_cat{1, [1,3,4], 3}], :, band_i_), ...
    lfp_mod_sac([site_mod_cat{2, [1,3,4], 4}], :, band_i_), ...
    lfp_mod_sac([site_mod_cat{2, [1,3,4], 3}], :, band_i_)}, ...
    linspace(-2,3,250), ...
    {'Ado mod', 'Ado non-mod', 'Adu mod', 'Adu non-mod'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]}, ...
    [], 12, [3.2, 3.2])
%     {'Adolescent modulated', 'Adolescent non-modulated', 'Adult modulated', 'Adult non-modulated'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title(title_st{band_i_})
% legend off
xlabel('Time from saccadic onset (s)')
fig_name = ['high_gamma_power_modulation_by_gamma_power_modulation_sac_draft'];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
savefig(gcf, fullfile(project_dir, fig_lib, fig_name));
%%  Gamma power by gamma power modulation (sac aligned)
band_i_ = 3;
mod_plot_sac({lfp_mod_sac([site_mod_cat{1, [1,3,4], 2}], :, band_i_), ...
    lfp_mod_sac([site_mod_cat{1, [1,3,4], 1}], :, band_i_), ...
    lfp_mod_sac([site_mod_cat{2, [1,3,4], 2}], :, band_i_), ...
    lfp_mod_sac([site_mod_cat{2, [1,3,4], 1}], :, band_i_)}, ...
    linspace(-2,3,250), ...
    {'Ado mod', 'Ado non-mod', 'Adu mod', 'Adu non-mod'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]}, ...
    [], 12, [3.2, 3.2])
%     {'Adolescent modulated', 'Adolescent non-modulated', 'Adult modulated', 'Adult non-modulated'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title(title_st{band_i_})
xlabel('Time from saccadic onset (s)')
fig_name = ['gamma_power_modulation_by_gamma_power_modulation_sac_draft'];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
savefig(gcf, fullfile(project_dir, fig_lib, fig_name));
%%  LFP power modulation by spiking tuning (cue aligned)
for band_i_ = 1:4
mod_plot({lfp_mod_cue([site_tuning_cat{1, [1,3,4], :, 2}], :, band_i_), ...
    lfp_mod_cue([site_tuning_cat{1, [1,3,4], :, 1}], :, band_i_), ...
    lfp_mod_cue([site_tuning_cat{2, [1,3,4], :, 2}], :, band_i_), ...
    lfp_mod_cue([site_tuning_cat{2, [1,3,4], :, 1}], :, band_i_)}, ...
    linspace(-1,4,250), ...
    {'Adolescent inf', 'Adolescent non-inf', 'Adult inf', 'Adult non-inf'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]}, ...
    [], 12, [3.2, 3.2])
%     {'Adolescent informative', 'Adolescent non-informative', 'Adult informative', 'Adult non-informative'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title(title_st{band_i_})
xlabel('Time from cue onset (s)')
if band_i_ ~= 2
    legend off
end
fig_name = ['lfp_pow_band_', num2str(band_i_), 'power_modulation_by_spiking_tuning', '_cue_draft'];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
savefig(gcf, fullfile(project_dir, fig_lib, fig_name));
end
%%  Spiking firing rate raw
psth_plot({best_psth_raw([neuron_tuning_cat{1, [1,3,4], :, :}], :), ...
    best_psth_raw([neuron_tuning_cat{2, [1,3,4], :, :}], :)}, ...
    linspace(-1,4,250), ...
    {'Adolescent', 'Adult'}, {'b','r'}, ...
    [], 12, [3.2, 3.2])
title('Absolute Firing Rate')
xlabel('Time from cue onset (s)')
fig_name = ['raw_spiking_firing_draft'];
legend off
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
% savefig(gcf, fullfile(project_dir, fig_lib, fig_name));
%%  Spiking firing rate (1000ms baseline)
psth_plot({best_psth_1s_baseline([neuron_tuning_cat{1, [1,3,4], :, :}], :), ...
    best_psth_1s_baseline([neuron_tuning_cat{2, [1,3,4], :, :}], :)}, ...
    linspace(-1,4,250), ...
    {'Adolescent', 'Adult'}, {'b','r'}, ...
        [], 12, [3.2, 3.2])
title('Evoked Firing Rate')
xlabel('Time from cue onset (s)')
fig_name = ['evoked_spiking_firing_draft'];
% legend off
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
savefig(gcf, fullfile(project_dir, fig_lib, fig_name));
%%  Spiking firing rate (1000ms baseline) by spiking tuning
psth_plot({best_psth_1s_baseline([neuron_tuning_cat{1, [1,3,4], :, 2}], :), ...
    best_psth_1s_baseline([neuron_tuning_cat{1, [1,3,4], :, 1}], :), ...
    best_psth_1s_baseline([neuron_tuning_cat{2, [1,3,4], :, 2}], :), ...
    best_psth_1s_baseline([neuron_tuning_cat{2, [1,3,4], :, 1}], :)}, ...
    linspace(-1,4,250), ...
    {'Adolescent in', 'Adolescent non-in', 'Adult in', 'Adult non-in'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]}, ...
    [], 12, [3.2, 3.2])
%     {'Adolescent informative', 'Adolescent non-informative', 'Adult informative', 'Adult non-informative'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title('Evoked Firing Rate')
xlabel('Time from cue onset (s)')
fig_name = ['spiking_firing_rate_by_spiking_tuning_draft'];
% legend off
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
% savefig(gcf, fullfile(project_dir, fig_lib, fig_name));
%%  Neuron PEV by spiking tuning
pev_plot({neuron_pev_cue([neuron_tuning_cat{1, [1,3,4], :, 2}], :), ...
    neuron_pev_cue([neuron_tuning_cat{1, [1,3,4], :, 1}], :), ...
    neuron_pev_cue([neuron_tuning_cat{2, [1,3,4], :, 2}], :), ...
    neuron_pev_cue([neuron_tuning_cat{2, [1,3,4], :, 1}], :)}, ...
    linspace(-1,4,250), ...
    {'Adolescent informative', 'Adolescent non-informative', 'Adult informative', 'Adult non-informative'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]}, ...
    [0, .26], 12, [3.2, 3.2])
title('Spiking PEV')
xlabel('Time from cue onset (s)')
fig_name = ['neuron_pev_by_spiking_tuning_draft'];
legend off
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
savefig(gcf, fullfile(project_dir, fig_lib, fig_name));%
%% LFP PEV (cue aligned)
n_bands = 4;
title_st = {'Alpha (8-16 Hz)', 'Beta (16-32 Hz)', 'Gamma (32-64 Hz)', 'High-Gamma (64-128 Hz)'};
for i = 1:n_bands
pev_plot({lfp_pev_cue([t{1, :}], :, i), ...
    lfp_pev_cue([t{2, :}], :, i)}, ...
    linspace(-1,4,250), {'Adolescent', 'Adult'}, {'b' 'r'}, [], 12, [3.2, 3.2])
title(title_st{i})
xlabel('Time from cue onset (s)')
if i ~= 3
    legend off
end
fig_name = ['lfp_pev_band_', num2str(i), '_cue_draft'];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
savefig(gcf, fullfile(project_dir, fig_lib, fig_name));
end
%% LFP PEV by spiking tuning (cue aligned)
n_bands = 4;
title_st = {'Alpha (8-16 Hz)', 'Beta (16-32 Hz)', 'Gamma (32-64 Hz)', 'High-Gamma (64-128 Hz)'};
for i = 1:n_bands
pev_plot({lfp_pev_cue([site_tuning_cat{1, [1,3,4], :, 2}], :, i), ...
    lfp_pev_cue([site_tuning_cat{1, [1,3,4], :, 1}], :, i), ...
    lfp_pev_cue([site_tuning_cat{2, [1,3,4], :, 2}], :, i), ...
    lfp_pev_cue([site_tuning_cat{2, [1,3,4], :, 1}], :, i)}, ...
    linspace(-1,4,250), ...
    {'Adolescent informative', 'Adolescent non-informative', 'Adult informative', 'Adult non-informative'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]}, ...
    [0, .26], 12, [3.2, 3.2])
title(title_st{i})
xlabel('Time from cue onset (s)')
if i ~= 3
    legend off
end
fig_name = ['lfp_pev_by_sp_info_band_', num2str(i), '_cue_draft'];
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
% savefig(gcf, fullfile(project_dir, fig_lib, fig_name));
end
%%
beh_for_plot = cell(2, 4);
for i = 1:2
    for j = 1:4
        current_t = [t{i, j}];
        [iu, ia,~] = unique({behavior_repo(current_t).Filename});
        numel(iu)
        numel(ia)        
        beh_for_plot{i, j} = current_t(ia);
    end
end
%%
mean_s_rate = zeros(size(beh_for_plot));
sem_s_rate  = zeros(size(beh_for_plot));
for i = 1:2
    for j = 1:4
mean_s_rate(i, j) = mean([behavior_repo([beh_for_plot{i,j}]).s_rate])*100;
sem_s_rate(i, j) = std([behavior_repo([beh_for_plot{i,j}]).s_rate])/numel([beh_for_plot{i,j}]);
sem_s_rate(i, j) = std([behavior_repo([beh_for_plot{i,j}]).s_rate])*100;
    end
end
%%
x = [1,5,9]' + [-.8, .8];
x1 = [1,5,9]' + [-.8, .8];

% x = [1:3]' + [0, 0];

figure('Unit', 'inches', 'Position', [2, 2, 3.2, 3.2]);
hold on
ytickformat('percentage')
bt(1) = bar(x(:, 1), mean_s_rate(1, [1, 3, 4])', 'BarWidth', .3, 'LineWidth', 1.2)
bt(2) = bar(x(:, 2), mean_s_rate(2, [1, 3, 4])', 'BarWidth', .3, 'LineWidth', 1.2)
mr = reshape(mean_s_rate(:, [1, 3, 4])', 1, numel(mean_s_rate(:, [1, 3, 4])'));
mx = reshape(x1, 1, numel(x1));
ms = reshape(sem_s_rate(:, [1, 3, 4])', 1, numel(sem_s_rate(:, [1, 3, 4])'));
errorbar(mx, mr, -ms, ms, '.k', 'LineWidth', 1.2)
bt(1).FaceColor = 'b';
bt(2).FaceColor = 'r';
ylim([0, 105]);
xticks([1, 5, 9]);
xticklabels({'I', 'L', 'K'})
set(gca, {'FontSize'}, {12})
fig_name = ['beh_draft'];
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
% savefig(gcf, fullfile(project_dir, fig_lib, fig_name));