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
%   Loads: 
fname_ = 'temp_cwt_3_2_2021.mat';
load(fullfile(project_dir, output_database, fname_), 'temp_baseline', 'temp_cwt_cue_b', 'temp_cwt_sac_b', 'temp_cwt_cue', 'temp_cwt_sac');
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
%
fname_ = 'tfr_modulation.mat';
load(fullfile(project_dir, output_database, fname_), 'tfr_total_modulation_flag', 'tfr_cue_modulation_flag', 'tfr_delay_modulation_flag', 'tfr_modulation_output_cue');
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
fname_ = 'neuron_site_categories.mat';
load(fullfile(project_dir, output_database, fname_), 'site_mod_cat', 'site_tuning_cat', 'neuron_mod_cat', 'neuron_tuning_cat', 'mapping_mat', 'responsive');
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
%%
fname_ = 'gamma_power_best_class';
save(fullfile(project_dir, output_database, fname_), 'best_gamma_power', 'gamma_power_best_class');
%%
% Computes the best class for delay period high gamma power
high_gamma_power_best_class = zeros([size(tfr_modulation_output_cue, 1), 1]);
best_high_gamma_power = zeros([size(tfr_modulation_output_cue, 1), 1]);
for i = 1:numel(high_gamma_power_best_class)
    if ~isempty(tfr_modulation_output_cue(i, 4).mean)
        [y_, i_] = max(tfr_modulation_output_cue(i, 4).mean(2, :));
        high_gamma_power_best_class(i) = i_;
        best_high_gamma_power(i) = y_;
    end
end
%%
fname_ = 'high_gamma_power_best_class';
save(fullfile(project_dir, output_database, fname_), 'best_high_gamma_power', 'high_gamma_power_best_class');

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
figure('Unit', 'inches', 'Position', [2, 2, 3, 2.5]);
imagesc([-1, 4], [2,128], (squeeze(nanmean(temp_cwt_cue_b([t{1, [1,3,4]}], :, :), 1)) - 1)*100, c_range*100);
set(gca, 'YDir','normal'); 
colormap(clm_)
xlabel('Time from cue onset (s)')
ylabel('Frequency (Hz)')
h = colorbar();
h.Ruler.TickLabelFormat='%g%%';
ylabel(h, '% Power of baseline')
set_plot_3_4_2021();
title('Adolescent')
fig_name = 'spectrogram_change_young_cue_draft_-1000_2500'
add_epochline([1 1 1])
xlim([-1, 2.5])
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
savefig(gcf, fullfile(project_dir, fig_lib, fig_name));
% 
figure('Unit', 'inches', 'Position', [2, 2, 3, 2.5]);
imagesc([-1, 4], [2,128], (squeeze(nanmean(temp_cwt_cue_b([t{2, [1,3,4]}], :, :), 1)) - 1)*100, c_range*100);
set(gca, 'YDir','normal'); 
colormap(clm_)
xlabel('Time from cue onset (s)')
ylabel('Frequency (Hz)')
h = colorbar();
h.Ruler.TickLabelFormat='%g%%';
ylabel(h, '% Power of baseline')
set_plot_3_4_2021();
title('Adult')
fig_name = 'spectrogram_change_adult_cue_draft_-1000_2500'
add_epochline([1 1 1])
xlim([-1, 2.5])
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
savefig(gcf, fullfile(project_dir, fig_lib, fig_name));
%%  Spectrogram for adolescent and adult respectively (saccade aligned)
clm_n = 100;
clm_ = [[zeros(1, clm_n);zeros(1, clm_n); linspace(1, 0, clm_n)]'; [linspace(0, 1, clm_n); zeros(1, clm_n);zeros(1, clm_n)]'];
c_range = [-.51, .51];
figure('Unit', 'inches', 'Position', [2, 2, 3, 2.5]);
imagesc([-2, 3], [2,128], (squeeze(nanmean(temp_cwt_sac_b([t{1, [1,3,4]}], :, :), 1)) - 1)*100, c_range*100);
set(gca, 'YDir','normal'); 
colormap(clm_)
xlabel('Time from saccadic onset (s)')
ylabel('Frequency (Hz)')
h = colorbar();
h.Ruler.TickLabelFormat='%g%%';
ylabel(h, '% Power of baseline')
set_plot_3_4_2021();
title('Adolescent')
fig_name = 'spectrogram_change_young_sac_draft_-1000_2500'
add_sacline([1 1 1])
xlim([-1, 2.5])
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
savefig(gcf, fullfile(project_dir, fig_lib, fig_name));
% 
figure('Unit', 'inches', 'Position', [2, 2, 3, 2.5]);
imagesc([-2, 3], [2,128], (squeeze(nanmean(temp_cwt_sac_b([t{2, [1,3,4]}], :, :), 1)) - 1)*100, c_range*100);
set(gca, 'YDir','normal'); 
colormap(clm_)
xlabel('Time from saccadic onset (s)')
ylabel('Frequency (Hz)')
h = colorbar();
h.Ruler.TickLabelFormat='%g%%';
ylabel(h, '% Power of baseline')
set_plot_3_4_2021();
title('Adult')
fig_name = 'spectrogram_change_adult_sac_draft_-1000_2500'
add_sacline([1 1 1])
xlim([-1, 2.5])
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
savefig(gcf, fullfile(project_dir, fig_lib, fig_name));
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
    linspace(-1,4,250), {'Adolescent', 'Adult'}, {'b' 'r'})
xlim([-1, 2.5])
title(title_st{i})
xlabel('Time from cue onset (s)')
if i ~= 3
    legend off
end
fig_name = ['lfp_pow_band_', num2str(i), '_cue_draft_-1000_2500'];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
savefig(gcf, fullfile(project_dir, fig_lib, fig_name));
end
%%
%% LFP modulation (saccade aligned)
n_bands = 4;
title_st = {'Alpha (8-16 Hz)', 'Beta (16-32 Hz)', 'Gamma (32-64 Hz)', 'High-Gamma (64-128 Hz)'};
for i = 1:n_bands
% mod_plot({squeeze(nanmean(temp_cwt([t13;t14], target_frs(i, 1):target_frs(i, 2), :)./temp_baseline([t13;t14], target_frs(i, 1):target_frs(i, 2)), 2)), ...
%     squeeze(nanmean(temp_cwt([t23;t24], target_frs(i, 1):target_frs(i, 2), :)./temp_baseline([t23;t24], target_frs(i, 1):target_frs(i, 2)), 2))}, ...
%     linspace(-1,4,250), {'Adolescent', 'Adult'}, {'b' 'r'})
mod_plot_sac({squeeze(nanmean(temp_cwt_sac([t{1, [1,3,4]}], target_frs(i, 1):target_frs(i, 2), :)./temp_baseline([t{1, [1,3,4]}], target_frs(i, 1):target_frs(i, 2)), 2)), ...
    squeeze(nanmean(temp_cwt_sac([t{2, [1,3,4]}], target_frs(i, 1):target_frs(i, 2), :)./temp_baseline([t{2, [1,3,4]}], target_frs(i, 1):target_frs(i, 2)), 2))}, ...
    linspace(-2,3,250), {'Adolescent', 'Adult'}, {'b' 'r'})
title(title_st{i})
xlabel('Time from saccadic onset (s)')
if i ~= 3
    legend off
end
fig_name = ['lfp_pow_band_', num2str(i), '_sac_draft'];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
end
%%  High Gamma power by high gamma power modulation (cue aligned)
band_i_ = 4;
mod_plot({lfp_mod_cue([site_mod_cat{1, [1,3,4], 4}], :, band_i_), ...
    lfp_mod_cue([site_mod_cat{1, [1,3,4], 3}], :, band_i_), ...
    lfp_mod_cue([site_mod_cat{2, [1,3,4], 4}], :, band_i_), ...
    lfp_mod_cue([site_mod_cat{2, [1,3,4], 3}], :, band_i_)}, ...
    linspace(-1,4,250), ...
    {'Ado mod', 'Ado non-mod', 'Adu mod', 'Adu non-mod'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
%     {'Adolescent modulated', 'Adolescent non-modulated', 'Adult modulated', 'Adult non-modulated'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
xlim([-1, 2.5]);
title(title_st{band_i_})
% legend off
xlabel('Time from cue onset (s)')
fig_name = ['high_gamma_power_modulation_by_gamma_power_modulation_cue_draft_-1000_2500'];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
savefig(gcf, fullfile(project_dir, fig_lib, fig_name));
%%  Gamma power by gamma power modulation (cue aligned)
band_i_ = 3;
mod_plot({lfp_mod_cue([site_mod_cat{1, [1,3,4], 2}], :, band_i_), ...
    lfp_mod_cue([site_mod_cat{1, [1,3,4], 1}], :, band_i_), ...
    lfp_mod_cue([site_mod_cat{2, [1,3,4], 2}], :, band_i_), ...
    lfp_mod_cue([site_mod_cat{2, [1,3,4], 1}], :, band_i_)}, ...
    linspace(-1,4,250), ...
    {'Ado mod', 'Ado non-mod', 'Adu mod', 'Adu non-mod'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
%     {'Adolescent modulated', 'Adolescent non-modulated', 'Adult modulated', 'Adult non-modulated'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
xlim([-1, 2.5]);
title(title_st{band_i_})
xlabel('Time from cue onset (s)')
fig_name = ['gamma_power_modulation_by_gamma_power_modulation_cue_draft_-1000_2500'];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
savefig(gcf, fullfile(project_dir, fig_lib, fig_name));
% 
%%  High Gamma power by high gamma power modulation (sac aligned)
band_i_ = 4;
mod_plot_sac({lfp_mod_sac([site_mod_cat{1, [1,3,4], 4}], :, band_i_), ...
    lfp_mod_sac([site_mod_cat{1, [1,3,4], 3}], :, band_i_), ...
    lfp_mod_sac([site_mod_cat{2, [1,3,4], 4}], :, band_i_), ...
    lfp_mod_sac([site_mod_cat{2, [1,3,4], 3}], :, band_i_)}, ...
    linspace(-2,3,250), ...
    {'Ado mod', 'Ado non-mod', 'Adu mod', 'Adu non-mod'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
%     {'Adolescent modulated', 'Adolescent non-modulated', 'Adult modulated', 'Adult non-modulated'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title(title_st{band_i_})
% legend off
xlabel('Time from saccadic onset (s)')
fig_name = ['high_gamma_power_modulation_by_gamma_power_modulation_sac_draft'];
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
%%  Gamma power by gamma power modulation (sac aligned)
band_i_ = 3;
mod_plot_sac({lfp_mod_sac([site_mod_cat{1, [1,3,4], 2}], :, band_i_), ...
    lfp_mod_sac([site_mod_cat{1, [1,3,4], 1}], :, band_i_), ...
    lfp_mod_sac([site_mod_cat{2, [1,3,4], 2}], :, band_i_), ...
    lfp_mod_sac([site_mod_cat{2, [1,3,4], 1}], :, band_i_)}, ...
    linspace(-2,3,250), ...
    {'Ado mod', 'Ado non-mod', 'Adu mod', 'Adu non-mod'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
%     {'Adolescent modulated', 'Adolescent non-modulated', 'Adult modulated', 'Adult non-modulated'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title(title_st{band_i_})
xlabel('Time from saccadic onset (s)')
fig_name = ['gamma_power_modulation_by_gamma_power_modulation_sac_draft'];
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
% 
%%  LFP power modulation by spiking tuning (cue aligned)
for band_i_ = 1:4
mod_plot({lfp_mod_cue([site_tuning_cat{1, [1,3,4], :, 2}], :, band_i_), ...
    lfp_mod_cue([site_tuning_cat{1, [1,3,4], :, 1}], :, band_i_), ...
    lfp_mod_cue([site_tuning_cat{2, [1,3,4], :, 2}], :, band_i_), ...
    lfp_mod_cue([site_tuning_cat{2, [1,3,4], :, 1}], :, band_i_)}, ...
    linspace(-1,4,250), ...
    {'Adolescent inf', 'Adolescent non-inf', 'Adult inf', 'Adult non-inf'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
%     {'Adolescent informative', 'Adolescent non-informative', 'Adult informative', 'Adult non-informative'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title(title_st{band_i_})
xlabel('Time from cue onset (s)')
xlim([-1, 2.5]);
if band_i_ ~= 2
    legend off
end
fig_name = ['lfp_pow_band_', num2str(band_i_), 'power_modulation_by_spiking_tuning', '_cue_draft_-1000_2500'];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
savefig(gcf, fullfile(project_dir, fig_lib, fig_name));
end
%%  LFP power modulation by spiking tuning (sac aligned)
for band_i_ = 1:4
mod_plot_sac({lfp_mod_sac([site_tuning_cat{1, [1,3,4], :, 2}], :, band_i_), ...
    lfp_mod_sac([site_tuning_cat{1, [1,3,4], :, 1}], :, band_i_), ...
    lfp_mod_sac([site_tuning_cat{2, [1,3,4], :, 2}], :, band_i_), ...
    lfp_mod_sac([site_tuning_cat{2, [1,3,4], :, 1}], :, band_i_)}, ...
    linspace(-2,3,250), ...
    {'Adolescent inf', 'Adolescent non-inf', 'Adult inf', 'Adult non-inf'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
%     {'Adolescent informative', 'Adolescent non-informative', 'Adult informative', 'Adult non-informative'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})title(title_st{band_i_})
title(title_st{band_i_})
xlabel('Time from saccadic onset (s)')
if band_i_ ~= 2
    legend off
end
fig_name = ['lfp_pow_band_', num2str(band_i_), 'power_modulation_by_spiking_tuning', '_sac_draft'];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
end
%%  High Gamma power modulation by spiking tuning
band_i_ = 4;
mod_plot({lfp_mod_cue([site_tuning_cat{1, [1,3,4], :, 2}], :, band_i_), ...
    lfp_mod_cue([site_tuning_cat{1, [1,3,4], :, 1}], :, band_i_), ...
    lfp_mod_cue([site_tuning_cat{2, [1,3,4], :, 2}], :, band_i_), ...
    lfp_mod_cue([site_tuning_cat{2, [1,3,4], :, 1}], :, band_i_)}, ...
    linspace(-1,4,250), ...
    {'Adolescent informative', 'Adolescent non-informative', 'Adult informative', 'Adult non-informative'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title('High Gamma Power')
xlabel('Time from cue onset (s)')
fig_name = ['high_gamma_power_modulation_by_spiking_tuning_draft'];
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
%%
%%  Spiking firing rate raw
psth_plot({best_psth_raw([neuron_tuning_cat{1, [1,3,4], :, :}], :), ...
    best_psth_raw([neuron_tuning_cat{2, [1,3,4], :, :}], :)}, ...
    linspace(-1,4,250), ...
    {'Adolescent', 'Adult'}, {'b','r'})
title('Absolute Firing Rate')
xlabel('Time from cue onset (s)')
xlim([-1. 2.5])
fig_name = ['raw_spiking_firing_draft_-1000_2500'];
legend off
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
savefig(gcf, fullfile(project_dir, fig_lib, fig_name))
%%  Spiking firing rate (1000ms baseline)
psth_plot({best_psth_1s_baseline([neuron_tuning_cat{1, [1,3,4], :, :}], :), ...
    best_psth_1s_baseline([neuron_tuning_cat{2, [1,3,4], :, :}], :)}, ...
    linspace(-1,4,250), ...
    {'Adolescent', 'Adult'}, {'b','r'})
title('Evoked Firing Rate')
xlabel('Time from cue onset (s)')
xlim([-1. 2.5])
fig_name = ['evoked_spiking_firing_draft_-1000_2500'];
legend off
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
savefig(gcf, fullfile(project_dir, fig_lib, fig_name))
%%
%%  Spiking firing rate (1000ms baseline) by spiking tuning
psth_plot({best_psth_1s_baseline([neuron_tuning_cat{1, [1,3,4], :, 2}], :), ...
    best_psth_1s_baseline([neuron_tuning_cat{1, [1,3,4], :, 1}], :), ...
    best_psth_1s_baseline([neuron_tuning_cat{2, [1,3,4], :, 2}], :), ...
    best_psth_1s_baseline([neuron_tuning_cat{2, [1,3,4], :, 1}], :)}, ...
    linspace(-1,4,250), ...
    {'Adolescent in', 'Adolescent non-in', 'Adult in', 'Adult non-in'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
%     {'Adolescent informative', 'Adolescent non-informative', 'Adult informative', 'Adult non-informative'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
xlim([-1. 2.5])
title('Evoked Firing Rate')
xlabel('Time from cue onset (s)')
fig_name = ['spiking_firing_rate_by_spiking_tuning_draft_-1000_2500'];
% legend off
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
% savefig(gcf, fullfile(project_dir, fig_lib, fig_name))

%%  Neuron PEV by spiking tuning
pev_plot({neuron_pev_cue([neuron_tuning_cat{1, [1,3,4], :, 2}], :), ...
    neuron_pev_cue([neuron_tuning_cat{1, [1,3,4], :, 1}], :), ...
    neuron_pev_cue([neuron_tuning_cat{2, [1,3,4], :, 2}], :), ...
    neuron_pev_cue([neuron_tuning_cat{2, [1,3,4], :, 1}], :)}, ...
    linspace(-1,4,250), ...
    {'Adolescent informative', 'Adolescent non-informative', 'Adult informative', 'Adult non-informative'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]}, ...
    [0, 0.25])
xlim([-1. 2.5])
title('Spiking PEV')
xlabel('Time from cue onset (s)')
fig_name = ['neuron_pev_by_spiking_tuning_draft_-1000_2500'];
legend off
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
% savefig(gcf, fullfile(project_dir, fig_lib, fig_name))
%
%%
set(gca, 'Position', [0.1731 0.1669 0.7319 0.7469])
fig_name = ['neuron_pev_by_spiking_tuning_draft_-1000_2500_axes_position_fixed'];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
savefig(gcf, fullfile(project_dir, fig_lib, fig_name))

%%  Neuron PEV for non-significant neurons
pev_plot({neuron_pev_cue([neuron_tuning_cat_nonsig{1, [1,3,4], :, :}], :), ...
    neuron_pev_cue([neuron_tuning_cat_nonsig{2, [1,3,4], :, :}], :)}, ...
    linspace(-1,4,250), ...
    {'Adolescent nonsig', 'Adult nonsig'}, {[0.3, 0.85, 0.85], [0.85, 0.5, 0.3]})
title('Spiking PEV')
xlabel('Time from cue onset (s)')
fig_name = ['neuron_pev_nonsig_draft'];
% legend off
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
% 
% 
% 
%%  Spiking firing rate raw (nonsig)
psth_plot({best_psth_raw([neuron_tuning_cat_nonsig{1, [1,3,4], :, :}], :), ...
    best_psth_raw([neuron_tuning_cat_nonsig{2, [1,3,4], :, :}], :)}, ...
    linspace(-1,3,200), ...
    {'Adolescent nonsig', 'Adult nonsig'},  {[0.3, 0.85, 0.85], [0.85, 0.5, 0.3]})
title('Absolute Firing Rate')
xlabel('Time from cue onset (s)')
% ylim([0, 10])
fig_name = ['raw_spiking_firing_nonsig_draft'];
legend off
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
%%  Neuron PEV by gamma power modulation (sig + nonsig)
pev_plot({neuron_pev_cue([neuron_mod_cat{1, [1,3,4], 2}, neuron_mod_cat_nonsig{1, [1,3,4], 2}], :), ...
    neuron_pev_cue([neuron_mod_cat{1, [1,3,4], 1}, neuron_mod_cat_nonsig{1, [1,3,4], 1}], :), ...
    neuron_pev_cue([neuron_mod_cat{2, [1,3,4], 2}, neuron_mod_cat_nonsig{2, [1,3,4], 2}], :), ...
    neuron_pev_cue([neuron_mod_cat{2, [1,3,4], 1}, neuron_mod_cat_nonsig{2, [1,3,4], 1}], :)}, ...
    linspace(-1,4, 250), ...
    {'Ado mod', 'Ado non-mod', 'Adu mod', 'Adu non-mod'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
%     {'Adolescent gamma-modulated', 'Adolescent non-gamma-modulated', 'Adult gamma-modulated', 'Adult non-gamma-modulated'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title('Neuronal tuning')
xlabel('Time from cue onset (s)')
fig_name = ['neuron_pev_by_gamma_modulation_sig+nonsig'];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
%%  Neuron PEV by gamma power modulation (sig only)
pev_plot({neuron_pev_cue([neuron_mod_cat{1, [1,3,4], 2}], :), ...
    neuron_pev_cue([neuron_mod_cat{1, [1,3,4], 1}], :), ...
    neuron_pev_cue([neuron_mod_cat{2, [1,3,4], 2}], :), ...
    neuron_pev_cue([neuron_mod_cat{2, [1,3,4], 1}], :)}, ...
    linspace(-1,4,250), ...
    {'Ado mod', 'Ado non-mod', 'Adu mod', 'Adu non-mod'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
%     {'Adolescent gamma-modulated', 'Adolescent non-gamma-modulated', 'Adult gamma-modulated', 'Adult non-gamma-modulated'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title('Neuronal tuning')
xlabel('Time from cue onset (s)')
fig_name = ['neuron_pev_by_gamma_modulation_sig'];
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
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
imagesc(min_maxed_pev_(trace_order, :), [0, 1]);
% colormap('jet')
h = colorbar;
ylabel(h, 'Norm. Spiking PEV');
hold on
% plot(sorted_peak_timepoint+75, 1:numel(target_subset), 'r');

% 
band_i = 3;
mod_to_plot = lfp_mod_cue(map_neuron_to_site(mapping_mat, target_subset(trace_order)), :, band_i) - 1; %    Converting to change from baseline
mod_to_plot = mod_to_plot./max(mod_to_plot, [], 2);
ylabel('Neuron #')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1);
fig_name = ['adult_pev_trace_draft'];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');


figure
imagesc(mod_to_plot, [-1, 1])
colormap(clm_)
h = colorbar;
ylabel(h, 'Norm. Gamma Power');
ylabel('Neuron #')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1);
fig_name = ['adult_gamma_mod_trace_draft'];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');

%
band_i = 4;
mod_to_plot = lfp_mod_cue(map_neuron_to_site(mapping_mat, target_subset(trace_order)), :, band_i) - 1; %    Converting to change from baseline
mod_to_plot = mod_to_plot./max(mod_to_plot, [], 2);
figure
imagesc(mod_to_plot, [-1, 1])
colormap(clm_)
h = colorbar;
ylabel(h, 'Norm. High Gamma Power');
ylabel('Neuron #')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1);
fig_name = ['adult_high_gamma_mod_trace_draft'];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
%%  Sort single trials by delay period power peaks
fname_ = 'neuron_site_categories.mat';
load(fullfile(project_dir, output_database, fname_));
fname_ = 'neuron_site_categories_nonsig.mat';
load(fullfile(project_dir, output_database, fname_));
fname_ = 'best_psth.mat';
load(fullfile(project_dir, output_database, fname_), 'neuron_best_class');
fname_ = 'neuron_tbl_w_pfc.mat';
load(fullfile(project_dir, output_database, fname_));
fname_ = 'neuron_epoch_delay_tuning_flag.mat';
load(fullfile(project_dir, output_database, fname_), 'neuron_epoch_delay_tuning_flag');
fname_ = 'tfr_modulation.mat';
load(fullfile(project_dir, output_database, fname_), 'tfr_cue_modulation_flag', 'tfr_delay_modulation_flag', 'tfr_modulation_output_cue');
%%
fname_ = 'complete_cwt_repo_3_2_2-21.mat';
load(fullfile(project_dir, output_database, fname_));
%   Reduce repo size by emptying unwanted trials
for i = 1:numel(cwt_repo)
    if ~ismember(i, [site_mod_cat{:}])
        cwt_repo(i).class = [];
    end
end
%%
fname_ = 'upsample_neuron_repo_3_11_2021.mat';
load(fullfile(project_dir, output_database, fname_));
%%  Make neuron-cwt matched single trial repo of only required neurons/sites
neuron_repo_by_trial = struct();
neuron_tbl_by_trial = table();
trial_counter = 0;
neuron_list_ = [neuron_tuning_cat_nonsig{:, [1,3,4], :}, neuron_tuning_cat{:, [1,3,4], :}];
mapping_mat_total = or(mapping_mat, mapping_mat_nonsig);
%   The following warning would appear when assigning values to a subset of
%   table variables using subscripts (instead of variable name)
warning('off', 'MATLAB:table:RowsAddedExistingVars');
for i = neuron_list_
    lfp_id_ = find(mapping_mat_total(:, i));
    class_id_ = neuron_best_class(i);
    for k = 1:size(neuron_repo(i).class(class_id_).psth_cue_upsampled, 1)
        trial_counter = trial_counter + 1;
        neuron_repo_by_trial(trial_counter).psth_cue = neuron_repo(i).class(class_id_).psth_cue_upsampled(k, :);
        neuron_repo_by_trial(trial_counter).psth_sac = neuron_repo(i).class(class_id_).psth_sac_upsampled(k, :);
%         neuron_repo_by_trial(trial_counter).psth_cue = neuron_repo(i).class(class_id_).psth_cue_upsampled_for_pev(k, :);
%         neuron_repo_by_trial(trial_counter).psth_sac = neuron_repo(i).class(class_id_).psth_sac_upsampled_for_pev(k, :);
        neuron_repo_by_trial(trial_counter).cue_cwt  = cwt_repo(lfp_id_).class(class_id_).cue_cwt(k, :, :);
        neuron_repo_by_trial(trial_counter).cue_cwt_site_mean  = nanmean(cwt_repo(lfp_id_).class(class_id_).cue_cwt, 1);
        neuron_repo_by_trial(trial_counter).sac_cwt  = cwt_repo(lfp_id_).class(class_id_).sac_cwt(k, :, :);
        %
        neuron_tbl_by_trial(trial_counter, 1:size(neuron_tbl, 2)) = neuron_tbl(i, :);
        neuron_tbl_by_trial.neuron_repo_id(trial_counter) = i;
        neuron_tbl_by_trial.lfp_repo_id(trial_counter)    = lfp_id_;
        neuron_tbl_by_trial.best_class(trial_counter)     = class_id_;
        neuron_tbl_by_trial.trial_n(trial_counter)        = k;
        neuron_tbl_by_trial.responsive(trial_counter)     = responsive(i);
        neuron_tbl_by_trial.informative(trial_counter)    = neuron_epoch_delay_tuning_flag(i);
        neuron_tbl_by_trial.delay_gamma_mod(trial_counter)= tfr_delay_modulation_flag(lfp_id_, 3, 1); 
        if mod(trial_counter, 500) == 0
            trial_counter
        end
    end
end
neuron_tbl_by_trial.Properties.VariableNames(1:size(neuron_tbl, 2)) = neuron_tbl.Properties.VariableNames;
%%
fname_ = 'neuron_repo_by_trial';
save(fullfile(project_dir, output_database, fname_), 'neuron_tbl_by_trial', 'neuron_repo_by_trial');
%%
psth_cue_ = cat(1, neuron_repo_by_trial.psth_cue);
cue_cwt_  = cat(1, neuron_repo_by_trial.cue_cwt);
cue_gamma_ = squeeze(sum(cue_cwt_(:, 16:32, :), 2));
cue_gamma_mod_ = cue_gamma_./mean(cue_gamma_(:, 26:50), 2) - 1;
cue_high_gamma_ = squeeze(sum(cue_cwt_(:, 32:64, :), 2));
cue_high_gamma_mod_ = cue_high_gamma_./mean(cue_high_gamma_(:, 26:50), 2) - 1;
cue_all_gamma_ = squeeze(sum(cue_cwt_(:, 16:64, :), 2));
% 
cue_cwt_site_mean_ = cat(1, neuron_repo_by_trial.cue_cwt_site_mean);
cue_gamma_site_mean_ = squeeze(sum(cue_cwt_site_mean_(:, 16:32, :), 2));
%%
keep_trial = ones(size(cue_gamma_, 1), 1);
for i = 1:size(cue_gamma_, 1)
    if any(isnan(cue_gamma_(i, :)))
        keep_trial(i) = 0;
    end
end
%%
young_target_ = find(neuron_tbl_by_trial.responsive.*(neuron_tbl_by_trial.stage == 1).*keep_trial);
adult_target_ = find(neuron_tbl_by_trial.responsive.*(neuron_tbl_by_trial.stage == 2).*keep_trial);
%%  Candidate metrics
sort_window = 76:150;
min_maxed_psth_       = (psth_cue_ - min(psth_cue_(:, sort_window), [], 2))./(max(psth_cue_(:, sort_window), [], 2) - min(psth_cue_(:, sort_window), [], 2));
evoked_psth_ = psth_cue_ - mean(psth_cue_(:, 1:50), 2) ;
min_maxed_evoked_psth_ = evoked_psth_./max(abs(evoked_psth_(:, sort_window)), [], 2);
min_maxed_gamma_      = (cue_gamma_ - min(cue_gamma_(:, sort_window), [], 2))./(max(cue_gamma_(:, sort_window), [], 2) - min(cue_gamma_(:, sort_window), [], 2));
min_maxed_cue_gamma_mod_ = cue_gamma_mod_./max(abs(cue_gamma_mod_(:, sort_window)), [], 2);
min_maxed_cue_high_gamma_mod_ = cue_high_gamma_mod_./max(abs(cue_high_gamma_mod_(:, sort_window)), [], 2);
min_maxed_high_gamma_ = (cue_high_gamma_ - min(cue_high_gamma_(:, sort_window), [], 2))./(max(cue_high_gamma_(:, sort_window), [], 2) - min(cue_high_gamma_(:, sort_window), [], 2));
min_maxed_all_gamma_  = (cue_all_gamma_ - min(cue_all_gamma_(:, sort_window), [], 2))./(max(cue_all_gamma_(:, sort_window), [], 2) - min(cue_all_gamma_(:, sort_window), [], 2));
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
%%  Gamma YOUNG
figure
imagesc(min_maxed_gamma_(young_target_(trace_order), :), [0, 1]);
% colormap('jet')
h = colorbar;
ylabel(h, 'Norm. gamma power');
hold on
% plot(sorted_peak_timepoint+75, 1:numel(target_subset), 'r');
% 
ylabel('Trial #')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1);
fig_name = ['young_gamma_trace_draft'];
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
%%  Gamma mod YOUNG
figure
imagesc(min_maxed_cue_gamma_mod_(young_target_(trace_order), :), [-1, 1]);
% imagesc(cue_gamma_mod_(young_target_(trace_order), :), [-.6, .6]);
% colormap('jet')
h = colorbar;
colormap(clm_);
ylabel(h, 'Norm. gamma (32-64 Hz) modulation');
hold on
% plot(sorted_peak_timepoint+75, 1:numel(target_subset), 'r');
% 
ylabel('Trial #')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1);
fig_name = ['young_gamma_mod_trace_gamma_sort_3_15_21'];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
%%  High gamma YOUNG
figure
imagesc(min_maxed_high_gamma_(young_target_(trace_order), :), [0, 1]);
% colormap('jet')
h = colorbar;
ylabel(h, 'Norm. high gamma power');
hold on
% plot(sorted_peak_timepoint+75, 1:numel(target_subset), 'r');
% 
ylabel('Trial #')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1);
fig_name = ['young_high_gamma_trace_draft'];
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
%%  All gamma YOUNG
figure
imagesc(min_maxed_all_gamma_(young_target_(trace_order), :), [0, 1]);
% colormap('jet')
h = colorbar;
ylabel(h, 'Norm. gamma (36 - 128 Hz) power');
hold on
% plot(sorted_peak_timepoint+75, 1:numel(target_subset), 'r');
% 
ylabel('Trial #')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1);
fig_name = ['young_all_gamma_trace_draft'];
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
%%  PSTH YOUNG
figure
% imagesc(min_maxed_psth_(young_target_(trace_order), :), [0, 1]);
% imagesc(evoked_psth_(young_target_(trace_order), :), [0, 15]);
imagesc(min_maxed_evoked_psth_(young_target_(trace_order), :), [0, 1]);
% colormap('jet')
h = colorbar;
ylabel(h, 'Norm. evoked firing rate');
hold on
% plot(sorted_peak_timepoint+75, 1:numel(target_subset), 'r');
% 
ylabel('Trial #')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1);
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

%%  Gamma ADULT
figure
imagesc(min_maxed_gamma_(adult_target_(trace_order), :), [0, 1]);
% colormap('jet')
h = colorbar;
ylabel(h, 'Norm. gamma power');
hold on
% plot(sorted_peak_timepoint+75, 1:numel(adult_target_), 'r');
% plot(peak_timepoint_adult_gamma_site_mean(trace_order)+75, 1:numel(adult_target_), '.r', 'MarkerSize', 1);
% 
ylabel('Trial #')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1);
fig_name = ['young_gamma_trace_draft'];
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
%%  Gamma mod ADULT
figure('Unit', 'inches', 'Position', [2, 2, 2.5, 2])
imagesc(min_maxed_cue_gamma_mod_(adult_target_(trace_order), :), [-1, 1]);
% imagesc(cue_gamma_mod_(adult_target_(trace_order), :), [-.6, .6]);
% colormap('jet')
h = colorbar;
colormap(clm_);
ylabel(h, 'Norm. gamma (32-64 Hz) modulation');
hold on
% plot(sorted_peak_timepoint+75, 1:numel(target_subset), 'r');
% 
ylabel('Trial #')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1);
fig_name = ['adult_gamma_mod_trace_gamma_sort_3_15_21'];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');

%%  Comapre peak times
figure
hold on
% imagesc(nan(size(min_maxed_gamma_(adult_target_(trace_order), :))));

plot(peak_timepoint_adult_gamma_site_mean(trace_order)+75, 1:numel(adult_target_), 'dr', 'MarkerSize', 1);
plot(sorted_peak_timepoint + 75, 1:numel(adult_target_), 'dk', 'MarkerSize', 1);
plot(peak_timepoint_adult_psth(trace_order)+75, 1:numel(adult_target_), 'db', 'MarkerSize', 1);
set(gca, 'YDir', 'reverse');
%%  Gamma site mean ADULT
figure
imagesc(min_maxed_gamma_site_mean_(adult_target_(trace_order), :), [0, 1]);
% colormap('jet')
h = colorbar;
ylabel(h, 'Norm. gamma power');
hold on
% plot(sorted_peak_timepoint+75, 1:numel(target_subset), 'r');
% 
ylabel('Trial #')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1);
fig_name = ['young_gamma_trace_draft'];
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');

%%  High gamma ADULT
figure
imagesc(min_maxed_high_gamma_(adult_target_(trace_order), :), [0, 1]);
% colormap('jet')
h = colorbar;
ylabel(h, 'Norm. high gamma power');
hold on
% plot(sorted_peak_timepoint+75, 1:numel(target_subset), 'r');
% 
ylabel('Trial #')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1);
fig_name = ['adult_high_gamma_trace_draft'];
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
%%  PSTH ADULT
figure('Unit', 'inches', 'Position', [2, 2, 2.5, 2])
imagesc([-1, 4], 1:numel(trace_order), min_maxed_psth_(adult_target_(trace_order), :), [0, 1]);
% imagesc(psth_cue_(adult_target_(trace_order), :), [0, 100]);
% imagesc(evoked_psth_(adult_target_(trace_order), :), [0, 15]);
% imagesc(min_maxed_evoked_psth_(adult_target_(trace_order), :), [0, 1]);
% ylabel(h, 'Norm. evoked firing rate');
% colormap('jet')
h = colorbar;
ylabel(h, 'Norm. firing rate');
xlabel('Time from cue onset (s)')
hold on
% plot(sorted_peak_timepoint+75, 1:numel(target_subset), 'r');
% 
ylabel('Trial #')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1);
fig_name = ['adult_psth_trace_gamma_sort_3_15_21'];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
%%
%%
clm_n = 100;
clm_ = [[zeros(1, clm_n);zeros(1, clm_n); linspace(1, 0, clm_n)]'; [linspace(0, 1, clm_n); zeros(1, clm_n);zeros(1, clm_n)]'];
clm_2nd = [[zeros(1, clm_n);linspace(150/200, 0, clm_n); linspace(1, 0, clm_n)]'; [linspace(0, 1, clm_n); linspace(0, 150/200, clm_n); zeros(1, clm_n)]'];
%%  Cluster-corrected t values (spectrogram)
for i = 1
    for_plot_pos = squeeze(stat{i}.posclusterslabelmat);
    for_plot_neg = squeeze(stat{i}.negclusterslabelmat);
    for_plot     = zeros(size(for_plot_pos));
    for j = 1:numel(for_plot_pos)
        if for_plot_pos(j)
            for_plot(j) = stat{i}.posclusters(for_plot_pos(j)).prob <0.05;
        end
        if for_plot_neg(j)
            for_plot(j) = -(stat{i}.negclusters(for_plot_neg(j)).prob <0.05);
        end
    end
    figure('Unit', 'inches', 'Position', [2, 2, 3, 2.5]);
    hold on
    imagesc(in_y.time, in_y.freq, for_plot, [-1,1]);
    pos_center = zw_find_cluster_center(for_plot == 1);
    if ~any(isnan(pos_center))
        text(in_y.time(pos_center(2)), in_y.freq(pos_center(1)), sprintf('%.3f', stat{i}.posclusters(1).prob), 'HorizontalAlignment', 'center', 'Color', 'w');
    end
    neg_center = zw_find_cluster_center(for_plot == -1);
    if ~any(isnan(neg_center))
        text(in_y.time(neg_center(2)), in_y.freq(neg_center(1)), sprintf('%.3f', stat{i}.negclusters(1).prob), 'HorizontalAlignment', 'center', 'Color', 'w');
    end
    plot([in_y.time(1), in_y.time(end)], [1, 1].*32, '--', 'Color', 'w')
    xlabel('Time from stim. onset (s)')
    ylabel('Frequency (Hz)')
    xlim([in_y.time(1), in_y.time(end)]);
    ylim([in_y.freq(1), in_y.freq(end)]);
    xticks([0:0.5:2])
    yt = yticks;
    yticks(sort([yt, 32]));
    title(title_strings{i})
    set(gca, {'FontSize', 'FontWeight'}, {11, 'bold'})
%     colormap([0, 0, .8; 0, 0, 0; .8, 0, 0]);
    colormap([0, 150/200, 1; 0, 0, 0; 1, 150/200, 0]);
    h = colorbar;
    h.Ticks = [-.7, 0, .7];
    h.TickLabels = {'Adolescent\newline<\newlineAdult', 'Non-sig.', 'Adolescent\newline>\newlineAdult'};
%     h.FontSize = 7;
    add_epochline([1 1 1])
set_plot_3_4_2021();
    fig_name = 'power_diff_cluster_standard_size';
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
% savefig(gcf, fullfile(project_dir, fig_lib, fig_name));
end
%%
% 
figure('Unit', 'inches', 'Position', [2, 2, 3, 2.5]);
imagesc(in_y.time, in_y.freq, squeeze(stat{1}.stat), [-10, 10]);
set(gca, 'YDir','normal'); 
colormap(clm_2nd)
xlabel('Time from cue onset (s)')
ylabel('Frequency (Hz)')
h = colorbar();
% h.Ruler.TickLabelFormat='%g%%';
ylabel(h, 'T score')
% xlim([-1, 2.5])
title('Adolscent - Adult')
fig_name = 'spectrogram_diff_t_standard_size'
add_epochline([1 1 1])
set_plot_3_4_2021();
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
savefig(gcf, fullfile(project_dir, fig_lib, fig_name));