responsive = any(table2array(neuron_tbl(:, 9:12)), 2);
%%  Converts maunally rejected session flag to 0 and 1
valid_lfp = zeros(size(lfp_tbl, 1), 1);
for i = 1:size(lfp_tbl, 1)
    if any(target_cwts(~logical(session_flag_list)) == i)
        valid_lfp(i) = 1;
    end
end     
%%  Base 2D matrix mapping LFP and neuron by electrodes, eliminating invalid sessions
mapping_mat = zeros(size(lfp_tbl, 1), size(neuron_tbl, 1));
for i = 1:size(lfp_tbl, 1)
%     i
    mapping_mat(i, :) = prod(table2array(lfp_tbl(i, 2:5)) == table2array(neuron_tbl(:, [6, 7, 8, 3])), 2).*neuron_tbl.PFC.*responsive;
end
% mapping_mat = mapping_mat.*valid_lfp.*(lfp_tbl.task_id == 1);
mapping_mat = mapping_mat.*valid_lfp.*(lfp_tbl.ps_file == 1);
% mapping_mat = mapping_mat.*(lfp_tbl.ps_file == 1);
%%  Indexing each stage-by-animal in 't' (LFP) and 'n' (neuron)
t = cell(2, 4);
n = cell(2, 4);
for i_stage = 1:2
    for i_monk = 1:4
        t{i_stage, i_monk} = find(any((lfp_tbl.stage == i_stage).*(lfp_tbl.monkey_id == i_monk).*mapping_mat, 2))';
        n{i_stage, i_monk} = find(any((lfp_tbl.stage == i_stage).*(lfp_tbl.monkey_id == i_monk).*mapping_mat, 1));
    end
end
%%
save(fullfile(project_dir, output_database, 'lfp_neuron_matching.mat'), 't', 'n', 'mapping_mat');
%%  Neuron tuning categorization (gamma X spiking).
temp_counts_ = {};
neuron_tuning_cat = cell(2, 4, 2, 2);
for i_stage = 1:2
    for i_monk = 1:4
        %               stage,   monkey, gamma tuning, spiking tuning
        neuron_tuning_cat{i_stage, i_monk, 1, 1} = intersect(n{i_stage, i_monk}, find(any(mapping_mat.*~tfr_delay_tuning_flag(:, 3).*~neuron_epoch_delay_tuning_flag, 1)));        
        neuron_tuning_cat{i_stage, i_monk, 1, 2} = intersect(n{i_stage, i_monk}, find(any(mapping_mat.*~tfr_delay_tuning_flag(:, 3).*neuron_epoch_delay_tuning_flag, 1)));
        neuron_tuning_cat{i_stage, i_monk, 2, 1} = intersect(n{i_stage, i_monk}, find(any(mapping_mat.*tfr_delay_tuning_flag(:, 3).*~neuron_epoch_delay_tuning_flag, 1)));
        neuron_tuning_cat{i_stage, i_monk, 2, 2} = intersect(n{i_stage, i_monk}, find(any(mapping_mat.*tfr_delay_tuning_flag(:, 3).*neuron_epoch_delay_tuning_flag, 1)));
        temp_counts_{i_stage, i_monk} = [numel(neuron_tuning_cat{i_stage, i_monk, 1, 1}), numel(neuron_tuning_cat{i_stage, i_monk, 1, 2}); ...
            numel(neuron_tuning_cat{i_stage, i_monk, 2, 1}), numel(neuron_tuning_cat{i_stage, i_monk, 2, 2})];
    end
end
%%  Site tuning categorization (gamma X spiking)
site_tuning_cat = cell(2, 4, 2, 2);
for i_stage = 1:2
    for i_monk = 1:4
        %               stage,   monkey, gamma tuning, spiking tuning
        site_tuning_cat{i_stage, i_monk, 1, 1} = intersect(t{i_stage, i_monk}, find(any(mapping_mat.*~tfr_delay_tuning_flag(:, 3).*~neuron_epoch_delay_tuning_flag, 2))');        
        site_tuning_cat{i_stage, i_monk, 1, 2} = intersect(t{i_stage, i_monk}, find(any(mapping_mat.*~tfr_delay_tuning_flag(:, 3).*neuron_epoch_delay_tuning_flag, 2))');
        site_tuning_cat{i_stage, i_monk, 2, 1} = intersect(t{i_stage, i_monk}, find(any(mapping_mat.*tfr_delay_tuning_flag(:, 3).*~neuron_epoch_delay_tuning_flag, 2))');
        site_tuning_cat{i_stage, i_monk, 2, 2} = intersect(t{i_stage, i_monk}, find(any(mapping_mat.*tfr_delay_tuning_flag(:, 3).*neuron_epoch_delay_tuning_flag, 2))');
    end
end

%%  Neuron Gamma power modulation categorization.
neuron_mod_cat = cell(2, 4, 2);
for i_stage = 1:2
    for i_monk = 1:4
        %               stage,   monkey,   delay gamma mod
        neuron_mod_cat{i_stage, i_monk, 1} = intersect(n{i_stage, i_monk}, find(any(mapping_mat.*~tfr_delay_modulation_flag(:, 3, 1), 1)));
        neuron_mod_cat{i_stage, i_monk, 2} = intersect(n{i_stage, i_monk}, find(any(mapping_mat.*tfr_delay_modulation_flag(:, 3, 1), 1)));
    end
end
%%  Site Gamma power modulation categorization.
site_mod_cat = cell(2, 4, 2);
for i_stage = 1:2
    for i_monk = 1:4
        %               stage,   monkey, delay gamma mod
        site_mod_cat{i_stage, i_monk, 1} = intersect(t{i_stage, i_monk}, find(any(mapping_mat.*~tfr_delay_modulation_flag(:, 3, 1), 2))');
        site_mod_cat{i_stage, i_monk, 2} = intersect(t{i_stage, i_monk}, find(any(mapping_mat.*tfr_delay_modulation_flag(:, 3, 1), 2))');
    end
end
%%  
% n_1_t = intersect([n{1,3:4}], find(any(mapping_mat.*tfr_delay_tuning_flag(:, 3), 1)));
% n_1_nt = intersect([n{1,3:4}], find(any(mapping_mat.*~tfr_delay_tuning_flag(:, 3), 1)));
% n_2_t = intersect([n{2,3:4}], find(any(mapping_mat.*tfr_delay_tuning_flag(:, 3), 1)));
% n_2_nt = intersect([n{2,3:4}], find(any(mapping_mat.*~tfr_delay_tuning_flag(:, 3), 1)));
n_1_t = intersect([n{1,[1,3,4]}], find(any(mapping_mat.*tfr_delay_tuning_flag(:, 3), 1)));
n_1_nt = intersect([n{1,[1,3,4]}], find(any(mapping_mat.*~tfr_delay_tuning_flag(:, 3), 1)));
n_2_t = intersect([n{2,[1,3,4]}], find(any(mapping_mat.*tfr_delay_tuning_flag(:, 3), 1)));
n_2_nt = intersect([n{2,[1,3,4]}], find(any(mapping_mat.*~tfr_delay_tuning_flag(:, 3), 1)));

%%
% t_1_t = intersect([t{1,3:4}], find(any(mapping_mat.*tfr_delay_tuning_flag(:, 3), 2)));
% t_1_nt = intersect([t{1,3:4}], find(any(mapping_mat.*~tfr_delay_tuning_flag(:, 3), 2)));
% t_2_t = intersect([t{2,3:4}], find(any(mapping_mat.*tfr_delay_tuning_flag(:, 3), 2)));
% t_2_nt = intersect([t{2,3:4}], find(any(mapping_mat.*~tfr_delay_tuning_flag(:, 3), 2)));
t_1_t = intersect([t{1,[1,3,4]}], find(any(mapping_mat.*tfr_delay_tuning_flag(:, 3), 2)));
t_1_nt = intersect([t{1,[1,3,4]}], find(any(mapping_mat.*~tfr_delay_tuning_flag(:, 3), 2)));
t_2_t = intersect([t{2,[1,3,4]}], find(any(mapping_mat.*tfr_delay_tuning_flag(:, 3), 2)));
t_2_nt = intersect([t{2,[1,3,4]}], find(any(mapping_mat.*~tfr_delay_tuning_flag(:, 3), 2)));
%%
n_bands = 4;
target_frs = [4,8; 8,16; 16, 32; 32, 64];
baseline_bin = [26, 50];
title_st = {'Alpha (8-16 Hz)', 'Beta (16-32 Hz)', 'Gamma (32-64 Hz)', 'High-Gamma (64-128 Hz)'};

%%  alpha_tuning
figure();
plot(nanmean(lfp_pev_cue([t13; t14], :, 1), 1))
hold on
plot(nanmean(lfp_pev_cue([t23; t24], :, 1), 1))
%%  LFP tuning
for i = 1:n_bands
pev_plot({lfp_pev_cue([t{1, [1,3,4]}], :, i), lfp_pev_cue([t{2, [1,3,4]}], :, i)}, linspace(-1,3,200), {'Adolescent', 'Adult'}, {'b' 'r'})
% pev_plot({lfp_pev_cue([t13; t14], :, i), lfp_pev_cue([t23; t24], :, i)}, linspace(-1,3,200), {'Adolescent', 'Adult'}, {'b' 'r'})
xlabel('Time from cue onset (s)')
title(title_st{i})
fig_name = ['lfp_pev_band_', num2str(i)];
yl = ylim;
l1 = line([0, 0.5], [yl(1) + 0.001, yl(1) + 0.001], 'LineWidth',2, 'Color', 'g');
l2 = line([2, 2], yl, 'LineStyle', '--', 'Color', 'k');
set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylim(yl);
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
end
%%  Neuron PEV
% pev_plot({neuron_pev_cue([n13; n14], :), neuron_pev_cue([n23; n24], :)}, linspace(-1,3,200), {'Adolescent', 'Adult'}, {'b' 'r'})
pev_plot({neuron_pev_cue([n{1, [1,3,4]}], :), neuron_pev_cue([n{2, [1,3,4]}], :)}, linspace(-1,3,200), {'Adolescent', 'Adult'}, {'b' 'r'})
title('Spiking')
xlabel('Time from cue onset (s)')
fig_name = ['spiking_pev'];
yl = ylim;
l1 = line([0, 0.5], [yl(1) + 0.01, yl(1) + 0.01], 'LineWidth',2, 'Color', 'g');
l2 = line([2, 2], yl, 'LineStyle', '--', 'Color', 'k');
set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylim(yl);
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
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
fig_name = ['lfp_pow_band_', num2str(i)];
yl = ylim;
l1 = line([0, 0.5], [yl(1) + 0.01, yl(1) + 0.01], 'LineWidth',2, 'Color', 'g');
l2 = line([2, 2], yl, 'LineStyle', '--', 'Color', 'k');
set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylim(yl);
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
end
%%  Spiking (1000ms baseline subtracted)
psth_plot({best_psth_1s_baseline(n_1_t, :), best_psth_1s_baseline(n_1_nt, :), best_psth_1s_baseline(n_2_t, :), best_psth_1s_baseline(n_2_nt, :)}, linspace(-1,3,200), ...
    {'Adolescent w/ gamma tuning', 'Adolescent w/o gamma tuning', 'Adult w/ gamma tuning', 'Adult w/o gamma tuning'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title('Spiking')
xlabel('Time from cue onset (s)')
fig_name = ['PSTH'];
yl = ylim;
l1 = line([0, 0.5], [yl(1) + 0.01, yl(1) + 0.01], 'LineWidth',2, 'Color', 'g');
l2 = line([2, 2], yl, 'LineStyle', '--', 'Color', 'k');
set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylim(yl);
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');

%%  Spiking (500ms baseline subtracted)
psth_plot({best_psth(n_1_t, :), best_psth(n_1_nt, :), best_psth(n_2_t, :), best_psth(n_2_nt, :)}, linspace(-1,3,200), ...
    {'Adolescent w/ gamma tuning', 'Adolescent w/o gamma tuning', 'Adult w/ gamma tuning', 'Adult w/o gamma tuning'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title('Spiking')
xlabel('Time from cue onset (s)')
fig_name = ['PSTH'];
yl = ylim;
l1 = line([0, 0.5], [yl(1) + 0.01, yl(1) + 0.01], 'LineWidth',2, 'Color', 'g');
l2 = line([2, 2], yl, 'LineStyle', '--', 'Color', 'k');
set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylim(yl);
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
%%  Spiking (baseline not subtracted)
psth_plot({best_psth_raw(n_1_t, :), best_psth_raw(n_1_nt, :), best_psth_raw(n_2_t, :), best_psth_raw(n_2_nt, :)}, linspace(-1,3,200), ...
    {'Adolescent w/ gamma tuning', 'Adolescent w/o gamma tuning', 'Adult w/ gamma tuning', 'Adult w/o gamma tuning'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title('Spiking')
xlabel('Time from cue onset (s)')
fig_name = ['PSTH'];
yl = ylim;
l1 = line([0, 0.5], [yl(1) + 0.01, yl(1) + 0.01], 'LineWidth',2, 'Color', 'g');
l2 = line([2, 2], yl, 'LineStyle', '--', 'Color', 'k');
set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylim(yl);
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');

%%  Neuron PEV by gamma tuning
pev_plot({neuron_pev_cue(n_1_t, :), neuron_pev_cue(n_1_nt, :), neuron_pev_cue(n_2_t, :),neuron_pev_cue(n_2_nt, :)}, linspace(-1,3,200), ...
    {'Adolescent w/ gamma tuning', 'Adolescent w/o gamma tuning', 'Adult w/ gamma tuning', 'Adult w/o gamma tuning'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title('Spiking')
xlabel('Time from cue onset (s)')
fig_name = ['spiking_pev_by_gamma'];
yl = ylim;
l1 = line([0, 0.5], [yl(1) + 0.01, yl(1) + 0.01], 'LineWidth',2, 'Color', 'g');
l2 = line([2, 2], yl, 'LineStyle', '--', 'Color', 'k');
set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylim([yl(1), yl(2)+0.04]);
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
%%  Alpha PEV by spiking tuning
band_i_ = 2;
pev_plot({lfp_pev_cue([site_tuning_cat{1, [1,3,4], :, 2}], :, band_i_), ...
    lfp_pev_cue([site_tuning_cat{1, [1,3,4], :, 1}], :, band_i_), ...
    lfp_pev_cue([site_tuning_cat{2, [1,3,4], :, 2}], :, band_i_), ...
    lfp_pev_cue([site_tuning_cat{2, [1,3,4], :, 1}], :, band_i_)}, ...
    linspace(-1,3,200), ...
    {'Adolescent informative', 'Adolescent non-informative', 'Adult informative', 'Adult non-informative'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title('Alpha Tuning')
xlabel('Time from cue onset (s)')
fig_name = ['Alpha_pev_by_spiking_tuning'];
yl = ylim;
l1 = line([0, 0.5], [yl(1) + 0.01, yl(1) + 0.01], 'LineWidth',2, 'Color', 'g');
l2 = line([2, 2], yl, 'LineStyle', '--', 'Color', 'k');
set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylim([yl(1), yl(2)+0.04]);
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
%%  Beta PEV by spiking tuning
band_i_ = 2;
pev_plot({lfp_pev_cue([site_tuning_cat{1, [1,3,4], :, 2}], :, band_i_), ...
    lfp_pev_cue([site_tuning_cat{1, [1,3,4], :, 1}], :, band_i_), ...
    lfp_pev_cue([site_tuning_cat{2, [1,3,4], :, 2}], :, band_i_), ...
    lfp_pev_cue([site_tuning_cat{2, [1,3,4], :, 1}], :, band_i_)}, ...
    linspace(-1,3,200), ...
    {'Adolescent informative', 'Adolescent non-informative', 'Adult informative', 'Adult non-informative'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title('Beta Tuning')
xlabel('Time from cue onset (s)')
fig_name = ['beta_pev_by_spiking_tuning'];
yl = ylim;
l1 = line([0, 0.5], [yl(1) + 0.01, yl(1) + 0.01], 'LineWidth',2, 'Color', 'g');
l2 = line([2, 2], yl, 'LineStyle', '--', 'Color', 'k');
set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylim([yl(1), yl(2)+0.04]);
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
%%  Gamma PEV by spiking tuning
band_i_ = 3;
pev_plot({lfp_pev_cue([site_tuning_cat{1, [1,3,4], :, 2}], :, band_i_), ...
    lfp_pev_cue([site_tuning_cat{1, [1,3,4], :, 1}], :, band_i_), ...
    lfp_pev_cue([site_tuning_cat{2, [1,3,4], :, 2}], :, band_i_), ...
    lfp_pev_cue([site_tuning_cat{2, [1,3,4], :, 1}], :, band_i_)}, ...
    linspace(-1,3,200), ...
    {'Adolescent informative', 'Adolescent non-informative', 'Adult informative', 'Adult non-informative'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title('Gamma Tuning')
xlabel('Time from cue onset (s)')
fig_name = ['gamma_pev_by_spiking_tuning'];
yl = ylim;
l1 = line([0, 0.5], [yl(1) + 0.01, yl(1) + 0.01], 'LineWidth',2, 'Color', 'g');
l2 = line([2, 2], yl, 'LineStyle', '--', 'Color', 'k');
set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylim([yl(1), yl(2)+0.04]);
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
%%  High gamma PEV by spiking tuning
band_i_ = 4;
pev_plot({lfp_pev_cue([site_tuning_cat{1, [1,3,4], :, 2}], :, band_i_), ...
    lfp_pev_cue([site_tuning_cat{1, [1,3,4], :, 1}], :, band_i_), ...
    lfp_pev_cue([site_tuning_cat{2, [1,3,4], :, 2}], :, band_i_), ...
    lfp_pev_cue([site_tuning_cat{2, [1,3,4], :, 1}], :, band_i_)}, ...
    linspace(-1,3,200), ...
    {'Adolescent informative', 'Adolescent non-informative', 'Adult informative', 'Adult non-informative'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title('High Gamma Tuning')
xlabel('Time from cue onset (s)')
fig_name = ['high_gamma_pev_by_spiking_tuning'];
yl = ylim;
l1 = line([0, 0.5], [yl(1) + 0.01, yl(1) + 0.01], 'LineWidth',2, 'Color', 'g');
l2 = line([2, 2], yl, 'LineStyle', '--', 'Color', 'k');
set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylim([yl(1), yl(2)+0.04]);
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');

%%  Alpha power modulation by spiking tuning
band_i_ = 1;
mod_plot({lfp_mod_cue([site_tuning_cat{1, [1,3,4], :, 2}], :, band_i_), ...
    lfp_mod_cue([site_tuning_cat{1, [1,3,4], :, 1}], :, band_i_), ...
    lfp_mod_cue([site_tuning_cat{2, [1,3,4], :, 2}], :, band_i_), ...
    lfp_mod_cue([site_tuning_cat{2, [1,3,4], :, 1}], :, band_i_)}, ...
    linspace(-1,3,200), ...
    {'Adolescent informative', 'Adolescent non-informative', 'Adult informative', 'Adult non-informative'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title('Alpha Power')
xlabel('Time from cue onset (s)')
fig_name = ['alpha_power_modulation_by_spiking_tuning'];
yl = ylim;
l1 = line([0, 0.5], [yl(1) + 0.01, yl(1) + 0.01], 'LineWidth',2, 'Color', 'g');
l2 = line([2, 2], yl, 'LineStyle', '--', 'Color', 'k');
set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylim(yl);
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
%%  Beta power modulation by spiking tuning
band_i_ = 2;
mod_plot({lfp_mod_cue([site_tuning_cat{1, [1,3,4], :, 2}], :, band_i_), ...
    lfp_mod_cue([site_tuning_cat{1, [1,3,4], :, 1}], :, band_i_), ...
    lfp_mod_cue([site_tuning_cat{2, [1,3,4], :, 2}], :, band_i_), ...
    lfp_mod_cue([site_tuning_cat{2, [1,3,4], :, 1}], :, band_i_)}, ...
    linspace(-1,3,200), ...
    {'Adolescent informative', 'Adolescent non-informative', 'Adult informative', 'Adult non-informative'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title('Beta Power')
xlabel('Time from cue onset (s)')
fig_name = ['beta_power_modulation_by_spiking_tuning'];
yl = ylim;
l1 = line([0, 0.5], [yl(1) + 0.01, yl(1) + 0.01], 'LineWidth',2, 'Color', 'g');
l2 = line([2, 2], yl, 'LineStyle', '--', 'Color', 'k');
set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylim(yl);
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');

%%  Gamma power modulation by spiking tuning
band_i_ = 3;
mod_plot({lfp_mod_cue([site_tuning_cat{1, [1,3,4], :, 2}], :, band_i_), ...
    lfp_mod_cue([site_tuning_cat{1, [1,3,4], :, 1}], :, band_i_), ...
    lfp_mod_cue([site_tuning_cat{2, [1,3,4], :, 2}], :, band_i_), ...
    lfp_mod_cue([site_tuning_cat{2, [1,3,4], :, 1}], :, band_i_)}, ...
    linspace(-1,3,200), ...
    {'Adolescent informative', 'Adolescent non-informative', 'Adult informative', 'Adult non-informative'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title('Gamma Power')
xlabel('Time from cue onset (s)')
fig_name = ['gamma_power_modulation_by_spiking_tuning'];
yl = ylim;
l1 = line([0, 0.5], [yl(1) + 0.01, yl(1) + 0.01], 'LineWidth',2, 'Color', 'g');
l2 = line([2, 2], yl, 'LineStyle', '--', 'Color', 'k');
set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylim(yl);
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
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
fig_name = ['high_gamma_power_modulation_by_spiking_tuning'];
yl = ylim;
l1 = line([0, 0.5], [yl(1) + 0.01, yl(1) + 0.01], 'LineWidth',2, 'Color', 'g');
l2 = line([2, 2], yl, 'LineStyle', '--', 'Color', 'k');
set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylim(yl);
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
%%  Spiking firing rate (baseline not subtracted) by spiking tuning
psth_plot({best_psth_raw([neuron_tuning_cat{1, [1,3,4], :, 2}], :), ...
    best_psth_raw([neuron_tuning_cat{1, [1,3,4], :, 1}], :), ...
    best_psth_raw([neuron_tuning_cat{2, [1,3,4], :, 2}], :), ...
    best_psth_raw([neuron_tuning_cat{2, [1,3,4], :, 1}], :)}, ...
    linspace(-1,3,200), ...
    {'Adolescent informative', 'Adolescent non-informative', 'Adult informative', 'Adult non-informative'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title('Neuronal Firing Rate')
xlabel('Time from cue onset (s)')
fig_name = ['spiking_firing_rate_by_spiking_tuning'];
yl = ylim;
l1 = line([0, 0.5], [yl(1) + 0.01, yl(1) + 0.01], 'LineWidth',2, 'Color', 'g');
l2 = line([2, 2], yl, 'LineStyle', '--', 'Color', 'k');
set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylim(yl);
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
%%  Spiking firing rate (1000ms baseline) by spiking tuning
psth_plot({best_psth_1s_baseline([neuron_tuning_cat{1, [1,3,4], :, 2}], :), ...
    best_psth_1s_baseline([neuron_tuning_cat{1, [1,3,4], :, 1}], :), ...
    best_psth_1s_baseline([neuron_tuning_cat{2, [1,3,4], :, 2}], :), ...
    best_psth_1s_baseline([neuron_tuning_cat{2, [1,3,4], :, 1}], :)}, ...
    linspace(-1,3,200), ...
    {'Adolescent informative', 'Adolescent non-informative', 'Adult informative', 'Adult non-informative'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title('Evoked Firing Rate')
xlabel('Time from cue onset (s)')
fig_name = ['spiking_firing_rate_by_spiking_tuning'];
yl = ylim;
l1 = line([0, 0.5], [yl(1) + 0.01, yl(1) + 0.01], 'LineWidth',2, 'Color', 'g');
l2 = line([2, 2], yl, 'LineStyle', '--', 'Color', 'k');
set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylim(yl);
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
%%  Neuron PEV by spiking tuning
pev_plot({neuron_pev_cue([neuron_tuning_cat{1, [1,3,4], :, 2}], :), ...
    neuron_pev_cue([neuron_tuning_cat{1, [1,3,4], :, 1}], :), ...
    neuron_pev_cue([neuron_tuning_cat{2, [1,3,4], :, 2}], :), ...
    neuron_pev_cue([neuron_tuning_cat{2, [1,3,4], :, 1}], :)}, ...
    linspace(-1,3,200), ...
    {'Adolescent informative', 'Adolescent non-informative', 'Adult informative', 'Adult non-informative'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title('Spiking PEV')
xlabel('Time from cue onset (s)')
fig_name = ['neuron_pev_by_spiking_tuning'];
yl = ylim;
l1 = line([0, 0.5], [yl(1) + 0.01, yl(1) + 0.01], 'LineWidth',2, 'Color', 'g');
l2 = line([2, 2], yl, 'LineStyle', '--', 'Color', 'k');
set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylim(yl);
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
%%  Spiking firing rate (baseline not subtracted) by gamma modulation
psth_plot({...
    best_psth_raw([neuron_mod_cat{1, [1,3,4], 2}], :), ...
    best_psth_raw([neuron_mod_cat{1, [1,3,4], 1}], :), ...
    best_psth_raw([neuron_mod_cat{2, [1,3,4], 2}], :), ...
    best_psth_raw([neuron_mod_cat{2, [1,3,4], 1}], :)}, ...
    linspace(-1,3,200), ...
    {'Adolescent gamma-modulated', 'Adolescent non-gamma-modulated', 'Adult gamma-modulated', 'Adult non-gamma-modulated'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title('Neuronal Firing Rate')
xlabel('Time from cue onset (s)')
fig_name = ['spiking_firing_rate_by_gamma_modulation'];
yl = ylim;
l1 = line([0, 0.5], [yl(1) + 0.01, yl(1) + 0.01], 'LineWidth',2, 'Color', 'g');
l2 = line([2, 2], yl, 'LineStyle', '--', 'Color', 'k');
set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylim(yl);
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
%%  Spiking firing rate (1000ms baseline subtracted) by gamma modulation
psth_plot({...
    best_psth_1s_baseline([neuron_mod_cat{1, [1,3,4], 2}], :), ...
    best_psth_1s_baseline([neuron_mod_cat{1, [1,3,4], 1}], :), ...
    best_psth_1s_baseline([neuron_mod_cat{2, [1,3,4], 2}], :), ...zhewan
    best_psth_1s_baseline([neuron_mod_cat{2, [1,3,4], 1}], :)}, ...
    linspace(-1,3,200), ...
    {'Adolescent gamma-modulated', 'Adolescent non-gamma-modulated', 'Adult gamma-modulated', 'Adult non-gamma-modulated'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title('Neuronal Firing Rate')
xlabel('Time from cue onset (s)')
fig_name = ['spiking_firing_rate_by_gamma_modulation'];
yl = ylim;
l1 = line([0, 0.5], [yl(1) + 0.01, yl(1) + 0.01], 'LineWidth',2, 'Color', 'g');
l2 = line([2, 2], yl, 'LineStyle', '--', 'Color', 'k');
set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylim(yl);
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
%%  Gamma power by gamma power modulation 
band_i_ = 3;
mod_plot({lfp_mod_cue([site_mod_cat{1, [1,3,4], 2}], :, band_i_), ...
    lfp_mod_cue([site_mod_cat{1, [1,3,4], 1}], :, band_i_), ...
    lfp_mod_cue([site_mod_cat{2, [1,3,4], 2}], :, band_i_), ...
    lfp_mod_cue([site_mod_cat{2, [1,3,4], 1}], :, band_i_)}, ...
    linspace(-1,3,200), ...
    {'Adolescent gamma-modulated', 'Adolescent non-gamma-modulated', 'Adult gamma-modulated', 'Adult non-gamma-modulated'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title('Gamma Power')
xlabel('Time from cue onset (s)')
fig_name = ['gamma_power_modulation_by_gamma_power_modulation'];
yl = ylim;
l1 = line([0, 0.5], [yl(1) + 0.01, yl(1) + 0.01], 'LineWidth',2, 'Color', 'g');
l2 = line([2, 2], yl, 'LineStyle', '--', 'Color', 'k');
set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylim(yl);
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
%%  Beta power by gamma power modulation 
band_i_ = 2;
mod_plot({lfp_mod_cue([site_mod_cat{1, [1,3,4], 2}], :, band_i_), ...
    lfp_mod_cue([site_mod_cat{1, [1,3,4], 1}], :, band_i_), ...
    lfp_mod_cue([site_mod_cat{2, [1,3,4], 2}], :, band_i_), ...
    lfp_mod_cue([site_mod_cat{2, [1,3,4], 1}], :, band_i_)}, ...
    linspace(-1,3,200), ...
    {'Adolescent gamma-modulated', 'Adolescent non-gamma-modulated', 'Adult gamma-modulated', 'Adult non-gamma-modulated'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title('Beta Power')
xlabel('Time from cue onset (s)')
fig_name = ['gamma_power_modulation_by_gamma_power_modulation'];
yl = ylim;
l1 = line([0, 0.5], [yl(1) + 0.01, yl(1) + 0.01], 'LineWidth',2, 'Color', 'g');
l2 = line([2, 2], yl, 'LineStyle', '--', 'Color', 'k');
set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylim(yl);
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
%%  Gamma PEV by gamma power modulation 
band_i_ = 3;
pev_plot({lfp_pev_cue([site_mod_cat{1, [1,3,4], 2}], :, band_i_), ...
    lfp_pev_cue([site_mod_cat{1, [1,3,4], 1}], :, band_i_), ...
    lfp_pev_cue([site_mod_cat{2, [1,3,4], 2}], :, band_i_), ...
    lfp_pev_cue([site_mod_cat{2, [1,3,4], 1}], :, band_i_)}, ...
    linspace(-1,3,200), ...
    {'Adolescent gamma-modulated', 'Adolescent non-gamma-modulated', 'Adult gamma-modulated', 'Adult non-gamma-modulated'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title('Gamma tuning')
xlabel('Time from cue onset (s)')
fig_name = ['gamma_pev_by_gamma_modulation'];
yl = ylim;
l1 = line([0, 0.5], [yl(1) + 0.01, yl(1) + 0.01], 'LineWidth',2, 'Color', 'g');
l2 = line([2, 2], yl, 'LineStyle', '--', 'Color', 'k');
set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylim([yl(1), yl(2)+0.04]);
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
%%  Neuron PEV by gamma power modulation
pev_plot({neuron_pev_cue([neuron_mod_cat{1, [1,3,4], 2}], :), ...
    neuron_pev_cue([neuron_mod_cat{1, [1,3,4], 1}], :), ...
    neuron_pev_cue([neuron_mod_cat{2, [1,3,4], 2}], :), ...
    neuron_pev_cue([neuron_mod_cat{2, [1,3,4], 1}], :)}, ...
    linspace(-1,3,200), ...
    {'Adolescent gamma-modulated', 'Adolescent non-gamma-modulated', 'Adult gamma-modulated', 'Adult non-gamma-modulated'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title('Neuronal tuning')
xlabel('Time from cue onset (s)')
fig_name = ['neuron_pev_by_gamma_modulation'];
yl = ylim;
l1 = line([0, 0.5], [yl(1) + 0.01, yl(1) + 0.01], 'LineWidth',2, 'Color', 'g');
l2 = line([2, 2], yl, 'LineStyle', '--', 'Color', 'k');
set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylim([yl(1), yl(2)+0.04]);
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
%%  Neuron PEV by gamma power modulation (Only at gamma-modulated site)
pev_plot({neuron_pev_cue([neuron_mod_cat{1, [1,3,4], 2}], :), ...
    neuron_pev_cue([neuron_mod_cat{2, [1,3,4], 2}], :)}, ...
    linspace(-1,3,200), ...
    {'Adolescent gamma-modulated', 'Adult gamma-modulated'}, {'b','r'})
title('Neuronal tuning')
xlabel('Time from cue onset (s)')
fig_name = ['neuron_pev_by_gamma_modulation'];
yl = ylim;
l1 = line([0, 0.5], [yl(1) + 0.01, yl(1) + 0.01], 'LineWidth',2, 'Color', 'g');
l2 = line([2, 2], yl, 'LineStyle', '--', 'Color', 'k');
set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylim([yl(1), yl(2)+0.04]);
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
%%  Best Gamma power modulation by spiking tuning
band_i_ = 3;
mod_plot({lfp_mod_cue([site_tuning_cat{1, [1,3,4], :, 2}], :, band_i_), ...
    lfp_mod_cue([site_tuning_cat{1, [1,3,4], :, 1}], :, band_i_), ...
    lfp_mod_cue([site_tuning_cat{2, [1,3,4], :, 2}], :, band_i_), ...
    lfp_mod_cue([site_tuning_cat{2, [1,3,4], :, 1}], :, band_i_)}, ...
    linspace(-1,3,200), ...
    {'Adolescent informative', 'Adolescent non-informative', 'Adult informative', 'Adult non-informative'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title('Gamma Power')
xlabel('Time from cue onset (s)')
fig_name = ['gamma_power_modulation_by_spiking_tuning'];
yl = ylim;
l1 = line([0, 0.5], [yl(1) + 0.01, yl(1) + 0.01], 'LineWidth',2, 'Color', 'g');
l2 = line([2, 2], yl, 'LineStyle', '--', 'Color', 'k');
set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylim(yl);
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
%%  Gamma PEV by gamma tuning 
band_i_ = 3;
pev_plot({lfp_pev_cue(t_1_t, :, band_i_), ...
    lfp_pev_cue(t_1_nt, :, band_i_), ...
    lfp_pev_cue(t_2_t, :, band_i_), ...
    lfp_pev_cue(t_2_nt, :, band_i_)}, ...
    linspace(-1,3,200), ...
    {'Adolescent gamma-tuned', 'Adolescent non-gamma-tuned', 'Adult gamma-tuned', 'Adult gamma-tuned'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
title('Gamma tuning')
xlabel('Time from cue onset (s)')
fig_name = ['gamma_pev_by_gamma_tuning'];
yl = ylim;
l1 = line([0, 0.5], [yl(1) + 0.01, yl(1) + 0.01], 'LineWidth',2, 'Color', 'g');
l2 = line([2, 2], yl, 'LineStyle', '--', 'Color', 'k');
set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylim([yl(1), yl(2)+0.04]);
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');

%%  beta_tuning
figure();
plot(nanmean(lfp_pev_cue([t13; t14], :, 2), 1))
hold on
plot(nanmean(lfp_pev_cue([t23; t24], :, 2), 1))
%%  gamma_tuning
figure();
plot(nanmean(lfp_pev_cue([t13; t14], :, 3), 1))
hold on
plot(nanmean(lfp_pev_cue([t23; t24], :, 3), 1))
%%  high_gamma_tuning
figure();
plot(nanmean(lfp_pev_cue([t13; t14], :, 4), 1))
hold on
plot(nanmean(lfp_pev_cue([t23; t24], :, 4), 1))
%%
average_neuron_pev_by_site = zeros(size(mapping_mat, 1), size(neuron_pev_cue, 2));
for i = 1:size(average_neuron_pev_by_site, 1)
    n_ = find(mapping_mat(i, :));
    if isempty(n_)
        average_neuron_pev_by_site(i, :) = nan;
        continue
    end
    average_neuron_pev_by_site(i, :) = nanmean(neuron_pev_cue(n_, :), 1);
end
%%
figure()
plot(nanmean(average_neuron_pev_by_site([t13;t14], 40:50), 2), nanmean(lfp_pev_cue([t13;t14], 100:125, 4), 2), '.')
hold on
plot(nanmean(average_neuron_pev_by_site([t23;t24], 40:50), 2), nanmean(lfp_pev_cue([t23;t24], 100:125, 4), 2), '.')
%%
