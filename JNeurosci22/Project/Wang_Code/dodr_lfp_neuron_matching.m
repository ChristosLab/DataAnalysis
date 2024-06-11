%   Load informative vs. non-informative neuron flags for the 1 - 2 second
%   epoch. Script for generating this file is not saved by mistake. Refer
%   to 'Whole_delay_period_neuron_tuning.m' for general pattern
fname_ = 'dodr_neuron_epoch_delay_tuning_flag.mat';
load(fullfile(project_dir, output_database, fname_), 'neuron_epoch_delay_tuning_flag', 'neuron_epoch_analysis_output_cue');
fname_ = 'dodr_valid_lfp.mat';
load(fullfile(project_dir, output_database, fname_), 'valid_lfp');
%%
responsive = any(table2array(dodr_neuron_tbl(:, 9:10)), 2);
dodr_mapping_mat = zeros(size(dodr_lfp_tbl, 1), size(dodr_neuron_tbl, 1));
%   Irregularities in data identified from manual inspection
% further_site_mask = ones(size(dodr_lfp_tbl, 1), 1);
% further_site_mask([19, 35, 106, 109, 118, 119, 130]) = 0;
for i = 1:size(dodr_lfp_tbl, 1)
%     i
    dodr_mapping_mat(i, :) = prod(table2array(dodr_lfp_tbl(i, 2:5)) == table2array(dodr_neuron_tbl(:, [6, 7, 8, 3])), 2).*responsive;
end
dodr_mapping_mat = dodr_mapping_mat.*valid_lfp;
%%  Indexing each stage-by-animal in 't' (LFP) and 'n' (neuron)
t = cell(2, 4);
n = cell(2, 4);
for i_stage = 1:2
    for i_monk = 1:4
        t{i_stage, i_monk} = find(any((dodr_lfp_tbl.stage == i_stage).*(dodr_lfp_tbl.monkey_id == i_monk).*dodr_mapping_mat, 2))';
        n{i_stage, i_monk} = find(any((dodr_lfp_tbl.stage == i_stage).*(dodr_lfp_tbl.monkey_id == i_monk).*dodr_mapping_mat, 1));
    end
end
%%  Neuron tuning categorization (spiking).
temp_counts_ = {};
dodr_neuron_tuning_cat = cell(2, 4, 2);
for i_stage = 1:2
    for i_monk = 1:4
        %               stage,   monkey, gamma tuning, spiking tuning
        dodr_neuron_tuning_cat{i_stage, i_monk, 1} = intersect(n{i_stage, i_monk}, find(any(dodr_mapping_mat.*~neuron_epoch_delay_tuning_flag, 1)));        
        dodr_neuron_tuning_cat{i_stage, i_monk, 2} = intersect(n{i_stage, i_monk}, find(any(dodr_mapping_mat.*neuron_epoch_delay_tuning_flag, 1)));
    end
end
%%  Site tuning categorization (spiking)
dodr_site_tuning_cat = cell(2, 4, 2);
for i_stage = 1:2
    for i_monk = 1:4
        %               stage,   monkey, gamma tuning, spiking tuning
        dodr_site_tuning_cat{i_stage, i_monk, 2} = intersect(t{i_stage, i_monk}, find(any(dodr_mapping_mat.*neuron_epoch_delay_tuning_flag, 2))');
        dodr_site_tuning_cat{i_stage, i_monk, 1} = setxor(t{i_stage, i_monk}, [dodr_site_tuning_cat{i_stage, i_monk, 2}]);        
%         dodr_site_tuning_cat{i_stage, i_monk, 1} = intersect(t{i_stage, i_monk}, find(any(dodr_mapping_mat.*~neuron_epoch_delay_tuning_flag, 2))');        
    end
end
%%
save(fullfile(project_dir, output_database, 'dodr_lfp_neuron_matching.mat'), 't', 'n', 'dodr_mapping_mat', 'dodr_site_tuning_cat', 'dodr_neuron_tuning_cat');
%%  Neuron Gamma power modulation categorization. (Not yet processed)
dodr_neuron_mod_cat = cell(2, 4, 2);
for i_stage = 1:2
    for i_monk = 1:4
        %               stage,   monkey,   delay gamma mod
        dodr_neuron_mod_cat{i_stage, i_monk, 1} = intersect(n{i_stage, i_monk}, find(any(dodr_mapping_mat.*~tfr_delay_modulation_flag(:, 3, 1), 1)));
        dodr_neuron_mod_cat{i_stage, i_monk, 2} = intersect(n{i_stage, i_monk}, find(any(dodr_mapping_mat.*tfr_delay_modulation_flag(:, 3, 1), 1)));
    end
end
%%  Site Gamma power modulation categorization.
dodr_site_mod_cat = cell(2, 4, 2);
for i_stage = 1:2
    for i_monk = 1:4
        %               stage,   monkey, delay gamma mod
        dodr_site_mod_cat{i_stage, i_monk, 1} = intersect(t{i_stage, i_monk}, find(any(dodr_mapping_mat.*~tfr_delay_modulation_flag(:, 3, 1), 2))');
        dodr_site_mod_cat{i_stage, i_monk, 2} = intersect(t{i_stage, i_monk}, find(any(dodr_mapping_mat.*tfr_delay_modulation_flag(:, 3, 1), 2))');
        dodr_site_mod_cat{i_stage, i_monk, 3} = intersect(t{i_stage, i_monk}, find(any(dodr_mapping_mat.*~tfr_delay_modulation_flag(:, 4, 1), 2))');
        dodr_site_mod_cat{i_stage, i_monk, 4} = intersect(t{i_stage, i_monk}, find(any(dodr_mapping_mat.*tfr_delay_modulation_flag(:, 4, 1), 2))');        
    end
end
%%  Compute session mean band power for plotting
for i = 1:size(target_frs, 1)
    lfp_mod_cue(:, :, i) = squeeze(nanmean(temp_cwt_cue_b(:, target_frs(i, 1):target_frs(i, 2), :), 2));
end
%%  Compute class mean band power for plotting
for i = 1:size(target_frs, 1)
    lfp_class_mod_cue(:, :, :, i) = squeeze(nanmean(temp_class_cwt_cue_b(:, :, target_frs(i, 1):target_frs(i, 2), :), 3));
end
%%  Band power according to neuron preferred class
title_st = {'Alpha (8-16 Hz)', 'Beta (16-32 Hz)', 'Gamma (32-64 Hz)', 'High-Gamma (64-128 Hz)'};
% n_y = [n{1, 3}];
% n_a = [n{2, 3}];
n_y = [n{1, :}];
n_a = [n{2, :}];
% n_y = [dodr_neuron_tuning_cat{1, :, 2}];
% n_a = [dodr_neuron_tuning_cat{2, :, 2}];
n_matched_site_y = zeros(size(n_y));
n_matched_site_a = zeros(size(n_a));
for i = 1:numel(n_y)
    n_matched_site_y(i) = find(dodr_mapping_mat(:, n_y(i)));
end
for i = 1:numel(n_a)
    n_matched_site_a(i) = find(dodr_mapping_mat(:, n_a(i)));
end
% dodr_neuron_best_class_1s_baseline_target
n_matched_lfp_mod_y    = zeros(numel(n_matched_site_y), size(lfp_class_mod_cue, 3), size(lfp_class_mod_cue, 4));
n_matched_lfp_mod_y_op = n_matched_lfp_mod_y;
n_matched_lfp_mod_a    = zeros(numel(n_matched_site_a), size(lfp_class_mod_cue, 3), size(lfp_class_mod_cue, 4));
n_matched_lfp_mod_a_op = n_matched_lfp_mod_a;

for i = 1:numel(n_matched_site_y)
    current_site     = n_matched_site_y(i);
    current_class    = dodr_neuron_best_class_1s_baseline_target(n_y(i));
    if current_class ~= 4
        current_class_op = mod(current_class + 4, 8);
    else
        current_class_op = current_class + 4;
    end
    n_matched_lfp_mod_y(i, :, :)    = lfp_class_mod_cue(current_site, current_class, :, :);
    n_matched_lfp_mod_y_op(i, :, :) = lfp_class_mod_cue(current_site, current_class_op, :, :);
end
% 
for i = 1:numel(n_matched_site_a)
    current_site     = n_matched_site_a(i);
    current_class    = dodr_neuron_best_class_1s_baseline_target(n_a(i));
    if current_class ~= 4
        current_class_op = mod(current_class + 4, 8);
    else
        current_class_op = current_class + 4;
    end
    n_matched_lfp_mod_a(i, :, :)    = lfp_class_mod_cue(current_site, current_class, :, :);
    n_matched_lfp_mod_a_op(i, :, :) = lfp_class_mod_cue(current_site, current_class_op, :, :);
end
%%  Power modulation between stages
for band_i_ = 1:4
    dodr_mod_plot({lfp_mod_cue([t{1, :, :}], :, band_i_), ...
        lfp_mod_cue([t{2, :, :}], :, band_i_)}, ...
        linspace(-1,4,250), ...
        {'Adolescent', 'Adult'}, {'b', 'r'})
    title(title_st{band_i_})
    xlabel('Time from cue onset (s)')
    xlim([-1, 3])
    % ylim([0.7, 2.4])
    if band_i_ ~= 2
        legend off
    end
    fig_name = ['dodr_lfp_pow_band_', num2str(band_i_), 'power_modulation_by_stage', '_draft'];
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',10,'FontWeight','Bold', 'LineWidth', 1);
    print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
end
%%  Power modulation to neuron's preferred and non-preferred class
for band_i_ = 1:4
    dodr_mod_plot({n_matched_lfp_mod_y(:, :, band_i_), ...
        n_matched_lfp_mod_y_op(:, :, band_i_), ...
        n_matched_lfp_mod_a(:, :, band_i_), ...
        n_matched_lfp_mod_a_op(:, :, band_i_)}, ...
        linspace(-1,4,250), ...
        {'Ado p', 'Ado d', 'Adu p', 'Adu d'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
    %     {'Adolescent preferred', 'Adolescent diametric', 'Adult preferred', 'Adult diametric'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
    title(title_st{band_i_})
    xlabel('Time from cue onset (s)')
    % ylim([0.7, 2.4])
    xlim([-1, 3])
    if band_i_ ~= 2
        legend off
    end
    fig_name = ['dodr_lfp_pow_band_', num2str(band_i_), 'power_modulation_by_neuron_preference', '_draft'];
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',10,'FontWeight','Bold', 'LineWidth', 1);
    print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
end
%%  Power modulation difference to neuron's preferred and non-preferred class
for band_i_ = 1:4
mod_plot({n_matched_lfp_mod_y(:, :, band_i_) - ...
    n_matched_lfp_mod_y_op(:, :, band_i_), ...
    n_matched_lfp_mod_a(:, :, band_i_) - ...
    n_matched_lfp_mod_a_op(:, :, band_i_)}, ...
    linspace(-1,4,250), ...
    {'Ado p - d', 'Adu p - d'}, {'b', 'r'}, [-.5, .5] ...
    )
title(title_st{band_i_})
ylabel('\Delta relative power')
xlabel('Time from cue onset (s)')
xlim([-1, 3])
p_ = plot([-1, 3], [0, 0], 'Color', [.3, .3, .3], 'LineWidth', 1);
set(get(get(p_,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
if band_i_ ~= 2
    legend off
end
fig_name = ['dodr_lfp_pow_band_', num2str(band_i_), 'power_modulation_neuron_preference_diff', '_draft'];
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',10,'FontWeight','Bold', 'LineWidth', 1);
    print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
end
%%  PSTH to neuron's preferred and non-preferred class (1-second baseline corrected)
psth_plot({best_psth_1s_baseline_target(n_y, :), ...
    best_psth_1s_baseline_opp(n_y, :), ...
    best_psth_1s_baseline_target(n_a, :), ...
    best_psth_1s_baseline_opp(n_a, :)}, ...
    linspace(-1,4,250), ...
    {'Adolescent p', 'Adolescent d', 'Adult p', 'Adult d'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
xlabel('Time from cue onset (s)')
xlim([-1, 3]);
fig_name = ['dodr_psth_draft'];
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',10,'FontWeight','Bold', 'LineWidth', 1);
    print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
%%  PSTH to neuron's preferred and non-preferred class (raw)
psth_plot({best_psth_raw_target(n_y, :), ...
    best_psth_raw_opp(n_y, :), ...
    best_psth_raw_target(n_a, :), ...
    best_psth_raw_opp(n_a, :)}, ...
    linspace(-1,4,250), ...
    {'Adolescent p', 'Adolescent d', 'Adult p', 'Adult d'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
xlabel('Time from cue onset (s)')
xlim([-1, 3]);
fig_name = ['dodr_psth_draft'];
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',10,'FontWeight','Bold', 'LineWidth', 1);
    print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
%%
figure
hold on
plot(mean(lfp_mod_cue([dodr_site_tuning_cat{1,:,1}], :, 3), 1))
plot(mean(lfp_mod_cue([dodr_site_tuning_cat{1,:,2}], :, 3), 1))
plot(mean(lfp_mod_cue([dodr_site_tuning_cat{2,:,1}], :, 3), 1))
plot(mean(lfp_mod_cue([dodr_site_tuning_cat{2,:,2}], :, 3), 1))
%%
figure
hold on
plot(mean(lfp_mod_cue([dodr_site_tuning_cat{1,:,1}], :, 4), 1))
plot(mean(lfp_mod_cue([dodr_site_tuning_cat{1,:,2}], :, 4), 1))
plot(mean(lfp_mod_cue([dodr_site_tuning_cat{2,:,1}], :, 4), 1))
plot(mean(lfp_mod_cue([dodr_site_tuning_cat{2,:,2}], :, 4), 1))
%%
figure
hold on
plot(mean(lfp_mod_cue([dodr_site_tuning_cat{1,:,1}], :, 2), 1))
plot(mean(lfp_mod_cue([dodr_site_tuning_cat{1,:,2}], :, 2), 1))
plot(mean(lfp_mod_cue([dodr_site_tuning_cat{2,:,1}], :, 2), 1))
plot(mean(lfp_mod_cue([dodr_site_tuning_cat{2,:,2}], :, 2), 1))
%%  Power modulation by spiking tuning
for band_i_ = 1:4
mod_plot({lfp_mod_cue([dodr_site_tuning_cat{1, :, 2}], :, band_i_), ...
    lfp_mod_cue([dodr_site_tuning_cat{1, :, 1}], :, band_i_), ...
    lfp_mod_cue([dodr_site_tuning_cat{2, :, 2}], :, band_i_), ...
    lfp_mod_cue([dodr_site_tuning_cat{2, :, 1}], :, band_i_)}, ...
    linspace(-1,4,250), ...
    {'Adolescent informative', 'Adolescent non-informative', 'Adult informative', 'Adult non-informative'}, {'b', [0.2, 0.2, 0.6],'r', [0.6, 0.2, 0.2]})
% title(title_st{band_i_})
xlabel('Time from cue onset (s)')
ylim([0.7, 2.4])
if band_i_ ~= 2
    legend off
end
fig_name = ['dodr_lfp_pow_band_', num2str(band_i_), 'power_modulation_by_spiking_tuning', '_draft'];
% print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
end
%%
for band_i_ = 1:4
mod_plot({lfp_mod_cue(setxor([t{1, :}], suspect_n_y), :, band_i_), ...
    lfp_mod_cue(setxor([t{2, :}], suspect_n_a), :, band_i_)}, ...
    linspace(-1,4,250), ...
     {'Adolescent', 'Adult'}, {'b','r'});
end
%%  Check single site relative spectrogram
for i = [t{2,:}]
    c_range = [-1, 1];
    figure
    imagesc([-1, 4], [2,128], squeeze(temp_cwt_cue_b(i, :, :) - 1)*100, c_range*100)
    set(gca, 'YDir','normal'); 
    colormap(clm_)
    xlabel('Time from cue onset (s)')
    ylabel('Frequency (Hz)')
    h = colorbar();
    h.Ruler.TickLabelFormat='%g%%';
    ylabel(h, '% Power of baseline')
    set_plot_12_4();
    title(num2str(i))
    pause;
    close gcf
end