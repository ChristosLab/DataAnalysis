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
%%  Plot spectrum
figure; hold on
plot(log(2:2:128), squeeze(mean(mean(log(temp_cwt_cue([site_tuning_cat{1, :, 1:2,:}], :, 26:50)), 3), 1)))
plot(log(2:2:128), squeeze(mean(mean(log(temp_cwt_cue([site_tuning_cat{2, :, 1:2,:}], :, 26:50)), 3), 1)))
%%
figure; hold on
plot(log(2:2:128), squeeze(mean(mean(log(temp_cwt_cue([site_tuning_cat{1, :, 1:2,:}], :, 51:75)), 3), 1)))
plot(log(2:2:128), squeeze(mean(mean(log(temp_cwt_cue([site_tuning_cat{2, :, 1:2,:}], :, 51:75)), 3), 1)))
%%
figure; hold on
plot(log(2:2:128), squeeze(mean(mean(log(temp_cwt_cue([site_tuning_cat{1, :, 1:2,:}], :, 76:150)), 3), 1)))
plot(log(2:2:128), squeeze(mean(mean(log(temp_cwt_cue([site_tuning_cat{2, :, 1:2,:}], :, 76:150)), 3), 1)))
%%  Correlation between delay period evoked firing and LFP power by frequency (Full set)
font_size = 10;
clc
f_range = 2:2:128;
% corr_method_ = 'Pearson';
corr_method_ = 'Spearman';
r_adolescent = zeros(1, numel(f_range));
r_adult      = zeros(1, numel(f_range));
p_adolescent = zeros(1, numel(f_range));
p_adult      = zeros(1, numel(f_range));
p_comp_corr  = zeros(1, numel(f_range));
for i = 1:numel(f_range)
[r_adolescent(i), p_adolescent(i)] = corr(mean(best_psth_1s_baseline([neuron_tuning_cat{1, [1,3,4], :, :}], 76:150), 2),  mean(temp_cwt_cue_b(map_neuron_to_site(mapping_mat, [neuron_tuning_cat{1, [1,3,4], :, :}]), i, 76:150), 3), 'type', corr_method_);
[r_adult(i), p_adult(i)] = corr(mean(best_psth_1s_baseline([neuron_tuning_cat{2, [1,3,4], :, :}], 76:150), 2),  mean(temp_cwt_cue_b(map_neuron_to_site(mapping_mat, [neuron_tuning_cat{2, [1,3,4], :, :}]), i, 76:150), 3), 'type', corr_method_);

% disp(newline)
p_comp_corr(i) = compare_correlation_coefficients(r_adolescent(i), r_adult(i), numel([neuron_tuning_cat{1, [1,3,4], :, :}]), numel([neuron_tuning_cat{2, [1,3,4], :, :}]));
% disp(newline)
end
%
figure('Unit', 'inches', 'Position', [2, 2, 3, 2.5]);
; hold on
plot(f_range, r_adolescent, 'LineWidth', 1.2, 'Color', 'b');
plot(f_range, r_adult, 'LineWidth', 1.2, 'Color', 'r');
[~, h_] = bonf_holm(p_adolescent, 0.05);
plot(f_range(h_), h_(h_) - 0.99, '.', 'MarkerSize', 5, 'Color', [0.2, 0.2, 0.6]);
[~, h_] = bonf_holm(p_adult, 0.05);
plot(f_range(h_), h_(h_) - 1.01, '.', 'MarkerSize', 5, 'Color', [0.6, 0.2, 0.2]);
[~, h_] = bonf_holm(p_comp_corr, 0.05);
plot(f_range(h_), h_(h_) - 0.4, '*k');%%
xlabel('Frequency (Hz)')
ylabel('Correlation coefficient')
xlim([2, 128])
ylim([-0.2, 0.6])
set_plot_poster
title('Firing rate vs. LFP power')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',font_size,'FontWeight','Bold', 'LineWidth', 1);
fig_name = 'n_by_n_corr_fr_lfp_full';
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
savefig(gcf, fullfile(project_dir, fig_lib, fig_name));
%%  Correlation between delay period neuronal PEV and LFP power by frequency (Full set)
clc
f_range = 2:2:128;
% corr_method_ = 'Pearson';
corr_method_ = 'Spearman';
r_adolescent = zeros(1, numel(f_range));
r_adult      = zeros(1, numel(f_range));
p_adolescent = zeros(1, numel(f_range));
p_adult      = zeros(1, numel(f_range));
p_comp_corr  = zeros(1, numel(f_range));
for i = 1:numel(f_range)
[r_adolescent(i), p_adolescent(i)] = corr(mean(neuron_pev_cue([neuron_tuning_cat{1, [1,3,4], :, :}], 76:150), 2),  mean(temp_cwt_cue_b(map_neuron_to_site(mapping_mat, [neuron_tuning_cat{1, [1,3,4], :, :}]), i, 76:150), 3), 'type', corr_method_);
[r_adult(i), p_adult(i)] = corr(mean(neuron_pev_cue([neuron_tuning_cat{2, [1,3,4], :, :}], 76:150), 2),  mean(temp_cwt_cue_b(map_neuron_to_site(mapping_mat, [neuron_tuning_cat{2, [1,3,4], :, :}]), i, 76:150), 3), 'type', corr_method_);

% disp(newline)
p_comp_corr(i) = compare_correlation_coefficients(r_adolescent(i), r_adult(i), numel([neuron_tuning_cat{1, [1,3,4], :, :}]), numel([neuron_tuning_cat{2, [1,3,4], :, :}]));
% disp(newline)
end
% 
figure('Unit', 'inches', 'Position', [2, 2, 3, 2.5]);
; hold on
plot(f_range, r_adolescent, 'LineWidth', 1.2, 'Color', 'b');
plot(f_range, r_adult, 'LineWidth', 1.2, 'Color', 'r');
[~, h_] = bonf_holm(p_adolescent, 0.05);
plot(f_range(h_), h_(h_) - 0.99, '.', 'MarkerSize', 5, 'Color', [0.2, 0.2, 0.6]);
[~, h_] = bonf_holm(p_adult, 0.05);
plot(f_range(h_), h_(h_) - 1.01, '.', 'MarkerSize', 5, 'Color', [0.6, 0.2, 0.2]);
[~, h_] = bonf_holm(p_comp_corr, 0.05);
plot(f_range(h_), h_(h_) - 0.4, '*k');%%
xlabel('Frequency (Hz)')
ylabel('Correlation coefficient')
xlim([2, 128])
ylim([-0.2, 0.6])
set_plot_poster
title('Neuron PEV vs. LFP power')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',font_size,'FontWeight','Bold', 'LineWidth', 1);
fig_name = 'n_by_n_corr_pev_lfp_full';
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
savefig(gcf, fullfile(project_dir, fig_lib, fig_name));
% 
% 
% 
% 
% 
%%  Correlation between delay period evoked firing and LFP power by frequency (Reduced set)
clc
f_range = 2:2:128;
% corr_method_ = 'Pearson';
corr_method_ = 'Spearman';
r_adolescent = zeros(1, numel(f_range));
r_adult      = zeros(1, numel(f_range));
p_adolescent = zeros(1, numel(f_range));
p_adult      = zeros(1, numel(f_range));
p_comp_corr  = zeros(1, numel(f_range));
for i = 1:numel(f_range)
[r_adolescent(i), p_adolescent(i)] = corr(mean(best_psth_1s_baseline(redueced_adolescent_set, 76:150), 2),  mean(temp_cwt_cue_b(map_neuron_to_site(mapping_mat, redueced_adolescent_set), i, 76:150), 3), 'type', corr_method_);
[r_adult(i), p_adult(i)] = corr(mean(best_psth_1s_baseline(redueced_adult_set, 76:150), 2),  mean(temp_cwt_cue_b(map_neuron_to_site(mapping_mat, redueced_adult_set), i, 76:150), 3), 'type', corr_method_);

% disp(newline)
p_comp_corr(i) = compare_correlation_coefficients(r_adolescent(i), r_adult(i), numel([neuron_tuning_cat{1, [1,3,4], :, :}]), numel([neuron_tuning_cat{2, [1,3,4], :, :}]));
% disp(newline)
end
%
figure('Unit', 'inches', 'Position', [2, 2, 3, 2.5]);
; hold on
plot(f_range, r_adolescent, 'LineWidth', 1.2, 'Color', 'b');
plot(f_range, r_adult, 'LineWidth', 1.2, 'Color', 'r');
[~, h_] = bonf_holm(p_adolescent, 0.05);
h_1 = h_;
plot(f_range(h_), h_(h_) - 0.99, '.', 'MarkerSize', 5, 'Color', [0.2, 0.2, 0.6]);
[~, h_] = bonf_holm(p_adult, 0.05);
h_2 = h_;
plot(f_range(h_), h_(h_) - 1.01, '.', 'MarkerSize', 5, 'Color', [0.6, 0.2, 0.2]);
[~, h_] = bonf_holm(p_comp_corr, 0.05);
h_ = and(h_, or(h_1, h_2));
plot(f_range(h_), h_(h_) - 0.4, '*k');%%
xlabel('Frequency (Hz)')
ylabel('Correlation coefficient')
xlim([2, 128])
ylim([-0.2, 0.6])
title('Firing rate vs. LFP power (reduced set)')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',font_size,'FontWeight','Bold', 'LineWidth', 1);

fig_name = 'n_by_n_corr_fr_lfp_reduced';
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
savefig(gcf, fullfile(project_dir, fig_lib, fig_name));%%  Correlation between delay period neuronal PEV and LFP power by frequency (Reduced set)
%%
clc
f_range = 2:2:128;
% corr_method_ = 'Pearson';
corr_method_ = 'Spearman';
r_adolescent = zeros(1, numel(f_range));
r_adult      = zeros(1, numel(f_range));
p_adolescent = zeros(1, numel(f_range));
p_adult      = zeros(1, numel(f_range));
p_comp_corr  = zeros(1, numel(f_range));
for i = 1:numel(f_range)
[r_adolescent(i), p_adolescent(i)] = corr(mean(neuron_pev_cue(redueced_adolescent_set, 76:150), 2),  mean(temp_cwt_cue_b(map_neuron_to_site(mapping_mat, redueced_adolescent_set), i, 76:150), 3), 'type', corr_method_);
[r_adult(i), p_adult(i)] = corr(mean(neuron_pev_cue(redueced_adult_set, 76:150), 2),  mean(temp_cwt_cue_b(map_neuron_to_site(mapping_mat, redueced_adult_set), i, 76:150), 3), 'type', corr_method_);

% disp(newline)
p_comp_corr(i) = compare_correlation_coefficients(r_adolescent(i), r_adult(i), numel([neuron_tuning_cat{1, [1,3,4], :, :}]), numel([neuron_tuning_cat{2, [1,3,4], :, :}]));
% disp(newline)
end
% 
figure('Unit', 'inches', 'Position', [2, 2, 3, 2.5]);
; hold on
plot(f_range, r_adolescent, 'LineWidth', 1.2, 'Color', 'b');
plot(f_range, r_adult, 'LineWidth', 1.2, 'Color', 'r');
[~, h_] = bonf_holm(p_adolescent, 0.05);
plot(f_range(h_), h_(h_) - 0.99, '.', 'MarkerSize', 5, 'Color', [0.2, 0.2, 0.6]);
[~, h_] = bonf_holm(p_adult, 0.05);
plot(f_range(h_), h_(h_) - 1.01, '.', 'MarkerSize', 5, 'Color', [0.6, 0.2, 0.2]);
[~, h_] = bonf_holm(p_comp_corr, 0.05);
plot(f_range(h_), h_(h_) - 0.4, '*k');%%
xlabel('Frequency (Hz)')
ylabel('Correlation coefficient')
xlim([2, 128])
ylim([-0.2, 0.6])
set_plot_poster
title('Neuron PEV vs. LFP power (reduced set)')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',font_size,'FontWeight','Bold', 'LineWidth', 1);
fig_name = 'n_by_n_corr_pev_lfp_reduced';
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpng', '-r400');
savefig(gcf, fullfile(project_dir, fig_lib, fig_name));