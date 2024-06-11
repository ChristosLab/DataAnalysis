%%  Adolescent
%  Sort by neuron PEV peaks, then plot 1) Norm. PEV; 2) Norm. gamma power; 3) Norm. high-gamma power
target_subset = [neuron_tuning_cat{1, [1,3,4],:, 2}];
min_maxed_pev_ = (neuron_pev_cue(target_subset, :) - min(neuron_pev_cue(target_subset, :), [], 2))./(max(neuron_pev_cue(target_subset, :), [], 2) - min(neuron_pev_cue(target_subset, :), [], 2));

peak_timepoint = zeros(numel(target_subset), 1);
for i = 1:numel(target_subset)
    [~, i_] = max(min_maxed_pev_(i, 76:200));
    peak_timepoint(i) = i_;
%         peak_timepoint(i) = zw_center_of_mass(min_maxed_pev_(i, 76:150));

end
[sorted_peak_timepoint, trace_order] = sort(peak_timepoint);


figure('Unit', 'inches', 'Position', [2, 2, 2.5, 2])
imagesc([-1, 3], 1:numel(trace_order), min_maxed_pev_(trace_order, :), [0, 1]);
% colormap('jet')
h = colorbar;
ylabel(h, 'Norm. Spiking PEV');
xlabel('Time from cue onset (s)')
hold on
% plot(sorted_peak_timepoint+75, 1:numel(target_subset), 'r');

% 
band_i = 3;
mod_to_plot = lfp_mod_cue(map_neuron_to_site(mapping_mat, target_subset(trace_order)), :, band_i) - 1; %    Converting to change from baseline
mod_to_plot = mod_to_plot./max(mod_to_plot, [], 2);
ylabel('Neuron #')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',10,'FontWeight','Bold', 'LineWidth', 1);
fig_name = ['adolescent_pev_trace_draft'];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');


figure('Unit', 'inches', 'Position', [2, 2, 2.5, 2])
imagesc([-1, 3], 1:numel(trace_order), mod_to_plot, [-1, 1])
colormap(clm_)
h = colorbar;
ylabel(h, 'Norm. Gamma Power');
ylabel('Neuron #')
xlabel('Time from cue onset (s)')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',10,'FontWeight','Bold', 'LineWidth', 1);
fig_name = ['adolescent_gamma_mod_trace_draft'];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');

%
band_i = 4;
mod_to_plot = lfp_mod_cue(map_neuron_to_site(mapping_mat, target_subset(trace_order)), :, band_i) - 1; %    Converting to change from baseline
mod_to_plot = mod_to_plot./max(mod_to_plot, [], 2);
figure('Unit', 'inches', 'Position', [2, 2, 2.5, 2])
imagesc([-1, 3], 1:numel(trace_order), mod_to_plot, [-1, 1])
colormap(clm_)
h = colorbar;
ylabel(h, 'Norm. High Gamma Power');
ylabel('Neuron #')
xlabel('Time from cue onset (s)')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',10,'FontWeight','Bold', 'LineWidth', 1);
fig_name = ['adolescent_high_gamma_mod_trace_draft'];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
%%  Adult
%  Sort by neuron PEV peaks, then plot 1) Norm. PEV; 2) Norm. gamma power; 3) Norm. high-gamma power
target_subset = [neuron_tuning_cat{2, [1,3,4],:, 2}];
min_maxed_pev_ = (neuron_pev_cue(target_subset, :) - min(neuron_pev_cue(target_subset, :), [], 2))./(max(neuron_pev_cue(target_subset, :), [], 2) - min(neuron_pev_cue(target_subset, :), [], 2));

peak_timepoint = zeros(numel(target_subset), 1);
for i = 1:numel(target_subset)
    [~, i_] = max(min_maxed_pev_(i, 76:200));
    peak_timepoint(i) = i_;
%         peak_timepoint(i) = zw_center_of_mass(min_maxed_pev_(i, 76:150));

end
[sorted_peak_timepoint, trace_order] = sort(peak_timepoint);


figure('Unit', 'inches', 'Position', [2, 2, 2.5, 2])
imagesc([-1, 3], 1:numel(trace_order), min_maxed_pev_(trace_order, :), [0, 1]);
% colormap('jet')
h = colorbar;
ylabel(h, 'Norm. Spiking PEV');
xlabel('Time from cue onset (s)')
hold on
% plot(sorted_peak_timepoint+75, 1:numel(target_subset), 'r');

% 
band_i = 3;
mod_to_plot = lfp_mod_cue(map_neuron_to_site(mapping_mat, target_subset(trace_order)), :, band_i) - 1; %    Converting to change from baseline
mod_to_plot = mod_to_plot./max(mod_to_plot, [], 2);
ylabel('Neuron #')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',10,'FontWeight','Bold', 'LineWidth', 1);
fig_name = ['adult_pev_trace_draft'];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');


figure('Unit', 'inches', 'Position', [2, 2, 2.5, 2])
imagesc([-1, 3], 1:numel(trace_order), mod_to_plot, [-1, 1])
colormap(clm_)
h = colorbar;
ylabel(h, 'Norm. Gamma Power');
ylabel('Neuron #')
xlabel('Time from cue onset (s)')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',10,'FontWeight','Bold', 'LineWidth', 1);
fig_name = ['adult_gamma_mod_trace_draft'];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');

%
band_i = 4;
mod_to_plot = lfp_mod_cue(map_neuron_to_site(mapping_mat, target_subset(trace_order)), :, band_i) - 1; %    Converting to change from baseline
mod_to_plot = mod_to_plot./max(mod_to_plot, [], 2);
figure('Unit', 'inches', 'Position', [2, 2, 2.5, 2])
imagesc([-1, 3], 1:numel(trace_order), mod_to_plot, [-1, 1])
colormap(clm_)
h = colorbar;
ylabel(h, 'Norm. High Gamma Power');
ylabel('Neuron #')
xlabel('Time from cue onset (s)')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',10,'FontWeight','Bold', 'LineWidth', 1);
fig_name = ['adult_high_gamma_mod_trace_draft'];
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');