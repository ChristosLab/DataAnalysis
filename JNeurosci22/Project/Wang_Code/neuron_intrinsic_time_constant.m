%   [tau, exit_code]
tau_repo = nan(numel(neuron_repo), 2);
rho_mean = nan(numel(neuron_repo), 12);
bin_width = 50;
max_lag = 600;
tar_dur = 1:20;
tic
for i = [neuron_mod_cat{:}]
    input_mat = cat(1, neuron_repo(i).class.psth_cue);
    [tau_repo(i, 1), tau_repo(i, 2), rho_mean(i, :)] = zw_spike_time_constant(input_mat(:, tar_dur), bin_width, max_lag, 1);
%     if ~isnan(tau_repo(i, 1))
%         pause
%     end
    toc
end
%%  Young, all traces
figure('Units', 'inches','Position', [2,2,4,3])
hold on
for i = [neuron_mod_cat{1,:,:}]
    if any(isnan(rho_mean(i, :)))
        continue
    end
    l_ = plot(50:50:600, rho_mean(i, :), '-k');
    l_.Color(4) = 0.2;
end
l_m_ = plot(50:50:600, nanmean(rho_mean([neuron_mod_cat{1,:,:}], :), 1), 'b-d', 'LineWidth', 1.3, 'DisplayName','Adolescent mean');
ylim([-0.2, 0.6])
ylabel('Autocorrelation')
xlabel('Lag (ms)')
legend(l_m_);
set(gca, {'FontSize', 'LineWidth'}, {10, 1.2})
fig_name = 'time_constant_all_traces_young';
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
%%  Adult, all traces
figure('Units', 'inches','Position', [2,2,4,3])
hold on
for i = [neuron_mod_cat{2,:,:}]
    if any(isnan(rho_mean(i, :)))
        continue
    end
    l_ = plot(50:50:600, rho_mean(i, :), '-k');
    l_.Color(4) = 0.2;
end
l_m_ = plot(50:50:600, nanmean(rho_mean([neuron_mod_cat{2,:,:}], :), 1), 'r-d', 'LineWidth', 1.3, 'DisplayName','Adult mean');
ylim([-0.2, 0.6])
ylabel('Autocorrelation')
xlabel('Lag (ms)')
legend(l_m_);
set(gca, {'FontSize', 'LineWidth'}, {10, 1.2})
fig_name = 'time_constant_all_traces_adult';
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
%%
n_y = sum(any(~isnan(rho_mean([neuron_mod_cat{1,:,:}], :)), 2));
n_a = sum(any(~isnan(rho_mean([neuron_mod_cat{2,:,:}], :)), 2));
figure('Units', 'inches','Position', [2,2,4,3])
hold on
errorbar(50:50:600, nanmean(rho_mean([neuron_mod_cat{1,:,:}], :), 1), ...
    nanstd(rho_mean([neuron_mod_cat{1,:,:}], :), 1)/sqrt(n_y), ...
    'b-o', 'MarkerSize', 4, 'LineWidth', 1.3, 'DisplayName', sprintf('Adolescent, N = %d', n_y))
errorbar(50:50:600, nanmean(rho_mean([neuron_mod_cat{2,:,:}], :), 1), ...
    nanstd(rho_mean([neuron_mod_cat{2,:,:}], :), 1)/sqrt(n_a), ...
    'r-o', 'MarkerSize', 4, 'LineWidth', 1.3, 'DisplayName', sprintf('Adult, N = %d', n_a))
ylim([0, 0.2])
ylabel('Autocorrelation')
xlabel('Lag (ms)')
legend
set(gca, {'FontSize', 'LineWidth'}, {10, 1.2})
fig_name = 'time_constant_mean_traces';
print(gcf, fullfile(project_dir, fig_lib, fig_name), '-dpdf');
%%
input_mat = [];
for i = [neuron_mod_cat{:}]
    input_mat_ = cat(1, neuron_repo(i).class.psth_cue);
    input_mat = [input_mat; input_mat_];
end
%%
input_mat_y = [];
for i = [neuron_mod_cat{1,:,:,:}]
    input_mat_ = cat(1, neuron_repo(i).class.psth_cue);
    input_mat_y = [input_mat_y; input_mat_];
end
%%
input_mat_a = [];
for i = [neuron_mod_cat{2,:,:,:}]
    input_mat_ = cat(1, neuron_repo(i).class.psth_cue);
    input_mat_a = [input_mat_a; input_mat_];
end
%%
[tau_, exit_code_] = zw_spike_time_constant(input_mat(:, tar_dur), bin_width, max_lag);
%%
[tau_y, exit_code_, rho_mean_y, ci_y] = zw_spike_time_constant(input_mat_y(:, tar_dur), bin_width, max_lag, 0);
[tau_a, exit_code_, rho_mean_a, ci_a] = zw_spike_time_constant(input_mat_a(:, tar_dur), bin_width, max_lag, 0);
