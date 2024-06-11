%%  RMS distribution
f1 = histogram(reshape(abs_pow(:,:,:,:,1), numel(abs_pow(:,:,:,:,1)), 1), 10.^(2:0.1:4.5));
hold on
f2 = histogram(reshape(abs_pow(:,:,:,:,2), numel(abs_pow(:,:,:,:,2)), 1), 10.^(2:0.1:4.5));
set(gca, 'XScale', 'log')
% set(gca, 'Xticks', [1])

xlabel('RMS')
ylabel('Session count')
legend('Baseline', 'Total')
set(gca, 'FontSize', 14);
set(gca, 'FontWeight', 'bold');
print(gcf, fullfile(project_dir, fig_lib, 'abs_pow_dist'), '-dpng', '-r400')
%%  Spectrogram plotting parameters
% t_stft = -1.25:0.02:3.25;
% t_cwt = -1.5:0.02:3.5;
t_stft = [-1.25, 3.25];
t_cwt = [-1.5, 3.5];
ps_segment_steps = [-1, 0, 1, 1.5, 3];
as_segment_Steps = [-1, 0, 1, 1.1];
t_lim = [-1.1, 3.1];
y_range_stft = [log(100), log(800)];
y_range_cwt = [log(100), log(800)];
%%
function do_plot_(in_cell, tf_method, task_id)
%   tf_method - 0:  stft
%               1:  cwt
if tf_method == 0
    t_for_plot = t_stft;
    y_range = y_range_stft;
elseif tf_method == 1
    t_for_plot = t_cwt;
    y_range = y_cwt;
end

if task_id == 1
    t_step = ps_segment_steps;
elseif task_id == 2
    t_step = as_segment_steps;
end

imagesc(mean(log(in_cell{1}), 3), y_range_);
hold on
xlim(t_lim);
set(gca, 'YDir','normal');
plot()
% colorbar operations
h = colorbar();
h.Label = 'Logarithmic power (a.u.)';
end