t_cwt = [-2.5, 2.5]; % Time for fixation onset
ps_segment_steps = [-1, 0, 1, 1.5, 3] - 1;
as_segment_Steps = [-1, 0, 1, 1.1] - 1;
t_lim = [-1.1, 3.1];

for subj_to_plot = 1:3
    figure();
    imagesc(t_cwt, [2, 100],mean(cwt_young_ps{subj_to_plot}(2:end, :, :), 3), cwt_y);set(gca, 'YDir','normal'); xlim(t_lim);colorbar()
    plot_lines(task_id, subj_to_plot);
    
    figure();
    imagesc(t_cwt, [2, 100],mean(cwt_adult_ps{subj_to_plot}(2:end, :, :), 3), cwt_y);set(gca, 'YDir','normal'); xlim(t_lim);colorbar()
    plot_lines(task_id, subj_to_plot);

    figure();
    imagesc(t_cwt, [2, 100],mean(cwt_young_as{subj_to_plot}(2:end, :, :), 3), cwt_y);set(gca, 'YDir','normal'); xlim(t_lim);colorbar()
    plot_lines(task_id, subj_to_plot);

    
    figure();
    imagesc(t_cwt, [2, 100],mean(cwt_adult_as{subj_to_plot}(2:end, :, :), 3), cwt_y);set(gca, 'YDir','normal'); xlim(t_lim);colorbar()
    plot_lines(task_id, subj_to_plot);
    
end
function plot_lines(task_i, subj_i)
if task_i == 1
    steps = ps_segment_steps;
else
    steps = as_segment_steps;
end
for i = 1:numel(steps)
    line([0,0] + steps(i), [2, 100], '--k');    
end
xlabel('Time from fixation onset (s)')
xlim(t_lim);
ylabel('Frequency (Hz)')
h = colorbar();
h.Label = 'Power (a.u.)';
% 
end
