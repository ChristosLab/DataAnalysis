% t_cwt_1 = [-2.49, -0.51]; % Time for fixation onset
% t_cwt_2 = [-1.11, 2.49]; % Time for fixation onset

%%
% for subj_to_plot = 1:3
%     figure();
%     imagesc(t_cwt_1, [2, 100],mean(cwt_young_ps{subj_to_plot, 1}(2:end, :, :), 3), cwt_y);set(gca, 'YDir','normal'); xlim(t_lim);colorbar()
%     plot_lines(task_id, subj_to_plot);
%     
%     figure();
%     imagesc(t_cwt_1, [2, 100],mean(cwt_adult_ps{subj_to_plot, 1}(2:end, :, :), 3), cwt_y);set(gca, 'YDir','normal'); xlim(t_lim);colorbar()
%     plot_lines(task_id, subj_to_plot);
% 
%     figure();
%     imagesc(t_cwt_1, [2, 100],mean(cwt_young_as{subj_to_plot, 1}(2:end, :, :), 3), cwt_y);set(gca, 'YDir','normal'); xlim(t_lim);colorbar()
%     plot_lines(task_id, subj_to_plot);
% 
%     
%     figure();
%     imagesc(t_cwt_1, [2, 100],mean(cwt_adult_as{subj_to_plot, 1}(2:end, :, :), 3), cwt_y);set(gca, 'YDir','normal'); xlim(t_lim);colorbar()
%     plot_lines(task_id, subj_to_plot);
%     
% end
%%  Subject to plot/Proportion of spec to plot
subj_to_plot = 2;
sample_spec = [0, 1];
%%  Average as
task_i = 2;
bridge_length = .2;
plot_cwt(...
    cwt_young_as, ...
    subj_to_plot, ...
    bridge_length, task_i, sample_spec...
    )
title('Young-AS');
print(gcf, fullfile(project_dir, fig_lib, sprintf('ya_subj_%1.0d_5_15', subj_to_plot)), '-dpng', '-r400')
%
plot_cwt(...
    cwt_adult_as, ...
    subj_to_plot, ...
    bridge_length, task_i, sample_spec...
    )
title('Adult-AS');
print(gcf, fullfile(project_dir, fig_lib, sprintf('aa_subj_%1.0d_5_15', subj_to_plot)), '-dpng', '-r400')

%%  Average ps
task_i = 1;
bridge_length = .2;
plot_cwt(...
    cwt_young_ps, ...
    subj_to_plot, ...
    bridge_length, task_i, sample_spec...
    )
title('Young-PS');
print(gcf, fullfile(project_dir, fig_lib, sprintf('yp_subj_%1.0d_5_15', subj_to_plot)), '-dpng', '-r400')

plot_cwt(...
    cwt_adult_ps, ...
    subj_to_plot, ...
    bridge_length, task_i, sample_spec...
    )
title('Adult-PS');
print(gcf, fullfile(project_dir, fig_lib, sprintf('ap_subj_%1.0d_5_15', subj_to_plot)), '-dpng', '-r400')

%%  Session by session examination
in_cwt = cwt_adult_as;
task_i = 2;
for i = 1:size(in_cwt{subj_to_plot, 1}, 3)
    sample_spec = [i, 1];
    plot_cwt(...
        in_cwt, ...
        subj_to_plot, ...
        bridge_length, task_i, sample_spec...
        )
    i
    pause;
    close all
end
%%  Session extraction
target_sessions = [12,19,32,93,134];
for i = target_sessions
    sample_spec = [i, 1];
    plot_cwt(...
        in_cwt, ...
        subj_to_plot, ...
        bridge_length, task_i, sample_spec...
        )
    print(gcf, fullfile(project_dir, fig_lib, sprintf('aa_subj_%1.0d_session_%1.0d_5_15', subj_to_plot, i)), '-dpng', '-r400')
end
%%  K-means clustering within condition
cluster_k = 5;
%
in_cwt = cwt_young_as;
task_i = 2;
% 

cluster_cwt = [in_cwt{subj_to_plot, 1}, in_cwt{subj_to_plot, 2}];
in_size = size(cluster_cwt);
cluster_cwt = reshape(cluster_cwt, [prod(in_size([1,2])), in_size(3)])'; % Each row is 1 observation
cluster_cwt = cluster_cwt./std(cluster_cwt, [], 2);
idx = kmeans(cluster_cwt, cluster_k);
%
for k = 1:cluster_k
    for i = find(idx == k)'
        sample_spec = [i, 1];
        plot_cwt(...
            in_cwt, ...
            subj_to_plot, ...
            bridge_length, task_i, sample_spec...
            )
        sprintf('Cluster No.: %d', k)
        i
        pause;
        close all
    end
end
%%
function kmeans_session(in_cwt, subj_to_plot)
end
%%
function plot_cwt(in_cwt, subj_to_plot, bridge_length, task_i, sample_spec)
cwt_1 = in_cwt{subj_to_plot, 1};
cwt_2 = in_cwt{subj_to_plot, 2};
if all(sample_spec) % 0 values results in the whole signal 
    cwt_1 = cwt_1(:, : , sample_spec(1):(sample_spec(2) + sample_spec(1) - 1));
    cwt_2 = cwt_2(:, : , sample_spec(1):(sample_spec(2) + sample_spec(1) - 1));
% if diff(sample_portion) ~= 1
%     cwt_1 = cwt_1(:, : , 1 + floor(end*sample_portion(1)):floor(end*sample_portion(2)));
%     cwt_2 = cwt_2(:, : , 1 + floor(end*sample_portion(1)):floor(end*sample_portion(2)));
end
cwt_y = [1000, 16000]./(500^2);
pad = 0.04;
fs = 50;
t_lim = [-1.99 - pad, 2.49 + bridge_length];
fig_size = [5 1 4.5 3];
figure('Units', 'inches','Position', fig_size);
hold on
t_cwt_1 = [-1.99 - pad, -1.01 + pad]; % Time for fixation onset
t_bridge = [-1.0 + pad, -1.0 + bridge_length - 0.02 - pad];
t_cwt_2 = [-1.01 - pad, 2.49] + bridge_length; % Time for fixation onset
range_1 = (26 - pad*fs):(75 + pad*fs);
range_2 = (6 - pad*fs):181;
to_plot_1 = nanmean(cwt_1(2:end, range_1, :), 3)./(500^2);
to_plot_2 = nanmean(cwt_2(2:end, range_2, :), 3)./(500^2);
imagesc(t_cwt_1, [2, 100], to_plot_1, cwt_y);
imagesc(t_cwt_2, [2, 100], to_plot_2, cwt_y);
fill([t_bridge, fliplr(t_bridge)], [2,2,100,100],[.5,.5,.5], 'LineStyle' , 'none')
% h = colorbar();
set(gca, 'YDir','normal'); 
xlim(t_lim);
ylim([2, 100]);

ps_segment_steps = [-2, -1, -1 + bridge_length, 0 + bridge_length, 0.5 + bridge_length, 2 + bridge_length];
ps_xticks = [...
    -2, -1.5, -1, -1 + bridge_length, -.5 + bridge_length, 0 + bridge_length, ...
    0.25 + bridge_length, 0.5 + bridge_length, 1 + bridge_length, ...
    1.25 + bridge_length, 1.5 + bridge_length, 2 + bridge_length...
    ];
ps_xlabels = {...
    '', 'baseline', '', '-1', 'prep', '0', ...
    'cue', '.5', '', ...
    'delay', '', '2'};
as_segment_steps = [-2, -1, -1 + bridge_length, 0 + bridge_length, 0.1 + bridge_length];
as_xticks = [...
    -2, -1.5, -1 , -1 + bridge_length, -.5 + bridge_length, 0 + bridge_length, ...
    0.1+ bridge_length...
    ];
as_xlabels = {...
    '', 'baseline', '', '-1', 'prep', '0', ''...
    };
if task_i == 1
    steps = ps_segment_steps;
    xt = ps_xticks;
    x_lab = ps_xlabels;
else
    steps = as_segment_steps;
    xt = as_xticks;
    x_lab = as_xlabels;
end
for i = 1:numel(steps)
    line([0,0] + steps(i), [1, 101], 'Color', 'k', 'LineStyle', '--');    
end
h = colorbar();
xticks(xt);
xticklabels(x_lab);
xlabel('Time from cue onset (s)')
ylabel('Frequency (Hz)')
ylabel(h, 'Power (a.u.)');
% 
set(gca, 'FontSize', 9);
set(gca, 'FontWeight', 'bold');
end