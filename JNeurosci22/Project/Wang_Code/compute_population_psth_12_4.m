close all
clear
zw_setpath
load(fullfile(project_dir, output_database, 'neuron_tbl_w_pfc.mat'), 'neuron_tbl');
load(fullfile(project_dir, output_database, 'complete_neuron_repo_w_psth.mat'), 'neuron_repo');
%%
temp_psth = zeros([numel(neuron_repo), 8, 80]); % [neuron, class, time]
temp_psth_raw = zeros([numel(neuron_repo), 8, 80]);
temp_psth_1s_baseline = zeros([numel(neuron_repo), 8, 80]);
for i = 1:numel(neuron_repo)
    if ~(neuron_tbl.task_id(i) == 1) || ~(neuron_tbl.PFC(i) == 1)
        temp_psth_raw(i, :) = nan;
        temp_psth(i, :) = nan;
        continue
    end
    for j = 1:8
        temp_psth(i, j, :) = nanmean(neuron_repo(i).class(j).psth_cue - mean(neuron_repo(i).class(j).psth_cue(:, 11:20), 2), 1);
        temp_psth_1s_baseline(i, j, :) = nanmean(neuron_repo(i).class(j).psth_cue - mean(neuron_repo(i).class(j).psth_cue(:, 1:20), 2), 1);
        temp_psth_raw(i, j, :) = nanmean(neuron_repo(i).class(j).psth_cue, 1);

    end
end
%%
neuron_best_class_raw = zeros(size(temp_psth_raw, 1), 1);
for i = 1:size(temp_psth_raw, 1)
    [y_, i_] = max(sum(temp_psth_raw(i, :, 21:60), 3));
    neuron_best_class_raw(i) = i_;
end
%%
best_psth_raw = zeros(size(temp_psth_raw, 1), size(temp_psth_raw, 3));
for i = 1:size(best_psth_raw, 1)
    best_psth_raw(i, :) = temp_psth_raw(i, neuron_best_class_raw(i), :);
end
%%
neuron_best_class = zeros(size(temp_psth, 1), 1);
for i = 1:size(temp_psth, 1)
    [y_, i_] = max(sum(temp_psth(i, :, 21:60), 3));
    neuron_best_class(i) = i_;
end
%%
best_psth = zeros(size(temp_psth, 1), size(temp_psth, 3));
for i = 1:size(best_psth, 1)
    best_psth(i, :) = temp_psth(i, neuron_best_class(i), :);
end
%%
neuron_best_class_1s_baseline = zeros(size(temp_psth, 1), 1);
for i = 1:size(temp_psth, 1)
    [y_, i_] = max(sum(temp_psth_1s_baseline(i, :, 21:60), 3));
    neuron_best_class_1s_baseline(i) = i_;
end
%%
best_psth_1s_baseline = zeros(size(temp_psth, 1), size(temp_psth, 3));
for i = 1:size(best_psth, 1)
    best_psth_1s_baseline(i, :) = temp_psth(i, neuron_best_class_1s_baseline(i), :);
end
%%
save(fullfile(project_dir, output_database, 'best_psth.mat'), ...
    'temp_psth', 'best_psth', 'neuron_best_class', ...
    'temp_psth_1s_baseline', 'best_psth_1s_baseline', 'neuron_best_class_1s_baseline', ...
    'temp_psth_raw', 'best_psth_raw', 'neuron_best_class_raw');