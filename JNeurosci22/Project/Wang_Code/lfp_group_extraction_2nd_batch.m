%%  Path setting
zw_setpath
%%  Load LFP_repo
tic
load(fullfile(project_dir, output_database, 'lfp_repo_2.mat'));
toc
%%
LFP_repo = lfp_repo_2;
clear lfp_repo_2;
%%  Define target labels and values for entry extraction
entry_labels = {'channel_id', 'task_id', 'session_id', 'monkey_id', 'class'}; 
entry_label_values = cell(size(entry_labels));
%%  Extract rows of LFP_repo according to entry_labels
LFP_repo_lite = [];
group_dims = zeros(1, numel(entry_labels));
for el_index = 1:numel(entry_labels)
    LFP_repo_lite = [LFP_repo_lite, [LFP_repo.(entry_labels{el_index})]'];
    entry_label_values{el_index} = unique([LFP_repo.(entry_labels{el_index})]);
    group_dims(1, el_index) = numel(entry_label_values{el_index});
end
groups = cell(group_dims);
%%
save(fullfile(project_dir, output_database, 'LFP_repo_2_lite.mat'), 'LFP_repo_lite');
%%
group_ids = ndgrid_w(group_dims);
%%  Group trial ids by the labels and save in cell arrays
tic
for group_index = 1:size(group_ids, 1)
    groups{sub2ind_w(group_dims, group_ids(group_index, :))} = ismember(...
        LFP_repo_lite, ...
        ref_arrays_in_cell(entry_label_values, group_ids(group_index, :)), ...
        'rows'...
        );
    if mod(group_index, 100) < 1
        sprintf('%d/%d', group_index, size(group_ids, 1))
        toc
    end
               
end
%%
save(fullfile(project_dir, output_database, 'groups_lfp_repo_2.mat'), 'groups', 'entry_labels', 'entry_label_values', 'LFP_repo_lite', 'group_ids');
%%
addpath(fullfile(project_dir, code_lib, 'External\fieldtrip\'));
ft_defaults
addpath(genpath(fullfile(project_dir, code_lib, 'External\chronux_2_12\')));
%%  CWT parameters
pre_dur = 0.5;
post_dur = 3.5;
normalizer = 3;
kernel_flag = 3;
down_sample = 10;
fs = 500;
f_range = 2:2:100;
avg_method = 1; %   0 - arithmetic; 1 - geometric
%%  Compute single trial post cwt
seg_size = 10000;
for i = 1:ceil(size(LFP_repo_lite, 1)/seg_size)
    tt = clock;
    seg_start = 1 + (i - 1)*seg_size;
    seg_end = min(seg_start + seg_size -1, size(LFP_repo_lite, 1));
    fprintf(...
        'Segment %d to %d started at %s\n', ...
        seg_start, seg_end, string(datetime)...
        );
    repo_cwt = zw_repo_cwt_no_pre(...
        LFP_repo(1, seg_start:seg_end), ...
        post_dur, normalizer, ...
        kernel_flag, down_sample, f_range, fs...
        );
    save(...
        fullfile(...
        project_dir, ...
        output_database, ...
        sprintf('repo_cwt_post_only_%d_to_%d.mat', seg_start, seg_end)...
        ), ...
        'repo_cwt'...
        );
    repo_cwt = [];
end
%%

%%
function output = ref_arrays_in_cell(input_cell, array_sub)
output = zeros(1, numel(input_cell));
for i = 1:numel(input_cell)
    output(i) = input_cell{i}(array_sub(i));
end
end
%
function out_grid = ndgrid_w(dims)
%NDGRID_W wraps the NDGRID function and outputs permutations of 1-based
%indexing in a single matrix
cols_ = cell(1, numel(dims)); % Output 
dims_ = cell(1, numel(dims)); % Input 
for i = 1:numel(dims)
    dims_{i} = 1:dims(i);
end
[cols_{:}] = ndgrid(dims_{:});
for i = 1:numel(dims)
    cols_{i} = reshape(cols_{i}, numel(cols_{i}), 1);
end
out_grid = [cols_{:}];
end
