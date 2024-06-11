%%  Path setting
zw_setpath
%%  Load LFP_repo
tic
load(fullfile(project_dir, output_database, 'LFP_repo_w_pretrial.mat'));
toc
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
save(fullfile(project_dir, output_database, 'LFP_repo_lite.mat'), 'LFP_repo_lite');
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
save(fullfile(project_dir, output_database, 'groups_9_14.mat'), 'groups', 'entry_labels', 'entry_label_values', 'LFP_repo_lite', 'group_ids');
%%
%
% 
% 
% 
% 
%%
sub2ind_w([4,2,3], [2,1,1])
sub2ind([4,2,3], 2,1,1)
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
