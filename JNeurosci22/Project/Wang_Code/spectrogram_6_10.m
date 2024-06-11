clear
%%  Path setting
project_dir = 'F:\CCLab\Projects\LFP_project_v0.2';
code_lib = 'Wang_Code';
fig_lib = 'Wang_Figure';
output_database = 'Wang_Database';
addpath(fullfile(project_dir, code_lib));
addpath(genpath(fullfile(project_dir, code_lib, '/External/chronux_2_12')));
disp('Paths added');
datetime
%%  Load trial id groups
load(fullfile(project_dir, output_database, 'LFP_repo_w_pretrial.mat'));
load(fullfile(project_dir, output_database, 'groups_5_5.mat'));
disp('Repos loaded ...')
datetime
%%  Define target labels and values for entry extraction
channels = unique([LFP_repo.channel_id]);
tasks = unique([LFP_repo.task_id]);
sessions = unique([LFP_repo.session_id]);
monkeys = unique([LFP_repo.monkey_id]);
entry_labels = {'channel_id', 'task_id', 'session_id', 'monkey_id'};
%%  Extract rows of LFP_repo according to entry_labels
% field_names = fieldnames(LFP_repo);
% LFP_repo_lite = rmfield(LFP_repo, field_names(~ismember(field_names, entry_labels)));

LFP_repo_lite = [];
for e_l = entry_labels
    LFP_repo_lite = [LFP_repo_lite, [LFP_repo.(e_l{1})]'];
end
%%  CWT parameters
pre_dur = 0.5;
post_dur = 3.5;
normalizer = 3;
kernel_flag = [2, 3];
down_sample = 10;
fs = 500;
f_range = 2:2:100;
avg_method = 1; %   0 - arithmetic; 1 - geometric
%%  CWT computation
for i = 1:2
    cwt_1{i} = zeros([size(groups), numel(f_range), fs * pre_dur/down_sample]);
    cwt_2{i} = zeros([size(groups), numel(f_range), fs * post_dur/down_sample]);
end

tic
for a = 1:size(groups, 1)
    for b = 1:size(groups, 2)
        for c = 1:size(groups, 3)
            c
            for d = 1:size(groups, 4)
                if ~any(groups{a, b, c ,d})
                    cwt_1{1}(a, b, c, d, :, :) = nan;
                    cwt_2{1}(a, b, c, d, :, :) = nan;
                    cwt_1{2}(a, b, c, d, :, :) = nan;
                    cwt_2{2}(a, b, c, d, :, :) = nan;
                    continue;
                end
                % 
                [temp_1, temp_2] = zw_group_cwt_ver_2(...
                    cwt_1{1}(a, b, c, d, :, :), cwt_2{1}(a, b, c, d, :, :), ...
                    LFP_repo(groups{a, b, c ,d}), ...
                    pre_dur, post_dur, normalizer, ...
                    kernel_flag(1), down_sample, f_range, fs, ...
                    avg_method...
                    );
                cwt_1{1}(a, b, c, d, :, :) = temp_1;
                cwt_2{1}(a, b, c, d, :, :) = temp_2;
                %  
                [temp_1, temp_2] = zw_group_cwt_ver_2(...
                    cwt_1{2}(a, b, c, d, :, :), cwt_2{2}(a, b, c, d, :, :), ...
                    LFP_repo(groups{a, b, c ,d}), ...
                    pre_dur, post_dur, normalizer, ...
                    kernel_flag(2), down_sample, f_range, fs, ...
                    avg_method...
                    );
                cwt_1{2}(a, b, c, d, :, :) = temp_1;
                cwt_2{2}(a, b, c, d, :, :) = temp_2; 
            end
            toc
        end
    end
end
%%  Save CWT
disp('Saving...')
datetime
fname_ = 'cwt_pre_cohen_chronux_v2.mat';
save(fullfile(project_dir, output_database, fname_), 'cwt_1');
%
fname_ = 'cwt_post_cohen_chronux_v2.mat';
save(fullfile(project_dir, output_database, fname_), 'cwt_2');
%
fname_ = 'global_params_v2.mat';
save(...
fullfile(project_dir, output_database, fname_), ...
'f_range', 'pre_dur', 'post_dur', 'kernel_flag', ...
'avg_method', 'normalizer', 'down_sample', 'sessions', ...
'channels', 'tasks', 'monkeys', 'entry_labels'...
);
