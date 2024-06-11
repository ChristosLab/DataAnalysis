%%  Path setting
project_dir = 'F://CCLab//';
code_lib = 'Wang_Code';
fig_lib = 'Wang_Figure';
output_database = 'Wang_Database';
antisac_database = fullfile(project_dir, '//LFPDatabychannel_Antisac//'); % Repository to antisaccade data.
prosac_database = fullfile(project_dir, '//LFPDatabychannel_Prosac//'); % Repository to prosaccade data.
file_identifier = '*CH*.mat';
antisac_files = dir(fullfile(antisac_database, file_identifier));
prosac_files = dir(fullfile(prosac_database, file_identifier));
all_files = [antisac_files; prosac_files];
addpath(fullfile(project_dir, code_lib));
%%  Define target labels and values for entry extraction
channels = unique([LFP_repo.channel_id]);
tasks = unique([LFP_repo.task_id]);
sessions = unique([LFP_repo.session_id]);
monkeys = unique([LFP_repo.monkey_id]);
groups = cell([numel(channels), numel(tasks), numel(sessions), numel(monkeys)]);
entry_labels = {'channel_id', 'task_id', 'session_id', 'monkey_id'};
%%  Load trial id groups
load(fullfile(project_dir, output_database, 'LFP_repo_w_pretrial.mat'));
load(fullfile(project_dir, output_database, 'groups_5_5.mat'));
%%  Extract rows of LFP_repo according to entry_labels
% field_names = fieldnames(LFP_repo);
% LFP_repo_lite = rmfield(LFP_repo, field_names(~ismember(field_names, entry_labels)));

LFP_repo_lite = [];
for e_l = entry_labels
    LFP_repo_lite = [LFP_repo_lite, [LFP_repo.(e_l{1})]'];
end

%%  CMorlet wavelet parameters
interval_post = [-1, 2];
interval_pre = [-1, 0];
f_range = 2:2:100;
fs = 500;
%%
spectrogram_change = zeros([size(groups), numel(f_range), diff(interval_post)*fs]);
tic
for a = 1:size(groups, 1)
    for b = 1:size(groups, 2)
        for c = 1:size(groups, 3)
            c
            for d = 1:size(groups, 4)
                if ~any(groups{a, b, c ,d})
                    spectrogram_change(a, b, c, d, :, :) = nan;
                    continue;
                end 
                spectrogram_change(a, b, c, d, :, :) = zw_group_power_change(...
                    spectrogram_change(a, b, c, d, :, :), LFP_repo(groups{a, b, c ,d}), ...
                    interval_post, interval_pre, f_range, fs...
                    );
            end
            toc
        end
    end
end
%%  Save spectogram_change file with pre and post intervals as name.
fname_ = sprintf('spectrogram_change_%d_%d_%d_%d.mat',interval_pre(1),interval_pre(2),interval_post(1),interval_post(2));
save(fullfile(project_dir, output_database, fname_), 'spectrogram_change');
%%
monkey_id = {'ind', 'jac', 'lem', 'ken'};

