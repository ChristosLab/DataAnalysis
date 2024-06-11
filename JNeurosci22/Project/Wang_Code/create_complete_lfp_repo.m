%% Create an LFP repo with all available LFP files
%%
clear
zw_setpath;
%%  Create LFP table
%%  Preallocate structure field names.
field_names = {'Filename', 'monkey_id', 'session_id', 'task_id', 'channel_id', 'stage', 'lfp_filename'};
%%
monkey_id = {'ind', 'jac', 'lem', 'ken'};
last_young_session = {20, 26, 66, 27};
stage_cutoff = struct('ind', 20,'jac', 26, 'lem', 66, 'ken', 27);
%  Loops through each MAT file.
excluding_flag = 'catch';
tic
n_file = numel(lfp_files);
lfp_repo_2 = cell2struct(cell(size(field_names)), field_names, 2);
n_valid = 0;
for ind = 1:n_file
    if regexp(lfp_files(ind).name, excluding_flag)
        continue
    end
    n_valid = n_valid + 1;
    lfp_repo_2(n_valid) = add_trials(lfp_repo_2(1), lfp_files(ind), stage_cutoff);
    toc
    sprintf('%d/%d\n', ind, n_file)
end
%%  Save lfp tbl
lfp_tbl = struct2table(lfp_repo_2);
fname_ = 'lfp_tbl.mat';
save(fullfile(project_dir, output_database, fname_), 'lfp_tbl');
%%  Load lfp tbl
fname_ = 'lfp_tbl.mat';
load(fullfile(project_dir, output_database, fname_));
%%  Create new lfp_repo
for i = 1:size(lfp_tbl, 1)
    load(fullfile(lfp_database, lfp_tbl.lfp_filename{i}));
    lfp_repo(i) = LFPData;
    clear LFPData
    i
end
%%
fname_ = 'lfp_repo_complete.mat';
save(fullfile(project_dir, output_database, fname_), 'lfp_repo');
%%  Functions
function in_struct = add_trials(in_struct, file_dir, stage_cutoff)
%ADD_TRIALs 
%     Utility function that adds each trial in the file
%     to the trial-by-trial dataset.
%     in_struct: the trial-by-trial dataset to be written.
%     file_dir: Single file directory structure created with DIR().

monkey_id = {'ind', 'jac', 'lem', 'ken'};

data_struct_name = 'LFPData';
file_ = fullfile(file_dir.folder, file_dir.name);
assert(... %    Check for variable naming in the .MAT file.
    ismember(who('-file', file_), data_struct_name), ...
    sprintf('Data strcuture name not found: %s', data_struct_name)...
    );
%   Adds file-specific info
in_struct = add_file_info(...
    in_struct, file_dir...
    );
%   Add stage infomation - 1 young -2 adult
if in_struct.session_id > stage_cutoff.(in_struct.monkey_id)
    in_struct.stage = 2;
else
    in_struct.stage = 1;
end
in_struct.Filename = file_dir.name(1:end-8);
in_struct.lfp_filename = file_dir.name;
in_struct.monkey_id = find(ismember(monkey_id, in_struct.monkey_id));
end
%%

%%
function in_struct = add_file_info(in_struct, file_dir)
%ADD_FILE_INFO 
%     Utility function that adds information in the file name
%     to the trial-by-trial dataset.
%     in_struct: the trial-by-trial dataset to be written.
%     file_dir: Single file directory structure created with DIR(). Feeds
%     into ADD_TRIAL_INFO

[field_names, field_values] = read_lfp_fname(file_dir.name); % 'name' field referenced.
% out_struct = in_struct;
for j = 1:numel(field_names)
    value_ = field_values{j};
    if ischar(value_)
        if all(ismember(value_, '1234567890'))
            value_ = str2double(value_);
        end
    end
    in_struct.(field_names{j})= value_;
end       
end
%%
function [field_names, field_values] = read_lfp_fname(fname)
%READ_LFP_FNAME
%     Utility function that deliminate file name strings in LFP saving 
%     conventions. Feeds into ADD_FILE_INFO
    
field_names = {'monkey_id', 'session_id', 'task_id', 'channel_id'};
field_values = textscan(fname, '%s%s%s', 'Delimiter', '_'); % Deliminate FNAME using "_".
ch_scan = textscan(field_values{3}{1}, '%2c%f');
field_values = [textscan(field_values{1}{1}, '%3c%f'), field_values{2}{1}, ch_scan{2}]; % Seperate SUBJECT_ID and SESSION_ID.
end